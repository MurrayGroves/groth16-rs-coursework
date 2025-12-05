use crate::helpers::ark_de;
use crate::helpers::ark_se;
use crate::polynomial::Polynomial;
use ark_ff::FftField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use log::debug;
use rand::Rng;
use rootcause::{Report, report};
use serde::{Deserialize, Serialize};
use std::iter::zip;

/// Represents a Rank 1 Constraint System. Should be created using `R1CS::new(...)`,
/// which lets you provide matrices with any type that can be converted into the Scalar type.
/// (E.g. to allow vec literals)
#[derive(Clone, Debug)]
pub struct R1CS<S: FftField> {
    /// Column-wise, i.e. a vec of columns
    pub L: Vec<Vec<S>>,
    /// Column-wise, i.e. a vec of columns
    pub R: Vec<Vec<S>>,
    /// Column-wise, i.e. a vec of columns
    pub O: Vec<Vec<S>>,
    pub public_witness: Vec<S>,
}

impl<S: FftField> R1CS<S> {
    /// Create a new R1CS from matrices in **column-major** order.
    pub fn new<T, W>(l: Vec<Vec<T>>, r: Vec<Vec<T>>, o: Vec<Vec<T>>, public_witness: Vec<W>) -> Self
    where
        S: From<T> + From<W>,
        T: Copy,
        W: Copy,
    {
        R1CS {
            L: l.iter()
                .map(|column| column.iter().map(|x| S::from(*x)).collect())
                .collect(),
            R: r.iter()
                .map(|column| column.iter().map(|x| S::from(*x)).collect())
                .collect(),
            O: o.iter()
                .map(|column| column.iter().map(|x| S::from(*x)).collect())
                .collect(),
            public_witness: public_witness.iter().map(|x| S::from(*x)).collect(),
        }
    }

    pub(crate) fn verify(&self, witness: &Vec<S>) -> Result<bool, Report> {
        let o = zip(&self.O, witness)
            .map(|(o, w)| o.iter().map(|x| *x * *w).collect::<Vec<_>>())
            .reduce(|a, b| zip(a, b).map(|(a_i, b_i)| a_i + b_i).collect())
            .ok_or(report!("Empty vec"))?;
        let l = zip(&self.L, witness)
            .map(|(o, w)| o.iter().map(|x| *x * *w).collect::<Vec<_>>())
            .reduce(|a, b| zip(a, b).map(|(a_i, b_i)| a_i + b_i).collect())
            .ok_or(report!("Empty vec"))?;
        let r = zip(&self.R, witness)
            .map(|(o, w)| o.iter().map(|x| *x * *w).collect::<Vec<_>>())
            .reduce(|a, b| zip(a, b).map(|(a_i, b_i)| a_i + b_i).collect())
            .ok_or(report!("Empty vec"))?;

        debug!("{:?} == {:?} * {:?}", o, l, r);
        let rhs = zip(l, r).map(|(a_i, b_i)| a_i * b_i).collect::<Vec<_>>();

        Ok(o == rhs)
    }
}

/// Represents a Quadratic Arithmetic Program. Cannot be instantiated directly, should instead be derived from a Rank 1 Constraint System using `QAP::from(r1cs)`
#[derive(
    Debug, PartialEq, Eq, Clone, Serialize, Deserialize, CanonicalDeserialize, CanonicalSerialize,
)]
pub struct QAP<S>
where
    S: FftField,
{
    /// Output
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub w: Vec<Polynomial<S>>,
    /// LHS of multiplication
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub u: Vec<Polynomial<S>>,
    /// RHS of multiplication
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub v: Vec<Polynomial<S>>,
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub public_witness: Vec<S>,
}

impl<S: FftField> QAP<S> {
    /// A QAP has degree `n` where `n` is the number of rows in the R1CS it was formed from
    pub fn degree(&self) -> usize {
        self.max_polynomial_degree() + 1
    }

    pub fn max_polynomial_degree(&self) -> usize {
        vec![
            self.u.iter().map(|x| x.degree()).max().unwrap_or(0),
            self.v.iter().map(|x| x.degree()).max().unwrap_or(0),
            self.w.iter().map(|x| x.degree()).max().unwrap_or(0),
        ]
        .iter()
        .max()
        .map(|x| *x)
        .unwrap_or(0)
    }

    pub(crate) fn verify(&self, witness: &Vec<S>) -> bool {
        if witness.len() != self.u.len()
            || witness.len() != self.v.len()
            || witness.len() != self.w.len()
        {
            return false;
        }
        let mut rng = rand::rng();
        let tau = S::from(rng.random_range(0..1000));

        let a: Polynomial<S> = zip(&self.u, witness).map(|(u_i, a_i)| u_i * *a_i).sum();

        let b: Polynomial<S> = zip(&self.v, witness).map(|(v_i, a_i)| v_i * *a_i).sum();

        let w: Polynomial<S> = zip(&self.w, witness).map(|(w_i, a_i)| w_i * *a_i).sum();

        let ht = &(&a * &b) - &w;

        a.evaluate(&tau) * b.evaluate(&tau) == w.evaluate(&tau) + ht.evaluate(&tau)
    }
}

impl<S: FftField> From<R1CS<S>> for QAP<S> {
    fn from(r1cs: R1CS<S>) -> Self {
        QAP {
            u: r1cs
                .L
                .iter()
                .map(Polynomial::interpolate_from_vector)
                .collect(),
            v: r1cs
                .R
                .iter()
                .map(Polynomial::interpolate_from_vector)
                .collect(),
            w: r1cs
                .O
                .iter()
                .map(Polynomial::interpolate_from_vector)
                .collect(),
            public_witness: r1cs.public_witness,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::circuits::{QAP, R1CS};
    use crate::polynomial::Polynomial;
    use ark_ff::{Fp64, MontBackend};
    use log::debug;
    use rand::Rng;
    use rootcause::Report;

    #[derive(ark_ff::MontConfig)]
    #[modulus = "641"]
    #[generator = "3"]
    struct FieldConfig;
    type Field = Fp64<MontBackend<FieldConfig, 1>>;
    #[test]
    fn r1cs_to_qap() {
        // Test case from https://risencrypto.github.io/R1CSQAP/
        let L = vec![
            vec![0, 0, 0, 5],
            vec![0, 0, 0, 0],
            vec![1, 0, 1, 0],
            vec![0, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![0, 0, 0, 1],
        ]
        .iter()
        .map(|col| {
            col.iter()
                .map(|element| Field::from(*element))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

        let R = vec![
            vec![0, 0, 1, 1],
            vec![0, 0, 0, 0],
            vec![1, 1, 0, 0],
            vec![0, 0, 0, 0],
            vec![0, 0, 0, 0],
            vec![0, 0, 0, 0],
        ]
        .iter()
        .map(|col| {
            col.iter()
                .map(|element| Field::from(*element))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

        let O = vec![
            vec![0, 0, 0, 0],
            vec![0, 0, 0, 1],
            vec![0, 0, 0, 0],
            vec![1, 0, 0, 0],
            vec![0, 1, 0, 0],
            vec![0, 0, 1, 0],
        ]
        .iter()
        .map(|col| {
            col.iter()
                .map(|element| Field::from(*element))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

        let r1cs = R1CS {
            L,
            R,
            O,
            public_witness: Vec::new(),
        };

        let qap = QAP::from(r1cs);

        let known_good: QAP<Field> = QAP {
            u: vec![
                Polynomial::from(vec![636, 116, 636, 535]),
                Polynomial::from(vec![0, 0, 0, 0]),
                Polynomial::from(vec![8, 416, 5, 213]),
                Polynomial::from(vec![635, 330, 637, 321]),
                Polynomial::from(vec![4, 634, 324, 320]),
                Polynomial::from(vec![640, 536, 640, 107]),
            ],
            v: vec![
                Polynomial::from(vec![3, 529, 323, 427]),
                Polynomial::from(vec![0, 0, 0, 0]),
                Polynomial::from(vec![639, 112, 318, 214]),
                Polynomial::from(vec![0, 0, 0, 0]),
                Polynomial::from(vec![0, 0, 0, 0]),
                Polynomial::from(vec![0, 0, 0, 0]),
            ],
            w: vec![
                Polynomial::from(vec![0, 0, 0, 0]),
                Polynomial::from(vec![640, 536, 640, 107]),
                Polynomial::from(vec![0, 0, 0, 0]),
                Polynomial::from(vec![4, 423, 322, 534]),
                Polynomial::from(vec![635, 330, 637, 321]),
                Polynomial::from(vec![4, 634, 324, 320]),
            ],
            public_witness: Vec::new(),
        };

        assert_eq!(qap, known_good)
    }

    #[test]
    fn r1cs_verification() -> Result<(), Report> {
        let l = vec![
            vec![0, 0, 0],
            vec![0, 0, 0],
            vec![1, 0, 0],
            vec![0, 0, 0],
            vec![0, 1, 0],
            vec![0, 0, 0],
            vec![0, 0, 1],
            vec![0, 0, 0],
        ];

        let r = vec![
            vec![0, 0, 0],
            vec![0, 0, 0],
            vec![0, 0, 0],
            vec![1, 0, 0],
            vec![0, 0, 0],
            vec![0, 1, 0],
            vec![0, 0, 0],
            vec![0, 0, 1],
        ];

        let o = vec![
            vec![0, 0, 0],
            vec![0, 0, 1],
            vec![0, 0, 0],
            vec![0, 0, 0],
            vec![0, 0, 0],
            vec![0, 0, 0],
            vec![1, 0, 0],
            vec![0, 1, 0],
        ];

        debug!("Matrices initialised");
        let r1cs: R1CS<Field> = R1CS::new(l, r, o, Vec::<i32>::new());
        let mut rng = rand::rng();
        let x = Field::from(rng.random_range(0..1000));
        let y = Field::from(rng.random_range(0..1000));
        let z = Field::from(rng.random_range(0..1000));
        let u = Field::from(rng.random_range(0..1000));
        let r = x * y * z * u;
        let v1 = x * y;
        let v2 = z * u;
        let w = vec![Field::from(1), r, x, y, z, u, v1, v2];

        debug!("Generating proof for witness {:?}", w);
        assert_eq!(r, x * y * z * u);
        assert!(r1cs.verify(&w)?);
        Ok(())
    }
}
