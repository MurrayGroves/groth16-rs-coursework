use crate::polynomial::Polynomial;
use ark_ff::{FftField, Field};
use clap::Parser;
use rootcause::hooks::report_creation::AttachmentCollectorHook;

pub struct R1CS<S: FftField> {
    /// Column-wise, i.e. a vec of columns
    L: Vec<Vec<S>>,
    /// Column-wise, i.e. a vec of columns
    R: Vec<Vec<S>>,
    /// Column-wise, i.e. a vec of columns
    O: Vec<Vec<S>>,
    pub public_witness: Vec<S>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct QAP<S>
where
    S: FftField,
{
    /// Output
    pub w: Vec<Polynomial<S>>,
    /// LHS of multiplication
    pub u: Vec<Polynomial<S>>,
    /// RHS of multiplication
    pub v: Vec<Polynomial<S>>,
    pub public_witness: Vec<S>,
}

impl<S: FftField> QAP<S> {
    pub fn degree(&self) -> usize {
        self.w.len()
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
                Polynomial::from_ints(vec![636, 116, 636, 535]),
                Polynomial::from_ints(vec![0, 0, 0, 0]),
                Polynomial::from_ints(vec![8, 416, 5, 213]),
                Polynomial::from_ints(vec![635, 330, 637, 321]),
                Polynomial::from_ints(vec![4, 634, 324, 320]),
                Polynomial::from_ints(vec![640, 536, 640, 107]),
            ],
            v: vec![
                Polynomial::from_ints(vec![3, 529, 323, 427]),
                Polynomial::from_ints(vec![0, 0, 0, 0]),
                Polynomial::from_ints(vec![639, 112, 318, 214]),
                Polynomial::from_ints(vec![0, 0, 0, 0]),
                Polynomial::from_ints(vec![0, 0, 0, 0]),
                Polynomial::from_ints(vec![0, 0, 0, 0]),
            ],
            w: vec![
                Polynomial::from_ints(vec![0, 0, 0, 0]),
                Polynomial::from_ints(vec![640, 536, 640, 107]),
                Polynomial::from_ints(vec![0, 0, 0, 0]),
                Polynomial::from_ints(vec![4, 423, 322, 534]),
                Polynomial::from_ints(vec![635, 330, 637, 321]),
                Polynomial::from_ints(vec![4, 634, 324, 320]),
            ],
            public_witness: Vec::new(),
        };

        assert_eq!(qap, known_good)
    }
}
