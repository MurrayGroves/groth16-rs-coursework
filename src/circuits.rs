use crate::polynomial::Polynomial;
use ark_ff::{FftField, Field};
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
