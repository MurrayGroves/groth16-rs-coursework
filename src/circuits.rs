use crate::polynomial::Polynomial;
use ark_ff::{FftField, Field};

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
