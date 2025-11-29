use ark_ec::CurveGroup;
use ark_ff::Field;
use rootcause::{Report, report};
use std::fmt::Debug;
use std::ops::{Mul, MulAssign};

pub struct Polynomial<F>
where
    F: Field,
{
    coefficients: Vec<F>,
}

impl<F: Field> Polynomial<F> {
    pub fn evaluate_over_srs<T>(&self, srs: &Vec<T>) -> Result<T, Report>
    where
        T: MulAssign<F> + CurveGroup + Debug,
    {
        if srs.len() < self.coefficients.len() {
            return Err(report!("SRS too small for polynomial")
                .attach(format!("SRS: {:?}", srs))
                .attach(format!("Polynomial: {:?}", self.coefficients)));
        }

        self.coefficients
            .iter()
            .enumerate()
            .map(|(degree, coefficient)| {
                let mut result = srs[degree];
                result *= *coefficient;
                result
            })
            .reduce(std::ops::Add::add)
            .ok_or(report!("Polynomial has no coefficients"))
    }

    pub fn degree(&self) -> usize {
        self.coefficients.len()
    }
}
