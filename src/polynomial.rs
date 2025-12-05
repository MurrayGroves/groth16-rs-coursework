use ark_ec::CurveGroup;
use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::iterable::Iterable;
use log::{debug, trace};
use rootcause::prelude::ResultExt;
use rootcause::{Report, bail, report};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::fmt::Debug;
use std::iter::{Sum, zip};
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Sub, SubAssign};

#[derive(
    Clone, PartialEq, Eq, Debug, Serialize, Deserialize, CanonicalSerialize, CanonicalDeserialize,
)]
pub struct Polynomial<F>
where
    F: Field,
{
    /// x^0, x^1, x^2, ...
    coefficients: Vec<F>,
}

impl<F: Field> Polynomial<F> {
    pub fn evaluate_over_srs<T>(&self, srs: &Vec<T>) -> Result<T, Report>
    where
        T: MulAssign<F> + CurveGroup + Debug,
    {
        if srs.len() < self.coefficients.len() {
            return Err(report!("SRS too small for polynomial")
                .attach(format!(
                    "SRS degree: {:?} (supports polynomial of degree {:?})",
                    srs.len(),
                    srs.len() - 1
                ))
                .attach(format!("Polynomial degree: {:?}", self.degree())));
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

    pub fn evaluate(&self, tau: &F) -> F {
        self.coefficients
            .iter()
            .enumerate()
            .map(|(degree, coefficient)| *coefficient * tau.pow([degree as u64]))
            .reduce(std::ops::Add::add)
            .unwrap_or(F::default())
    }
    pub fn degree(&self) -> usize {
        self.coefficients.len().checked_sub(1).unwrap_or(0)
    }

    pub fn is_zero(&self) -> bool {
        self.coefficients.len() == 0 || self.coefficients.iter().all(|x| *x == F::default())
    }

    pub fn lead(&self) -> Polynomial<F> {
        Polynomial {
            coefficients: self
                .coefficients
                .iter()
                .enumerate()
                .map(|(i, coefficient)| {
                    if i == self.coefficients.len() - 1 {
                        coefficient.clone()
                    } else {
                        F::default()
                    }
                })
                .collect(),
        }
    }

    pub fn is_lead(&self) -> bool {
        self.coefficients
            .iter()
            .rev()
            .skip(1)
            .all(|x| *x == F::default())
    }

    pub fn new(vec: Vec<F>) -> Self {
        Polynomial { coefficients: vec }
    }

    /// Find a polynomial by doing Lagrange interpolation over a vector
    pub fn interpolate_from_vector(vec: &Vec<F>) -> Self {
        vec.iter()
            .enumerate()
            .map(|(x, y)| {
                let x = F::from((x + 1) as u128);
                &(&(1..vec.len() + 1)
                    .filter_map(|x_i| {
                        if F::from(x_i as u128) == x {
                            None
                        } else {
                            Some(Polynomial {
                                coefficients: vec![-F::from(x_i as u128), F::from(1)],
                            })
                        }
                    })
                    .reduce(std::ops::Mul::mul)
                    .unwrap_or(Polynomial::from(vec![1]))
                    / (1..vec.len() + 1)
                        .filter_map(|x_i| {
                            if F::from(x_i as u128) == x {
                                None
                            } else {
                                Some(x - F::from(x_i as u128))
                            }
                        })
                        .reduce(std::ops::Mul::mul)
                        .unwrap_or(F::from(1)))
                    * *y
            })
            .sum()
    }
}

impl<T, F: Field> From<Vec<T>> for Polynomial<F>
where
    T: Copy,
    F: From<T>,
{
    fn from(vec: Vec<T>) -> Self {
        Polynomial {
            coefficients: vec.iter().map(|x| F::from(*x)).collect(),
        }
    }
}

impl<F: Field> Add for Polynomial<F> {
    type Output = Polynomial<F>;

    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl<F: Field> Add for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn add(self, rhs: Self) -> Self::Output {
        let mut a = self.clone();
        let mut b = rhs.clone();

        match a.coefficients.len().cmp(&b.coefficients.len()) {
            Ordering::Less => {
                // Pad a
                let new_elems = b.coefficients.len() - a.coefficients.len();
                let padding = vec![F::default(); new_elems];
                a.coefficients = [a.coefficients, padding].concat()
            }
            Ordering::Greater => {
                // Pad b
                let new_elems = a.coefficients.len() - b.coefficients.len();
                let padding = vec![F::default(); new_elems];
                b.coefficients = [b.coefficients, padding].concat();
            }
            _ => {}
        }

        a.coefficients = zip(a.coefficients, b.coefficients)
            .map(|(a_i, b_i)| a_i + b_i)
            .collect();
        a
    }
}

impl<F: Field> AddAssign for Polynomial<F> {
    fn add_assign(&mut self, rhs: Self) {
        *self = &*self + &rhs;
    }
}

impl<F: Field> Sub for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut a = self.clone();
        let mut b = (*rhs).clone();

        match a.degree().cmp(&b.degree()) {
            Ordering::Less => {
                // Pad a
                let new_elems = b.degree() - a.degree();
                let padding = vec![F::default(); new_elems];
                a.coefficients = [a.coefficients, padding].concat()
            }
            Ordering::Greater => {
                // Pad b
                let new_elems = a.degree() - b.degree();
                let padding = vec![F::default(); new_elems];
                b.coefficients = [b.coefficients, padding].concat();
            }
            _ => {}
        }

        a.coefficients = zip(a.coefficients, b.coefficients)
            .map(|(a_i, b_i)| a_i - b_i)
            .collect();
        a
    }
}

impl<F: Field> SubAssign for Polynomial<F> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = &*self - &rhs;
    }
}

impl<F: Field> Mul for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn mul(self, rhs: &Polynomial<F>) -> Self::Output {
        let mut a = self.clone();
        let mut b = rhs.clone();

        match a.degree().cmp(&b.degree()) {
            Ordering::Less => {
                // Pad a
                let new_elems = b.degree() - a.degree();
                let padding = vec![F::default(); new_elems];
                a.coefficients = [a.coefficients, padding].concat()
            }
            Ordering::Greater => {
                // Pad b
                let new_elems = a.degree() - b.degree();
                let padding = vec![F::default(); new_elems];
                b.coefficients = [b.coefficients, padding].concat();
            }
            _ => {}
        }

        // n^2 approach from https://home.cse.ust.hk/~dekai/271/notes/L03/L03.pdf page 4.
        // TODO - Use n log n algorithm with FFT
        let mut out = vec![F::default(); 1 + a.degree() + b.degree()]; // +1 for constant
        for k in 0..out.len() {
            let mut coefficient = F::default();
            for i in 0..k + 1 {
                coefficient += *a.coefficients.get(i).unwrap_or(&F::default())
                    * *b.coefficients.get(k - i).unwrap_or(&F::default())
            }
            out[k] = coefficient
        }

        // Truncate trailing zeroes
        if let Some(pos) = out.iter().rposition(|x| *x != F::default()) {
            out.truncate(pos + 1)
        } else {
            out = Vec::new()
        }
        Polynomial { coefficients: out }
    }
}

impl<F: Field> Mul for Polynomial<F> {
    type Output = Polynomial<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
    }
}

impl<F: Field> Div for Polynomial<F> {
    type Output = Result<Polynomial<F>, Report>;

    fn div(self, rhs: Self) -> Self::Output {
        if rhs.is_zero() {
            bail!("Divisor is zero")
        } else if self.is_zero() {
            return Ok(Polynomial::new(vec![]));
        }

        if self.is_lead() && rhs.is_lead() {
            let degree = self.degree() - rhs.degree();

            let mut coefficients = vec![F::default(); degree + 1];
            coefficients[degree] =
                *self.coefficients.last().unwrap() / rhs.coefficients.last().unwrap();

            if coefficients.iter().all(|x| *x == F::default()) {
                return Ok(Polynomial::new(vec![]));
            }

            let out = Polynomial::from(coefficients);
            trace!("{:?}/{:?} == {:?}", self, rhs, out);
            return Ok(out);
        }

        trace!("Dividing {:?}/{:?}", self, rhs);
        let mut quotient = Polynomial::new(vec![]);
        let mut remainder = self.clone();

        while !remainder.is_zero() && remainder.degree() >= rhs.degree() {
            let tmp = (remainder.lead().clone() / rhs.lead())
                .context("Dividing lead")
                .attach(format!("LHS: {:?}", self.lead()))
                .attach(format!("RHS: {:?}", rhs.lead()))?;
            quotient += tmp.clone();
            remainder -= &tmp * &rhs;
            trace!(
                "Quotient: {:?}\nRemainder: {:?}\nTmp: {:?}",
                quotient, remainder, tmp
            );
        }

        if !remainder.is_zero() {
            Err(report!("Non zero remainder").attach(format!("Remainder: {:?}", remainder)))
        } else {
            Ok(quotient)
        }
    }
}

impl<F: Field> Mul<F> for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn mul(self, rhs: F) -> Self::Output {
        Polynomial {
            coefficients: self.coefficients.iter().map(|x| *x * rhs).collect(),
        }
    }
}
impl<F: Field> Div<F> for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn div(self, rhs: F) -> Self::Output {
        Polynomial {
            coefficients: self.coefficients.iter().map(|x| *x / rhs).collect(),
        }
    }
}

impl<F: Field> Sum for Polynomial<F> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(std::ops::Add::add).unwrap_or(Polynomial {
            coefficients: Vec::new(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use log::debug;
    use rand::{Rng, RngCore};

    type Field = ark_mnt6_753::Fr;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn polynomial_add() {
        let a = Polynomial {
            coefficients: vec![Field::from(1), Field::from(3), Field::from(5)],
        };
        let b = Polynomial {
            coefficients: vec![Field::from(0), Field::from(10), Field::from(2)],
        };
        let c = Polynomial {
            coefficients: vec![Field::from(1), Field::from(13), Field::from(7)],
        };
        assert_eq!(a + b, c);
    }

    #[test]
    fn polynomial_add_different_order() {
        init();
        let a = Polynomial {
            coefficients: vec![Field::from(1), Field::from(3), Field::from(5)],
        };
        let b = Polynomial {
            coefficients: vec![Field::from(10), Field::from(2)],
        };
        let c = Polynomial {
            coefficients: vec![Field::from(11), Field::from(5), Field::from(5)],
        };
        assert_eq!(a + b, c);

        let a = Polynomial::<Field>::new(vec![]);
        let b = Polynomial::from(vec![0, 0, 10]);
        assert_eq!(&a + &b, b);
    }

    #[test]
    fn polynomial_mult() {
        let a = Polynomial {
            coefficients: vec![Field::from(1), Field::from(3), Field::from(5)],
        };
        let b = Polynomial {
            coefficients: vec![Field::from(10), Field::from(2)],
        };
        let c = Polynomial {
            coefficients: vec![
                Field::from(10),
                Field::from(32),
                Field::from(56),
                Field::from(10),
            ],
        };
        assert_eq!(a * b, c);

        let a: Polynomial<Field> = Polynomial::from(vec![-3, 1]);
        let b: Polynomial<Field> = Polynomial::from(vec![-4, 1]);
        let c: Polynomial<Field> = Polynomial::from(vec![12, -7, 1]);
        assert_eq!(a * b, c);
    }

    #[test]
    fn polynomial_interpolation() {
        let mut rng = rand::rng();
        let length: usize = rng.random_range(1..30);
        let vec: Vec<Field> = (1..length + 1)
            .map(|_| {
                let mut bytes = [0u8; 32];
                rng.fill_bytes(&mut bytes);
                ark_ff::Field::from_random_bytes(&bytes).unwrap()
            })
            .collect();
        assert_eq!(vec.len(), length);
        let poly = Polynomial::interpolate_from_vector(&vec);
        assert_eq!(poly.degree(), length - 1);
        let out: Vec<Field> = (1..length + 1)
            .map(|x| poly.evaluate(&Field::from(x as u128)))
            .collect();

        debug!("Poly: {:?}", poly);
        assert_eq!(vec, out);
    }

    #[test]
    fn polynomial_evaluation() {
        let poly = Polynomial::<Field>::from(vec![3, 2, 4]);
        let out = poly.evaluate(&Field::from(3));
        assert_eq!(out, Field::from(45))
    }

    #[test]
    fn polynomial_division() -> Result<(), Report> {
        init();
        let a: Polynomial<Field> = Polynomial::from(vec![3, 5, 10]);
        let b = Polynomial::from(vec![3, 5, 10]);
        let c = &a * &b;
        assert_eq!((c / b)?, a);
        Ok(())
    }
}
