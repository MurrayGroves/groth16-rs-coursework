use ark_ec::CurveGroup;
use ark_ff::Field;
use ark_std::iterable::Iterable;
use rootcause::{Report, report};
use std::cmp::Ordering;
use std::fmt::Debug;
use std::iter::{Sum, zip};
use std::ops::{Add, Div, Mul, MulAssign, Sub};

#[derive(Clone, PartialEq, Eq, Debug)]
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

    pub fn evaluate(&self, tau: F) -> F {
        self.coefficients
            .iter()
            .enumerate()
            .map(|(degree, coefficient)| *coefficient * tau.pow([degree as u64]))
            .reduce(std::ops::Add::add)
            .unwrap_or(F::default())
    }
    pub fn degree(&self) -> usize {
        self.coefficients.len()
    }

    /// Find a polynomial by doing Lagrange interpolation over a vector
    pub fn interpolate_from_vector(vec: &Vec<F>) -> Self {
        vec.iter()
            .enumerate()
            .map(|(x, y)| {
                let x = F::from(x as u128);
                &(&(0..vec.len())
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
                    .unwrap_or(Polynomial {
                        coefficients: Vec::new(),
                    })
                    / (0..vec.len())
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

impl<F: Field> Add for Polynomial<F> {
    type Output = Polynomial<F>;

    fn add(self, rhs: Self) -> Self::Output {
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

        a.coefficients = zip(a.coefficients, b.coefficients)
            .map(|(a_i, b_i)| a_i + b_i)
            .collect();
        a
    }
}
impl<F: Field> Sub for Polynomial<F> {
    type Output = Polynomial<F>;

    fn sub(self, rhs: Self) -> Self::Output {
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

        a.coefficients = zip(a.coefficients, b.coefficients)
            .map(|(a_i, b_i)| a_i - b_i)
            .collect();
        a
    }
}

impl<F: Field> Mul for Polynomial<F> {
    type Output = Polynomial<F>;

    fn mul(self, rhs: Polynomial<F>) -> Self::Output {
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
        let mut out = vec![F::default(); a.degree() + b.degree()];
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
    use rand::{Rng, RngCore, random};

    type Field = ark_mnt6_753::Fr;

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
    }

    #[test]
    fn polynomial_interpolation() {
        let mut rng = rand::rng();
        let length: usize = rng.random_range(0..30);
        let vec = (0..length)
            .map(|_| {
                let mut bytes = [0u8; 32];
                rng.fill_bytes(&mut bytes);
                ark_ff::Field::from_random_bytes(&bytes).unwrap()
            })
            .collect();
        let poly = Polynomial::interpolate_from_vector(&vec);
        let out: Vec<Field> = (0..vec.len())
            .map(|x| poly.evaluate(Field::from(x as u128)))
            .collect();

        assert_eq!(vec, out)
    }
}
