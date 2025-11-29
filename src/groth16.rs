use crate::circuits::QAP;
use crate::helpers::rand_scalar;
use ark_ec::mnt6::{MNT6, MNT6Config};
use ark_ec::pairing::Pairing;
use ark_ec::{CurveGroup, PrimeGroup, ScalarMul, mnt6};
use ark_ff::fields::Field;
use ark_ff::{FftField, Fp, MontBackend};
use ark_mnt6_753::{Fq, Fr, FrConfig, G1Projective, G2Projective, MNT6_753};
use rand::{RngCore, SeedableRng};
use rootcause::{Report, report};
use std::borrow::Borrow;
use std::iter::zip;
use std::ops::Mul;

// type G1 = G1Projective;
// type G2 = G2Projective;
//
// pub(crate) type Scalar = Fr;

struct Proof<C: Pairing> {
    a: C::G1,
    b: C::G2,
    c: C::G1,
}

impl<C: Pairing> Proof<C> {
    fn verify(
        &self,
        trusted_setup: TrustedSetupOutput<C>,
        public_witness: Vec<C::ScalarField>,
    ) -> bool {
        let lhs = C::pairing(self.a, self.b);
        let alpha_beta = C::pairing(trusted_setup.alpha, trusted_setup.beta_2);
        let x1 = public_witness
            .iter()
            .enumerate()
            .map(|(i, a_i)| trusted_setup.psi_polynomials[i] * a_i)
            .reduce(|a, b| a + b)
            .unwrap_or(C::G1::generator());
        let x1_gamma = C::pairing(x1, trusted_setup.gamma);
        let c_delta = C::pairing(self.c, trusted_setup.delta_2);

        lhs == alpha_beta + x1_gamma + c_delta
    }
}

struct TrustedSetupOutput<C: Pairing> {
    qap: QAP<C::ScalarField>,
    alpha: C::G1,
    beta_1: C::G1,
    beta_2: C::G2,
    gamma: C::G2,
    delta_1: C::G1,
    delta_2: C::G2,
    group_1_srs: Vec<C::G1>,
    group_2_srs: Vec<C::G2>,
    zero_polynomial_srs: Vec<C::G1>,
    psi_polynomials: Vec<C::G1>,
}

impl<C: Pairing> TrustedSetupOutput<C> {
    fn new(qap: QAP<C::ScalarField>) -> TrustedSetupOutput<C> {
        let mut rng = rand::rngs::StdRng::from_os_rng();

        let alpha: C::ScalarField = rand_scalar(&mut rng);
        let beta: C::ScalarField = rand_scalar(&mut rng);
        let tau: C::ScalarField = rand_scalar(&mut rng);
        let gamma: C::ScalarField = rand_scalar(&mut rng);
        let delta: C::ScalarField = rand_scalar(&mut rng);

        let group_1_srs = (0..qap.degree())
            .rev()
            .map(|i| C::G1::generator() * tau.pow([i as u64]))
            .collect();
        let group_2_srs = (0..qap.degree())
            .rev()
            .map(|i| C::G2::generator() * tau.pow([i as u64]))
            .collect();

        TrustedSetupOutput {
            qap,
            alpha: C::G1::generator() * alpha,
            beta_1: C::G1::generator() * beta,
            beta_2: C::G2::generator() * beta,
            gamma: C::G2::generator() * gamma,
            delta_1: C::G1::generator() * delta,
            delta_2: C::G2::generator() * delta,
            group_1_srs,
            group_2_srs,
            zero_polynomial_srs: Vec::new(),
            psi_polynomials: Vec::new(),
        }
    }

    fn prove(&self, witness: Vec<C::ScalarField>) -> Result<Proof<C>, Report> {
        let mut rng = rand::rngs::StdRng::from_os_rng();

        let r: C::ScalarField = rand_scalar(&mut rng);
        let s: C::ScalarField = rand_scalar(&mut rng);

        let evaluated_u = self
            .qap
            .u
            .iter()
            .map(|x| x.evaluate_over_srs(&self.group_1_srs))
            .collect::<Result<Vec<_>, Report>>()?;

        let a = self.alpha
            + zip(evaluated_u, &witness)
                .map(|(p, a_i)| p * a_i)
                .collect::<Vec<_>>()
                .into_iter()
                .reduce(std::ops::Add::add)
                .ok_or(report!("Empty witness"))?
            + (self.delta_1 * r);

        let evaluated_v_2 = self
            .qap
            .v
            .iter()
            .map(|x| x.evaluate_over_srs(&self.group_2_srs))
            .collect::<Result<Vec<_>, Report>>()?;

        let b_2 = self.beta_2
            + zip(evaluated_v_2, &witness)
                .map(|(p, a_i)| p * a_i)
                .collect::<Vec<_>>()
                .into_iter()
                .reduce(std::ops::Add::add)
                .ok_or(report!("Empty witness"))?
            + (self.delta_2 * s);

        let evaluated_v_1 = self
            .qap
            .v
            .iter()
            .map(|x| x.evaluate_over_srs(&self.group_1_srs))
            .collect::<Result<Vec<_>, Report>>()?;

        let b_1 = self.beta_1
            + zip(evaluated_v_1, &witness)
                .map(|(p, a_i)| p * a_i)
                .collect::<Vec<_>>()
                .into_iter()
                .reduce(std::ops::Add::add)
                .ok_or(report!("Empty witness"))?
            + (self.delta_1 * s);

        let c = zip(&self.psi_polynomials, &witness)
            .skip(self.qap.public_witness.len())
            .map(|(psi, a_i)| *psi * a_i)
            .reduce(std::ops::Add::add)
            .ok_or(report!("Empty witness"))?
            + (a * s)
            + (b_1 * r)
            - (self.delta_1 * (r * s));

        Ok(Proof { a, b: b_2, c })
    }
}
