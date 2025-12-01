use crate::circuits::QAP;
use crate::helpers::rand_scalar;
use crate::polynomial::Polynomial;
use ark_ec::PrimeGroup;
use ark_ec::pairing::Pairing;
use ark_ff::fields::Field;
use itertools::izip;
use log::debug;
use rand::SeedableRng;
use rootcause::prelude::ResultExt;
use rootcause::{Report, bail, report};
use std::iter::zip;

struct Proof<C: Pairing> {
    a: C::G1,
    b: C::G2,
    c: C::G1,
}

impl<C: Pairing> Proof<C> {
    pub fn verify(
        &self,
        trusted_setup: TrustedSetupOutput<C>,
        public_witness: &Vec<C::ScalarField>,
    ) -> bool {
        debug!("Verifying with public witness: {:?}", public_witness);
        let lhs = C::pairing(self.a, self.b);
        let alpha_beta = C::pairing(trusted_setup.alpha, trusted_setup.beta_2);
        let x1 = public_witness
            .iter()
            .enumerate()
            .map(|(i, a_i)| trusted_setup.psi_polynomials[i] * a_i)
            .reduce(|a, b| a + b);
        let x1_gamma = if let Some(x1) = x1 {
            Some(C::pairing(x1, trusted_setup.gamma))
        } else {
            None
        };
        let c_delta = C::pairing(self.c, trusted_setup.delta_2);

        debug!("{} == {} + {:?} + {}", lhs, alpha_beta, x1_gamma, c_delta);
        let rhs = if let Some(x1_gamma) = x1_gamma {
            alpha_beta + x1_gamma + c_delta
        } else {
            alpha_beta + c_delta
        };
        lhs.0 == rhs.0
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
    fn group_1_srs(length: usize, tau: C::ScalarField) -> Vec<C::G1> {
        (0..length)
            .map(|i| C::G1::generator() * tau.pow([i as u64]))
            .collect()
    }

    /// Get zero polynomial (x - 1)(x -2)(...)(x - n)
    fn t(num_roots: usize) -> Result<Polynomial<C::ScalarField>, Report> {
        Ok((1..num_roots + 1)
            .map(|x| {
                Polynomial::new(vec![
                    -C::ScalarField::from(x as u128),
                    C::ScalarField::from(1),
                ])
            })
            .reduce(std::ops::Mul::mul)
            .ok_or(report!("QAP has degree zero"))?)
    }

    /// Generate SRS for the zero polynomial of form [t(tau)/delta, tau * t(tau)/delta, tau^2 * t(tau)/delta, ...]
    ///
    /// # Arguments
    ///
    /// * `num_evaluation_points`: Roots of zero polynomial, i.e. `[1,2,...].len()`, i.e. the degree of the QAP
    /// * `srs_length`: The length of the SRS is 1 more than the degree of polynomial it needs to support
    /// * `delta`: Secret scalar used to ensure separation of public/private witness
    /// * `group_1_srs`: SRS for G1
    ///
    /// returns: Result<Vec<<C as Pairing>::G1, Global>, Report<dyn Any, Mutable, SendSync>>
    fn zero_polynomial_srs(
        num_evaluation_points: usize,
        srs_length: usize,
        delta: C::ScalarField,
        group_1_srs: &Vec<C::G1>,
    ) -> Result<Vec<C::G1>, Report> {
        let t_tau = Self::t(num_evaluation_points)?;

        debug!("Generated t(tau)");

        Ok((0..srs_length)
            .map(|i| {
                let mut Zx = vec![0; i + 1];
                Zx[i] = 1;
                let poly = Polynomial::from(Zx);
                (&(&t_tau * &poly) / delta)
                    .evaluate_over_srs(&group_1_srs)
                    .context("Evaluating over group_1_srs")
                    .attach(format!("t_tau: {:?}", t_tau))
                    .attach(format!("poly: {:?}", poly))
            })
            .collect::<Result<_, _>>()?)
    }

    fn psi_polynomials(
        qap: &QAP<C::ScalarField>,
        group_1_srs: &Vec<C::G1>,
        alpha: C::ScalarField,
        beta: C::ScalarField,
        gamma: C::ScalarField,
        delta: C::ScalarField,
    ) -> Result<Vec<C::G1>, Report> {
        izip!(&qap.u, &qap.v, &qap.w)
            .enumerate()
            .map(|(i, (u_i, v_i, w_i))| {
                let divisor = if i < qap.public_witness.len() {
                    gamma
                } else {
                    delta
                };

                Ok(
                    (&((v_i * alpha) + (u_i * beta)) + w_i).evaluate_over_srs(&group_1_srs)?
                        * (C::ScalarField::from(1) / divisor),
                )
            })
            .collect::<Result<Vec<_>, Report>>()
    }

    pub fn new(qap: QAP<C::ScalarField>) -> Result<TrustedSetupOutput<C>, Report> {
        debug!("Starting trusted setup");
        let mut rng = rand::rngs::StdRng::from_os_rng();
        debug!("Got RNG");

        let alpha: C::ScalarField = rand_scalar(&mut rng);
        let beta: C::ScalarField = rand_scalar(&mut rng);
        let tau: C::ScalarField = rand_scalar(&mut rng);
        let gamma: C::ScalarField = rand_scalar(&mut rng);
        let delta: C::ScalarField = rand_scalar(&mut rng);

        debug!("Generated random scalars");

        let group_1_srs = Self::group_1_srs((2 * qap.degree()) - 1, tau);

        debug!("Generated Group 1 SRS");

        let group_2_srs = (0..qap.degree())
            .map(|i| C::G2::generator() * tau.pow([i as u64]))
            .collect();

        debug!("Generated Group 2 SRS");

        let zero_polynomial_srs =
            Self::zero_polynomial_srs(qap.degree(), qap.degree() - 1, delta, &group_1_srs)
                .context("Calculating zero polynomial SRS")?;

        debug!("Generated zero polynomial srs");

        let psi_polynomials = Self::psi_polynomials(&qap, &group_1_srs, alpha, beta, gamma, delta)?;

        debug!("Generated psi polynomials");

        Ok(TrustedSetupOutput {
            qap,
            alpha: C::G1::generator() * alpha,
            beta_1: C::G1::generator() * beta,
            beta_2: C::G2::generator() * beta,
            gamma: C::G2::generator() * gamma,
            delta_1: C::G1::generator() * delta,
            delta_2: C::G2::generator() * delta,
            group_1_srs,
            group_2_srs,
            zero_polynomial_srs,
            psi_polynomials,
        })
    }

    fn calculate_zero_polynomial(
        &self,
        witness: &Vec<C::ScalarField>,
    ) -> Result<Polynomial<C::ScalarField>, Report> {
        let au_sum: Polynomial<C::ScalarField> =
            zip(&self.qap.u, witness).map(|(u_i, a_i)| u_i * *a_i).sum();
        let av_sum: Polynomial<C::ScalarField> =
            zip(&self.qap.v, witness).map(|(v_i, a_i)| v_i * *a_i).sum();
        let aw_sum: Polynomial<C::ScalarField> =
            zip(&self.qap.w, witness).map(|(w_i, a_i)| w_i * *a_i).sum();
        Ok(((&(au_sum * av_sum) - &aw_sum) / Self::t(self.qap.degree())?)?)
    }

    fn evaluate_u(&self, witness: &Vec<C::ScalarField>) -> Result<C::G1, Report> {
        if witness.len() != self.qap.u.len() {
            bail!("Witness wrong size for QAP")
        }
        let evaluated_u = self
            .qap
            .u
            .iter()
            .map(|x| x.evaluate_over_srs(&self.group_1_srs))
            .collect::<Result<Vec<_>, Report>>()?;

        Ok(zip(evaluated_u, witness)
            .map(|(p, a_i)| p * a_i)
            .collect::<Vec<_>>()
            .into_iter()
            .reduce(std::ops::Add::add)
            .ok_or(report!("Empty witness"))?)
    }

    fn evaluate_v(&self, witness: &Vec<C::ScalarField>) -> Result<C::G2, Report> {
        if witness.len() != self.qap.v.len() {
            bail!("Witness wrong size for QAP")
        }
        let evaluated_v = self
            .qap
            .v
            .iter()
            .map(|x| x.evaluate_over_srs(&self.group_2_srs))
            .collect::<Result<Vec<_>, Report>>()?;

        Ok(zip(evaluated_v, witness)
            .map(|(p, a_i)| p * a_i)
            .collect::<Vec<_>>()
            .into_iter()
            .reduce(std::ops::Add::add)
            .ok_or(report!("Empty witness"))?)
    }

    fn evaluate_v_1(&self, witness: &Vec<C::ScalarField>) -> Result<C::G1, Report> {
        if witness.len() != self.qap.v.len() {
            bail!("Witness wrong size for QAP")
        }
        let evaluated_v = self
            .qap
            .v
            .iter()
            .map(|x| x.evaluate_over_srs(&self.group_1_srs))
            .collect::<Result<Vec<_>, Report>>()?;

        Ok(zip(evaluated_v, witness)
            .map(|(p, a_i)| p * a_i)
            .collect::<Vec<_>>()
            .into_iter()
            .reduce(std::ops::Add::add)
            .ok_or(report!("Empty witness"))?)
    }
    fn evaluate_w(&self, witness: &Vec<C::ScalarField>) -> Result<C::G1, Report> {
        if witness.len() != self.qap.w.len() {
            bail!("Witness wrong size for QAP")
        }
        let evaluated_w = self
            .qap
            .w
            .iter()
            .map(|x| x.evaluate_over_srs(&self.group_1_srs))
            .collect::<Result<Vec<_>, Report>>()?;

        Ok(zip(evaluated_w, witness)
            .map(|(p, a_i)| p * a_i)
            .collect::<Vec<_>>()
            .into_iter()
            .reduce(std::ops::Add::add)
            .ok_or(report!("Empty witness"))?)
    }

    fn prove(&self, witness: &Vec<C::ScalarField>) -> Result<Proof<C>, Report> {
        let mut rng = rand::rngs::StdRng::from_os_rng();

        let r: C::ScalarField = rand_scalar(&mut rng);
        let s: C::ScalarField = rand_scalar(&mut rng);

        let a = self.alpha + self.evaluate_u(witness)? + (self.delta_1 * r);

        let b_2 = self.beta_2 + self.evaluate_v(witness)? + (self.delta_2 * s);

        let b_1 = self.beta_1 + self.evaluate_v_1(witness)? + (self.delta_1 * s);

        let ht = self.calculate_zero_polynomial(witness)?;
        let ht_tau = ht.evaluate_over_srs(&self.zero_polynomial_srs)?;

        let c = zip(&self.psi_polynomials, witness)
            .skip(self.qap.public_witness.len())
            .map(|(psi, a_i)| *psi * a_i)
            .reduce(std::ops::Add::add)
            .ok_or(report!("Empty witness"))?
            + ht_tau
            + (a * s)
            + (b_1 * r)
            - (self.delta_1 * (r * s));
        Ok(Proof { a, b: b_2, c })
    }
}

#[cfg(test)]
mod tests {
    use crate::circuits::{QAP, R1CS};
    use crate::groth16::TrustedSetupOutput;
    use crate::helpers::rand_scalar;
    use crate::polynomial::Polynomial;
    use ark_ec::PrimeGroup;
    use ark_ec::pairing::{MillerLoopOutput, Pairing, PairingOutput};
    use ark_ff::{MontConfig, PrimeField};
    use ark_mnt6_753::MNT6_753;
    use log::debug;
    use rand::Rng;
    use rootcause::prelude::ResultExt;
    use rootcause::{Report, report};
    use std::iter::zip;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    type Field = ark_mnt6_753::Fr;
    #[test]
    fn groth16() -> Result<(), Report> {
        init();

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
        let r1cs: R1CS<Field> = R1CS::new(l, r, o, Vec::new());

        debug!("R1CS initialised: {:?}", r1cs);
        let qap = QAP::from(r1cs.clone());

        debug!("QAP derived");
        let trusted_setup: TrustedSetupOutput<ark_mnt6_753::MNT6_753> =
            TrustedSetupOutput::new(qap.clone())?;

        debug!("Trusted Setup complete");
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
        assert!(qap.verify(&w)?);
        let proof = trusted_setup.prove(&w)?;

        debug!("Proof generated");
        assert!(proof.verify(trusted_setup, &qap.public_witness));
        Ok(())
    }

    #[test]
    fn group_1_srs() -> Result<(), Report> {
        let mut rng = rand::rng();
        let tau: ark_mnt6_753::Fr = rand_scalar(&mut rng);

        let srs = TrustedSetupOutput::<MNT6_753>::group_1_srs(16, tau);

        let polynomial: Polynomial<ark_mnt6_753::Fr> = Polynomial::from(vec![3, 5, 10, 20]);

        assert_eq!(
            ark_mnt6_753::G1Projective::generator() * polynomial.evaluate(&tau),
            polynomial.evaluate_over_srs(&srs)?
        );
        Ok(())
    }

    #[test]
    fn zero_polynomial_is_zero() -> Result<(), Report> {
        init();
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
        let r1cs: R1CS<Field> = R1CS::new(l, r, o, Vec::new());

        debug!("R1CS initialised: {:?}", r1cs);
        let qap = QAP::from(r1cs.clone());

        debug!("QAP derived");
        let trusted_setup: TrustedSetupOutput<ark_mnt6_753::MNT6_753> =
            TrustedSetupOutput::new(qap.clone())?;

        debug!("Trusted Setup complete");

        let mut rng = rand::rng();
        let x = Field::from(rng.random_range(0..1000));
        let y = Field::from(rng.random_range(0..1000));
        let z = Field::from(rng.random_range(0..1000));
        let u = Field::from(rng.random_range(0..1000));
        let r = x * y * z * u;
        let v1 = x * y;
        let v2 = z * u;
        let w = vec![Field::from(1), r, x, y, z, u, v1, v2];

        let zero_polynomial = TrustedSetupOutput::<MNT6_753>::t(qap.degree())?;

        debug!("QAP has degree {}", qap.max_polynomial_degree());
        for i in 1..qap.max_polynomial_degree() + 1 {
            debug!("Running on X={}", i);
            assert_eq!(
                zero_polynomial.evaluate(&Field::from(i as u128)),
                Field::default()
            )
        }

        Ok(())
    }

    #[test]
    fn naive() -> Result<(), Report> {
        init();
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
        let r1cs: R1CS<Field> = R1CS::new(l, r, o, Vec::new());

        debug!("R1CS initialised: {:?}", r1cs);
        let qap = QAP::from(r1cs.clone());
        assert_eq!(qap.degree(), 3);

        debug!("QAP derived");
        let trusted_setup: TrustedSetupOutput<ark_mnt6_753::MNT6_753> =
            TrustedSetupOutput::new(qap.clone())?;

        debug!("Trusted Setup complete");

        let mut rng = rand::rng();
        let x = Field::from(rng.random_range(0..1000));
        let y = Field::from(rng.random_range(0..1000));
        let z = Field::from(rng.random_range(0..1000));
        let u = Field::from(rng.random_range(0..1000));
        let r = x * y * z * u;
        let v1 = x * y;
        let v2 = z * u;
        let w = vec![Field::from(1), r, x, y, z, u, v1, v2];

        debug!(
            "QAP Max Polynomial Degree: {:?}",
            qap.max_polynomial_degree()
        );

        let zero_polynomial_srs: Vec<<MNT6_753 as Pairing>::G1> =
            TrustedSetupOutput::<MNT6_753>::zero_polynomial_srs(
                qap.degree(),
                (qap.degree()) - 1,
                <MNT6_753 as Pairing>::ScalarField::from(1),
                &trusted_setup.group_1_srs,
            )?;

        assert_eq!(zero_polynomial_srs.len(), qap.degree() - 1);

        let ht = trusted_setup.calculate_zero_polynomial(&w)?;
        debug!("zero polynomial: {:?}", ht);
        let ht_tau = ht
            .evaluate_over_srs(&zero_polynomial_srs)
            .context("Evaluating ht_tau")?;

        let a = trusted_setup.evaluate_u(&w).context("Evaluating u")?;
        let b = trusted_setup.evaluate_v(&w).context("Evaluating u")?;
        let c = trusted_setup.evaluate_w(&w).context("Evaluating w")? + ht_tau;

        let lhs = MNT6_753::pairing(a, b).0;
        let rhs = MNT6_753::pairing(c, <MNT6_753 as Pairing>::G2::generator()).0;

        debug!("{:?} == {:?}", lhs, rhs);

        assert_eq!(lhs, rhs);

        Ok(())
    }
    #[test]
    fn naive_plus_alpha_beta() -> Result<(), Report> {
        init();
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
        let r1cs: R1CS<Field> = R1CS::new(l, r, o, Vec::new());

        debug!("R1CS initialised: {:?}", r1cs);
        let qap = QAP::from(r1cs.clone());
        assert_eq!(qap.degree(), 3);

        debug!("QAP derived");
        let trusted_setup: TrustedSetupOutput<ark_mnt6_753::MNT6_753> =
            TrustedSetupOutput::new(qap.clone())?;

        debug!("Trusted Setup complete");

        let mut rng = rand::rng();
        let x = Field::from(rng.random_range(0..1000));
        let y = Field::from(rng.random_range(0..1000));
        let z = Field::from(rng.random_range(0..1000));
        let u = Field::from(rng.random_range(0..1000));
        let r = x * y * z * u;
        let v1 = x * y;
        let v2 = z * u;
        let w = vec![Field::from(1), r, x, y, z, u, v1, v2];

        debug!(
            "QAP Max Polynomial Degree: {:?}",
            qap.max_polynomial_degree()
        );

        let zero_polynomial_srs: Vec<<MNT6_753 as Pairing>::G1> =
            TrustedSetupOutput::<MNT6_753>::zero_polynomial_srs(
                qap.degree(),
                (qap.degree()) - 1,
                <MNT6_753 as Pairing>::ScalarField::from(1),
                &trusted_setup.group_1_srs,
            )?;

        assert_eq!(zero_polynomial_srs.len(), qap.degree() - 1);

        let ht = trusted_setup.calculate_zero_polynomial(&w)?;
        debug!("zero polynomial: {:?}", ht);
        let ht_tau = ht
            .evaluate_over_srs(&zero_polynomial_srs)
            .context("Evaluating ht_tau")?;

        let alpha = rand_scalar(&mut rng);
        let beta = rand_scalar(&mut rng);
        let alpha_1 = <MNT6_753 as Pairing>::G1::generator() * alpha;
        let beta_2 = <MNT6_753 as Pairing>::G2::generator() * beta;
        let psi_polynomials = TrustedSetupOutput::<MNT6_753>::psi_polynomials(
            &qap,
            &trusted_setup.group_1_srs,
            alpha,
            beta,
            <MNT6_753 as Pairing>::ScalarField::from(1),
            <MNT6_753 as Pairing>::ScalarField::from(1),
        )?;

        assert_eq!(psi_polynomials.len(), w.len());

        let a = alpha_1 + trusted_setup.evaluate_u(&w).context("Evaluating u")?;
        let b = beta_2 + trusted_setup.evaluate_v(&w).context("Evaluating v")?;
        let c = zip(&psi_polynomials, &w)
            .skip(qap.public_witness.len())
            .map(|(psi, a_i)| *psi * a_i)
            .reduce(std::ops::Add::add)
            .ok_or(report!("Empty witness"))?
            + ht_tau;

        let lhs = MNT6_753::pairing(a, b).0;
        let rhs = (MNT6_753::pairing(alpha_1, beta_2)
            + MNT6_753::pairing(c, <MNT6_753 as Pairing>::G2::generator()))
        .0;

        assert_eq!(lhs, rhs);

        Ok(())
    }
    #[test]
    fn naive_plus_gamma_delta() -> Result<(), Report> {
        init();
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
        let r1cs: R1CS<Field> = R1CS::new(l, r, o, Vec::new());

        debug!("R1CS initialised: {:?}", r1cs);
        let qap = QAP::from(r1cs.clone());
        assert_eq!(qap.degree(), 3);

        debug!("QAP derived");
        let trusted_setup: TrustedSetupOutput<ark_mnt6_753::MNT6_753> =
            TrustedSetupOutput::new(qap.clone())?;

        debug!("Trusted Setup complete");

        let mut rng = rand::rng();
        let x = Field::from(rng.random_range(0..1000));
        let y = Field::from(rng.random_range(0..1000));
        let z = Field::from(rng.random_range(0..1000));
        let u = Field::from(rng.random_range(0..1000));
        let r = x * y * z * u;
        let v1 = x * y;
        let v2 = z * u;
        let w = vec![Field::from(1), r, x, y, z, u, v1, v2];

        debug!(
            "QAP Max Polynomial Degree: {:?}",
            qap.max_polynomial_degree()
        );

        let delta = rand_scalar(&mut rng);
        let zero_polynomial_srs: Vec<<MNT6_753 as Pairing>::G1> =
            TrustedSetupOutput::<MNT6_753>::zero_polynomial_srs(
                qap.degree(),
                (qap.degree()) - 1,
                delta,
                &trusted_setup.group_1_srs,
            )?;

        assert_eq!(zero_polynomial_srs.len(), qap.degree() - 1);

        let ht = trusted_setup.calculate_zero_polynomial(&w)?;
        debug!("zero polynomial: {:?}", ht);
        let ht_tau = ht
            .evaluate_over_srs(&zero_polynomial_srs)
            .context("Evaluating ht_tau")?;

        let alpha = rand_scalar(&mut rng);
        let beta = rand_scalar(&mut rng);
        let alpha_1 = <MNT6_753 as Pairing>::G1::generator() * alpha;
        let beta_2 = <MNT6_753 as Pairing>::G2::generator() * beta;
        let psi_polynomials = TrustedSetupOutput::<MNT6_753>::psi_polynomials(
            &qap,
            &trusted_setup.group_1_srs,
            alpha,
            beta,
            <MNT6_753 as Pairing>::ScalarField::from(1),
            delta,
        )?;

        assert_eq!(psi_polynomials.len(), w.len());

        let a = alpha_1 + trusted_setup.evaluate_u(&w).context("Evaluating u")?;
        let b = beta_2 + trusted_setup.evaluate_v(&w).context("Evaluating v")?;
        let c = zip(&psi_polynomials, &w)
            .skip(qap.public_witness.len())
            .map(|(psi, a_i)| *psi * a_i)
            .reduce(std::ops::Add::add)
            .ok_or(report!("Empty witness"))?
            + ht_tau;

        let delta_2 = <MNT6_753 as Pairing>::G2::generator() * delta;
        let lhs = MNT6_753::pairing(a, b).0;
        let rhs = (MNT6_753::pairing(alpha_1, beta_2)
            // + MNT6_753::pairing(
            //     <MNT6_753 as Pairing>::G1::generator(),
            //     <MNT6_753 as Pairing>::G2::generator()
            //         * <MNT6_753 as Pairing>::ScalarField::from(1),
            // )
            + MNT6_753::pairing(c, delta_2))
        .0;

        assert_eq!(lhs, rhs);

        Ok(())
    }
}
