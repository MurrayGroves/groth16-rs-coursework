use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
use rand::Rng;

pub(crate) fn rand_scalar<T, S>(rng: &mut T) -> S
where
    T: Rng,
    S: Field,
{
    let mut bytes = [0; 256];
    let mut out = None;
    while out.is_none() {
        rng.fill_bytes(&mut bytes);
        out = S::from_random_bytes(&bytes);
    }

    out.unwrap()
}

pub(crate) fn ark_se<S, A: CanonicalSerialize>(a: &A, s: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    // ark_se and ark_de taken from https://github.com/arkworks-rs/algebra/issues/178#issuecomment-1413219278
    let mut bytes = vec![];
    a.serialize_with_mode(&mut bytes, Compress::Yes)
        .map_err(serde::ser::Error::custom)?;
    s.serialize_bytes(&bytes)
}

pub(crate) fn ark_de<'de, D, A: CanonicalDeserialize>(data: D) -> Result<A, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    let s: Vec<u8> = serde::de::Deserialize::deserialize(data)?;
    let a = A::deserialize_with_mode(s.as_slice(), Compress::Yes, Validate::Yes);
    a.map_err(serde::de::Error::custom)
}

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
    fn r1cs_matrices() -> (Vec<Vec<i32>>, Vec<Vec<i32>>, Vec<Vec<i32>>) {
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

        (l, r, o)
    }

    #[test]
    fn serialisation_matches() -> Result<(), Report> {
        init();
        let (l, r, o) = r1cs_matrices();

        debug!("Matrices initialised");
        let r1cs: R1CS<Field> = R1CS::new(l, r, o, Vec::<i32>::new());

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
        assert!(qap.verify(&w));
        let proof = trusted_setup.prove(&w)?;

        debug!("Proof generated");
        assert!(proof.verify(&trusted_setup, &qap.public_witness));

        let trusted_setup_serialized = serde_json::to_string(&trusted_setup)?;
        let trusted_setup_deserialized = serde_json::from_str(&trusted_setup_serialized)?;
        assert_eq!(trusted_setup, trusted_setup_deserialized);

        let proof_serialized = serde_json::to_string(&proof)?;
        let proof_deserialized = serde_json::from_str(&proof_serialized)?;
        assert_eq!(proof, proof_deserialized);
        Ok(())
    }
}
