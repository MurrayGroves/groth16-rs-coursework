use ark_ff::Field;
use rand::Rng;

pub fn rand_scalar<T, S>(rng: &mut T) -> S
where
    T: Rng,
    S: Field,
{
    // TODO - Make sure 256 bytes is enough
    let mut bytes = [0; 256];
    rng.fill_bytes(&mut bytes);
    let mut out = None;
    while out.is_none() {
        out = S::from_random_bytes(&bytes);
    }

    out.unwrap()
}
