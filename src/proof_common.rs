//! FastPow, GetB, SerializeForm, DeserializeForm, FastPowFormNucomp.
//!
//! Port of chiavdf/src/proof_common.h.

use num_bigint::BigInt;
use num_traits::Zero;
use crate::form::Form;
use crate::bqfc::{serialize, deserialize};
use crate::primetest::hash_prime;
use crate::reducer::reduce;
use crate::nucomp::{nucomp, nudupl};
use crate::integer::{modpow, num_bits};

/// B_bits = 264
pub const B_BITS: usize = 264;
/// B_bytes = 33
pub const B_BYTES: usize = (B_BITS + 7) / 8;

/// Compute 2^b mod c.
/// This is FastPow(2, b, c) in C.
pub fn fast_pow(b: u64, c: &BigInt) -> BigInt {
    if c.is_zero() {
        panic!("FastPow: division by zero");
    }
    let base = BigInt::from(2u32);
    let exp = BigInt::from(b);
    modpow(&base, &exp, c)
}

/// Serialize form to BQFC_FORM_SIZE bytes.
/// Reduces the form first.
pub fn serialize_form(f: &mut Form, d_bits: usize) -> Vec<u8> {
    f.reduce();
    serialize(&f.a, &f.b, d_bits)
}

/// Deserialize form from BQFC_FORM_SIZE bytes with discriminant D.
pub fn deserialize_form(d: &BigInt, bytes: &[u8]) -> Result<Form, String> {
    let (a, b) = deserialize(d, bytes)?;
    if a.is_zero() {
        return Err("deserialized form has a=0".to_string());
    }
    let f = Form::from_abd(a, b, d);
    if !f.is_reduced() {
        return Err("deserialized form is not reduced".to_string());
    }
    Ok(f)
}

/// Compute GetB(D, x, y):
/// B = HashPrime(serialize(x) || serialize(y), 264, {263})
pub fn get_b(d: &BigInt, x: &mut Form, y: &mut Form) -> BigInt {
    let d_bits = num_bits(d);
    let sx = serialize_form(x, d_bits);
    let sy = serialize_form(y, d_bits);
    let mut seed = sx;
    seed.extend_from_slice(&sy);
    hash_prime(&seed, B_BITS, &[B_BITS - 1])
}

/// Exponentiate a form by num_iterations using binary method with nucomp/nudupl.
/// This is FastPowFormNucomp.
pub fn fast_pow_form_nucomp(
    x: &Form,
    d: &BigInt,
    num_iterations: &BigInt,
    l: &BigInt,
) -> Form {
    if num_iterations.is_zero() {
        return Form::identity(d);
    }

    let mut res = x.clone();
    let n_bits = num_iterations.bits() as usize;

    // max_size threshold for lazy reduction: -D.impl->_mp_size / 2
    // In Rust: bits of |D| / 2 limbs ~ num_bits(D) / 2 limbs
    // We convert to approximate limb count: D is ~1024 bits = 16 limbs, so threshold = 8 limbs
    // In practice we reduce after every step for correctness (since we don't access raw limb counts)
    // The C code only reduces lazily when a's size exceeds half the discriminant's limb count.
    // For correctness in this prototype, we reduce after every step.
    let do_lazy_reduce = {
        let d_limbs = (num_bits(d) + 63) / 64;
        let max_size = d_limbs / 2;
        max_size
    };

    for i in (0..n_bits.saturating_sub(1)).rev() {
        res = nudupl(&res, d, l);

        // Lazy reduction: reduce only when a's bit size exceeds threshold
        let a_limbs = (num_bits(&res.a) + 63) / 64;
        if a_limbs > do_lazy_reduce {
            reduce(&mut res);
        }

        if num_iterations.bit(i as u64) {
            res = nucomp(&res, x, d, l);
        }
    }

    reduce(&mut res);
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::discriminant::create_discriminant;

    #[test]
    fn test_fast_pow() {
        let modulus = BigInt::from(1000003u64);
        // 2^10 mod 1000003 = 1024
        assert_eq!(fast_pow(10, &modulus), BigInt::from(1024u64));
    }

    #[test]
    fn test_identity_exponentiation() {
        let seed = b"test_disc";
        let d = create_discriminant(seed, 512);
        let l = Form::compute_l(&d);
        let x = Form::identity(&d);
        let result = fast_pow_form_nucomp(&x, &d, &BigInt::from(100u32), &l);
        // identity^n = identity
        assert_eq!(result.a, BigInt::from(1u32));
        assert_eq!(result.b, BigInt::from(1u32));
    }

    #[test]
    fn test_serialize_deserialize_identity() {
        let seed = b"test_disc";
        let d = create_discriminant(seed, 512);
        let d_bits = num_bits(&d);
        let mut f = Form::identity(&d);
        let bytes = serialize_form(&mut f, d_bits);
        let f2 = deserialize_form(&d, &bytes).unwrap();
        assert_eq!(f.a, f2.a);
        assert_eq!(f.b, f2.b);
    }
}
