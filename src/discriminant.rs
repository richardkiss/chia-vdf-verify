//! CreateDiscriminant: generate a negative discriminant from a seed.
//!
//! Port of chiavdf/src/create_discriminant.h.

use crate::primetest::hash_prime;
use num_bigint::BigInt;

/// Create a discriminant D from a seed and bit length.
/// D = -HashPrime(seed, length, {0, 1, 2, length-1})
/// D ≡ 7 (mod 8), so D ≡ 1 (mod 8) after negation... wait.
/// Actually D = -p where p ≡ 7 (mod 8), so D ≡ 1 (mod 8).
/// The bitmask {0, 1, 2, length-1} sets:
///   - bit 0: forces odd
///   - bit 1, 2: forces 7 mod 8 (bits 1 and 2 both set → value & 7 has bits 1,2,0 set = 7)
///   - bit length-1: sets the MSB to ensure exact bit length
pub fn create_discriminant(seed: &[u8], length: usize) -> BigInt {
    assert!(
        length > 0 && length % 8 == 0,
        "length must be positive multiple of 8"
    );
    let bitmask = vec![0usize, 1, 2, length - 1];
    let p = hash_prime(seed, length, &bitmask);
    -p
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primetest::is_prime_bpsw;
    use num_traits::Signed;

    #[test]
    fn test_discriminant_is_negative() {
        let seed = b"test_seed";
        let d = create_discriminant(seed, 512);
        assert!(d.is_negative(), "discriminant should be negative");
    }

    #[test]
    fn test_discriminant_mod8() {
        let seed = b"test_seed";
        let d = create_discriminant(seed, 512);
        // D = -p, p ≡ 7 mod 8 → D ≡ -7 ≡ 1 mod 8
        use num_bigint::BigInt;
        use num_integer::Integer;
        let r = d.mod_floor(&BigInt::from(8));
        assert_eq!(r, BigInt::from(1), "discriminant should be ≡ 1 mod 8");
    }

    #[test]
    fn test_discriminant_is_prime_magnitude() {
        let seed = b"small_test";
        let d = create_discriminant(seed, 256);
        assert!(
            is_prime_bpsw(&(-d)),
            "discriminant magnitude should be prime"
        );
    }
}
