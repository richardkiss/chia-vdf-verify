//! BigInt wrapper with GMP-compatible operations.

use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

pub use num_bigint::BigInt as Int;

/// Number of bits in a machine limb (for Lehmer extraction).
pub const LIMB_BITS: usize = 64;

/// Compute floor division (matching GMP mpz_fdiv_q behavior).
/// GMP floors toward -infinity; Rust truncates toward zero.
pub fn fdiv_q(a: &BigInt, b: &BigInt) -> BigInt {
    a.div_floor(b)
}

/// Compute floor remainder (matching GMP mpz_fdiv_r).
pub fn fdiv_r(a: &BigInt, b: &BigInt) -> BigInt {
    a.mod_floor(b)
}

/// Compute truncating division (Rust default `/`).
pub fn tdiv_q(a: &BigInt, b: &BigInt) -> BigInt {
    a / b
}

/// Compute truncating remainder (Rust default `%`).
pub fn tdiv_r(a: &BigInt, b: &BigInt) -> BigInt {
    a % b
}

/// Exact division — panics (debug) or gives wrong result (release) if remainder != 0.
pub fn divexact(a: &BigInt, b: &BigInt) -> BigInt {
    let (q, r) = a.div_rem(b);
    debug_assert!(r.is_zero(), "divexact: remainder is not zero");
    q
}

/// Integer square root (floor).
pub fn isqrt(n: &BigInt) -> BigInt {
    assert!(!n.is_negative(), "isqrt: negative argument");
    BigInt::from(n.magnitude().sqrt())
}

/// nth root (floor).
pub fn nth_root(n: &BigInt, k: u32) -> BigInt {
    assert!(!n.is_negative(), "nth_root: negative argument");
    BigInt::from(n.magnitude().nth_root(k))
}

/// Jacobi symbol (a/n). n must be odd and positive.
pub fn jacobi(a: &BigInt, n: &BigInt) -> i32 {
    use num_bigint::BigUint;

    debug_assert!(n > &BigInt::one(), "n must be > 1, got {}", n);
    debug_assert!(n.is_odd(), "n must be odd, got {}", n);

    // Work with non-negative a via a mod n
    let mut a: BigUint = {
        let a_mod = a.mod_floor(n);
        // a_mod is in [0, n), and n > 0, so result is non-negative
        a_mod.magnitude().clone()
    };
    let mut n: BigUint = n.magnitude().clone();
    let mut result = 1i32;

    loop {
        // Reduce a mod n
        a = &a % &n;

        if a.is_zero() {
            return 0;
        }

        // Factor out powers of 2 from a
        let two_exp = a.trailing_zeros().unwrap_or(0);
        if two_exp > 0 {
            a >>= two_exp;
            // (2/n) = (-1)^((n^2-1)/8): negative iff n ≡ 3 or 5 (mod 8)
            if two_exp % 2 == 1 {
                let n_mod8 = (&n % 8u32).to_u32_digits().into_iter().next().unwrap_or(0);
                if n_mod8 == 3 || n_mod8 == 5 {
                    result = -result;
                }
            }
        }

        if a == BigUint::one() {
            return result;
        }

        // Apply quadratic reciprocity: (a/n)(n/a) = (-1)^((a-1)(n-1)/4)
        // Which negates result iff both a ≡ 3 (mod 4) and n ≡ 3 (mod 4)
        let a_mod4 = (&a % 4u32).to_u32_digits().into_iter().next().unwrap_or(0);
        let n_mod4 = (&n % 4u32).to_u32_digits().into_iter().next().unwrap_or(0);
        if a_mod4 == 3 && n_mod4 == 3 {
            result = -result;
        }

        // Swap: (a/n) -> (n mod a / a)
        let tmp = &n % &a;
        n = a;
        a = tmp;

        if n.is_zero() {
            return 0;
        }
    }
}

/// modpow: a^b mod m (all positive).
pub fn modpow(base: &BigInt, exp: &BigInt, modulus: &BigInt) -> BigInt {
    base.modpow(exp, modulus)
}

/// Import big-endian bytes as a non-negative integer.
pub fn from_bytes_be(bytes: &[u8]) -> BigInt {
    BigInt::from_bytes_be(Sign::Plus, bytes)
}

/// Export as big-endian bytes with the given byte length (zero-padded on left).
pub fn to_bytes_be_padded(n: &BigInt, len: usize) -> Vec<u8> {
    let (_, bytes) = n.to_bytes_be();
    if bytes.len() >= len {
        bytes[bytes.len() - len..].to_vec()
    } else {
        let mut out = vec![0u8; len];
        out[len - bytes.len()..].copy_from_slice(&bytes);
        out
    }
}

/// Import little-endian bytes as a non-negative integer.
pub fn from_bytes_le(bytes: &[u8]) -> BigInt {
    BigInt::from_bytes_le(Sign::Plus, bytes)
}

/// Export as little-endian bytes with the given byte length (zero-padded on right).
pub fn to_bytes_le_padded(n: &BigInt, len: usize) -> Vec<u8> {
    let (_, bytes) = n.to_bytes_le();
    if bytes.len() >= len {
        bytes[..len].to_vec()
    } else {
        let mut out = bytes;
        out.resize(len, 0);
        out
    }
}

/// Number of bits needed to represent the absolute value.
pub fn num_bits(n: &BigInt) -> usize {
    if n.is_zero() {
        return 1;
    }
    n.bits() as usize
}

/// Extract a 64-bit signed integer and exponent approximation of n,
/// mirroring the C `mpz_get_si_2exp` used in the Reducer.
/// Returns (mantissa, exponent) where n ≈ mantissa * 2^(exponent - 63).
pub fn get_si_2exp(n: &BigInt) -> (i64, i64) {
    if n.is_zero() {
        return (0, 0);
    }
    let bits = num_bits(n) as i64;
    // Get the top 64 bits
    let shift = if bits > 64 { bits - 64 } else { 0 };
    let (_, mag_bytes) = n.to_bytes_be();
    // Extract via shifting
    let shifted = n.magnitude() >> shift as usize;
    let (_, top_bytes) = BigInt::from(shifted).to_bytes_be();
    let mut top = 0u64;
    for b in top_bytes.iter().take(8) {
        top = (top << 8) | (*b as u64);
    }
    // Normalize to 63 bits for signed representation
    let lg2 = 64 - top.leading_zeros() as i64;
    let mantissa_shift = 63 - lg2;
    let mantissa = if mantissa_shift >= 0 {
        (top << mantissa_shift) as i64
    } else {
        (top >> (-mantissa_shift)) as i64
    };
    let exp = bits;
    let mantissa = if n.is_negative() { -mantissa } else { mantissa };
    (mantissa, exp)
}

/// Extract the low word of (|x| >> shift_bits), where x is non-negative.
pub fn extract_uword_from_shift_nonneg(x: &BigInt, shift_bits: i64) -> i64 {
    if shift_bits <= 0 {
        let (_, digits) = x.to_u64_digits();
        return digits.first().copied().unwrap_or(0) as i64;
    }
    let shifted = x.magnitude() >> shift_bits as usize;
    let digits = BigUint::from(shifted).to_u64_digits();
    digits.first().copied().unwrap_or(0) as i64
}

/// Get bit length of the absolute value (matching chiavdf_mpz_bitlen_nonneg).
pub fn bitlen_nonneg(x: &BigInt) -> i64 {
    if x.is_zero() {
        return 1;
    }
    x.bits() as i64
}

/// Trailing zeros in the absolute value (number of factors of 2).
pub fn trailing_zeros(n: &BigInt) -> u64 {
    n.trailing_zeros().unwrap_or(0)
}

/// Lehmer-accelerated full extended GCD.
///
/// Returns (gcd, x, y) such that gcd = x * a + y * b.
/// Roughly 4–5x faster than `num_integer::extended_gcd` for multi-limb inputs.
pub fn fast_extended_gcd(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    let a_neg = a.is_negative();
    let b_neg = b.is_negative();
    let mut r0 = a.abs();
    let mut r1 = b.abs();

    // Maintain: r0 = s0 * |a| + t0 * |b|
    //           r1 = s1 * |a| + t1 * |b|
    let mut s0 = BigInt::one();  let mut t0 = BigInt::zero();
    let mut s1 = BigInt::zero(); let mut t1 = BigInt::one();

    while !r1.is_zero() {
        let (p, q, r, s, steps) = lehmer_inner_loop(&r0, &r1);
        if steps == 0 {
            let qq = &r0 / &r1;
            let rem = &r0 - &qq * &r1;
            r0 = std::mem::replace(&mut r1, rem);
            let ns0 = s0 - &qq * &s1; let old = std::mem::replace(&mut s1, ns0); s0 = old;
            let nt0 = t0 - &qq * &t1; let old = std::mem::replace(&mut t1, nt0); t0 = old;
        } else {
            let nr0 = BigInt::from(p) * &r0 + BigInt::from(q) * &r1;
            let nr1 = BigInt::from(r) * &r0 + BigInt::from(s) * &r1;
            let ns0 = BigInt::from(p) * &s0 + BigInt::from(q) * &s1;
            let ns1 = BigInt::from(r) * &s0 + BigInt::from(s) * &s1;
            let nt0 = BigInt::from(p) * &t0 + BigInt::from(q) * &t1;
            let nt1 = BigInt::from(r) * &t0 + BigInt::from(s) * &t1;
            // Fix signs: keep residues non-negative
            if nr1.is_negative() {
                r1 = -nr1; s1 = -ns1; t1 = -nt1;
            } else {
                r1 = nr1; s1 = ns1; t1 = nt1;
            }
            if nr0.is_negative() {
                r0 = -nr0; s0 = -ns0; t0 = -nt0;
            } else {
                r0 = nr0; s0 = ns0; t0 = nt0;
            }
        }
    }

    let x = if a_neg { -s0 } else { s0 };
    let y = if b_neg { -t0 } else { t0 };
    (r0, x, y)
}

/// Half extended GCD: returns (gcd, y) where gcd ≡ y * b (mod a).
/// Only tracks the coefficient of b — about 30% faster than `fast_extended_gcd`
/// when you only need one Bezout coefficient.
pub fn fast_gcd_coeff_b(a: &BigInt, b: &BigInt) -> (BigInt, BigInt) {
    let b_neg = b.is_negative();
    let mut r0 = a.abs();
    let mut r1 = b.abs();

    // Only track: r0 = ??? + t0 * |b|,  r1 = ??? + t1 * |b|
    let mut t0 = BigInt::zero();
    let mut t1 = BigInt::one();

    while !r1.is_zero() {
        let (p, q, r, s, steps) = lehmer_inner_loop(&r0, &r1);
        if steps == 0 {
            let qq = &r0 / &r1;
            let rem = &r0 - &qq * &r1;
            r0 = std::mem::replace(&mut r1, rem);
            let nt0 = t0 - &qq * &t1; let old = std::mem::replace(&mut t1, nt0); t0 = old;
        } else {
            let nr0 = BigInt::from(p) * &r0 + BigInt::from(q) * &r1;
            let nr1 = BigInt::from(r) * &r0 + BigInt::from(s) * &r1;
            let nt0 = BigInt::from(p) * &t0 + BigInt::from(q) * &t1;
            let nt1 = BigInt::from(r) * &t0 + BigInt::from(s) * &t1;
            if nr1.is_negative() { r1 = -nr1; t1 = -nt1; } else { r1 = nr1; t1 = nt1; }
            if nr0.is_negative() { r0 = -nr0; t0 = -nt0; } else { r0 = nr0; t0 = nt0; }
        }
    }

    let y = if b_neg { -t0 } else { t0 };
    (r0, y)
}

/// Lehmer inner loop: compute 2×2 matrix of small scalars [[p,q],[r,s]].
/// Invariant: r0_new ≈ p*r0+q*r1,  r1_new ≈ r*r0+s*r1 (using top-word approximation).
#[inline]
fn lehmer_inner_loop(r0: &BigInt, r1: &BigInt) -> (i64, i64, i64, i64, usize) {
    let bits0 = r0.bits() as i64;
    let bits1 = r1.bits() as i64;
    let shift = std::cmp::max(std::cmp::max(bits0, bits1) - LIMB_BITS as i64 + 1, 0) as usize;
    let mut rr0 = extract_word_unsigned(r0, shift);
    let mut rr1 = extract_word_unsigned(r1, shift);
    let mut p: i64 = 1; let mut q: i64 = 0;
    let mut r: i64 = 0; let mut s: i64 = 1;
    let mut i = 0usize;
    loop {
        if rr1 == 0 { break; }
        let qq = rr0 / rr1;
        let t1_ = rr0 - qq * rr1;
        let tp = p - qq * r;
        let tq = q - qq * s;
        if i & 1 == 0 {
            if t1_ < -tq || rr1 - t1_ < tp - p { break; }
        } else {
            if t1_ < -tp || rr1 - t1_ < tq - q { break; }
        }
        rr0 = rr1; rr1 = t1_;
        p = r; q = s; r = tp; s = tq;
        i += 1;
    }
    (p, q, r, s, i)
}

/// Extract the low 64 bits of (n >> shift) treating n as non-negative.
fn extract_word_unsigned(n: &BigInt, shift: usize) -> i64 {
    if shift == 0 {
        let (_, digits) = n.to_u64_digits();
        return digits.first().copied().unwrap_or(0) as i64;
    }
    let shifted = n.magnitude() >> shift;
    let digits = BigUint::from(shifted).to_u64_digits();
    digits.first().copied().unwrap_or(0) as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fast_extended_gcd() {
        let cases: Vec<(i64, i64, i64)> = vec![
            (5, 3, 1), (6, 4, 2), (100, 37, 1), (15, 10, 5),
            (7, 5, 1), (35, 20, 5), (1000000, 999983, 1),
            (0, 7, 7), (7, 0, 7), (1, 1, 1),
        ];
        for (a, b, expected_gcd) in cases {
            let ba = BigInt::from(a);
            let bb = BigInt::from(b);
            let (gcd, x, y) = fast_extended_gcd(&ba, &bb);
            assert_eq!(gcd, BigInt::from(expected_gcd),
                "gcd({},{}) wrong: got {}", a, b, gcd);
            assert_eq!(&x * &ba + &y * &bb, gcd,
                "Bezout identity failed for ({},{}): {}*{}+{}*{}≠{}",
                a, b, x, a, y, b, gcd);
        }
        // Test with negative inputs
        let (gcd, x, y) = fast_extended_gcd(&BigInt::from(-15i64), &BigInt::from(10i64));
        assert_eq!(gcd, BigInt::from(5));
        assert_eq!(&x * BigInt::from(-15i64) + &y * BigInt::from(10i64), gcd);
    }

    #[test]
    fn test_jacobi() {
        // Known values
        assert_eq!(jacobi(&BigInt::from(2), &BigInt::from(7)), 1);
        assert_eq!(jacobi(&BigInt::from(3), &BigInt::from(7)), -1);
        assert_eq!(jacobi(&BigInt::from(5), &BigInt::from(9)), 1); // 9 is composite but jacobi still defined
        assert_eq!(jacobi(&BigInt::from(0), &BigInt::from(7)), 0);
        assert_eq!(jacobi(&BigInt::from(1), &BigInt::from(7)), 1);
    }

    #[test]
    fn test_fdiv() {
        // -7 / 2 = -4 (floor) vs -3 (truncating)
        let a = BigInt::from(-7i64);
        let b = BigInt::from(2i64);
        assert_eq!(fdiv_q(&a, &b), BigInt::from(-4i64));
        assert_eq!(fdiv_r(&a, &b), BigInt::from(1i64));
    }

    #[test]
    fn test_isqrt() {
        assert_eq!(isqrt(&BigInt::from(16)), BigInt::from(4));
        assert_eq!(isqrt(&BigInt::from(15)), BigInt::from(3));
        assert_eq!(isqrt(&BigInt::from(0)), BigInt::from(0));
    }
}
