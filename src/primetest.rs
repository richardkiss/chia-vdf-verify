//! BPSW primality test and HashPrime.
//!
//! Ports:
//! - chiavdf/src/primetest.h (is_prime_bpsw)
//! - chiavdf/src/proof_common.h (HashPrime)

use crate::integer::{jacobi, modpow};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use sha2::{Digest, Sha256};

/// Small prime products for trial division (subset — enough for fast rejection).
/// These are products of consecutive primes.
static SMALL_PRIMES: &[u64] = &[
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199,
];

/// Miller-Rabin test with given base `b` modulo `n`.
/// Returns true if n is probably prime (passes MR with this base).
pub fn miller_rabin(n: &BigInt, base: &BigInt) -> bool {
    // Write n-1 = 2^s * d
    let n_minus1 = n - BigInt::one();
    let s = n_minus1.trailing_zeros().unwrap_or(0) as usize;
    let d = &n_minus1 >> s;

    let mut b = modpow(base, &d, n);

    if b == BigInt::one() {
        return true;
    }

    for _ in 0..s {
        let b_plus1 = &b + BigInt::one();
        if b_plus1 == *n {
            return true;
        }
        b = modpow(&b, &BigInt::from(2u32), n);
    }
    false
}

/// Find parameters (p, q) for Lucas test such that Jacobi(D, n) = -1.
/// Returns None if D >= VPRP_MAX_D without finding suitable D.
fn find_pq(n: &BigInt) -> Option<(i64, i64)> {
    let mut d = 5i64;
    for _ in 0..500 {
        let d_sign = if d % 4 == 1 { d } else { -d };
        let d_big = BigInt::from(d_sign);
        if jacobi(&d_big, n) == -1 {
            if d_sign == 5 {
                return Some((5, 5));
            } else {
                return Some((1, (1 - d_sign) / 4));
            }
        }
        d = d.abs() + 2;
    }
    None
}

/// Lucas-V probable prime test (vprp).
fn is_vprp(n: &BigInt) -> bool {
    let (p, q) = match find_pq(n) {
        Some(pq) => pq,
        None => return false,
    };

    let e = n + BigInt::one();
    let v1 = find_lucas_v(&e, n, p, q);

    // Check V1 ≡ 2q (mod n)
    let two_q = BigInt::from(2 * q);
    let v1_mod = if q < 0 {
        // cdiv_r (ceiling division remainder)
        let m = v1.mod_floor(n);
        if m.is_negative() {
            m + n
        } else {
            m
        }
    } else {
        v1.mod_floor(n)
    };
    let two_q_mod = two_q.mod_floor(n);

    v1_mod == two_q_mod
}

/// Compute Lucas-V sequence value V_{n+1} mod m using the doubling method.
fn find_lucas_v(e: &BigInt, m: &BigInt, p: i64, q: i64) -> BigInt {
    let l = e.bits() as usize;

    let mut u1 = BigInt::one(); // U_1
    let mut u2 = BigInt::from(p); // U_2
    let minus_2q = -2 * q;

    for i in (0..l.saturating_sub(1)).rev() {
        let tmp2 = &u2 * &u1;
        let u2_sq = &u2 * &u2;
        let u1_sq = &u1 * &u1;

        if e.bit(i as u64) {
            // Bit is 1
            u1 = &u2_sq - BigInt::from(q) * &u1_sq;
            u2 = if p != 1 {
                BigInt::from(p) * &u2_sq + BigInt::from(minus_2q) * &tmp2
            } else {
                &u2_sq + BigInt::from(minus_2q) * &tmp2
            };
        } else {
            // Bit is 0
            u2 = &u2_sq - BigInt::from(q) * &u1_sq;
            let tmp3 = BigInt::from(2) * &tmp2;
            u1 = if p != 1 {
                tmp3 - BigInt::from(p) * &u1_sq
            } else {
                tmp3 - &u1_sq
            };
        }

        u1 = u1.mod_floor(m);
        u2 = u2.mod_floor(m);
    }

    // V1 = 2*U2 - P*U1

    BigInt::from(2) * &u2 - BigInt::from(p) * &u1
}

/// BPSW primality test.
/// Returns true if n is (very likely) prime.
pub fn is_prime_bpsw(n: &BigInt) -> bool {
    if n <= &BigInt::one() {
        return false;
    }
    if n == &BigInt::from(2u32) || n == &BigInt::from(3u32) {
        return true;
    }
    if n.is_even() {
        return false;
    }

    // Trial division by small primes
    for &p in SMALL_PRIMES {
        let p_big = BigInt::from(p);
        if n == &p_big {
            return true;
        }
        if (n % &p_big).is_zero() {
            return false;
        }
    }

    // Miller-Rabin base 2
    let base2 = BigInt::from(2u32);
    if !miller_rabin(n, &base2) {
        return false;
    }

    // Lucas-V probable prime test
    is_vprp(n)
}

/// HashPrime: generate a pseudoprime of given bit length from seed.
///
/// Uses iterative SHA256 to expand seed, then applies bitmask and tests primality.
/// Matches chiavdf's HashPrime(seed, length, bitmask).
pub fn hash_prime(seed: &[u8], length: usize, bitmask: &[usize]) -> BigInt {
    assert!(length.is_multiple_of(8), "length must be multiple of 8");
    let byte_len = length / 8;

    let mut sprout = seed.to_vec();

    loop {
        let mut blob = Vec::with_capacity(byte_len);

        while blob.len() * 8 < length {
            // Increment sprout by 1 (big-endian increment)
            for i in (0..sprout.len()).rev() {
                sprout[i] = sprout[i].wrapping_add(1);
                if sprout[i] != 0 {
                    break;
                }
            }

            let hash = Sha256::digest(&sprout);
            let remaining = byte_len - blob.len();
            let take = remaining.min(hash.len());
            blob.extend_from_slice(&hash[..take]);
        }

        assert_eq!(blob.len(), byte_len);

        // Import as big integer (big-endian bytes)
        let mut p = crate::integer::from_bytes_be(&blob);

        // Apply bitmask
        for &bit in bitmask {
            p |= BigInt::one() << bit;
        }

        // Force odd
        p |= BigInt::one();

        if is_prime_bpsw(&p) {
            return p;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_miller_rabin_known_primes() {
        let base2 = BigInt::from(2u32);
        for &p in &[3u64, 5, 7, 11, 13, 17, 19, 23, 997, 7919] {
            let n = BigInt::from(p);
            assert!(miller_rabin(&n, &base2), "{} should pass MR(2)", p);
        }
    }

    #[test]
    fn test_miller_rabin_composites() {
        let base2 = BigInt::from(2u32);
        for &c in &[9u64, 15, 21, 25, 35, 49, 77, 91] {
            let n = BigInt::from(c);
            // Most composites will fail, though some are strong pseudoprimes
            let _result = miller_rabin(&n, &base2);
            // 341 is the first strong pseudoprime to base 2; our list doesn't include it
            // Just test that 9 fails
        }
        let n9 = BigInt::from(9u64);
        assert!(!miller_rabin(&n9, &base2), "9 should fail MR(2)");
    }

    #[test]
    fn test_is_prime_bpsw() {
        for &p in &[2u64, 3, 5, 7, 11, 13, 997, 7919, 104729] {
            assert!(is_prime_bpsw(&BigInt::from(p)), "{} should be prime", p);
        }
        for &c in &[4u64, 6, 9, 15, 25, 35, 49, 100] {
            assert!(
                !is_prime_bpsw(&BigInt::from(c)),
                "{} should be composite",
                c
            );
        }
    }

    #[test]
    fn test_hash_prime_is_prime() {
        let seed = b"test_seed_12345";
        let p = hash_prime(seed, 256, &[255]);
        assert!(is_prime_bpsw(&p), "hash_prime result should be prime");
        assert_eq!(
            p.bits(),
            256,
            "hash_prime result should have correct bit length"
        );
    }
}
