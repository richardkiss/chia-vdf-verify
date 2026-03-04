//! BQFC compressed form serialization/deserialization.
//!
//! Port of chiavdf/src/bqfc.c and bqfc.h.
//!
//! Format (100 bytes for 1024-bit discriminant):
//!   Byte 0: flags (b_sign, t_sign, is_identity, is_generator)
//!   Byte 1: g_size (size of 'g' in bytes minus 1)
//!   d/16 - g_size bytes: a' = a/g
//!   d/32 - g_size bytes: t' = t/g
//!   g_size+1 bytes: g
//!   g_size+1 bytes: b0

use crate::integer::{divexact, fdiv_r, from_bytes_le, isqrt, tdiv_q, to_bytes_le_padded};
use crate::xgcd_partial::xgcd_partial;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, Zero};

/// Size of the serialized form (100 bytes for 1024-bit max discriminant).
pub const BQFC_FORM_SIZE: usize = (1024 + 31) / 32 * 3 + 4;

/// Flag bits
const BQFC_B_SIGN: u8 = 1 << 0;
const BQFC_T_SIGN: u8 = 1 << 1;
const BQFC_IS_1: u8 = 1 << 2;
const BQFC_IS_GEN: u8 = 1 << 3;

/// Compressed form intermediate representation.
struct QfbC {
    a: BigInt,
    t: BigInt,
    g: BigInt,
    b0: BigInt,
    b_sign: bool,
}

/// Compress (a, b) to intermediate representation.
fn bqfc_compr(a: &BigInt, b: &BigInt) -> QfbC {
    if a == b {
        return QfbC {
            a: a.clone(),
            t: BigInt::zero(),
            g: BigInt::zero(),
            b0: BigInt::zero(),
            b_sign: false,
        };
    }

    let sign = b.is_negative();
    let a_sqrt = isqrt(a);

    let mut a_copy = a.clone();
    let mut b_copy = if sign { -b } else { b.clone() };

    let mut dummy = BigInt::zero();
    let mut t = BigInt::zero();

    xgcd_partial(&mut dummy, &mut t, &mut a_copy, &mut b_copy, &a_sqrt);
    // xgcd_partial sets the opposite sign for t
    t = -t;

    let g = a.gcd(&t);

    let (a_out, b0) = if g == BigInt::from(1u32) {
        (a.clone(), BigInt::zero())
    } else {
        let a_new = divexact(a, &g);
        let t_new = divexact(&t, &g);
        let b0 = tdiv_q(b, &a_new);
        let b0 = if sign { -b0 } else { b0 };
        t = t_new;
        (a_new, b0)
    };

    QfbC {
        a: a_out,
        t,
        g,
        b_sign: sign,
        b0,
    }
}

/// Decompress intermediate representation to (a, b) given discriminant D.
fn bqfc_decompr(d: &BigInt, c: &QfbC) -> Result<(BigInt, BigInt), String> {
    if c.t.is_zero() {
        return Ok((c.a.clone(), c.a.clone()));
    }

    let t = if c.t.is_negative() {
        &c.t + &c.a
    } else {
        c.t.clone()
    };

    if c.a.is_zero() {
        return Err("bqfc_decompr: a is zero".to_string());
    }

    // t_inv = modular inverse of t mod a
    let egcd = t.extended_gcd(&c.a);
    if egcd.gcd != BigInt::from(1u32) {
        return Err(format!("bqfc_decompr: gcd(t, a) = {} != 1", egcd.gcd));
    }
    let t_inv = if egcd.x.is_negative() {
        egcd.x + &c.a
    } else {
        egcd.x
    };

    // d_mod_a = D mod a
    let d_mod_a = fdiv_r(d, &c.a);

    // tmp = sqrt(t^2 * d_mod_a mod a)
    let t_sq_mod = c.t.modpow(&BigInt::from(2u32), &c.a);
    let tmp_prod = &t_sq_mod * &d_mod_a;
    let tmp_mod = &tmp_prod % &c.a; // tdiv_r

    // Check perfect square
    let tmp_sqrt = isqrt(&tmp_mod);
    if &tmp_sqrt * &tmp_sqrt != tmp_mod {
        return Err("bqfc_decompr: not a perfect square".to_string());
    }

    // out_b = tmp_sqrt * t_inv mod a
    let out_b = (&tmp_sqrt * &t_inv) % &c.a; // tdiv_r

    let out_a = if c.g > BigInt::from(1u32) {
        &c.a * &c.g
    } else {
        c.a.clone()
    };

    let out_b = if c.b0 > BigInt::zero() {
        out_b + &c.a * &c.b0
    } else {
        out_b
    };

    let out_b = if c.b_sign { -out_b } else { out_b };

    Ok((out_a, out_b))
}

/// Export n as little-endian bytes of exactly `size` bytes.
fn export_le(out: &mut Vec<u8>, n: &BigInt, size: usize) {
    let bytes = to_bytes_le_padded(n, size);
    out.extend_from_slice(&bytes);
}

/// Serialize (a, b) with discriminant bit length d_bits to 100-byte output.
/// This corresponds to bqfc_serialize.
pub fn serialize(a: &BigInt, b: &BigInt, d_bits: usize) -> Vec<u8> {
    let mut out = vec![0u8; BQFC_FORM_SIZE];

    // Special cases: identity (1,1) and generator (2,1)
    if b == &BigInt::from(1u32) && a <= &BigInt::from(2u32) {
        out[0] = if a == &BigInt::from(2u32) {
            BQFC_IS_GEN
        } else {
            BQFC_IS_1
        };
        return out;
    }

    let d_bits_rounded = (d_bits + 31) & !31usize;
    let c = bqfc_compr(a, b);
    let valid_size = bqfc_get_compr_size(d_bits);

    let mut buf = Vec::with_capacity(valid_size);
    let mut flags = 0u8;
    if c.b_sign {
        flags |= BQFC_B_SIGN;
    }
    if c.t.is_negative() {
        flags |= BQFC_T_SIGN;
    }

    let g_size = if c.g.is_zero() {
        0usize
    } else {
        let bits = c.g.bits() as usize;
        (bits + 7) / 8
    };
    let g_size = if g_size == 0 { 0 } else { g_size - 1 };

    buf.push(flags);
    buf.push(g_size as u8);

    export_le(&mut buf, &c.a, d_bits_rounded / 16 - g_size);
    let t_abs = if c.t.is_negative() {
        -c.t.clone()
    } else {
        c.t.clone()
    };
    export_le(&mut buf, &t_abs, d_bits_rounded / 32 - g_size);
    export_le(&mut buf, &c.g, g_size + 1);
    export_le(&mut buf, &c.b0, g_size + 1);

    // Copy into fixed-size output
    let copy_len = buf.len().min(BQFC_FORM_SIZE);
    out[..copy_len].copy_from_slice(&buf[..copy_len]);
    out
}

/// Deserialize a form from 100-byte input with discriminant D.
/// Returns (a, b) or error string.
pub fn deserialize(d: &BigInt, data: &[u8]) -> Result<(BigInt, BigInt), String> {
    if data.len() != BQFC_FORM_SIZE {
        return Err(format!(
            "expected {} bytes, got {}",
            BQFC_FORM_SIZE,
            data.len()
        ));
    }

    // Check special forms
    if data[0] & (BQFC_IS_1 | BQFC_IS_GEN) != 0 {
        let a = if data[0] & BQFC_IS_GEN != 0 {
            BigInt::from(2u32)
        } else {
            BigInt::from(1u32)
        };
        return Ok((a, BigInt::from(1u32)));
    }

    let d_bits = crate::integer::num_bits(d);
    let d_bits_rounded = (d_bits + 31) & !31usize;

    let g_size = data[1] as usize;
    if g_size >= d_bits_rounded / 32 {
        return Err("g_size out of range".to_string());
    }

    let mut offset = 2usize;

    let a_bytes = d_bits_rounded / 16 - g_size;
    let a = from_bytes_le(&data[offset..offset + a_bytes]);
    offset += a_bytes;

    let t_bytes = d_bits_rounded / 32 - g_size;
    let t_raw = from_bytes_le(&data[offset..offset + t_bytes]);
    offset += t_bytes;

    let g = from_bytes_le(&data[offset..offset + g_size + 1]);
    offset += g_size + 1;

    let b0 = from_bytes_le(&data[offset..offset + g_size + 1]);

    let b_sign = data[0] & BQFC_B_SIGN != 0;
    let t_sign = data[0] & BQFC_T_SIGN != 0;
    let t = if t_sign { -t_raw } else { t_raw };

    let c = QfbC {
        a,
        t,
        g,
        b0,
        b_sign,
    };

    let (out_a, out_b) = bqfc_decompr(d, &c)?;

    // Verify canonical serialization
    let canon = serialize(&out_a, &out_b, d_bits);
    if canon != data {
        return Err("non-canonical serialization".to_string());
    }

    Ok((out_a, out_b))
}

/// Compute the serialization size for a given discriminant bit length.
pub fn bqfc_get_compr_size(d_bits: usize) -> usize {
    let d_bits_rounded = (d_bits + 31) & !31usize;
    d_bits_rounded / 32 * 3 + 4
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::discriminant::create_discriminant;
    use crate::form::Form;
    use crate::reducer::reduce;

    #[test]
    fn test_bqfc_form_size() {
        // For 1024-bit discriminant: (1024+31)/32 = 32, *3 + 4 = 100
        assert_eq!(bqfc_get_compr_size(1024), 100);
        assert_eq!(BQFC_FORM_SIZE, 100);
    }

    #[test]
    fn test_serialize_identity() {
        let d = BigInt::from(-47i64);
        let a = BigInt::from(1u32);
        let b = BigInt::from(1u32);
        let data = serialize(&a, &b, crate::integer::num_bits(&d));
        assert_eq!(data[0], BQFC_IS_1);
    }

    #[test]
    fn test_serialize_generator() {
        let d = BigInt::from(-47i64);
        let a = BigInt::from(2u32);
        let b = BigInt::from(1u32);
        let data = serialize(&a, &b, crate::integer::num_bits(&d));
        assert_eq!(data[0], BQFC_IS_GEN);
    }

    #[test]
    fn test_bqfc_roundtrip_identity() {
        let d = BigInt::from(-47i64);
        let d_bits = crate::integer::num_bits(&d);
        let a = BigInt::from(1u32);
        let b = BigInt::from(1u32);
        let data = serialize(&a, &b, d_bits);
        let (ra, rb) = deserialize(&d, &data).unwrap();
        assert_eq!(ra, a);
        assert_eq!(rb, b);
    }
}
