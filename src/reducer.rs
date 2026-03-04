//! Pulmark form reducer.
//!
//! Port of chiavdf/src/Reducer.h.
//! Reduces a quadratic form (a, b, c) to its canonical reduced representative.

use crate::form::Form;
use crate::integer::get_si_2exp;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::Signed;

const THRESH: i64 = 1i64 << 31;
const EXP_THRESH: i64 = 31;

/// Reduce form f in place.
pub fn reduce(f: &mut Form) {
    while !is_reduced(f) {
        let (a_val, a_exp) = get_si_2exp(&f.a);
        let (b_val, b_exp) = get_si_2exp(&f.b);
        let (c_val, c_exp) = get_si_2exp(&f.c);

        let max_exp = *[a_exp, b_exp, c_exp].iter().max().unwrap() + 1;
        let min_exp = *[a_exp, b_exp, c_exp].iter().min().unwrap();

        if max_exp - min_exp > EXP_THRESH {
            // Fall back to simple step
            reducer_simple(f);
            continue;
        }

        let a_sh = max_exp - a_exp;
        let b_sh = max_exp - b_exp;
        let c_sh = max_exp - c_exp;

        let a = a_val >> a_sh;
        let b = b_val >> b_sh;
        let c = c_val >> c_sh;

        let (u, v, w, x) = calc_uvwx(a, b, c);

        // Apply the (u,v,w,x) matrix to (a, b, c)
        let new_a = f.a.clone() * (u * u) + f.b.clone() * (u * w) + f.c.clone() * (w * w);
        let new_b =
            f.a.clone() * (2 * u * v) + f.b.clone() * (u * x + v * w) + f.c.clone() * (2 * w * x);
        let new_c = f.a.clone() * (v * v) + f.b.clone() * (v * x) + f.c.clone() * (x * x);

        f.a = new_a;
        f.b = new_b;
        f.c = new_c;
    }
}

/// Simple reducer step (used when exponent spread is large).
fn reducer_simple(f: &mut Form) {
    // s = floor((b + c) / (2c))  — this is the `mpz_mdiv` step
    // mdiv is "round away from zero", equivalent to "ceil(b/c)" when b/c could be negative
    // The C code: s = (c + b)/(2c) using mpz_mdiv which is ceiling division
    // Actually: mpz_mdiv(r, b, c) = ceiling(b/c)
    // Then: s = (r + 1) / 2
    let r = ceildiv(&f.b, &f.c);
    let r_plus1 = r + BigInt::from(1u32);
    let s = r_plus1 >> 1usize;

    let cs = &f.c * &s;
    let cs2 = &cs << 1usize;

    // m = cs - b
    let m = &cs - &f.b;

    // new_b = -b + 2cs = 2cs - b
    let new_b = &cs2 - &f.b;

    // new_c = old_a
    let old_a = f.a.clone();

    // new_a = old_c
    f.a = f.c.clone();

    // new_c = old_a + cs^2 - bs = old_a + s*m
    f.c = old_a + &s * &m;
    f.b = new_b;
}

/// Ceiling division: ceil(a/b).
fn ceildiv(a: &BigInt, b: &BigInt) -> BigInt {
    // ceil(a/b) = floor((a + b - 1) / b) for b > 0
    // But GMP's mpz_mdiv is "magnitude division", i.e., it rounds toward zero for the quotient
    // magnitude and then adjusts sign.
    // Actually from the C code context: mpz_mdiv(ctx.r, ctx.b, ctx.c) where c > 0 during
    // reduction (c is always positive). For c > 0:
    // mdiv(a, b) = ceil(a / b) when a/b >= 0, = floor(a/b) when a/b < 0
    // Which is actually the same as truncating division toward zero, then adding 1 if remainder > 0.
    // Let's use: mdiv = sign(a/b) * ceil(|a/b|) = truncating div with round-up magnitude.
    // For GMP mpz_mdiv: quotient * divisor >= dividend (rounds towards 0 or away from 0?).
    // Looking at GMP docs: mpz_mdiv truncates toward zero (same as cdiv_q for positive divisor).
    // For b > 0: cdiv_q rounds toward +inf.
    // Let's use: (a + b - 1) / b for a >= 0, and a / b (truncating) for a < 0 with b > 0.
    if b.is_positive() {
        a.div_ceil(b)
    } else {
        // This shouldn't happen in normal reduction (c is always positive)
        a / b
    }
}

/// Check if the form is reduced and normalize if needed.
/// Returns true if already reduced (but may have swapped a/c or negated b).
fn is_reduced(f: &mut Form) -> bool {
    use num_traits::Signed;
    let _abs_b = f.b.abs();

    let a_cmpabs_b = f.a.magnitude().cmp(f.b.magnitude());
    let c_cmpabs_b = f.c.magnitude().cmp(f.b.magnitude());

    if a_cmpabs_b == std::cmp::Ordering::Less || c_cmpabs_b == std::cmp::Ordering::Less {
        return false;
    }

    // a >= |b| and c >= |b|, so it might be reduced
    let a_cmp_c = f.a.cmp(&f.c);
    if a_cmp_c == std::cmp::Ordering::Greater {
        std::mem::swap(&mut f.a, &mut f.c);
        f.b = -f.b.clone();
    } else if a_cmp_c == std::cmp::Ordering::Equal && f.b.is_negative() {
        f.b = -f.b.clone();
    }
    true
}

/// Lehmer acceleration step: compute (u, v, w, x) 2x2 matrix.
fn calc_uvwx(mut a: i64, mut b: i64, mut c: i64) -> (i64, i64, i64, i64) {
    let mut u_ = 1i64;
    let mut v_ = 0i64;
    let mut w_ = 0i64;
    let mut x_ = 1i64;

    let mut u;
    let mut v;
    let mut w;
    let mut x;

    loop {
        u = u_;
        v = v_;
        w = w_;
        x = x_;

        if c == 0 {
            break;
        }

        let s = if b >= 0 {
            (b + c) / (c << 1)
        } else {
            -(-b + c) / (c << 1)
        };

        let a_ = a;
        let b_ = b;

        a = c;
        b = -b + (c.wrapping_mul(s) << 1);
        c = a_ - s * (b_ - c.wrapping_mul(s));

        u_ = v;
        v_ = -u + s.wrapping_mul(v);
        w_ = x;
        x_ = -w + s.wrapping_mul(x);

        let below_threshold = (v_.abs() | x_.abs()) <= THRESH;
        if !(below_threshold && a > c && c > 0) {
            if below_threshold {
                u = u_;
                v = v_;
                w = w_;
                x = x_;
            }
            break;
        }
    }

    (u, v, w, x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::form::Form;
    use num_bigint::BigInt;

    fn disc_check(f: &Form, d: &BigInt) -> bool {
        let disc = &f.b * &f.b - BigInt::from(4) * &f.a * &f.c;
        &disc == d
    }

    #[test]
    fn test_reduce_preserves_discriminant() {
        // D = -47
        let d = BigInt::from(-47i64);
        // Form (3, 1, 4) with discriminant 1 - 48 = -47
        let mut f = Form::new(BigInt::from(3), BigInt::from(1), BigInt::from(4));
        assert!(disc_check(&f, &d));
        reduce(&mut f);
        assert!(disc_check(&f, &d), "discriminant changed after reduction");
        assert!(
            f.is_reduced(),
            "form not reduced: a={}, b={}, c={}",
            f.a,
            f.b,
            f.c
        );
    }

    #[test]
    fn test_reduce_idempotent() {
        let d = BigInt::from(-47i64);
        let f = Form::new(BigInt::from(5), BigInt::from(3), BigInt::from(3));
        let _disc = &f.b * &f.b - BigInt::from(4) * &f.a * &f.c;
        // 9 - 60 = -51, not -47. Let's use a form that's actually valid
        // Form (5, 1, c) where c = (1+47)/20 = 48/20 is not integer...
        // Use (2, 1, 6): disc = 1 - 48 = -47
        let mut f = Form::new(BigInt::from(2), BigInt::from(1), BigInt::from(6));
        assert!(disc_check(&f, &d));
        reduce(&mut f);
        let f2 = f.clone();
        reduce(&mut f);
        assert_eq!(f, f2, "reduction should be idempotent");
    }
}
