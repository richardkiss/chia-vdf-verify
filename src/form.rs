//! Quadratic form (a, b, c) with discriminant D, where b^2 - 4ac = D.
//! Always in reduced form: |b| <= a <= c, with b >= 0 when a == c or |b| == a.

use crate::integer::{divexact, num_bits};
use crate::reducer::reduce;
use num_bigint::BigInt;
use num_traits::{One, Zero};

/// A binary quadratic form with coefficients (a, b, c) and discriminant D.
#[derive(Clone, Debug)]
pub struct Form {
    pub a: BigInt,
    pub b: BigInt,
    pub c: BigInt,
}

impl Form {
    pub fn new(a: BigInt, b: BigInt, c: BigInt) -> Self {
        Form { a, b, c }
    }

    /// Construct form from (a, b) and discriminant D, computing c = (b^2 - D) / (4a).
    pub fn from_abd(a: BigInt, b: BigInt, d: &BigInt) -> Self {
        // c = (b^2 - D) / (4a)
        let b2 = &b * &b;
        let num = b2 - d;
        let denom = BigInt::from(4) * &a;
        let c = divexact(&num, &denom);
        Form { a, b, c }
    }

    /// Identity form: (1, 1, (1-D)/4).
    /// The C code uses a=1, b=1 as the identity.
    pub fn identity(d: &BigInt) -> Self {
        let a = BigInt::one();
        let b = BigInt::one();
        // c = (b^2 - D) / (4a) = (1 - D) / 4
        let num = BigInt::one() - d;
        let c = divexact(&num, &BigInt::from(4));
        Form { a, b, c }
    }

    /// Generator form: (2, 1, (1-D)/8) — only valid when D ≡ 1 mod 8.
    pub fn generator(d: &BigInt) -> Self {
        let a = BigInt::from(2);
        let b = BigInt::one();
        let num = BigInt::one() - d;
        let c = divexact(&num, &BigInt::from(8));
        Form { a, b, c }
    }

    /// Check if this is the identity form (a=1, b=1).
    pub fn is_identity(&self) -> bool {
        self.a == BigInt::one() && self.b == BigInt::one()
    }

    /// Check if this is the generator form (a=2, b=1).
    pub fn is_generator(&self) -> bool {
        self.a == BigInt::from(2) && self.b == BigInt::one()
    }

    /// Check if this form is reduced: |b| <= a <= c, with b >= 0 when a == c or |b| == a.
    pub fn is_reduced(&self) -> bool {
        use num_traits::Signed;
        let abs_b = self.b.abs();
        if abs_b > self.a {
            return false;
        }
        if self.a > self.c {
            return false;
        }
        if self.a == self.c && self.b < BigInt::zero() {
            return false;
        }
        if abs_b == self.a && self.b < BigInt::zero() {
            return false;
        }
        true
    }

    /// Reduce this form in place using the Pulmark reducer.
    pub fn reduce(&mut self) {
        reduce(self);
    }

    /// The half-max size parameter L = floor((-D)^(1/4)).
    /// Used as a threshold in nucomp.
    pub fn compute_l(d: &BigInt) -> BigInt {
        use num_traits::Signed;
        let neg_d = d.abs();
        crate::integer::nth_root(&neg_d, 4)
    }

    /// Discriminant size in bits.
    pub fn d_bits(d: &BigInt) -> usize {
        num_bits(d)
    }
}

impl PartialEq for Form {
    fn eq(&self, other: &Self) -> bool {
        self.a == other.a && self.b == other.b && self.c == other.c
    }
}

impl Eq for Form {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity_form() {
        // Discriminant -47 ≡ 1 mod 8? -47 mod 8 = 1 (since -47 = -6*8 + 1)
        let d = BigInt::from(-47i64);
        let f = Form::identity(&d);
        // Check discriminant: b^2 - 4ac = D
        let disc = &f.b * &f.b - BigInt::from(4) * &f.a * &f.c;
        assert_eq!(disc, d, "identity form discriminant check");
    }
}
