use std::fmt;
use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};

pub trait LinearSpace:
    Sized
    + Neg<Output = Self>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + Mul<Scalar, Output = Self>
    + Mul<f64, Output = Self>
    + Div<Scalar, Output = Self>
    + Div<f64, Output = Self>
    + fmt::Display
    + fmt::Debug
{
    type Value;

    fn from_self(u: Self) -> Self;
    fn get_value(&self) -> Self::Value;
    fn zero() -> Self;
    fn size(&self) -> usize;
    fn scale(&self, other: Scalar) -> Self;
    fn scale_f64(&self, other: f64) -> Self;
}

pub trait ScalarSpace:
    Sized + LinearSpace + Mul<Vector2, Output = Vector2> + PartialOrd<Self>
{
    fn abs(&self) -> Self;
    fn sqrt(&self) -> Self;
    fn sin(&self) -> Self;
    fn cos(&self) -> Self;
}

pub trait VectorSpace: Sized + LinearSpace + Index<usize> + IndexMut<usize> {
    fn get(&self, i: usize) -> Scalar;
    fn magnitude(&self) -> Scalar;
    fn magnitude_square(&self) -> Scalar;
    fn normalize(&self) -> Self;
}

pub trait InnerProduct: VectorSpace {
    fn inner_product(&self, other: Self) -> Scalar;
}

pub trait OuterProduct: VectorSpace {
    fn outer_product(&self, other: Self) -> Self;
}

#[derive(Clone, Copy, PartialEq, PartialOrd)]
pub struct Scalar(pub f64);
#[derive(Clone, Copy, PartialEq)]
pub struct Vector2 {
    pub data: [Scalar; 2],
}

impl Scalar {
    pub fn new(num: f64) -> Self {
        Self(num)
    }
}

impl Vector2 {
    pub fn new(e0: f64, e1: f64) -> Self {
        Self {
            data: [Scalar::new(e0), Scalar::new(e1)],
        }
    }
}

impl LinearSpace for Scalar {
    type Value = f64;

    fn from_self(u: Self) -> Self {
        u
    }
    fn get_value(&self) -> Self::Value {
        self.0
    }
    fn zero() -> Self {
        Self(0.0)
    }
    fn size(&self) -> usize {
        1
    }
    fn scale(&self, other: Scalar) -> Self {
        Self(self.0 * other.0)
    }
    fn scale_f64(&self, other: f64) -> Self {
        Self(self.0 * other)
    }
}

impl ScalarSpace for Scalar {
    fn abs(&self) -> Self {
        Self(self.0.abs())
    }
    fn sqrt(&self) -> Self {
        Self(self.0.sqrt())
    }
    fn sin(&self) -> Self {
        Self(self.0.sin())
    }
    fn cos(&self) -> Self {
        Self(self.0.cos())
    }
}

impl LinearSpace for Vector2 {
    type Value = [f64; 2];

    fn from_self(u: Self) -> Self {
        u
    }
    fn get_value(&self) -> Self::Value {
        [self.data[0].get_value(), self.data[1].get_value()]
    }
    fn zero() -> Self {
        Self::new(0.0, 0.0)
    }
    fn size(&self) -> usize {
        2
    }
    fn scale(&self, other: Scalar) -> Self {
        (*self) * other
    }
    fn scale_f64(&self, other: f64) -> Self {
        (*self) * other
    }
}

impl VectorSpace for Vector2 {
    fn get(&self, i: usize) -> Scalar {
        self.data[i]
    }
    fn magnitude_square(&self) -> Scalar {
        (self.data[0] * self.data[0]) + (self.data[1] * self.data[1])
    }
    fn magnitude(&self) -> Scalar {
        self.magnitude_square().sqrt()
    }
    fn normalize(&self) -> Self {
        match self.magnitude() {
            Scalar(0.0) => Self::zero(),
            _ => *self / self.magnitude(),
        }
    }
}

impl Index<usize> for Vector2 {
    type Output = Scalar;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl IndexMut<usize> for Vector2 {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl InnerProduct for Vector2 {
    fn inner_product(&self, other: Self) -> Scalar {
        (self.data[0] * other.data[0]) + (self.data[1] * other.data[1])
    }
}

impl fmt::Display for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Scalar({:.6})", self.0)
    }
}

impl fmt::Debug for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Scalar({:.6})", self.0)
    }
}

impl fmt::Display for Vector2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Vector2({:.6}, {:.6})", self.data[0], self.data[1])
    }
}

impl fmt::Debug for Vector2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Vector2({:.6}, {:.6})", self.data[0], self.data[1])
    }
}

impl Neg for Scalar {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl Neg for Vector2 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            data: [-self.data[0], -self.data[1]],
        }
    }
}

impl Add for Scalar {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl Add<f64> for Scalar {
    type Output = Self;

    fn add(self, rhs: f64) -> Self::Output {
        Self(self.0 + rhs)
    }
}

impl Add<Scalar> for f64 {
    type Output = Scalar;

    fn add(self, rhs: Scalar) -> Self::Output {
        Scalar(self + rhs.0)
    }
}

impl Sub for Scalar {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl Sub<f64> for Scalar {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self::Output {
        Self(self.0 - rhs)
    }
}

impl Sub<Scalar> for f64 {
    type Output = Scalar;

    fn sub(self, rhs: Scalar) -> Self::Output {
        Scalar(self - rhs.0)
    }
}

impl Mul<f64> for Scalar {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self(self.0 * rhs)
    }
}

impl Mul<Scalar> for f64 {
    type Output = Scalar;

    fn mul(self, rhs: Scalar) -> Self::Output {
        Scalar(self * rhs.0)
    }
}

impl Div for Scalar {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        Self(self.0 / rhs.0)
    }
}

impl Div<f64> for Scalar {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self(self.0 / rhs)
    }
}

impl Div<Scalar> for f64 {
    type Output = Scalar;

    fn div(self, rhs: Scalar) -> Self::Output {
        Scalar(self / rhs.0)
    }
}

impl Add for Vector2 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            data: [self.data[0] + rhs.data[0], self.data[1] + rhs.data[1]],
        }
    }
}

impl Sub for Vector2 {
    type Output = Self;

    fn sub(self, rhs: Vector2) -> Self::Output {
        Self {
            data: [self.data[0] - rhs.data[0], self.data[1] - rhs.data[1]],
        }
    }
}

impl Mul for Vector2 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            data: [rhs.data[0] * self.data[0], rhs.data[1] * self.data[1]],
        }
    }
}

impl Mul<Scalar> for Vector2 {
    type Output = Self;

    fn mul(self, rhs: Scalar) -> Self::Output {
        Self {
            data: [rhs * self.data[0], rhs * self.data[1]],
        }
    }
}

impl Mul<f64> for Vector2 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            data: [rhs * self.data[0], rhs * self.data[1]],
        }
    }
}

impl<L: LinearSpace> Mul<L> for Scalar {
    type Output = L;

    fn mul(self, rhs: L) -> Self::Output {
        rhs.scale(self)
    }
}

impl Mul<Vector2> for f64 {
    type Output = Vector2;

    fn mul(self, rhs: Vector2) -> Self::Output {
        Vector2 {
            data: [self * rhs.data[0], self * rhs.data[1]],
        }
    }
}

impl Div<Scalar> for Vector2 {
    type Output = Self;

    fn div(self, rhs: Scalar) -> Self::Output {
        Self {
            data: [self.data[0] / rhs, self.data[1] / rhs],
        }
    }
}

impl Div<f64> for Vector2 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        Self {
            data: [self.data[0] / rhs, self.data[1] / rhs],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    fn assert_scalar_eq(a: Scalar, b: Scalar) {
        let diff = (a - b).abs();
        assert!(
            diff < Scalar::new(EPS),
            "Assertion failed: {:?} != {:?} (diff: {})",
            a,
            b,
            diff
        );
    }

    fn assert_vector_eq(a: Vector2, b: Vector2) {
        let diff_e0 = (a - b)[0].abs();
        let diff_e1 = (a - b)[1].abs();
        assert!(
            diff_e0 < Scalar::new(EPS) && diff_e1 < Scalar::new(EPS),
            "Assertion failed: {:?} != {:?}",
            a,
            b
        )
    }

    #[test]
    fn test_scalar_op() {
        let a = Scalar::new(2.0);
        let b = Scalar::new(3.0);

        // Scalar <op> Scalar
        assert_scalar_eq(a + b, Scalar::new(5.0));
        assert_scalar_eq(a - b, Scalar::new(-1.0));
        assert_scalar_eq(a * b, Scalar::new(6.0));
        assert_scalar_eq(a / b, Scalar::new(2.0 / 3.0));

        // Neg
        assert_scalar_eq(-a, Scalar::new(-2.0));
    }

    #[test]
    fn test_scalar_f64_op() {
        let s = Scalar::new(2.0);
        let f = 3.0;

        // Scalar <op> f64
        assert_scalar_eq(s + f, Scalar::new(5.0));
        assert_scalar_eq(s - f, Scalar::new(-1.0));
        assert_scalar_eq(s * f, Scalar::new(6.0));
        assert_scalar_eq(s / f, Scalar::new(2.0 / 3.0));

        // f64 <op> Scalar
        assert_scalar_eq(f + s, Scalar::new(5.0));
        assert_scalar_eq(f - s, Scalar::new(1.0));
        assert_scalar_eq(f * s, Scalar::new(6.0));
        assert_scalar_eq(f / s, Scalar::new(1.5));
    }

    #[test]
    fn test_vector2_op() {
        let v1 = Vector2::new(1.0, 2.0);
        let v2 = Vector2::new(3.0, 4.0);

        // Vector2 <op> Vector2
        assert_vector_eq(v1 + v2, Vector2::new(4.0, 6.0));
        assert_vector_eq(v1 - v2, Vector2::new(-2.0, -2.0));
        assert_vector_eq(v1 * v2, Vector2::new(3.0, 8.0));

        // Neg
        assert_vector_eq(-v1, Vector2::new(-1.0, -2.0));
    }

    #[test]
    fn test_vector_scalar_op() {
        let v = Vector2::new(1.0, 2.0);
        let s = Scalar::new(2.0);
        let f = 3.0;

        // Mul
        assert_vector_eq(v * s, Vector2::new(2.0, 4.0));
        assert_vector_eq(s * v, Vector2::new(2.0, 4.0));
        assert_vector_eq(v * f, Vector2::new(3.0, 6.0));
        assert_vector_eq(f * v, Vector2::new(3.0, 6.0));

        // Div
        assert_vector_eq(v / s, Vector2::new(1.0 / 2.0, 1.0));
        assert_vector_eq(v / f, Vector2::new(1.0 / 3.0, 2.0 / 3.0));
    }

    #[test]
    fn test_vector_magnitude() {
        let v = Vector2::new(3.0, 4.0);

        assert_scalar_eq(v.magnitude(), Scalar::new(5.0));
        assert_scalar_eq(v.magnitude_square(), Scalar::new(25.0));
    }

    #[test]
    fn test_vector_normalize() {
        let v1 = Vector2::new(3.0, 4.0);
        let v2 = Vector2::zero();

        // Normalize
        let v1_normalized = v1.normalize();
        assert_scalar_eq(v1_normalized.magnitude(), Scalar::new(1.0));
        assert_vector_eq(v1_normalized, Vector2::new(0.6, 0.8));

        let v2_normalized = v2.normalize();
        assert_scalar_eq(v2_normalized.magnitude(), Scalar::new(0.0));
        assert_vector_eq(v2_normalized, Vector2::new(0.0, 0.0));
    }

    #[test]
    fn test_inner_product() {
        let v1 = Vector2::new(1.0, 0.0);
        let v2 = Vector2::new(0.0, 1.0);
        let v3 = Vector2::new(2.0, 2.0);

        // Orthogonal
        assert_scalar_eq(v1.inner_product(v2), Scalar::new(0.0));

        // Parallel
        assert_scalar_eq(v1.inner_product(v3), Scalar::new(2.0));

        // Self inner product
        assert_scalar_eq(v3.inner_product(v3), v3.magnitude_square());
    }
}
