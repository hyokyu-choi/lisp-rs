use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

pub trait LinearSpace:
    Sized +
    Neg<Output = Self> +
    Add<Self, Output = Self> +
    Sub<Self, Output = Self> +
    Mul<Scalar, Output = Self> +
    Mul<f64, Output = Self> +
    Div<Scalar, Output = Self> +
    Div<f64, Output = Self> +
    fmt::Display +
{
    fn zero() -> Self;
    fn size(&self) -> usize;
}

pub trait VectorSpace :
    LinearSpace +
    Mul<Self, Output = Scalar> +
{
    fn magnitude(&self) -> Scalar;
    fn magnitude_square(&self) -> Scalar;
    fn normalize(&self) -> Self;
}

#[derive(Clone, Copy)]
pub struct Scalar(pub f64);
#[derive(Clone, Copy)]
pub struct Vector2(pub Scalar, pub Scalar);

impl Scalar {
    pub fn new(num: f64) -> Self {
        Self(num)
    }
    fn get_value(&self) -> f64 {
        self.0
    }
    fn abs(&self) -> Self {
        Self(self.0.abs())
    }
    fn sqrt(&self) -> Self {
        Self(self.0.sqrt())
    }

}

impl Vector2 {
    pub fn new(e0: f64, e1: f64) -> Self {
        Self(Scalar(e0), Scalar(e1))
    }
    fn get_value(&self) -> Vec<f64> {
        vec![self.0.get_value(), self.1.get_value()]
    }
}

impl LinearSpace for Scalar {
    fn zero() -> Self {
        Self(0.0)
    }
    fn size(&self) -> usize {
        1
    }
}

impl LinearSpace for Vector2 {
    fn zero() -> Self {
        Self::new(0.0, 0.0)
    }
    fn size(&self) -> usize{
        2
    }
}

impl VectorSpace for Vector2 {
    fn magnitude_square(&self) -> Scalar {
        self.0 * self.0 + self.1 * self.1
    }
    fn magnitude(&self) -> Scalar {
        self.magnitude_square().sqrt()
    }
    fn normalize(&self) -> Self {
        *self/self.magnitude()
    }
}

impl fmt::Display for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Scalar({:.6})", self.0)
    }
}

impl fmt::Display for Vector2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Vector2({:.6}, {:.6})", self.0.get_value(), self.1.get_value())
    }
}

impl Neg for Scalar {
    type Output = Self;

    fn neg(self) -> Self {
        Self(-self.0)
    }
}

impl Neg for Vector2 {
    type Output = Self;

    fn neg(self) -> Self {
        Self(-self.0, -self.1)
    }
}

impl Add for Scalar {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl Add<f64> for Scalar {
    type Output = Self;

    fn add(self, rhs: f64) -> Self {
        Self(self.0 + rhs)
    }
}

impl Add<Scalar> for f64 {
    type Output = Scalar;

    fn add(self, rhs: Scalar) -> Scalar {
        Scalar(self + rhs.0)
    }
}

impl Sub for Scalar {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl Sub<f64> for Scalar {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self {
        Self(self.0 - rhs)
    }
}

impl Sub<Scalar> for f64 {
    type Output = Scalar;

    fn sub(self, rhs: Scalar) -> Scalar {
        Scalar(self - rhs.0)
    }
}

impl Mul for Scalar {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self(self.0 * rhs.0)
    }
}

impl Mul<f64> for Scalar {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self {
        Self(self.0 * rhs)
    }
}

impl Mul<Scalar> for f64 {
    type Output = Scalar;

    fn mul(self, rhs: Scalar) -> Scalar {
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

    fn div(self, rhs: f64) -> Self {
        Self(self.0 / rhs)
    }
}

impl Div<Scalar> for f64 {
    type Output = Scalar;

    fn div(self, rhs: Scalar) -> Scalar {
        Scalar(self / rhs.0)
    }
}

impl Add for Vector2 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl Sub for Vector2 {
    type Output = Self;

    fn sub(self, rhs: Vector2) -> Vector2 {
        Self(self.0 - rhs.0, self.1 - rhs.1)
    }
}

impl Mul<Scalar> for Vector2 {
    type Output = Self;

    fn mul(self, num: Scalar) -> Self {
        Self(num * self.0, num * self.1)
    }
}

impl Mul<f64> for Vector2 {
    type Output = Self;

    fn mul(self, num: f64) -> Self {
        Self(num * self.0, num * self.1)
    }
}

impl Mul<Vector2> for Scalar {
    type Output = Vector2;

    fn mul(self, rhs: Vector2) -> Vector2 {
        rhs * self
    }
}

impl Mul<Vector2> for f64 {
    type Output = Vector2;

    fn mul(self, rhs: Vector2) -> Vector2 {
        Vector2(self * rhs.0, self * rhs.1)
    }
}

impl Mul for Vector2 {
    type Output = Scalar;

    fn mul(self, rhs: Self) -> Scalar {
        self.0 * rhs.0 + self.1 * rhs.1
    }
}

impl Div<Scalar> for Vector2 {
    type Output = Self;

    fn div(self, rhs: Scalar) -> Self {
        Self(self.0 / rhs, self.1 / rhs)
    }
}

impl Div<f64> for Vector2 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        Self(self.0 / rhs, self.1 / rhs)
    }
}
