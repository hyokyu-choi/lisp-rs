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
        *self / self.magnitude()
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
