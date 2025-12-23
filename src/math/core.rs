use std::fmt;
use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};

pub trait LinearSpace:
    Sized
    + Neg<Output = Self>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
    + fmt::Display
    + fmt::Debug
{
    type Data;

    fn new(data: Self::Data) -> Self;
    fn zero() -> Self;
    fn size(&self) -> usize;
    fn get_data(&self) -> Self::Data;
}

pub trait ScalarSpace: LinearSpace {
    fn abs(&self) -> Self;
    fn sqrt(&self) -> Self;
    fn sin(&self) -> Self;
    fn cos(&self) -> Self;
}

pub trait InnerProduct: LinearSpace {
    fn dot(&self, other: Self) -> f64;
}

pub trait VectorSpace<const N: usize>: InnerProduct + Index<usize> + IndexMut<usize> {
    type Valule;

    fn get(&self, index: usize) -> Self::Valule;

    fn magnitude(&self) -> f64;
    fn magnitude_square(&self) -> f64;
    fn normalize(&self) -> Self;
}

pub trait CrossProduct<const N: usize>: VectorSpace<N> {
    fn cross(&self, other: Self) -> Self;
}

#[derive(Clone, Copy, PartialEq)]
pub struct Vector<const N: usize> {
    data: [f64; N],
}

impl LinearSpace for f64 {
    type Data = f64;

    fn new(data: Self::Data) -> Self {
        data
    }
    fn zero() -> Self {
        0.0
    }
    fn get_data(&self) -> Self::Data {
        *self
    }
    fn size(&self) -> usize {
        1
    }
}

impl ScalarSpace for f64 {
    fn abs(&self) -> Self {
        <f64>::abs(*self)
    }
    fn sqrt(&self) -> Self {
        <f64>::sqrt(*self)
    }
    fn sin(&self) -> Self {
        <f64>::sin(*self)
    }
    fn cos(&self) -> Self {
        <f64>::cos(*self)
    }
}

impl<const N: usize> LinearSpace for Vector<N> {
    type Data = [f64; N];

    fn new(data: Self::Data) -> Self {
        Self { data: data }
    }
    fn zero() -> Self {
        Self::new([0.0; N])
    }
    fn size(&self) -> usize {
        N
    }
    fn get_data(&self) -> Self::Data {
        self.data
    }
}

impl<const N: usize> VectorSpace<N> for Vector<N> {
    type Valule = f64;

    fn get(&self, index: usize) -> Self::Valule {
        self.data[index]
    }

    fn magnitude_square(&self) -> f64 {
        self.data.iter().map(|e| e * e).sum()
    }
    fn magnitude(&self) -> f64 {
        self.magnitude_square().sqrt()
    }
    fn normalize(&self) -> Self {
        match self.magnitude_square() {
            0.0 => Self::zero(),
            _ => *self / self.magnitude(),
        }
    }
}

impl<const N: usize> Index<usize> for Vector<N> {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<const N: usize> IndexMut<usize> for Vector<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl<const N: usize> InnerProduct for Vector<N> {
    fn dot(&self, other: Self) -> f64 {
        self.data
            .iter()
            .zip(other.data.iter())
            .map(|(e1, e2)| e1 * e2)
            .sum()
    }
}

impl CrossProduct<3> for Vector<3> {
    fn cross(&self, other: Self) -> Self {
        Self {
            data: [
                self.data[1] * other.data[2] - self.data[2] * self.data[1],
                self.data[2] * other.data[0] - self.data[0] * self.data[2],
                self.data[0] * other.data[1] - self.data[1] * self.data[0],
            ],
        }
    }
}

impl<const N: usize> fmt::Display for Vector<N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Vector<{N}>{:?}", self.data)
    }
}

impl<const N: usize> fmt::Debug for Vector<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Vector<{N}>{:?}", self.data)
    }
}

impl<const N: usize> Neg for Vector<N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            data: self.data.map(|e| -e),
        }
    }
}

impl<const N: usize> Add for Vector<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            data: std::array::from_fn(|i| self.data[i] + rhs.data[i]),
        }
    }
}

impl<const N: usize> Sub for Vector<N> {
    type Output = Self;

    fn sub(self, rhs: Vector<N>) -> Self::Output {
        Self {
            data: std::array::from_fn(|i| self.data[i] - rhs.data[i]),
        }
    }
}

impl<const N: usize> Mul<f64> for Vector<N> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            data: std::array::from_fn(|i| self.data[i] * rhs),
        }
    }
}

impl<const N: usize> Div<f64> for Vector<N> {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            data: std::array::from_fn(|i| self.data[i] / rhs),
        }
    }
}

impl<const N: usize> Mul<Vector<N>> for f64 {
    type Output = Vector<N>;

    fn mul(self, rhs: Vector<N>) -> Self::Output {
        Self::Output {
            data: std::array::from_fn(|i| self * rhs.data[i]),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vector_op() {
        let v1 = Vector::new([1.0, 2.0, 3.0]);
        let v2 = Vector::new([3.0, 4.0, 5.0]);
        let f = 3.0;

        assert_eq!(v1 + v2, Vector::new([4.0, 6.0, 8.0]), "Vector + Vector");
        assert_eq!(v1 - v2, Vector::new([-2.0, -2.0, -2.0]), "Vector + Vector");

        assert_eq!(-v1, Vector::new([-1.0, -2.0, -3.0]), "Vector");

        assert_eq!(v1.dot(v2), 26.0, "Vector * Vector (Inner Product)");
        assert_eq!(v1 * f, Vector::new([3.0, 6.0, 9.0]), "Vector * f64");
        assert_eq!(f * v1, Vector::new([3.0, 6.0, 9.0]), "f64 * Vector");

        assert_eq!(
            v1 / f,
            Vector::new([1.0 / 3.0, 2.0 / 3.0, 1.0]),
            "Vector / f64"
        );
    }

    #[test]
    fn test_vector_magnitude() {
        let v = Vector::new([3.0, 4.0]);

        assert_eq!(v.magnitude_square(), 25.0, "Vector.magnitude_square()");
        assert_eq!(v.magnitude(), 5.0, "Vector.magnitude()");
    }

    #[test]
    fn test_vector_normalize() {
        let v1 = Vector::new([3.0, 4.0]);
        let v2 = Vector::zero();

        // Normalize
        let v1_normalized = v1.normalize();
        assert_eq!(v1_normalized.magnitude(), 1.0, "Vector.normalize()");
        assert_eq!(v1_normalized, Vector::new([0.6, 0.8]), "Vector.normalize()");

        let v2_normalized = v2.normalize();
        assert_eq!(v2_normalized.magnitude(), 0.0, "Vector::zero().normalize()");
        assert_eq!(
            v2_normalized,
            Vector::new([0.0, 0.0]),
            "Vector::zero().normalize()"
        );
    }

    #[test]
    fn test_inner_product() {
        let v1 = Vector::new([1.0, 0.0, 0.0]);
        let v2 = Vector::new([0.0, 1.0, 0.0]);
        let v3 = Vector::new([0.0, 0.0, 1.0]);
        let v4 = Vector::new([2.0, 2.0, 0.0]);

        assert_eq!(v1.dot(v2), 0.0, "Vector Inner Product");
        assert_eq!(v2.dot(v3), 0.0, "Vector Inner Product");
        assert_eq!(v3.dot(v1), 0.0, "Vector Inner Product");

        assert_eq!(v1.dot(v4), 2.0, "Vector Inner Product");

        assert_eq!(
            v3.dot(v3),
            v3.magnitude_square(),
            "Vector Self Inner Product"
        );
    }

    #[test]
    fn test_outer_product() {
        let e0 = Vector::new([1.0, 0.0, 0.0]);
        let e1 = Vector::new([0.0, 1.0, 0.0]);
        let e2 = Vector::new([0.0, 0.0, 1.0]);
        let v = Vector::new([2.0, 2.0, 0.0]);

        assert_eq!(
            e0.cross(e1),
            Vector::new([0.0, 0.0, 1.0]),
            "Vector Outer Product"
        );
        assert_eq!(
            e1.cross(e2),
            Vector::new([1.0, 0.0, 0.0]),
            "Vector Outer Product"
        );
        assert_eq!(
            e2.cross(e0),
            Vector::new([0.0, 1.0, 0.0]),
            "Vector Outer Product"
        );

        assert_eq!(
            e0.cross(v),
            Vector::new([0.0, 0.0, 2.0]),
            "Vector Outer Product"
        );

        assert_eq!(e2.cross(e2), Vector::zero(), "Vector Self Outer Product");
    }
}
