use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::math::core::{LinearSpace, ScalarSpace, Vector, VectorSpace};

pub trait MatrixSpace<S: ScalarSpace, const N: usize, const M: usize>:
    LinearSpace + MatMul<Vector<S, M>, Output = Vector<S, N>> + Mul<Vector<S, M>, Output = Vector<S, N>>
{
    type Transpose: MatrixSpace<S, M, N>;

    fn rows(&self) -> usize {
        N
    }
    fn cols(&self) -> usize {
        M
    }

    fn get(&self, row: usize, col: usize) -> S;
    fn get_rows(&self, row: usize) -> Vector<S, M>;
    fn get_cols(&self, col: usize) -> Vector<S, N>;

    fn transpose(&self) -> Self::Transpose;
    fn adjoint(&self) -> Self::Transpose;
}

pub trait SquareMatrixSpace<S: ScalarSpace, const N: usize>:
    MatrixSpace<S, N, N> + Mul<Self, Output = Self> + MatMul<Self, Output = Self>
{
    fn identity() -> Self;
    fn is_invertible(&self) -> bool;
    fn invert(&self) -> Option<Self>;
    fn trace(&self) -> S;
    fn determinant(&self) -> S;
}

pub trait MatMul<RHS = Self> {
    type Output;

    fn matmul(&self, rhs: RHS) -> Self::Output;
}

#[derive(Clone, Copy, PartialEq)]
pub struct Matrix<S: ScalarSpace, const N: usize, const M: usize> {
    data: [[S; M]; N],
}
#[derive(Clone, Copy, PartialEq)]
pub struct SquareMatrix<S: ScalarSpace, const N: usize> {
    data: [[S; N]; N],
}

impl<S: ScalarSpace, const N: usize, const M: usize> LinearSpace for Matrix<S, N, M> {
    type Data = [[S; M]; N];

    fn new(data: Self::Data) -> Self {
        Self { data: data }
    }
    fn zero() -> Self {
        Self {
            data: [[S::zero(); M]; N],
        }
    }
    fn size(&self) -> usize {
        N * M
    }
    fn get_data(&self) -> Self::Data {
        self.data
    }
}

impl<S: ScalarSpace, const N: usize> LinearSpace for SquareMatrix<S, N> {
    type Data = [[S; N]; N];

    fn new(data: Self::Data) -> Self {
        Self {
            data: data,
        }
    }
    fn zero() -> Self {
        Self {
            data: [[S::zero(); N]; N],
        }
    }
    fn size(&self) -> usize {
        N * N
    }
    fn get_data(&self) -> Self::Data {
        self.data
    }
}

impl<S: ScalarSpace, const N: usize, const M: usize, const L: usize> MatMul<Matrix<S, M, L>>
    for Matrix<S, N, M>
{
    type Output = Matrix<S, N, L>;

    fn matmul(&self, rhs: Matrix<S, M, L>) -> Self::Output {
        Self::Output {
            data: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    (0..M)
                        .map(|k| self.data[i][k] * rhs.data[k][j])
                        .fold(S::zero(), |acc, val| acc + val)
                })
            }),
        }
    }
}
impl<S: ScalarSpace, const N: usize, const M: usize, const L: usize> Mul<Matrix<S, M, L>> for Matrix<S, N, M> {
    type Output = Matrix<S, N, L>;

    fn mul(self, rhs: Matrix<S, M, L>) -> Self::Output {
        self.matmul(rhs)
    }
}

impl<S: ScalarSpace, const N: usize, const M: usize> MatMul<Vector<S, M>> for Matrix<S, N, M> {
    type Output = Vector<S, N>;

    fn matmul(&self, rhs: Vector<S, M>) -> Self::Output {
        Self::Output::new(std::array::from_fn(|i| {
            self.data[i]
                .iter()
                .zip(rhs.as_array().iter())
                .map(|(a, b)| *a * *b)
                .fold(S::zero(), |acc, val| acc + val)
        }))
    }
}

impl<S: ScalarSpace, const N: usize, const M: usize> Mul<Vector<S, M>> for Matrix<S, N, M> {
    type Output = Vector<S, N>;

    fn mul(self, rhs: Vector<S, M>) -> Self::Output {
        self.matmul(rhs)
    }
}

impl<S: ScalarSpace, const N: usize> MatMul for SquareMatrix<S, N> {
    type Output = Self;

    fn matmul(&self, rhs: Self) -> Self::Output {
        Self::Output {
            data: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    (0..N)
                        .map(|k| self.data[i][k] * rhs.data[k][j])
                        .fold(S::zero(), |acc, val| acc + val)
                })
            }),
        }
    }
}
impl<S: ScalarSpace, const N: usize> Mul for SquareMatrix<S, N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.matmul(rhs)
    }
}

impl<S: ScalarSpace, const N: usize> MatMul<Vector<S, N>> for SquareMatrix<S, N> {
    type Output = Vector<S, N>;

    fn matmul(&self, rhs: Vector<S, N>) -> Self::Output {
        Self::Output::new(std::array::from_fn(|i| {
            self.data[i]
                .iter()
                .zip(rhs.as_array().iter())
                .map(|(a, b)| *a * *b)
                .fold(S::zero(), |acc, val| acc + val)
        }))
    }
}

impl<S: ScalarSpace, const N: usize> Mul<Vector<S, N>> for SquareMatrix<S, N> {
    type Output = Vector<S, N>;

    fn mul(self, rhs: Vector<S, N>) -> Self::Output {
        self.matmul(rhs)
    }
}

impl<S: ScalarSpace, const N: usize, const M: usize> MatrixSpace<S, N, M> for Matrix<S, N, M> {
    type Transpose = Matrix<S, M, N>;

    fn get(&self, row: usize, col: usize) -> S {
        self.data[row][col]
    }
    fn get_rows(&self, row: usize) -> Vector<S, M> {
        Vector::new(self.data[row])
    }
    fn get_cols(&self, col: usize) -> Vector<S, N> {
        Vector::new(std::array::from_fn(|i| self.data[i][col]))
    }
    fn transpose(&self) -> Self::Transpose {
        Matrix {
            data: std::array::from_fn(|i| std::array::from_fn(|j| self.data[j][i])),
        }
    }
    fn adjoint(&self) -> Self::Transpose {
        Matrix {
            data: std::array::from_fn(|i| std::array::from_fn(|j| self.data[j][i].conj())),
        }
    }
}

impl<S: ScalarSpace, const N: usize> MatrixSpace<S, N, N> for SquareMatrix<S, N> {
    type Transpose = Self;

    fn get(&self, row: usize, col: usize) -> S {
        self.data[row][col]
    }
    fn get_rows(&self, row: usize) -> Vector<S, N> {
        Vector::new(self.data[row])
    }
    fn get_cols(&self, col: usize) -> Vector<S, N> {
        Vector::new(std::array::from_fn(|i| self.data[i][col]))
    }
    fn transpose(&self) -> Self::Transpose {
        Self {
            data: std::array::from_fn(|i| std::array::from_fn(|j| self.data[j][i])),
        }
    }
    fn adjoint(&self) -> Self::Transpose {
        Self {
            data: std::array::from_fn(|i| std::array::from_fn(|j| self.data[j][i].conj())),
        }
    }
}

impl<S: ScalarSpace, const N: usize> SquareMatrixSpace<S, N> for SquareMatrix<S, N> {
    fn identity() -> Self {
        Self {
            data: std::array::from_fn(|i| {
                std::array::from_fn(|j| match j {
                    k if i == k => S::one(),
                    _ => S::zero(),
                })
            }),
        }
    }
    /// TODO: Implement
    fn is_invertible(&self) -> bool {
        false
    }

    /// TODO: Implement
    fn invert(&self) -> Option<Self> {
        Option::None
    }
    fn trace(&self) -> S {
        (0..N)
            .map(|i| self.data[i][i])
            .fold(S::zero(), |acc, var| acc + var)
    }
    /// TODO: Implement
    fn determinant(&self) -> S {
        S::zero()
    }
}

impl<S: ScalarSpace, const N: usize, const M: usize> fmt::Display for Matrix<S, N, M> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Matrix<{N}, {M}>{:?}", self.data)
    }
}

impl<S: ScalarSpace, const N: usize, const M: usize> fmt::Debug for Matrix<S, N, M> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Matrix<{N}, {M}>{:?}", self.data)
    }
}

impl<S: ScalarSpace, const N: usize> fmt::Display for SquareMatrix<S, N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "SquareMatrix<{N}>{:?}", self.data)
    }
}

impl<S: ScalarSpace, const N: usize> fmt::Debug for SquareMatrix<S, N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SquareMatrix<{N}>{:?}", self.data)
    }
}

impl<S: ScalarSpace, const N: usize, const M: usize> Neg for Matrix<S, N, M> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::Output {
            data: std::array::from_fn(|i| std::array::from_fn(|j| -self.data[i][j])),
        }
    }
}

impl<S: ScalarSpace, const N: usize, const M: usize> Add for Matrix<S, N, M> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::Output {
            data: std::array::from_fn(|i| {
                std::array::from_fn(|j| self.data[i][j] + rhs.data[i][j])
            }),
        }
    }
}

impl<S: ScalarSpace, const N: usize, const M: usize> Sub for Matrix<S, N, M> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            data: std::array::from_fn(|i| {
                std::array::from_fn(|j| self.data[i][j] - rhs.data[i][j])
            }),
        }
    }
}

impl<S: ScalarSpace, const N: usize, const M: usize> Mul<f64> for Matrix<S, N, M> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::Output {
            data: std::array::from_fn(|i| std::array::from_fn(|j| self.data[i][j] * rhs)),
        }
    }
}

impl<S: ScalarSpace, const N: usize, const M: usize> Div<f64> for Matrix<S, N, M> {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self::Output {
            data: std::array::from_fn(|i| std::array::from_fn(|j| self.data[i][j] / rhs)),
        }
    }
}

impl<S: ScalarSpace, const N: usize> Neg for SquareMatrix<S, N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::Output {
            data: std::array::from_fn(|i| std::array::from_fn(|j| -self.data[i][j])),
        }
    }
}

impl<S: ScalarSpace, const N: usize> Add for SquareMatrix<S, N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::Output {
            data: std::array::from_fn(|i| {
                std::array::from_fn(|j| self.data[i][j] + rhs.data[i][j])
            }),
        }
    }
}

impl<S: ScalarSpace, const N: usize> Sub for SquareMatrix<S, N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            data: std::array::from_fn(|i| {
                std::array::from_fn(|j| self.data[i][j] - rhs.data[i][j])
            }),
        }
    }
}

impl<S: ScalarSpace, const N: usize> Mul<f64> for SquareMatrix<S, N> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::Output {
            data: std::array::from_fn(|i| std::array::from_fn(|j| self.data[i][j] * rhs)),
        }
    }
}

impl<S: ScalarSpace, const N: usize> Div<f64> for SquareMatrix<S, N> {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self::Output {
            data: std::array::from_fn(|i| std::array::from_fn(|j| self.data[i][j] / rhs)),
        }
    }
}
