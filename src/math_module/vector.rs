use std::fmt;
use std::ops;

#[derive(Clone, Copy)]
pub struct Vector2 {
    x: f64,
    y: f64,
}

impl Vector2 {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x: x, y: y }
    }
    pub fn x(&self) -> f64 {
        self.x
    }
    pub fn y(&self) -> f64 {
        self.y
    }
    pub fn size(&self) -> usize {
        2
    }
    pub fn inner_product(&self, other: Vector2) -> f64 {
        self.x * other.x + self.y * other.y
    }
    pub fn magnitude_square(&self) -> f64 {
        self.x * self.x + self.y * self.y
    }
    pub fn magnitude(&self) -> f64 {
        self.magnitude_square().sqrt()
    }
    pub fn normalize(&self) -> Vector2 {
        *self / self.magnitude()
    }
}

impl fmt::Display for Vector2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Vector2({:.6}, {:.6})", self.x, self.y)
    }
}

impl ops::Add for Vector2 {
    type Output = Self;

    fn add(self, other: Vector2) -> Vector2 {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl ops::Sub for Vector2 {
    type Output = Self;

    fn sub(self, other: Vector2) -> Vector2 {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl ops::Mul<f64> for Vector2 {
    type Output = Self;

    fn mul(self, num: f64) -> Self {
        Self {
            x: num * self.x,
            y: num * self.y,
        }
    }
}

impl ops::Mul<Vector2> for f64 {
    type Output = Vector2;

    fn mul(self, other: Vector2) -> Vector2 {
        other * self
    }
}

impl ops::Mul for Vector2 {
    type Output = f64;

    fn mul(self, other: Vector2) -> f64 {
        self.x * other.x + self.y * other.y
    }
}

impl ops::Div<f64> for Vector2 {
    type Output = Self;

    fn div(self, num: f64) -> Vector2 {
        Self {
            x: self.x / num,
            y: self.y / num,
        }
    }
}
