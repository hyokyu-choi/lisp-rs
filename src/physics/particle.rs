use crate::math::core::{LinearSpace, Vector};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug)]
pub struct Particle {
    mass: f64,
    position: Vector<f64, 3>,
    velocity: Vector<f64, 3>,
}

impl Particle {
    pub fn new(m: f64, x: Vector<f64, 3>, v: Vector<f64, 3>) -> Self {
        Self {
            mass: m,
            position: x,
            velocity: v,
        }
    }
    pub fn mass(&self) -> f64 {
        self.mass
    }
    pub fn position(&self) -> Vector<f64, 3> {
        self.position
    }
    pub fn velocity(&self) -> Vector<f64, 3> {
        self.velocity
    }
    pub fn momentum(&self) -> Vector<f64, 3> {
        self.velocity * self.mass
    }
    pub fn set_position(&mut self, x: Vector<f64, 3>) {
        self.position = x;
    }
    pub fn set_velocity(&mut self, v: Vector<f64, 3>) {
        self.velocity = v;
    }
}

impl fmt::Display for Particle {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Particle(mass: {}, x: {:?}, v: {:?})",
            self.mass, self.position, self.velocity
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NBodyState<const N: usize> {
    pub components: Vector<Vector<f64, 3>, N>,
}

impl<const N: usize> NBodyState<N> {
    pub fn new(components: Vector<Vector<f64, 3>, N>) -> Self {
        Self { components }
    }

    pub fn from_particles_pos(ps: &[Particle; N]) -> Self {
        let mut components = Vector::zero();
        for i in 0..N {
            components[i] = ps[i].position;
        }
        Self { components }
    }

    pub fn from_particles_vel(ps: &[Particle; N]) -> Self {
        let mut components = Vector::zero();
        for i in 0..N {
            components[i] = ps[i].velocity;
        }
        Self { components }
    }
}

impl<const N: usize> LinearSpace for NBodyState<N> {
    type Data = [[f64; 3]; N];

    fn new(data: Self::Data) -> Self {
        let mut components = Vector::zero();
        for i in 0..N {
            components[i] = Vector::new(data[i]);
        }
        Self { components }
    }

    fn zero() -> Self {
        Self {
            components: Vector::zero(),
        }
    }

    fn size(&self) -> usize {
        N * 3
    }

    fn get_data(&self) -> Self::Data {
        std::array::from_fn(|i| self.components[i].get_data())
    }
}

impl<const N: usize> Neg for NBodyState<N> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            components: -self.components,
        }
    }
}

impl<const N: usize> Add for NBodyState<N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            components: self.components + rhs.components,
        }
    }
}

impl<const N: usize> Sub for NBodyState<N> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            components: self.components - rhs.components,
        }
    }
}

impl<const N: usize> Mul<f64> for NBodyState<N> {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            components: self.components * rhs,
        }
    }
}

impl<const N: usize> Div<f64> for NBodyState<N> {
    type Output = Self;
    fn div(self, rhs: f64) -> Self::Output {
        Self {
            components: self.components / rhs,
        }
    }
}

impl<const N: usize> fmt::Display for NBodyState<N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "NBodyState({:?})", self.components)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_particle_3d() {
        let p = Particle::new(
            1.0,
            Vector::new([1.0, 2.0, 3.0]),
            Vector::new([0.1, 0.2, 0.3]),
        );
        assert_eq!(p.position()[0], 1.0);
        assert_eq!(p.position()[1], 2.0);
        assert_eq!(p.position()[2], 3.0);
        assert_eq!(p.velocity()[0], 0.1);
        assert_eq!(p.velocity()[1], 0.2);
        assert_eq!(p.velocity()[2], 0.3);
        assert_eq!(p.mass(), 1.0);
    }

    #[test]
    fn test_nbody_state_linear_space() {
        // Test size 2 N-body system
        let s1 = NBodyState::<2>::new(Vector::new([
            Vector::new([1.0, 2.0, 3.0]),
            Vector::new([4.0, 5.0, 6.0]),
        ]));

        // Scalar multiplication
        let s2 = s1 * 2.0;
        let d2 = s2.get_data();
        assert_eq!(d2, [[2.0, 4.0, 6.0], [8.0, 10.0, 12.0]]);

        // Addition
        let s3 = s1 + s2;
        assert_eq!(s3.get_data(), [[3.0, 6.0, 9.0], [12.0, 15.0, 18.0]]);

        // Subtraction
        let s4 = s3 - s1;
        assert_eq!(s4, s2);

        // Negation
        let s5 = -s1;
        assert_eq!(s5.get_data(), [[-1.0, -2.0, -3.0], [-4.0, -5.0, -6.0]]);

        // Division
        let s6 = s2 / 2.0;
        assert_eq!(s6, s1);

        // Zero
        let sz = NBodyState::<2>::zero();
        assert_eq!(sz.get_data(), [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
    }

    #[test]
    fn test_nbody_state_conversions() {
        let p1 = Particle::new(
            1.0,
            Vector::new([1.0, 0.0, 0.0]),
            Vector::new([0.0, 1.0, 0.0]),
        );
        let p2 = Particle::new(
            2.0,
            Vector::new([0.0, 1.0, 0.0]),
            Vector::new([1.0, 0.0, 0.0]),
        );
        let ps = [p1, p2];

        let pos_state = NBodyState::from_particles_pos(&ps);
        assert_eq!(pos_state.get_data(), [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);

        let vel_state = NBodyState::from_particles_vel(&ps);
        assert_eq!(vel_state.get_data(), [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]);
    }

    #[test]
    fn test_nbody_state_data_roundtrip() {
        let data = [[1.1, 2.2, 3.3], [4.4, 5.5, 6.6], [7.7, 8.8, 9.9]];
        // Use LinearSpace::new to create from raw data
        let state = <NBodyState<3> as LinearSpace>::new(data);
        assert_eq!(state.get_data(), data);
    }
}
