use std::fmt;
use crate::math_module::core::{Scalar};

#[derive(Clone, Copy)]
pub struct Particle {
    mass: Scalar,
    position: Scalar,
    velocity: Scalar,
}

impl Particle {
    pub fn new(m: f64, x: f64, v: f64) -> Self {
        Self {
            mass: Scalar::new(m),
            position: Scalar::new(x),
            velocity: Scalar::new(v),
        }
    }
    pub fn mass(&self) -> Scalar {
        self.mass
    }
    pub fn position(&self) -> Scalar {
        self.position
    }
    pub fn velocity(&self) -> Scalar {
        self.velocity
    }
    pub fn momentum(&self) -> Scalar {
        self.mass * self.velocity
    }
    pub fn set_position(&mut self, x: Scalar) {
        self.position = x;
    }
    pub fn set_velocity(&mut self, x: Scalar) {
        self.velocity = x;
    }
}

impl fmt::Display for Particle {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Particle(mass: {}, x: {}, v: {})",
            self.mass, self.position, self.velocity
        )
    }
}
