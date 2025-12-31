use std::fmt;

#[derive(Clone, Copy)]
pub struct Particle {
    mass: f64,
    position: f64,
    velocity: f64,
}

impl Particle {
    pub fn new(m: f64, x: f64, v: f64) -> Self {
        Self {
            mass: m,
            position: x,
            velocity: v,
        }
    }
    pub fn mass(&self) -> f64 {
        self.mass
    }
    pub fn position(&self) -> f64 {
        self.position
    }
    pub fn velocity(&self) -> f64 {
        self.velocity
    }
    pub fn momentum(&self) -> f64 {
        self.mass * self.velocity
    }
    pub fn set_position(&mut self, x: f64) {
        self.position = x;
    }
    pub fn set_velocity(&mut self, x: f64) {
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
