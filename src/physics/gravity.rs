use crate::{math::{complex::Complex, core::LinearSpace, field::{Field, Field3D}}, physics::particle::Particle};
use std::f64::consts::PI;
use crate::physics::constants::G;

pub struct GravitationalPotential<const N: usize> {
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    z_min: f64,
    z_max: f64,
    field : Field3D<N>,
}

impl<const N: usize> GravitationalPotential<N> {
    pub fn new(
        x_min: f64,
        x_max: f64,
        y_min: f64,
        y_max: f64,
        z_min: f64,
        z_max: f64,
    ) -> Self {
        Self { 
            x_min,
            x_max,
            y_min,
            y_max,
            z_min,
            z_max,
            field: Field3D::zero(),
        }
    }
    fn init(&mut self) {
        for x in 0..N {
            for y in 0..N {
                for z in 0..N {
                    self.field[x][y][z] = Complex::zero();
                }
            }
        }
    }
    fn set_mass(&mut self, p: Particle) {
        // Get mass density distribution
        // Maybe object need mass density and radius (mass and volume)
        todo!()
    }
    fn solve_poisson_eq(&mut self) {
        self.field.fft();
        for x in 0..N {
            let kx = if x <= N/2 { x as f64 } else { (x as f64) - (N as f64)/2.0 };
            for y in 0..N {
                let ky = if y <= N/2 { y as f64 } else { (y as f64) - (N as f64)/2.0 };
                for z in 0..N {
                    let kz = if z <= N/2 { z as f64 } else { (z as f64) - (N as f64)/2.0 };
                    let k_sq = 1.0/((kx * kx + ky * ky + kz * kz) * PI);
                    self.field[x][y][z] = self.field[x][y][z] * G * -k_sq;
                }
            }
        }
        self.field.ifft();
    }
    pub fn step(&mut self, ps: Vec<Particle>) {
        self.init();
        for p in ps {
            self.set_mass(p);
        }
        // ps.iter().map(|p| self.set_mass(*p));
        self.solve_poisson_eq();
    }
    pub fn get_gravitational_field(&self, x: f64, y: f64, z: f64) {
        // 3D Interpolation and calculate gradient
        todo!()
    }
}
