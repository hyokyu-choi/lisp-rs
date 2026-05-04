use crate::{
    math::{
        complex::{Complex, ComplexSpace},
        core::{LinearSpace, Vector},
        field::{Field, Field3D},
        integrate::System,
    },
    physics::particle::{NBodyState, Particle},
};
use crate::physics::constants::G;
use std::cell::RefCell;

/// $O(N \log N)$ Particle-Mesh Gravity System
pub struct PMSystem<const N: usize, const GRID_SIZE: usize> {
    pub masses: [f64; N],
    pub potential_solver: RefCell<GravitationalPotential<GRID_SIZE>>,
}

impl<const N: usize, const GRID_SIZE: usize> PMSystem<N, GRID_SIZE> {
    pub fn new(masses: [f64; N], solver: GravitationalPotential<GRID_SIZE>) -> Self {
        Self {
            masses,
            potential_solver: RefCell::new(solver),
        }
    }

    /// Helper to rebuild Particle objects from state for the solver
    fn get_particles(&self, y: NBodyState<N>, y_prime: NBodyState<N>) -> Vec<Particle> {
        (0..N)
            .map(|i| {
                Particle::new(
                    self.masses[i],
                    y.components[i],
                    y_prime.components[i],
                )
            })
            .collect()
    }
}

impl<const N: usize, const GRID_SIZE: usize> System for PMSystem<N, GRID_SIZE> {
    type Vector = NBodyState<N>;

    fn derivative(&self, _t: f64, y: Self::Vector, y_prime: Self::Vector) -> Self::Vector {
        // 1. Reconstruct particles
        let particles = self.get_particles(y, y_prime);
        
        // 2. Compute gravitational field using PM method
        let mut solver = self.potential_solver.borrow_mut();
        solver.update(&particles);

        // 3. Get acceleration for each particle
        let mut acc_components = Vector::zero();
        for i in 0..N {
            acc_components[i] = solver.get_acceleration(&particles[i]);
        }
        
        NBodyState::new(acc_components)
    }
}

pub struct GravitationalPotential<const N: usize> {
    pub x_min: f64, pub _x_max: f64,
    pub y_min: f64, pub _y_max: f64,
    pub z_min: f64, pub _z_max: f64,
    pub dx: f64, pub dy: f64, pub dz: f64,
    pub density_field: Box<Field3D<N>>,
    pub potential_field: Box<Field3D<N>>,
    pub green_function: Box<Field3D<N>>,
}

impl<const N: usize> GravitationalPotential<N> {
    pub fn new(
        x_min: f64, x_max: f64,
        y_min: f64, y_max: f64,
        z_min: f64, z_max: f64,
    ) -> Self {
        let lx = x_max - x_min;
        let ly = y_max - y_min;
        let lz = z_max - z_min;

        let dx = lx / (N as f64 / 2.0);
        let dy = ly / (N as f64 / 2.0);
        let dz = lz / (N as f64 / 2.0);

        let density_field = Box::new(Field3D::zero());
        let potential_field = Box::new(Field3D::zero());
        let mut green_function = Box::new(Field3D::zero());

        // Precompute Isolated Green's Function: G(r) = -G / |r|
        for i in 0..N {
            let rx = if i <= N / 2 { i as f64 } else { (i as f64) - (N as f64) };
            let rx = rx * dx;
            for j in 0..N {
                let ry = if j <= N / 2 { j as f64 } else { (j as f64) - (N as f64) };
                let ry = ry * dy;
                for k in 0..N {
                    let rz = if k <= N / 2 { k as f64 } else { (k as f64) - (N as f64) };
                    let rz = rz * dz;

                    let r_sq = rx * rx + ry * ry + rz * rz;
                    if r_sq > 0.0 {
                        green_function[i][j][k] = Complex::new(-G / r_sq.sqrt(), 0.0);
                    } else {
                        // Softening at the origin to avoid infinity
                        green_function[i][j][k] = Complex::new(-G / (0.5 * (dx*dx + dy*dy + dz*dz).sqrt()), 0.0);
                    }
                }
            }
        }

        green_function.fft();

        Self { 
            x_min, _x_max: x_max,
            y_min, _y_max: y_max,
            z_min, _z_max: z_max,
            dx, dy, dz,
            density_field,
            potential_field,
            green_function,
        }
    }

    /// Step 1: Assign particle masses to the grid using CIC (Cloud-In-Cell)
    fn assign_mass_to_grid(&mut self, ps: &[Particle]) {
        *self.density_field = Field3D::zero();
        for p in ps {
            self.map_particle_to_grid(p);
        }
    }

    fn map_particle_to_grid(&mut self, p: &Particle) {
        let pos = p.position();
        
        let offset = (N as f64) / 4.0;
        let u = (pos[0] - self.x_min) / self.dx + offset;
        let v = (pos[1] - self.y_min) / self.dy + offset;
        let w = (pos[2] - self.z_min) / self.dz + offset;

        if u < 0.0 || u >= (N as f64) - 1.0 || v < 0.0 || v >= (N as f64) - 1.0 || w < 0.0 || w >= (N as f64) - 1.0 {
            return;
        }

        let i = u.floor() as usize;
        let j = v.floor() as usize;
        let k = w.floor() as usize;

        let du = u - i as f64;
        let dv = v - j as f64;
        let dw = w - k as f64;

        let m = p.mass();
        let w000 = (1.0 - du) * (1.0 - dv) * (1.0 - dw);
        let w100 = du * (1.0 - dv) * (1.0 - dw);
        let w010 = (1.0 - du) * dv * (1.0 - dw);
        let w001 = (1.0 - du) * (1.0 - dv) * dw;
        let w110 = du * dv * (1.0 - dw);
        let w101 = du * (1.0 - dv) * dw;
        let w011 = (1.0 - du) * dv * dw;
        let w111 = du * dv * dw;

        self.density_field[i][j][k] = self.density_field[i][j][k] + Complex::new(m * w000, 0.0);
        self.density_field[i+1][j][k] = self.density_field[i+1][j][k] + Complex::new(m * w100, 0.0);
        self.density_field[i][j+1][k] = self.density_field[i][j+1][k] + Complex::new(m * w010, 0.0);
        self.density_field[i][j][k+1] = self.density_field[i][j][k+1] + Complex::new(m * w001, 0.0);
        self.density_field[i+1][j+1][k] = self.density_field[i+1][j+1][k] + Complex::new(m * w110, 0.0);
        self.density_field[i+1][j][k+1] = self.density_field[i+1][j][k+1] + Complex::new(m * w101, 0.0);
        self.density_field[i][j+1][k+1] = self.density_field[i][j+1][k+1] + Complex::new(m * w011, 0.0);
        self.density_field[i+1][j+1][k+1] = self.density_field[i+1][j+1][k+1] + Complex::new(m * w111, 0.0);
    }

    /// Step 2: Solve Poisson equation in Fourier space
    fn solve_poisson_eq(&mut self) {
        self.density_field.fft();

        let volume = self.dx * self.dy * self.dz;
        for i in 0..N {
            for j in 0..N {
                for k in 0..N {
                    let hat_rho = self.density_field[i][j][k];
                    let hat_g = self.green_function[i][j][k];
                    self.potential_field[i][j][k] = hat_rho * hat_g * volume;
                }
            }
        }

        self.potential_field.ifft();
    }

    /// Step 3: Interpolate acceleration from potential field gradient
    pub fn get_acceleration(&self, p: &Particle) -> Vector<f64, 3> {
        let pos = p.position();
        
        let offset = (N as f64) / 4.0;
        let u = (pos[0] - self.x_min) / self.dx + offset;
        let v = (pos[1] - self.y_min) / self.dy + offset;
        let w = (pos[2] - self.z_min) / self.dz + offset;

        if u < 1.0 || u >= (N as f64) - 2.0 || v < 1.0 || v >= (N as f64) - 2.0 || w < 1.0 || w >= (N as f64) - 2.0 {
            return Vector::zero();
        }

        let i = u.floor() as usize;
        let j = v.floor() as usize;
        let k = w.floor() as usize;

        let du = u - i as f64;
        let dv = v - j as f64;
        let dw = w - k as f64;

        let get_g = |i: usize, j: usize, k: usize| {
            let gx = -(self.potential_field[i+1][j][k].re() - self.potential_field[i-1][j][k].re()) / (2.0 * self.dx);
            let gy = -(self.potential_field[i][j+1][k].re() - self.potential_field[i][j-1][k].re()) / (2.0 * self.dy);
            let gz = -(self.potential_field[i][j][k+1].re() - self.potential_field[i][j][k-1].re()) / (2.0 * self.dz);
            Vector::new([gx, gy, gz])
        };

        let g000 = get_g(i, j, k);
        let g100 = get_g(i+1, j, k);
        let g010 = get_g(i, j+1, k);
        let g001 = get_g(i, j, k+1);
        let g110 = get_g(i+1, j+1, k);
        let g101 = get_g(i+1, j, k+1);
        let g011 = get_g(i, j+1, k+1);
        let g111 = get_g(i+1, j+1, k+1);

        let w000 = (1.0 - du) * (1.0 - dv) * (1.0 - dw);
        let w100 = du * (1.0 - dv) * (1.0 - dw);
        let w010 = (1.0 - du) * dv * (1.0 - dw);
        let w001 = (1.0 - du) * (1.0 - dv) * dw;
        let w110 = du * dv * (1.0 - dw);
        let w101 = du * (1.0 - dv) * dw;
        let w011 = (1.0 - du) * dv * dw;
        let w111 = du * dv * dw;

        g000 * w000 + g100 * w100 + g010 * w010 + g001 * w001 + 
        g110 * w110 + g101 * w101 + g011 * w011 + g111 * w111
    }

    pub fn update(&mut self, ps: &[Particle]) {
        self.assign_mass_to_grid(ps);
        self.solve_poisson_eq();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::complex::ComplexSpace;
    use crate::math::core::VectorSpace;

    #[test]
    fn test_potential_solver_init() {
        let solver = GravitationalPotential::<8>::new(0.0, 4.0, 0.0, 4.0, 0.0, 4.0);
        assert_eq!(solver.dx, 1.0);
    }

    #[test]
    fn test_pm_system_init() {
        let masses = [1.0, 2.0];
        let solver = GravitationalPotential::<8>::new(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);
        let system = PMSystem::new(masses, solver);
        assert_eq!(system.masses[0], 1.0);
    }

    #[test]
    fn test_cic_mass_assignment_exact_point() {
        let mut solver = GravitationalPotential::<8>::new(0.0, 4.0, 0.0, 4.0, 0.0, 4.0);
        let p = Particle::new(10.0, Vector::new([1.0, 1.0, 1.0]), Vector::zero());
        solver.assign_mass_to_grid(&[p]);
        assert_eq!(solver.density_field[3][3][3].re(), 10.0);
    }

    #[test]
    fn test_boundary_particles() {
        let mut solver = GravitationalPotential::<8>::new(0.0, 4.0, 0.0, 4.0, 0.0, 4.0);
        let p_out = Particle::new(10.0, Vector::new([-5.0, 2.0, 2.0]), Vector::zero());
        let p_in = Particle::new(10.0, Vector::new([2.0, 2.0, 2.0]), Vector::zero());
        solver.assign_mass_to_grid(&[p_out, p_in]);
        let total_mass: f64 = (0..8).map(|x| {
            (0..8).map(|y| {
                (0..8).map(|z| solver.density_field[x][y][z].re()).sum::<f64>()
            }).sum::<f64>()
        }).sum();
        assert!((total_mass - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_superposition_principle() {
        let mut solver = GravitationalPotential::<8>::new(0.0, 4.0, 0.0, 4.0, 0.0, 4.0);
        let p_test = Particle::new(1.0, Vector::new([1.0, 2.0, 2.0]), Vector::zero());
        let p_a = Particle::new(50.0, Vector::new([2.0, 2.0, 2.0]), Vector::zero());
        solver.update(&[p_a]);
        let acc_a = solver.get_acceleration(&p_test);
        let p_b = Particle::new(50.0, Vector::new([2.0, 3.0, 2.0]), Vector::zero());
        solver.update(&[p_b]);
        let acc_b = solver.get_acceleration(&p_test);
        solver.update(&[p_a, p_b]);
        let acc_both = solver.get_acceleration(&p_test);
        let diff = acc_both - (acc_a + acc_b);
        assert!(diff.norm() < 0.1);
    }

    #[test]
    fn test_inverse_square_scaling() {
        // Use N=16 for better resolution, but keep function scope small for stack safety.
        let mut solver = GravitationalPotential::<16>::new(0.0, 8.0, 0.0, 8.0, 0.0, 8.0);
        let p_central = Particle::new(1000.0, Vector::new([4.0, 4.0, 4.0]), Vector::zero());
        
        // Test at r=2.0 and r=4.0. dx=1.0, so these are 2 and 4 cells away.
        let p_r2 = Particle::new(1.0, Vector::new([6.0, 4.0, 4.0]), Vector::zero());
        let p_r4 = Particle::new(1.0, Vector::new([8.0, 4.0, 4.0]), Vector::zero());
        
        solver.update(&[p_central]);
        let acc_r2 = solver.get_acceleration(&p_r2).norm();
        let acc_r4 = solver.get_acceleration(&p_r4).norm();
        
        let ratio = acc_r2 / acc_r4;
        println!("Acceleration ratio r=2.0 vs r=4.0 (expected ~4.0): {}", ratio);
        // In the far field, ratio should approach (4/2)^2 = 4.0
        assert!(ratio > 2.0, "Force must significantly decrease in the far field. Ratio: {}", ratio);
    }

    #[test]
    fn test_empty_field() {
        let mut solver = GravitationalPotential::<8>::new(0.0, 4.0, 0.0, 4.0, 0.0, 4.0);
        solver.update(&[]);
        let p_test = Particle::new(4.0, Vector::new([2.0, 2.0, 2.0]), Vector::zero());
        let acc = solver.get_acceleration(&p_test);
        assert_eq!(acc, Vector::zero());
    }
}
