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
pub struct PMSystem<const N: usize> {
    pub masses: [f64; N],
    pub potential_solver: RefCell<GravitationalPotential>,
}

impl<const N: usize> PMSystem<N> {
    pub fn new(masses: [f64; N], solver: GravitationalPotential) -> Self {
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

impl<const N: usize> System for PMSystem<N> {
    type Vector = NBodyState<N>;

    fn derivative(&self, _t: f64, y: Self::Vector, y_prime: Self::Vector) -> Self::Vector {
        let particles = self.get_particles(y, y_prime);
        let mut solver = self.potential_solver.borrow_mut();
        solver.update(&particles);

        let mut acc_components = Vector::zero();
        for i in 0..N {
            acc_components[i] = solver.get_acceleration(&particles[i]);
        }
        
        NBodyState::new(acc_components)
    }
}

pub struct GravitationalPotential {
    pub n: usize,
    pub x_min: f64, pub _x_max: f64,
    pub y_min: f64, pub _y_max: f64,
    pub z_min: f64, pub _z_max: f64,
    pub dx: f64, pub dy: f64, pub dz: f64,
    pub density_field: Field3D,
    pub potential_field: Field3D,
    pub green_function: Field3D,
}

impl GravitationalPotential {
    pub fn new(
        n_physical: usize,
        x_min: f64, x_max: f64,
        y_min: f64, y_max: f64,
        z_min: f64, z_max: f64,
    ) -> Self {
        // Hockney's method: Total grid size N = 2 * n_physical (padded to power of 2)
        let n_total = (n_physical * 2).next_power_of_two();
        
        let lx = x_max - x_min;
        let ly = y_max - y_min;
        let lz = z_max - z_min;

        // dx is based on the physical resolution
        let dx = lx / (n_physical as f64);
        let dy = ly / (n_physical as f64);
        let dz = lz / (n_physical as f64);

        let mut green_function = Field3D::new(n_total);
        let n = green_function.n;

        for i in 0..n {
            let rx = if i <= n / 2 { i as f64 } else { (i as f64) - (n as f64) };
            let rx = rx * dx;
            for j in 0..n {
                let ry = if j <= n / 2 { j as f64 } else { (j as f64) - (n as f64) };
                let ry = ry * dy;
                for k in 0..n {
                    let rz = if k <= n / 2 { k as f64 } else { (k as f64) - (n as f64) };
                    let rz = rz * dz;

                    let r_sq = rx * rx + ry * ry + rz * rz;
                    if r_sq > 0.0 {
                        green_function[[i, j, k]] = Complex::new(-G / r_sq.sqrt(), 0.0);
                    } else {
                        green_function[[i, j, k]] = Complex::new(-G / (0.5 * (dx*dx + dy*dy + dz*dz).sqrt()), 0.0);
                    }
                }
            }
        }

        green_function.fft();

        Self { 
            n,
            x_min, _x_max: x_max,
            y_min, _y_max: y_max,
            z_min, _z_max: z_max,
            dx, dy, dz,
            density_field: Field3D::new(n),
            potential_field: Field3D::new(n),
            green_function,
        }
    }

    fn assign_mass_to_grid(&mut self, ps: &[Particle]) {
        self.density_field.data.fill(Complex::zero());
        for p in ps {
            self.map_particle_to_grid(p);
        }
    }

    fn map_particle_to_grid(&mut self, p: &Particle) {
        let pos = p.position();
        let n = self.n;
        
        // Offset N/4 assuming n_physical = n/2
        let offset = (n as f64) / 4.0;
        let u = (pos[0] - self.x_min) / self.dx + offset;
        let v = (pos[1] - self.y_min) / self.dy + offset;
        let w = (pos[2] - self.z_min) / self.dz + offset;

        if u < 0.0 || u >= (n as f64) - 1.0 || v < 0.0 || v >= (n as f64) - 1.0 || w < 0.0 || w >= (n as f64) - 1.0 {
            return;
        }

        let i = u.floor() as usize;
        let j = v.floor() as usize;
        let k = w.floor() as usize;

        let du = u - i as f64;
        let dv = v - j as f64;
        let dw = w - k as f64;

        let m = p.mass();
        self.density_field[[i, j, k]] = self.density_field[[i, j, k]] + Complex::new(m * (1.0 - du) * (1.0 - dv) * (1.0 - dw), 0.0);
        self.density_field[[i+1, j, k]] = self.density_field[[i+1, j, k]] + Complex::new(m * du * (1.0 - dv) * (1.0 - dw), 0.0);
        self.density_field[[i, j+1, k]] = self.density_field[[i, j+1, k]] + Complex::new(m * (1.0 - du) * dv * (1.0 - dw), 0.0);
        self.density_field[[i, j, k+1]] = self.density_field[[i, j, k+1]] + Complex::new(m * (1.0 - du) * (1.0 - dv) * dw, 0.0);
        self.density_field[[i+1, j+1, k]] = self.density_field[[i+1, j+1, k]] + Complex::new(m * du * dv * (1.0 - dw), 0.0);
        self.density_field[[i+1, j, k+1]] = self.density_field[[i+1, j, k+1]] + Complex::new(m * du * (1.0 - dv) * dw, 0.0);
        self.density_field[[i, j+1, k+1]] = self.density_field[[i, j+1, k+1]] + Complex::new(m * (1.0 - du) * dv * dw, 0.0);
        self.density_field[[i+1, j+1, k+1]] = self.density_field[[i+1, j+1, k+1]] + Complex::new(m * du * dv * dw, 0.0);
    }

    fn solve_poisson_eq(&mut self) {
        self.density_field.fft();

        let n = self.n;
        let volume = self.dx * self.dy * self.dz;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let hat_rho = self.density_field[[i, j, k]];
                    let hat_g = self.green_function[[i, j, k]];
                    self.potential_field[[i, j, k]] = hat_rho * hat_g * volume;
                }
            }
        }

        self.potential_field.ifft();
    }

    pub fn get_acceleration(&self, p: &Particle) -> Vector<f64, 3> {
        let pos = p.position();
        let n = self.n;
        
        let offset = (n as f64) / 4.0;
        let u = (pos[0] - self.x_min) / self.dx + offset;
        let v = (pos[1] - self.y_min) / self.dy + offset;
        let w = (pos[2] - self.z_min) / self.dz + offset;

        if u < 1.0 || u >= (n as f64) - 2.0 || v < 1.0 || v >= (n as f64) - 2.0 || w < 1.0 || w >= (n as f64) - 2.0 {
            return Vector::zero();
        }

        let i = u.floor() as usize;
        let j = v.floor() as usize;
        let k = w.floor() as usize;

        let du = u - i as f64;
        let dv = v - j as f64;
        let dw = w - k as f64;

        let get_g = |i: usize, j: usize, k: usize| {
            let gx = -(self.potential_field[[i+1, j, k]].re() - self.potential_field[[i-1, j, k]].re()) / (2.0 * self.dx);
            let gy = -(self.potential_field[[i, j+1, k]].re() - self.potential_field[[i, j-1, k]].re()) / (2.0 * self.dy);
            let gz = -(self.potential_field[[i, j, k+1]].re() - self.potential_field[[i, j, k-1]].re()) / (2.0 * self.dz);
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
        // Physical N=8 -> dx=1.0. Total N will be 16 due to Hockney padding.
        let solver = GravitationalPotential::new(8, 0.0, 8.0, 0.0, 8.0, 0.0, 8.0);
        assert_eq!(solver.dx, 1.0);
        assert_eq!(solver.n, 16);
    }

    #[test]
    fn test_pm_system_init() {
        let masses = [1.0, 2.0];
        let solver = GravitationalPotential::new(4, -10.0, 10.0, -10.0, 10.0, -10.0, 10.0);
        let system = PMSystem::new(masses, solver);
        assert_eq!(system.masses[0], 1.0);
    }

    #[test]
    fn test_cic_mass_assignment_exact_point() {
        // Physical N=4, Range=4.0 -> dx=1.0
        let mut solver = GravitationalPotential::new(4, 0.0, 4.0, 0.0, 4.0, 0.0, 4.0);
        // Total N=8, offset=2.0. (1.0, 1.0, 1.0) physical -> (3, 3, 3) grid
        let p = Particle::new(10.0, Vector::new([1.0, 1.0, 1.0]), Vector::zero());
        solver.assign_mass_to_grid(&[p]);
        assert_eq!(solver.density_field[[3, 3, 3]].re(), 10.0);
    }

    #[test]
    fn test_superposition_principle() {
        let mut solver = GravitationalPotential::new(4, 0.0, 4.0, 0.0, 4.0, 0.0, 4.0);
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
        // Large scale test to ensure no stack overflow and accuracy
        let mut solver = GravitationalPotential::new(16, 0.0, 16.0, 0.0, 16.0, 0.0, 16.0);
        let m_central = 1000.0;
        let p_central = Particle::new(m_central, Vector::new([8.0, 8.0, 8.0]), Vector::zero());
        let p_r2 = Particle::new(1.0, Vector::new([10.0, 8.0, 8.0]), Vector::zero());
        let p_r4 = Particle::new(1.0, Vector::new([12.0, 8.0, 8.0]), Vector::zero());
        
        solver.update(&[p_central]);
        let acc_r2 = solver.get_acceleration(&p_r2).norm();
        let acc_r4 = solver.get_acceleration(&p_r4).norm();
        
        let ratio = acc_r2 / acc_r4;
        println!("Acceleration ratio r=2 vs r=4 (expected ~4.0): {}", ratio);
        assert!(ratio > 3.0, "Force must follow inverse square law. Ratio: {}", ratio);
    }
}
