mod utils {
    pub mod plot;
}

mod physics {
    pub mod harmonic_oscillator;
}

mod math {
    pub mod core;
    pub mod integrate;
}

use crate::math::core::{Scalar, ScalarSpace};
use crate::math::integrate::{RK4Method, Solver};
use crate::physics::harmonic_oscillator::DampedHarmonicOscillator;
use crate::utils::plot::plot_one;

use crate::physics::harmonic_oscillator::DrivenHarmonicOscillator;

fn main() {
    let dho_ode: DrivenHarmonicOscillator = DrivenHarmonicOscillator {
        k: Scalar::new(16.0),
        b: Scalar::new(1.0),
        f0: Scalar::new(5.0),
        omega: Scalar::new(2.0),
    };

    let y0 = Scalar::new(4.0);
    let y0_prime = Scalar::new(0.0);
    let h = Scalar::new(0.01);
    let steps = 2000;

    let mut sho_solver = Solver::new(RK4Method, dho_ode, y0, y0_prime);
    sho_solver.run(h, steps);
    let (ts, ys, ys_prime): (Vec<f64>, Vec<f64>, Vec<f64>) = sho_solver.get_results_f64();

    let _ = plot_one(
        String::from("Driven Harmonic Oscillation"),
        ts,
        vec![ys, ys_prime],
        [0.0, 20.0],
        [-5.0, 5.0],
        vec![String::from("position"), String::from("velocity")],
        String::from("plotters-doc-data/2.png"),
    );
}
