mod utils {
    pub mod plot;
}

mod physics {
    pub mod harmonic_oscillator;
}

mod math {
    pub mod core;
    pub mod integrate;
    pub mod matrix;
}

use crate::math::integrate::{RK4Method, Solver};
use crate::utils::plot::plot_one;

use crate::physics::harmonic_oscillator::DrivenHarmonicOscillator;

fn main() {
    let dho_ode = DrivenHarmonicOscillator {
        k: 16.0,
        b: 1.0,
        f0: 5.0,
        omega: 2.0,
    };

    let y0 = 4.0;
    let y0_prime = 0.0;
    let h = 0.01;
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
