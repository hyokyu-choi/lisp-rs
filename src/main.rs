mod utils;

mod physics_module {
    pub mod constants;
    pub mod particle;
}

mod math_module {
    pub mod numerical_integration;
    pub mod vector;
}

use physics_module::constants::PI;
use math_module::numerical_integration::{test, test_2};

fn main() {
    test_2(0.001);
}
