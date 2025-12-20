mod utils;

mod physics_module {
    pub mod constants;
}

mod math_module {
    pub mod numerical_integration;
    pub mod vector;
}

use math_module::numerical_integration::{test, test_2};

fn main() {
    test_2(0.1);
}
