pub type IntegrationMethod = fn(DerivativeEquation, f64, f64, f64) -> f64;
pub type DerivativeEquation = fn(f64, f64) -> f64;

/// Euler Method
/// 초기값 t0, y0
/// (1) dy/dt = f(t, y)
/// (2) y(t_i+1) = y(t_i) + h*y'(t_i) + O(h^2)
/// (1)과 (2)를 연립 -> * y_i+1 = y_i + h*f(t, y)
///
/// euler_method(f, t, y, h) -> y_i+1
pub fn euler_method(f: DerivativeEquation, t: f64, y: f64, h: f64) -> f64 {
    y + h * f(t, y)
}

/// 4th Runge-Kutta
/// 초기값 t0, y0
/// dy/dt = f(t, y)
/// k1 = f(t_i, y_i)
/// k2 = f(t_i + h/2, y_i + h*k1/2)
/// k3 = f(t_i + h/2, y_i + h*k2/2)
/// k4 = f(t_i + h, y_i + h*k3)
/// y_i+1 = y_i + h/6 * (k1 + 2*k2 + 2*k3 + k4)
///
/// RK4(f, t, y, h) -> y_i+1
pub fn runge_kutta_4th(f: DerivativeEquation, t: f64, y: f64, h: f64) -> f64 {
    let k1 = f(t, y);
    let k2 = f(t + h / 2.0, y + h * k1 / 2.0);
    let k3 = f(t + h / 2.0, y + h * k2 / 2.0);
    let k4 = f(t + h, y + h * k3);
    y + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}

pub fn print_results(ts: Vec<f64>, ys: Vec<f64>) {
    println!("   t :            y");
    for (t, y) in ts.iter().zip(ys.iter()) {
        println!("{t:.2} : {y:.10}")
    }
}

pub fn integrate_step(
    f: DerivativeEquation,
    integrate_method: IntegrationMethod,
    y0: f64,
    t_start: f64,
    t_end: f64,
    h: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut t = t_start; // 초기값
    let mut y = y0; // 초기값
    let mut ts = vec![t];
    let mut ys = vec![y];
    while t <= t_end {
        y = integrate_method(f, t, y, h);
        t += h;
        ts.push(t);
        ys.push(y);
    }
    (ts, ys)
}
