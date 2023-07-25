use crate::consts::*;
use std::f64::consts::PI;
struct ComputeEnv {
    F0: f64,
    n: f64,
}
fn f2(dxdt: f64) -> f64 {
    dxdt
}
impl ComputeEnv {
    fn phi_n(&self, t: f64) -> f64 {
        f64::sqrt(self.n / PI) * f64::exp(-self.n * t * t)
    }
    //d^2y/dx^2=
    fn f1(&self, t: f64, x: f64, dxdt: f64) -> f64 {
        -G * A / (h_b * m) * x + self.F0 / m * self.phi_n(t)
    }

    pub fn rk4(&self, t: f64, x: f64, dxdt: f64, h: f64) -> (f64, f64) {
        let k1 = h * f2(dxdt);
        let l1 = h * self.f1(t, x, dxdt);
        let k2 = h * f2(dxdt + 0.5 * l1);
        let l2 = h * self.f1(t + h / 2., x + 0.5 * k1, dxdt + 0.5 * l1);
        let k3 = h * f2(dxdt + 0.5 * l2);
        let l3 = h * self.f1(t + h / 2., x + 0.5 * k2, dxdt + 0.5 * l2);
        let k4 = h * f2(dxdt + l3);
        let l4 = h * self.f1(t + h, x + k3, dxdt + l3);

        (
            x + (k1 + 2. * (k2 + k3) + k4) / 6.,
            dxdt + (l1 + 2. * (l2 + l3) + l4) / 6.,
        )
    }
}

pub fn solve_rk4(tx: f64, n: i32, phi_n_n: f64, F0: f64) -> () {
    let h = tx / n as f64;
    let mut x = 0.;
    let mut dxdt = 0.;
    let mut sign_positive = false;
    let mut cnt = 0;
    let mut ts = Vec::new();
    let mut r = Vec::new();
    for step in -n..(2*n) {
        let env = ComputeEnv { F0, n: phi_n_n };
        let t = h * f64::from(step);
        (x, dxdt) = env.rk4(t, x, dxdt, h);
        if t > 0.  {
            if dxdt.is_sign_positive() != sign_positive {
                ts.push(t);
                r.push(x);
            }
            sign_positive = dxdt.is_sign_positive();
        }
        //println!("{}\t{}", h * f64::from(step), x);
    }
    println!(
        "{}",
        r.windows(2)
            .map(|e| (e[1] - e[0].abs()))
            .fold(0. / 0., f64::max)/2.
    );
}
