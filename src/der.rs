use crate::consts::*;
use std::f64::consts::PI;
struct ComputeEnv {
    F0: f64,
    n: f64,
}
fn f2(z: f64) -> f64 {
    z
}
impl ComputeEnv {
    fn phi_n(&self, t: f64) -> f64 {
        f64::sqrt(self.n / PI) * f64::exp(-self.n * t * t)
    }
    //d^2y/dx^2=
    fn f1(&self, x: f64, y: f64, z: f64) -> f64 {
        -G * A / (h_b * m) * y + self.F0 /m* self.phi_n(x)
    }

    pub fn rk4(&self, x: f64, y: f64, z: f64, h: f64) -> (f64, f64) {
        let k1 = h * f2(z);
        let l1 = h * self.f1(x, y, z);
        let k2 = h * f2(z + 0.5 * l1);
        let l2 = h * self.f1(x + h / 2., y + 0.5 * k1, z + 0.5 * l1);
        let k3 = h * f2(z + 0.5 * l2);
        let l3 = h * self.f1(x + h / 2., y + 0.5 * k2, z + 0.5 * l2);
        let k4 = h * f2(z + l3);
        let l4 = h * self.f1(x + h , y + k3, z + l3);

        (
            y + (k1 + 2. * (k2 + k3) + k4) / 6.,
            z + (l1 + 2. * (l2 + l3) + l4) / 6.,
        )
    }
}

pub fn solve_rk4(tx: f64, n: i32, phi_n_n: f64, F0: f64) -> Vec<(f64, f64)> {
    let h = tx / n as f64;
    let mut y = 0.;
    let mut z = 0.;
    let mut r = Vec::new();
    for step in -n..n {
        let env = ComputeEnv { F0, n: phi_n_n };
        (y, z) = env.rk4(h * f64::from(step), y, z, h);
        r.push((h * f64::from(step), y));
    }
    r
}
