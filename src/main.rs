use std::collections::HashMap;
use std::f64::consts::PI;
const dt: f64 = 0.0000001;
const G: f64 = 79E9;
const d_d: f64 = 0.15;
const A: f64 = PI / 4.0 * d_d * d_d;
const m: f64 = 6.0;
const h_b: f64 = 0.05;
const h: f64 = h_b;
const l_a: f64 = 0.6;
const l_b: f64 = 0.93;
const d_B: f64 = 0.5;
const l: f64 = d_B / 2.0 + (l_a + l_b) - 0.05;
const M_B: f64 = 70.0;
const d_a: f64 = 0.1;
const M_a: f64 = 4.0;
const I_Ga: f64 = (d_a * d_a / 16.0 + l_a * l_a / 12.0) * M_a;
const I_GB: f64 = 1.0 / 8.0 * M_B * d_B * d_B;
const M_b: f64 = 1.0;
const I_Gb: f64 = 1.0 / 12.0 * M_b * (h_b * h_b + l_b * l_b);
const I: f64 = I_GB
    + (I_Ga + 1.0 / 4.0 * M_a * (d_B + l_a) * (d_B + l_a))
    + I_Gb
    + 1.0 / 4.0 * M_b * (d_B + 2.0 * l_a + l_b) * (d_B + 2.0 * l_a + l_b);

struct ComputeEnv {
    F0: f64,
    n: f64,
    cache_x: HashMap<i32, f64>,
    cache_dxdt: HashMap<i32, f64>,
}
impl ComputeEnv {
    fn x(&mut self, step: i32) -> f64 {
        let r = match self.cache_x.get(&step) {
            Some(x) => *x,
            None => {
                let t = f64::from(step) * dt;
                if t < -0.3 {
                    return 0.0;
                }
                let x =
                    self.x(step - 1) + self.dxdt(step - 1) * dt + self.dxdt2(step - 1) * dt * dt;
                self.cache_x.insert(step, x);
                // println!("x@{step}:{x}");
                x
            }
        };

        r
    }

    fn dxdt2(&mut self, step: i32) -> f64 {
        let t = f64::from(step) * dt;
        self.F0 / m * self.phi_n(t) - G * A / (h * m) * self.x(step)
    }

    fn dxdt3(&mut self, step: i32) -> f64 {
        let t = f64::from(step) * dt;
        self.F0 / m * self.dphin_ndt(t) - G * A / (h * m) * self.dxdt(step)
    }

    fn dxdt(&mut self, step: i32) -> f64 {
        let t = f64::from(step) * dt;
        let r = match self.cache_dxdt.get(&step) {
            Some(dxdt) => *dxdt,
            None => {
                if t < -0.3 {
                    return 0.0;
                }

                let dxdt = self.dxdt(step - 1)
                    + self.dxdt2(step - 1) * dt
                    + self.dxdt3(step - 1) * dt * dt;
                self.cache_dxdt.insert(step, dxdt);
                // println!("dx@{step}:{dxdt}");
                dxdt
            }
        };

        r
    }
    fn phi_n(&self, t: f64) -> f64 {
        f64::sqrt(self.n / PI) * f64::exp(-self.n * t * t)
    }
    fn dphin_ndt(&self, t: f64) -> f64 {
        2.0 * t * (-self.n) * f64::sqrt(self.n / PI) * f64::exp(-self.n * t * t)
    }
}
fn compute_f0(omega: f64) -> (f64, f64) {
    let mut F0_min = 0.0;
    let mut F0_max = 1000.0;
    let mut F0 = 0.0;
    let mut d = 0.0;
    const target_step: i32 = 3000000;
    for it in 0..100 {
        eprintln!("{omega}:{it}");
        F0 = F0_max / 2.0 + F0_min / 2.0;
        let mut env = ComputeEnv {
            F0,
            n: 100.0,
            cache_dxdt: HashMap::new(),
            cache_x: HashMap::new(),
        };
        for step in -target_step..target_step{
            env.dxdt(step);
            
        }
        let left = env.dxdt(target_step);

        let right = l * (omega - l * F0 / I);
        if left > right {
            F0_max = F0;
        } else {
            F0_min = F0;
        }
        d = (left - right).abs();
        if d < 1.0 {
            break;
        }
    }
    return (F0, d);
}

fn sub_main() {
    //2.0 * PI * 33.0;

    for s in 0..6 {
        let omega = PI * f64::from(s)*11.0;
        let (F0, d) = compute_f0(omega);

        println!("{},{},{}", omega, F0, d);
    }
}
fn main() {
    const STACK_SIZE: usize = 512 * 1024 * 1024;
    std::thread::Builder::new()
        .stack_size(STACK_SIZE)
        .spawn(sub_main)
        .unwrap()
        .join()
        .unwrap();
}
