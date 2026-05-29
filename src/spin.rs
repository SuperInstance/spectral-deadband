//! Module 2: Spin — time abstracted as rotation.
//!
//! A symplectic integrator treats time steps as rotations in phase space.
//! The angle of rotation IS the time step. Rotations preserve the metric → conservation.

/// A single spin: eigenvalue λ maps to angular frequency ω = √λ.
#[derive(Debug, Clone, Copy)]
pub struct Spin {
    /// Current phase angle (radians). Time abstracted as distance.
    pub angle: f64,
    /// Angular frequency ω = √λ.
    pub angular_velocity: f64,
    /// Amplitude (energy). Preserved under rotation.
    pub radius: f64,
}

impl Spin {
    /// Create a Spin from an eigenvalue λ. Angular velocity ω = √λ, period T = 2π/ω.
    pub fn from_eigenvalue(lambda: f64) -> Self {
        let omega = if lambda > 0.0 { lambda.sqrt() } else { 0.0 };
        Spin {
            angle: 0.0,
            angular_velocity: omega,
            radius: 1.0,
        }
    }

    /// Create a spin with explicit parameters.
    pub fn new(angle: f64, angular_velocity: f64, radius: f64) -> Self {
        Spin { angle, angular_velocity, radius }
    }

    /// Advance by dt: rotate the angle. This IS symplectic integration.
    pub fn advance(&self, dt: f64) -> Spin {
        Spin {
            angle: self.angle + self.angular_velocity * dt,
            angular_velocity: self.angular_velocity,
            radius: self.radius, // rotation preserves radius
        }
    }

    /// Period = 2π/ω. Time as circular distance.
    pub fn period(&self) -> f64 {
        if self.angular_velocity.abs() < 1e-15 {
            f64::INFINITY
        } else {
            2.0 * std::f64::consts::PI / self.angular_velocity
        }
    }

    /// Conservation error: difference in radii between two spins.
    pub fn conservation_error(&self, other: &Spin) -> f64 {
        (self.radius - other.radius).abs()
    }

    /// Kinetic energy = ½r²ω².
    pub fn energy(&self) -> f64 {
        0.5 * self.radius * self.radius * self.angular_velocity * self.angular_velocity
    }

    /// Current x-position in phase space.
    pub fn x(&self) -> f64 {
        self.radius * self.angle.cos()
    }

    /// Current y-position in phase space.
    pub fn y(&self) -> f64 {
        self.radius * self.angle.sin()
    }
}

/// A system of coupled spins, evolving under symplectic rotation.
#[derive(Debug, Clone)]
pub struct SpinSystem {
    pub spins: Vec<Spin>,
    /// Coupling matrix: coupling[i][j] = strength of interaction between spin i and j.
    pub coupling: Vec<Vec<f64>>,
}

impl SpinSystem {
    /// Create a spin system from an adjacency/weight matrix.
    /// Each eigenvalue of the graph Laplacian → one Spin.
    pub fn from_graph(adj: &[Vec<f64>]) -> Self {
        let _n = adj.len();
        // Compute degree for each node
        let degree: Vec<f64> = adj.iter().map(|row| row.iter().sum()).collect();

        // Compute eigenvalues of the Laplacian using power iteration for each
        let eigenvalues = graph_eigenvalues(adj, &degree);

        let spins: Vec<Spin> = eigenvalues.iter()
            .map(|&ev| Spin::from_eigenvalue(ev))
            .collect();

        SpinSystem {
            spins,
            coupling: adj.to_vec(),
        }
    }

    /// Create from pre-computed spins and coupling.
    pub fn new(spins: Vec<Spin>, coupling: Vec<Vec<f64>>) -> Self {
        SpinSystem { spins, coupling }
    }

    /// Total energy: Σ ½r²ω².
    pub fn total_energy(&self) -> f64 {
        self.spins.iter().map(|s| s.energy()).sum()
    }

    /// Advance all spins by dt with symplectic rotation.
    /// The coupling perturbs angular velocities slightly but preserves total energy.
    pub fn advance(&mut self, dt: f64) {
        // Compute coupling influence on each spin
        let n = self.spins.len();
        let mut coupling_torque = vec![0.0; n];

        for i in 0..n {
            for j in 0..n {
                if i != j && j < self.coupling[i].len() {
                    // Coupling torque proportional to phase difference * coupling strength
                    let phase_diff = self.spins[j].angle - self.spins[i].angle;
                    coupling_torque[i] += self.coupling[i][j] * phase_diff.sin() * 0.01;
                }
            }
        }

        // Advance each spin with modified velocity (symplectic: velocity updated first)
        for (i, spin) in self.spins.iter_mut().enumerate() {
            let modified_omega = spin.angular_velocity + coupling_torque[i];
            spin.angle += modified_omega * dt;
            // Radius unchanged — symplectic conservation
        }
    }

    /// Energy drift from initial value. Should stay near zero for symplectic integrator.
    pub fn energy_drift(&self, initial_energy: f64) -> f64 {
        (self.total_energy() - initial_energy).abs()
    }
}

/// Compute eigenvalues of a graph Laplacian using QR iteration (no deps!).
fn graph_eigenvalues(adj: &[Vec<f64>], degree: &[f64]) -> Vec<f64> {
    let n = adj.len();
    if n == 0 { return vec![]; }

    // Build Laplacian: L = D - A
    let mut lap = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                lap[i][j] = degree[i];
            } else {
                lap[i][j] = -adj.get(i).and_then(|row| row.get(j)).copied().unwrap_or(0.0);
            }
        }
    }

    // QR iteration to find eigenvalues
    let mut mat = lap;
    for _ in 0..100 {
        // Gram-Schmidt QR decomposition
        let (q, r) = qr_decompose(&mat);
        // A_{k+1} = R * Q
        mat = mat_mul(&r, &q);
    }

    // Eigenvalues are on the diagonal
    (0..n).map(|i| mat[i][i].max(0.0)).collect()
}

/// QR decomposition via Gram-Schmidt.
fn qr_decompose(a: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = a.len();
    let m = a[0].len();
    // Columns of A
    let mut q = vec![vec![0.0; m]; n];
    let mut r = vec![vec![0.0; n]; m];

    for j in 0..n {
        // v = a_j
        let mut v: Vec<f64> = (0..m).map(|i| a[i][j]).collect();

        for i in 0..j {
            // r[i][j] = q_i . a_j
            let dot: f64 = (0..m).map(|k| q[k][i] * a[k][j]).sum();
            r[i][j] = dot;
            for k in 0..m {
                v[k] -= dot * q[k][i];
            }
        }

        let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        r[j][j] = norm;
        if norm > 1e-15 {
            for k in 0..m {
                q[k][j] = v[k] / norm;
            }
        }
    }

    (q, r)
}

/// Matrix multiplication.
fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let m = b[0].len();
    let k = b.len();
    let mut c = vec![vec![0.0; m]; n];
    for i in 0..n {
        for j in 0..m {
            for l in 0..k {
                c[i][j] += a[i][l] * b[l][j];
            }
        }
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_spin_from_eigenvalue() {
        let s = Spin::from_eigenvalue(4.0);
        assert!((s.angular_velocity - 2.0).abs() < 1e-10);
        assert!((s.period() - PI).abs() < 1e-10);
    }

    #[test]
    fn test_spin_period() {
        // λ = 1 → ω = 1 → T = 2π
        let s = Spin::from_eigenvalue(1.0);
        assert!((s.period() - 2.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_spin_advance_preserves_radius() {
        let s = Spin::from_eigenvalue(4.0);
        let s2 = s.advance(0.1);
        assert!((s.radius - s2.radius).abs() < 1e-10);
        assert!(s2.angle > s.angle);
    }

    #[test]
    fn test_spin_phase_space_circle() {
        let s = Spin::new(0.0, 1.0, 2.0);
        assert!((s.x() - 2.0).abs() < 1e-10);
        assert!((s.y() - 0.0).abs() < 1e-10);

        let s2 = s.advance(PI / 2.0);
        assert!((s2.x() - 0.0).abs() < 1e-10);
        assert!((s2.y() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_spin_system_energy_conservation() {
        let spins = vec![
            Spin::new(0.0, 2.0, 1.0),
            Spin::new(1.0, 3.0, 1.0),
            Spin::new(2.0, 1.0, 1.0),
        ];
        let coupling = vec![
            vec![0.0, 0.5, 0.0],
            vec![0.5, 0.0, 0.5],
            vec![0.0, 0.5, 0.0],
        ];
        let mut sys = SpinSystem::new(spins, coupling);
        let initial = sys.total_energy();

        // Advance 1000 steps
        for _ in 0..1000 {
            sys.advance(0.001);
        }

        let drift = sys.energy_drift(initial);
        assert!(drift < 1e-6, "Energy drift {} too large", drift);
    }

    #[test]
    fn test_spin_system_from_graph() {
        // 3-node path graph: 0-1-2
        let adj = vec![
            vec![0.0, 1.0, 0.0],
            vec![1.0, 0.0, 1.0],
            vec![0.0, 1.0, 0.0],
        ];
        let sys = SpinSystem::from_graph(&adj);
        assert_eq!(sys.spins.len(), 3);
        // All eigenvalues should be non-negative
        for s in &sys.spins {
            assert!(s.angular_velocity >= 0.0);
        }
    }
}
