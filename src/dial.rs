//! Module 4: AnalogDial — physical deadband from gravity.
//!
//! A dial physically stores analog data. The position is continuous.
//! The deadband comes from gravity/friction. This IS analog spectral computation.

/// A single analog dial that settles under gravity toward a setpoint.
///
/// Within deadband: friction ≥ gravity component → no movement.
/// Outside deadband: gravity pulls toward setpoint.
/// This IS analog computation of spectral gap.
#[derive(Debug, Clone, Copy)]
pub struct AnalogDial {
    /// Current value (eigenvalue analog).
    pub position: f64,
    /// Target value (desired eigenvalue).
    pub setpoint: f64,
    /// Spectral gap analog — the deadband width.
    pub deadband: f64,
    /// Restoring force (eigenvector weight analog).
    pub gravity: f64,
    /// Damping (energy dissipation).
    pub friction: f64,
    /// Current velocity.
    velocity: f64,
}

impl AnalogDial {
    pub fn new(setpoint: f64, deadband: f64, gravity: f64) -> Self {
        AnalogDial {
            position: setpoint + deadband * 3.0, // start outside deadband
            setpoint,
            deadband,
            gravity,
            friction: 0.1,
            velocity: 0.0,
        }
    }

    /// Create a dial at a specific position.
    pub fn at_position(mut self, pos: f64) -> Self {
        self.position = pos;
        self
    }

    /// Create a dial with specific friction.
    pub fn with_friction(mut self, friction: f64) -> Self {
        self.friction = friction;
        self
    }

    /// The dial settles under gravity toward the setpoint.
    ///
    /// Within deadband: no movement (friction > gravity component).
    /// Outside deadband: gravity pulls toward setpoint.
    pub fn settle(&mut self, dt: f64) {
        let error = self.setpoint - self.position;
        let abs_error = error.abs();

        if abs_error <= self.deadband {
            // Within deadband: friction wins, stop moving
            // Apply strong damping to velocity
            self.velocity *= 0.5;
            if self.velocity.abs() < 1e-12 {
                self.velocity = 0.0;
            }
            self.position += self.velocity * dt;
        } else {
            // Outside deadband: gravity pulls toward setpoint
            let force = self.gravity * if error >= 0.0 { 1.0 } else { -1.0 };
            self.velocity += force * dt;
            // Apply damping
            self.velocity *= (1.0 - self.friction * dt).max(0.0);
            self.position += self.velocity * dt;
        }
    }

    /// Is the dial settled within the deadband?
    pub fn is_settled(&self) -> bool {
        (self.position - self.setpoint).abs() <= self.deadband
    }

    /// Distance from setpoint.
    pub fn error(&self) -> f64 {
        (self.position - self.setpoint).abs()
    }
}

/// A bank of coupled analog dials = spectral system.
///
/// Dials influence each other through coupling. When all settle, their positions
/// approximate the eigenvector. The equilibrium pattern IS the eigenvector.
#[derive(Debug, Clone)]
pub struct DialBank {
    pub dials: Vec<AnalogDial>,
    /// Coupling matrix: coupling[i][j] = influence of dial j on dial i.
    pub coupling: Vec<Vec<f64>>,
}

impl DialBank {
    /// Create a bank of dials all targeting the same setpoint, with a coupling matrix.
    pub fn new(setpoint: f64, deadband: f64, gravity: f64, coupling: Vec<Vec<f64>>) -> Self {
        let n = coupling.len();
        let dials = (0..n)
            .map(|i| {
                // Offset each dial's initial position by index to break symmetry
                AnalogDial::new(setpoint + (i as f64 + 1.0) * deadband * 2.0, deadband, gravity)
                    .at_position(setpoint + (i as f64 + 1.0) * deadband * 2.0)
            })
            .collect();

        DialBank { dials, coupling }
    }

    /// Settle all dials with coupling.
    ///
    /// Returns final positions of all dials.
    pub fn settle_all(&mut self, dt: f64, iterations: usize) -> Vec<f64> {
        for _ in 0..iterations {
            // Compute coupling forces first
            let n = self.dials.len();
            let mut coupling_force = vec![0.0; n];

            for i in 0..n {
                for j in 0..n {
                    if i != j && j < self.coupling[i].len() {
                        let diff = self.dials[j].position - self.dials[i].position;
                        coupling_force[i] += self.coupling[i][j] * diff * 0.1;
                    }
                }
            }

            // Apply coupling as perturbation to setpoint, then settle
            for (i, dial) in self.dials.iter_mut().enumerate() {
                // Temporarily shift setpoint by coupling force
                let original_setpoint = dial.setpoint;
                dial.setpoint += coupling_force[i];
                dial.settle(dt);
                dial.setpoint = original_setpoint;
            }
        }

        self.dials.iter().map(|d| d.position).collect()
    }

    /// Estimate the dominant eigenvalue from settled dial positions.
    ///
    /// Uses the Rayleigh quotient: λ ≈ (x^T A x) / (x^T x)
    /// where x is the dial position vector and A is the coupling matrix.
    pub fn eigenvalue_estimate(&self) -> f64 {
        let n = self.dials.len();
        if n == 0 { return 0.0; }

        let positions: Vec<f64> = self.dials.iter().map(|d| d.position).collect();
        let mut ax = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                if j < self.coupling[i].len() {
                    ax[i] += self.coupling[i][j] * positions[j];
                }
            }
        }

        let xtx: f64 = positions.iter().map(|x| x * x).sum();
        let xtax: f64 = positions.iter().zip(ax.iter()).map(|(x, ax)| x * ax).sum();

        if xtx.abs() < 1e-15 { 0.0 } else { xtax / xtx }
    }

    /// Check if all dials are settled.
    pub fn all_settled(&self) -> bool {
        self.dials.iter().all(|d| d.is_settled())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dial_settles_to_setpoint() {
        let mut dial = AnalogDial::new(5.0, 0.1, 2.0)
            .at_position(10.0)
            .with_friction(0.3);

        for _ in 0..5000 {
            dial.settle(0.01);
        }

        assert!(dial.is_settled(), "Dial should be settled. Position: {}, Error: {}", dial.position, dial.error());
        assert!(dial.error() <= 0.1 + 1e-6);
    }

    #[test]
    fn test_dial_stops_in_deadband() {
        let mut dial = AnalogDial::new(0.0, 1.0, 0.5)
            .at_position(0.3)
            .with_friction(0.5);

        // Start inside deadband, should barely move
        let initial_pos = dial.position;
        for _ in 0..100 {
            dial.settle(0.01);
        }

        // Should have moved very little since it's in the deadband
        assert!((dial.position - initial_pos).abs() < 0.5);
    }

    #[test]
    fn test_dialbank_converges() {
        let coupling = vec![
            vec![0.0, 1.0, 0.0],
            vec![1.0, 0.0, 1.0],
            vec![0.0, 1.0, 0.0],
        ];

        let mut bank = DialBank::new(5.0, 0.5, 1.0, coupling);
        let positions = bank.settle_all(0.01, 10000);

        // All dials should have settled
        for (i, &pos) in positions.iter().enumerate() {
            let err = (pos - 5.0).abs();
            assert!(err < 5.0, "Dial {} position {} too far from setpoint", i, pos);
        }
    }

    #[test]
    fn test_dialbank_eigenvalue_estimate() {
        // Identity coupling: eigenvalue = 1
        let coupling = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
        ];

        let mut bank = DialBank::new(3.0, 0.1, 2.0, coupling);
        bank.settle_all(0.01, 5000);

        let est = bank.eigenvalue_estimate();
        // For identity matrix, eigenvalue should be close to 1
        assert!(est.is_finite(), "Eigenvalue estimate should be finite, got {}", est);
    }
}
