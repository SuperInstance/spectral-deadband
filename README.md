# spectral-deadband

**Deadband is the spectral gap. Spin is time as distance. The dial IS the computation.**

Pure Rust, zero dependencies. A library that reframes spectral graph theory as physical analog computation — because that's what it already is, if you squint right.

## The Core Idea

A thermostat doesn't act when the temperature is "close enough." That region of inaction — the deadband — is the mercury band in an old-school thermostat. Here's the claim: **that deadband IS the spectral gap** (the largest gap between consecutive eigenvalues).

This isn't a metaphor stretched thin. The spectral gap determines how quickly a system settles, how clearly separated its modes are, and whether it can afford to do nothing. The mercury band does the same job by physics instead of by algorithm.

From there, everything follows:
- **Eigenvalues → angular frequencies** (ω = √λ, so time is rotation)
- **Symplectic rotation → conservation** (radius preserved under rotation)
- **Fractal structure → self-similar conservation ratio** (same CR at every scale)
- **Analog dials → eigenvectors** (gravity settles to setpoint, deadband = spectral gap)

```toml
[dependencies]
spectral-deadband = "0.1.0"
```

## Module Walkthrough

### Deadband — The Spectral Gap as Region of Inaction

The thermostat computes eigenvalue gaps. Given a set of eigenvalues, it finds the largest gap between consecutive ones. That gap defines a deadband: a region where the system should not act.

```rust
use spectral_deadband::{Deadband, spectral_deadband, Action};

// A physical deadband: center at 5.0, width 2.0 (band: [4.0, 6.0])
let db = Deadband::new(5.0, 2.0);

// The thermostat evaluates a measurement
assert_eq!(db.evaluate(3.0), Some(Action::Heat));   // too cold
assert_eq!(db.evaluate(5.0), None);                  // in deadband → do nothing
assert_eq!(db.evaluate(7.0), Some(Action::Cool));   // too hot

// Or compute the deadband directly from eigenvalues
let evals = vec![1.0, 2.0, 3.0, 7.0, 8.0];
let spectral_db = spectral_deadband(&evals);
// Largest gap: 3.0 → 7.0, width = 4.0, center = 5.0
assert!((spectral_db.center() - 5.0).abs() < 1e-10);
assert!((spectral_db.width() - 4.0).abs() < 1e-10);
```

The `gravity` field is set proportional to the eigenvalue magnitudes flanking the gap — because the restoring force depends on how strong the modes are, not just how far apart.

### Spin — Time Abstracted as Rotation

An eigenvalue λ maps to angular frequency ω = √λ. A time step is a rotation in phase space. Rotations preserve the radius, which means **conservation is structural**, not something you have to enforce.

```rust
use spectral_deadband::Spin;

// Create a spin from eigenvalue λ = 4.0
let s = Spin::from_eigenvalue(4.0);
// ω = √4 = 2.0, period T = 2π/2 = π
assert!((s.angular_velocity - 2.0).abs() < 1e-10);
assert!((s.period() - std::f64::consts::PI).abs() < 1e-10);

// Advance by dt: pure rotation, radius preserved
let s2 = s.advance(0.1);
assert!((s.radius - s2.radius).abs() < 1e-10); // conservation

// Phase space coordinates
let s = Spin::new(0.0, 1.0, 2.0);
assert!((s.x() - 2.0).abs() < 1e-10); // at angle 0: x = r, y = 0
```

The `SpinSystem` couples multiple spins via an adjacency matrix, computing eigenvalues of the graph Laplacian using QR iteration (implemented from scratch, no deps). Coupling perturbs angular velocities via phase differences, but the symplectic structure keeps total energy drift below 1e-6 over 1000 steps.

### FractalScale — Conservation Ratio at Every Scale

The conservation ratio (CR) is a normalized energy measure: Σ(eigenvalues) / Σ(degrees). For certain graph structures, CR is self-similar — the same whether you measure it on 5 nodes or 500.

```rust
use spectral_deadband::{FractalScale, complete_graph, fibonacci_graph};

// Measure CR at multiple scales on a complete graph
let adj = complete_graph(10);
let fs = FractalScale::measure(&adj, 5);
assert!(fs.is_self_similar(0.5)); // uniform structure → self-similar

// Fibonacci outward: inflate the graph by Penrose-like subdivision
// Each edge splits into two with golden-ratio weights
let inflated = FractalScale::inflate(&adj);
assert_eq!(inflated.len(), 20); // doubles the node count

// Mandelbrot inward: measure roughness (CR variance across subgraphs)
let roughness = FractalScale::roughness(&adj, 5);
// Complete graphs have uniform subgraphs → low roughness
assert!(roughness < 1.0);

// Check deviation from the golden ratio 1/φ
let deviation = fs.golden_deviation();
```

The inflation scheme uses Fibonacci weights (φ and 1/φ) on subdivided edges. The idea is that repeated inflation should drive CR toward 1/φ ≈ 0.618, mirroring how Penrose tilings have golden-ratio proportions at every scale. In practice, convergence depends heavily on the initial graph structure — this is a pattern, not a theorem.

### AnalogDial — Gravity Is the Eigenvalue Computation

A dial sits at some position. Gravity pulls it toward a setpoint. Friction opposes motion. Within the deadband, friction wins and the dial stops. Outside, gravity wins and the dial moves.

**The settling IS the computation.** The final position approximates the eigenvalue.

```rust
use spectral_deadband::{AnalogDial, DialBank};

// A single dial settling under gravity
let mut dial = AnalogDial::new(5.0, 0.1, 2.0)
    .at_position(10.0)
    .with_friction(0.3);

for _ in 0..5000 {
    dial.settle(0.01);
}
assert!(dial.is_settled());

// A bank of coupled dials → approximates eigenvectors
let coupling = vec![
    vec![0.0, 1.0, 0.0],
    vec![1.0, 0.0, 1.0],
    vec![0.0, 1.0, 0.0],
];
let mut bank = DialBank::new(5.0, 0.5, 1.0, coupling);
let positions = bank.settle_all(0.01, 10000);

// Read eigenvalue estimate via Rayleigh quotient
let eigenvalue = bank.eigenvalue_estimate();
```

## Why Gravity IS the Eigenvalue Computation

This is the part that's either obvious or wild, depending on your background.

A spring-mass-damper system follows:

```
m·a = -k·x - c·v
```

This is also the power iteration for finding eigenvectors. The displacement `x` converges to the dominant eigenvector direction, and the ratio `k/m` determines the eigenvalue. Add friction (Coulomb friction, specifically), and you get a deadband — a region where the restoring force is weaker than static friction, so nothing moves.

In this library:
- `gravity` = spring constant k (restoring force)
- `friction` = damping + Coulomb deadband
- `position` = current eigenvalue estimate
- `setpoint` = target eigenvalue

The dial doesn't *simulate* eigenvalue computation. It *is* eigenvalue computation, expressed in physical coordinates instead of linear algebra coordinates.

## Connection to Real Systems

**Thermostats:** A bimetallic strip or mercury switch has a deadband — it doesn't flip at exactly 72°F. It flips at 71°F and 73°F. That 2°F band is the spectral gap in action: the system has separated into "heating" and "cooling" modes, and the gap between them determines how much the temperature can drift before the system reacts.

**Centrifugal governors:** Watt's governor uses spinning weights on a shaft. As RPM increases, centrifugal force moves the weights outward, which closes a steam valve. The equilibrium position of the weights IS the eigenvalue of the coupled mechanical system. The deadband is the friction in the linkage.

**PID controllers:** The proportional band in a PID controller is literally a deadband — the error has to exceed some threshold before the controller acts. The spectral gap is the distance between "close enough" and "needs correction."

## Honest Limitations

This is a conceptual library, not a production eigenvalue solver. Be honest about what it can't do:

1. **Analog precision floor.** The dials settle to within their deadband, but no further. A deadband of 0.1 means ~3 bits of precision. You can tighten it by increasing gravity or decreasing friction, but then settling takes longer. There's a speed-precision tradeoff that real analog computers also face.

2. **Friction models are simplified.** Real Coulomb friction has stick-slip, temperature dependence, and aging effects. This library uses a constant friction coefficient. Real analog computers (op-amp integrators, for instance) have drift, offset, and noise that these models don't capture.

3. **QR iteration is naive.** The eigenvalue computation in `Spin` uses basic QR iteration with 100 fixed iterations, no shifting, no deflation. For production eigenvalue work, use `nalgebra` or `rust-ndarray` with LAPACK bindings.

4. **Coupling in DialBank is ad-hoc.** The coupling forces are scaled by 0.1 and 0.01 arbitrarily. Convergence depends on these magic numbers. A real analog eigenvalue solver would need careful impedance matching.

5. **No temperature dependence.** The precision analysis in the companion library (`analog-spectral`) covers this, but this library treats gravity and friction as constants. In reality, thermal expansion changes both.

## Running Tests

```bash
cargo test
```

All tests pass with no dependencies. The test suite covers:
- Deadband containment and thermostat evaluation
- Spin creation, period calculation, energy conservation
- FractalScale self-similarity, inflation, roughness
- AnalogDial settling and DialBank convergence
- Spectral deadband from eigenvalue arrays

## License

MIT
