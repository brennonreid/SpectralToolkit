SpectralToolkit
Spectral Toolkit is a computational framework for analyzing the geometry, stability, and failure modes of the Riemann zeta function through raw, non-smoothed energy landscapes.

The toolkit discards smoothing, squaring, mollifying, and linearization to work directly with the raw magnitude:

E(σ, T) = |ζ(σ + iT)|
It targets the non-smooth geometry, induced gradient flow, and constrained dynamics of the function at high precision. This is a direct interrogation of the zeta function's magnitude as a physical terrain.

Capabilities
High-Precision Evaluation: Full GMP / MPFR / MPC / Arb implementation for arbitrary precision.

Direct Magnitude Sampling: Operates on |ζ(σ+iT)| without phase correction or rotation.

Geometric Detection:

Cusp minima at nontrivial zeros

Critical-line channeling (σ → 1/2)

Rigid basin separation (Lehmer pairs)

Integer-barrier wall formation

Catastrophic collapse under adversarial perturbation

Raw Signal Processing:

No |ζ|², no log |ζ|, no mollifiers.

No smoothing, artificial locking, snapping, or rounding toward expected results.

The system treats the zeta magnitude as a non-smooth energy functional.

Repository Layout
C++ sources: High-precision zeta evaluation, energy sampling, and sweep logic.

Experiments: Lehmer pair sweeps, integer barrier tests, and cone/ridge detection.

Sabotage: Adversarial stress-tests to verify landscape stability.

Artifacts: Raw numeric output for inspection.

VS Code: .vscode/tasks.json included for reproducible builds under MSYS2 UCRT64.

Build Environment (Windows)
Platform: Windows / MSYS2 (UCRT64) / g++ Dependencies: GMP, MPFR, MPC, FLINT, Arb

Installation:

Bash

pacman -S mingw-w64-ucrt-x86_64-gcc \
          mingw-w64-ucrt-x86_64-gmp \
          mingw-w64-ucrt-x86_64-mpfr \
          mingw-w64-ucrt-x86_64-mpc \
          mingw-w64-ucrt-x86_64-flint \
          mingw-w64-ucrt-x86_64-arb
VS Code Setup
The repository includes a configured .vscode/ directory.

Builds execute through UCRT64.

Tasks link directly against the shared zeta implementation.

Compilation flags are explicit; no abstraction layers.

Design Philosophy
Raw Geometry: Any transformation that removes cusp structure is rejected. The raw landscape is the ground truth.

Sabotage Verification: Failure is the primary metric. The landscape must survive adversarial data input.

Signal in Non-Smoothness: Computers can handle what analysis avoids. The non-smooth cusps are the signal, not an error condition.

Absolute Reproducibility: Numeric output is preserved verbatim.

Operational Evidence for Hilbert–Pólya
This toolkit demonstrates that the raw zeta magnitude behaves as a spectral object when smoothing assumptions are dropped. The computational evidence shows:

Zeros appear as ground states.

The critical line acts as an attractor.

Lehmer pairs remain separated by rigid geometric barriers.

Constrained domains form walls.

Adversarial perturbations cause topological collapse.

The software provides explicit, inspectable evidence of these spectral properties in the raw energy landscape.

Author
Brennon Reid

All code, experiments, and numerical results in this repository are my own.