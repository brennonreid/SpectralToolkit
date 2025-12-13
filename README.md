# SpectralToolkit

This repository contains the **Spectral Toolkit** — a computational framework for exploring the geometry, stability, and failure modes of the Riemann zeta function through *raw, non-smoothed energy landscapes*.

The core idea is simple but deliberate:  
instead of smoothing, squaring, mollifying, or linearizing the zeta function, this toolkit works directly with

```
E(σ, T) = |ζ(σ + iT)|
```

and studies its **non-smooth geometry**, induced gradient flow, and constrained dynamics at high precision.

This project is computational and geometric in nature. It does **not** introduce new analytic identities for ζ(s). All mathematics used here is classical. What is new is the *way the function is interrogated*.

---

## What this toolkit does

The code in this repository supports:

- High-precision evaluation of ζ(s) using **GMP / MPFR / MPC / Arb**
- Direct sampling of the raw magnitude `|ζ(σ+iT)|` without smoothing
- Detection of:
  - cusp minima at nontrivial zeros
  - critical-line channeling (`σ → 1/2`)
  - rigid basin separation (including Lehmer pairs)
  - integer-barrier wall formation
  - catastrophic collapse under adversarial perturbation
- Deterministic, reproducible numerical experiments (no heuristic snapping)

The system treats the zeta magnitude as a **non-smooth energy functional**, not as a squared norm or log-potential.

---

## What this toolkit intentionally does *not* do

- No `|ζ|²`, no `log |ζ|`, no mollifiers
- No smoothing to make things “nicer”
- No projection onto known zero tables
- No artificial locking, snapping, or rounding toward expected results
- No claim of a formal proof of the Riemann Hypothesis

If smoothing makes your result look better, this toolkit considers that a failure mode — not a feature.

---

## Repository layout (high-level)

The repo is intentionally minimal and direct.

- **C++ sources**
  - High-precision zeta evaluation
  - Energy sampling and sweep logic
  - Cone / ridge / wall detection experiments
- **VS Code integration**
  - `.vscode/tasks.json` included for reproducible builds under MSYS2 UCRT64
- **Experiments**
  - Lehmer pair sweeps
  - Integer barrier tests
  - Adversarial sabotage runs
- **Artifacts**
  - Raw numeric output intended to be inspected, not post-processed away

Nothing here relies on editor magic or hidden build steps.

---

## Build environment (Windows)

This repository is built and tested under:

- **Windows**
- **MSYS2 (UCRT64)**
- **g++ (UCRT64 toolchain)**

Required libraries:

- GMP
- MPFR
- MPC
- FLINT
- Arb

Install via MSYS2:

```bash
pacman -S mingw-w64-ucrt-x86_64-gcc \
          mingw-w64-ucrt-x86_64-gmp \
          mingw-w64-ucrt-x86_64-mpfr \
          mingw-w64-ucrt-x86_64-mpc \
          mingw-w64-ucrt-x86_64-flint \
          mingw-w64-ucrt-x86_64-arb
```

---

## VS Code setup

This repo includes a `.vscode/` directory with working build tasks.

Important points:

- Builds are executed through **UCRT64**
- The active-file build task links against the shared zeta implementation
- The flags reflect exactly how the code is compiled — no abstraction layers

If it builds in VS Code, it builds on the command line the same way.

---

## Design philosophy

This project is guided by a few strict rules:

1. **Raw geometry over analytic convenience**  
   If a transformation removes cusp structure, it is rejected.

2. **Failure matters more than success**  
   Sabotage tests are first-class citizens. If the landscape doesn’t collapse under bad data, something is wrong.

3. **Computers are allowed to do what analysis avoids**  
   Non-smoothness is not an error condition. It is the signal.

4. **Reproducibility over polish**  
   Numeric output is preserved verbatim. Nothing is massaged to “look right.”

---

## Relationship to Hilbert–Pólya

This toolkit does not attempt to construct a linear self-adjoint operator.

Instead, it explores whether the **raw zeta magnitude itself already behaves like a spectral object**, once smoothing assumptions are dropped:

- zeros appear as ground states
- the critical line acts as an attractor
- Lehmer pairs remain separated
- constrained domains form walls
- adversarial perturbations cause topological collapse

Whether this ultimately supports or undermines Hilbert–Pólya is left open — but the computational evidence is explicit and inspectable.

---

## Status

Active research code.

This is not a library.  
This is not a finished theory.  
This is not a toy.

It is a working computational apparatus designed to answer one question honestly:

> *What does the Riemann zeta function actually do when you stop smoothing it?*

---

## Author

**Brennon Reid**

All code, experiments, and numerical results in this repository are my own unless explicitly stated otherwise.
