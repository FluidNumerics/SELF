# CLAUDE.md
Claude Code Operating Rules for This Repository

This file defines **hard constraints** and **preferred behaviors** for Claude
working in this repository. These rules exist to protect numerical correctness,
performance, and long-term maintainability of a high-order spectral element
codebase.

Follow all rules below unless explicitly instructed otherwise by a human
maintainer.

---

## 1. Project Purpose

This project implements **high-order spectral element methods** for solving
systems of conservation laws (e.g., Euler, shallow water, Maxwell, MHD) on
structured and unstructured meshes.

Primary goals:
- Numerical correctness and stability
- Conservation properties
- High performance on CPUs and GPUs
- MPI scalability

Changes must NOT degrade:
- Accuracy order
- Stability properties
- Parallel scaling
- Memory behavior

---

## 2. Fortran Standards & Toolchain

### Language Standard
- Target: **Fortran 2008**
- Code must remain compatible with:
  - gfortran ≥ 11
  - ifx
  - nvfortran
  - amdflang

### Prohibited Language Features
Do NOT introduce:
- Coarrays
- Fortran 2018+ features
- Compiler-specific extensions
- Automatic polymorphism in performance-critical paths

### Formatting & Conventions
- Free-form source
- `implicit none` required in all program units
- Explicit `intent(in|out|inout)` on all dummy arguments
- Lowercase keywords preferred
- Line length ≤ 100 characters

---

## 3. Numerical & Algorithmic Constraints

### Hard Rules
- Do NOT change the mathematical formulation without approval
- Do NOT change discretization order
- Do NOT change basis, quadrature, or nodal ordering
- Do NOT reorder floating-point reductions
- Do NOT alter time integration schemes

### Floating-Point Behavior
- Bitwise reproducibility may be required
- Preserve operation ordering in loops
- Avoid algebraic “simplifications” unless mathematically justified

### Array Semantics
- Do NOT replace explicit loops with array syntax unless equivalence is proven
- Avoid implicit temporaries

---

## 4. Performance Rules (Critical)

This is an HPC codebase. Performance regressions are unacceptable.

### Memory
- Avoid temporary allocations in hot paths
- No automatic arrays in tight loops
- No hidden allocations via array slicing

### Loops
- Preserve loop ordering for cache locality
- Do NOT replace DO loops with WHERE / FORALL
- Vectorization-friendly structure must be preserved

### Abstraction
- Do NOT introduce runtime polymorphism in kernels
- Avoid excessive modularization inside hot loops

---

## 5. Parallel Programming

### MPI
- MPI calls must remain explicit
- Do NOT introduce blocking collectives inside time-stepping loops
- Do NOT change communicator usage
- Preserve rank-local data ownership

### OpenMP / GPU
- Preserve OpenMP semantics
- Do NOT move data regions without explicit instruction
- GPU kernels must preserve memory access patterns
- No implicit host/device transfers

---

## 6. Code Organization & APIs

### File & Module Structure
- Do NOT rename modules
- Do NOT move files between directories
- Do NOT change public interfaces without approval
- Preserve module dependency order

### Public APIs
- Public procedures are considered **stable**
- Backward compatibility is required unless stated otherwise

---

## 7. Testing & Validation

### Required
- All existing regression tests must pass
- Do NOT modify reference output files
- Numerical differences must be justified

### New Code
- New features require:
  - A test case
  - Clear validation criteria
- MPI tests must work on ≥ 2 ranks

---

## 8. Documentation & Comments

### Preserve Scientific Meaning
- Do NOT remove comments describing:
  - Equations
  - Algorithms
  - Numerical assumptions

### New Routines
Must include:
- Mathematical description
- Variable meaning and units
- Expected input ranges

---

## 9. Code Formatting

All Fortran source must be formatted with [`fprettify`](https://pypi.org/project/fprettify/) using the project's `fprettify.config`. PRs are checked for formatting before any other tests run.

To format all source files manually:

```shell
fprettify './src/'      --config-file ./fprettify.config --recursive --case 1 1 1 1
fprettify './test/'     --config-file ./fprettify.config --recursive --case 1 1 1 1
fprettify './examples/' --config-file ./fprettify.config --recursive --case 1 1 1 1
```

Alternatively, install the provided `pre-commit` hook to apply formatting automatically on each commit:

```shell
pip install pre-commit fprettify
pre-commit install   # run from repository root
```

When editing Fortran files, apply `fprettify` before committing. Do NOT manually reformat code by hand in ways that deviate from `fprettify` output.

---

## 10. Prohibited Actions (Explicit)

Do NOT:
- Rewrite code in another language
- Convert procedural code to OO Fortran
- Replace MPI with coarrays
- Introduce external dependencies
- “Modernize” syntax without benefit
- Delete legacy code without explanation

---

## 11. Domain-Specific Assumptions

- Grid indexing follows project conventions (do NOT reorder indices)
- Jacobians and metric terms are precomputed
- Flux routines assume nodal basis ordering
- Element-local operations must remain element-local
- Halo exchange patterns are fixed

---

## 12. Preferred Behavior

DO:
- Ask before changing algorithms
- Explain numerical and performance implications
- Provide minimal diffs
- Reference existing patterns in the codebase
- Flag any uncertainty explicitly

DO NOT:
- Make large refactors unless requested
- Assume intent beyond the explicit request

---

## 13. When in Doubt

If a change could affect:
- Numerical accuracy
- Stability
- Performance
- Parallel behavior

STOP and ask for clarification.

---

End of CLAUDE.md

