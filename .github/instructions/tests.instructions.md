---
applyTo: "test/**/*.f90"
---

# Test conventions

The Fortran source conventions in `.github/instructions/fortran.instructions.md`
apply here too. What follows is specific to the test suite.

- New behaviour requires a test, and the test states its validation criteria: the
  quantity measured and the tolerance it must meet.
- A bug fix requires a test that fails against the unfixed code. A bug that
  reached `main` is evidence that the suite had a gap, and closing that gap is
  part of the fix. Say in the pull request description why the existing tests did
  not catch it.
- A test asserts a numerical property, not an implementation detail. Prefer a
  convergence rate, a conserved quantity, an entropy budget, or an exact solution
  comparison over a golden value.
- Reference output files are not modified. A numerical difference against a
  reference is a result to explain, not a file to update.
- An MPI test must work on two or more ranks, and must not assume a particular
  decomposition.
- A test exercises one thing and names it. `test/` file names describe the case
  being validated, as in
  `test/ec_advection_2d_entropy_conservation.f90`.
- Keep the runtime short. The suite runs on every pull request across four
  compilers.
