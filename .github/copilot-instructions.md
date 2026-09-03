# Copilot instructions for SELF

SELF implements high-order spectral element methods for systems of conservation
laws. Numerical correctness, conservation, and performance on CPUs and GPUs are
the properties the project protects; everything below follows from that.

`CLAUDE.md` holds the full set of hard constraints and
`docs/Contributing/StyleGuide.md` holds the coding and documentation style.
Cite those rather than restating them.

## Stay inside the pull request

Review is scoped to the change in front of you.

- Comment only on lines this pull request adds or modifies. Do not comment on
  unchanged lines, on files absent from the diff, or on conditions that already
  existed on `main`.
- Read the pull request description before reviewing. It has four required
  sections, and two of them bound your review:
  - **Scope** states what the change is meant to do. A comment asking for
    something outside it is a feature request, not review feedback. Say so in
    one line if it matters, and do not block on it.
  - **Out of scope** lists problems the author already found and deliberately
    deferred. Do not raise them. They are known.
- Do not ask for unrelated cleanup, opportunistic refactoring, or improvements to
  neighbouring code that the author did not touch.
- One comment per defect. Do not restate what the diff already shows.

## Leave these to the automated checks

Commenting on any of the following is noise, because a job already reports it:

- Whitespace, indentation, alignment, keyword case, and line breaking. These are
  enforced by `fprettify-lint` against `fprettify.config`.
- Comment characters, markdown inside comments, emoji, em dashes, missing
  `implicit none`, and documentation prose shape. These are enforced by
  `style-check` through `.github/scripts/style_check.py` and
  `.github/scripts/docs_style_check.py`.
- Whether the pull request description is complete. That is `pr-template-check`.

## Do not propose these changes

Each of these is prohibited by `CLAUDE.md`, so suggesting one costs the author
time and cannot be accepted:

- Renaming a module, derived type, or public procedure, or moving a file between
  directories. Public interfaces are stable and backward compatibility is
  required.
- Rewriting an explicit `do` loop as array syntax, `where`, or `forall`, or
  reordering a loop nest. Loop ordering carries cache locality and vectorization
  behaviour that the code depends on.
- Reordering or otherwise rearranging a floating-point reduction, or applying an
  algebraic simplification. Bitwise reproducibility may be required, so operation
  order is part of the specification, not an implementation detail.
- Changing the mathematical formulation, discretization order, basis, quadrature,
  nodal ordering, or time integration scheme.
- Introducing runtime polymorphism, extra layers of abstraction, or automatic
  arrays inside a kernel or other hot path.
- Replacing MPI calls with coarrays, changing communicator usage, altering the
  halo exchange pattern, or adding a blocking collective inside a time stepping
  loop.
- Moving an OpenMP or GPU data region, or introducing an implicit host to device
  transfer.
- Modernizing to Fortran 2018 or later, or using a compiler specific extension.
  The target is Fortran 2008 across gfortran, ifx, nvfortran, and amdflang.
- Adding an external dependency.

## What is worth a comment

- A defect in the numerics: a wrong index, a metric term applied on the wrong
  face, a flux that does not match the stated formulation, a missing `_prec` on a
  real literal.
- A correctness problem in parallel code: a missing halo exchange, a race, an
  incorrect rank-local assumption.
- An allocation, temporary, or copy introduced into a hot path.
- A missing `intent` on a dummy argument.
- New behaviour with no accompanying test, or a bug fix whose test would still
  pass against the old code.
- A public interface change that the description does not mention.

When a change affects accuracy, stability, performance, or parallel behaviour and
you are unsure, say what you are unsure about and ask. Do not assert a numerical
claim you have not checked against the surrounding code.
