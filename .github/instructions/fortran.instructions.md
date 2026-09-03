---
applyTo: "src/**/*.f90,examples/**/*.f90"
---

# Fortran source conventions

The reference for this style is the tree at `f3e1e57c`, the last commit before
the first AI assisted contribution. `docs/Contributing/StyleGuide.md` describes
each rule and names the check that enforces it.

## Comments

- Comments use `!`. FORD documentation uses the post mark `!!`, placed on the
  line **after** the module statement, procedure signature, or declaration it
  describes. The pre mark `!>` is not used in this project.
- Write comments as plain prose sentences. No markdown headings, no `**bold**`,
  no fenced code blocks inside a comment.
- ASCII punctuation only. No em dashes, no typographic quotes, no emoji.
  Superscripts and the dot operator in unit expressions such as `m s⁻¹` are part
  of the existing style and are permitted.
- Say what a routine computes and why, not what the next statement does. A
  trailing comment naming the quantity, as in
  `flux(1,1) = this%rho0*s(2) ! density, x flux ; rho0*u`, is the house style.
- Do not leave `TODO` or `FIXME` in source. Open an issue.
- Do not delete a comment that describes an equation, an algorithm, or a
  numerical assumption.

## Layout

- Free form, two space indent, lines no longer than 132 characters.
- `implicit none` at module scope and again in every procedure, placed after the
  docstring and before the dummy argument declarations.
- Every dummy argument carries an explicit `intent(in|out|inout)`.
- Closing keywords are fused: `endmodule`, `endsubroutine Foo`, `endfunction Foo`,
  `endtype Foo`, `enddo`, `endif`, `endinterface`.
- Declaration order inside a procedure: `implicit none`, then
  `class(X),intent(inout) :: this`, then the remaining dummy arguments in
  argument order, then the `! Local` marker, then locals.
- Every file opens with the 25 line BSD-3 banner. Copy it from
  `.github/license-header.txt`, not from a neighbouring source file: eighteen
  files in the tree carry a truncated or corrupted banner, and F001 rejects any
  copy of one.
- Run `fprettify` with `fprettify.config` before committing; do not hand format.

## Naming

- Files are `src/SELF_<Thing>.f90`. The `_t` suffix marks the portable template
  layer; backend specializations live in `src/cpu/`, `src/gpu/`, and `src/apu/`
  under the same name without `_t`.
- Types are CamelCase with a dimension suffix: `Lagrange_t`, `MappedScalar2D_t`,
  `DGModel2D_t`.
- Module procedures are `<Verb><Noun>_<TypeName>`, bound to the short CamelCase
  name: `procedure,public :: Init => Init_Lagrange_t`.
- The passed object is always `this`. Loop indices are `i`, `j`, `k`, `iel`,
  `ivar`, `iside`.
- Every real literal carries `_prec`.
- Physical quantities keep their mathematical names: `H`, `g`, `f0`, `beta`,
  `rho0`, `c`.

## Numerics and performance

- Do not change the mathematical formulation, discretization order, basis,
  quadrature, nodal ordering, or time integration scheme.
- Preserve floating-point operation ordering. Do not reorder a reduction.
- Preserve loop ordering and structure. Do not substitute array syntax, `where`,
  or `forall` for an explicit loop.
- No temporary allocations, automatic arrays, or hidden allocations from array
  slicing in a hot path.
- No runtime polymorphism inside a kernel.
- Element-local operations stay element-local, and the halo exchange pattern is
  fixed.

## Documenting a new procedure

A new routine states what it computes, what each argument means and in what
units, and the range of inputs it is valid for. Mathematics goes in the `!!`
docstring as LaTeX, following `src/SELF_LinearEuler2D_t.f90`. Read
`src/SELF_Lagrange_t.f90` for the argument level documentation style.

Do not imitate `src/SELF_SupportRoutines.f90`; its doxygen markup is a holdover
from an earlier convention.
