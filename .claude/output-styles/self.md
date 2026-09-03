---
name: SELF
description: Write Fortran, comments, and documentation in the voice of the original SELF source
---

You are working in SELF, a high-order spectral element code for systems of
conservation laws. Write as the original authors of this code write.

The reference for this voice is the tree at commit `f3e1e57c`, the last commit
before the first AI assisted contribution. When you are unsure how something
should read, open one of the exemplar files listed below and match it. The hard
constraints in `CLAUDE.md` and the rule descriptions in
`docs/Contributing/StyleGuide.md` take precedence over anything here.

## Prose, everywhere

This applies to comments, documentation, commit messages, and pull request
descriptions alike.

Write plain declarative sentences in paragraphs. State what a thing computes and
why it is done that way. Do not write a summary of what you just wrote, do not
open with "Note that" or "It is important to", and do not close with a
restatement.

Use ASCII punctuation. No em dashes, no typographic quotes, no ellipsis
character, no emoji. Where a sentence wants an em dash, use a comma, a colon, or
two sentences.

Do not structure explanation as bullets, bold labels, or headings when a
paragraph will carry it. Emphasis is rare in this codebase; the sentence carries
the emphasis. Reserve lists for genuine enumerations such as a parameter list
with units.

Mathematical notation is welcome and expected. In documentation and `!!`
docstrings, write it as LaTeX. In body comments, write it in the ASCII form the
source already uses, as in `du/dt = -(2/J) sum_n D_split[n,i] * ...`.

## Fortran

Target Fortran 2008, compiling under gfortran, ifx, nvfortran, and amdflang.

Comments use `!`. FORD docstrings use the post mark `!!` on the line **after**
the module statement, procedure signature, or declaration being documented. The
pre mark `!>` is not used here. Never put markdown inside a comment.

Every file opens with the 25 line BSD-3 license banner. Copy it from
`.github/license-header.txt`, which is the canonical text. Do not copy it from a
neighbouring source file: eighteen files in the tree carry a banner that has been
truncated or corrupted, and copying one of those propagates the damage.

`implicit none` goes at module scope and again in every procedure, after the
docstring and before the dummy declarations. Every dummy argument carries an
explicit `intent`. Declarations run in argument order, then the `! Local` marker,
then locals.

Closing keywords are fused: `endmodule`, `endsubroutine Foo`, `endfunction Foo`,
`endtype Foo`, `enddo`, `endif`.

Two space indent, 132 column limit. Run `fprettify` with `fprettify.config`
rather than formatting by hand.

Name module procedures `<Verb><Noun>_<TypeName>` and bind them to the short
CamelCase name. The passed object is `this`. Loop indices are `i`, `j`, `k`,
`iel`, `ivar`, `iside`. Every real literal carries `_prec`. Physical quantities
keep their mathematical names.

A new routine documents what it computes, what each argument means and in what
units, and the valid input range.

## What not to do to this code

These are hard constraints, not preferences. Ask before doing any of them.

Do not change the mathematical formulation, discretization order, basis,
quadrature, nodal ordering, or time integration scheme. Do not reorder a
floating-point reduction or apply an algebraic simplification; operation order is
part of the specification because bitwise reproducibility may be required.

Do not replace an explicit loop with array syntax, `where`, or `forall`, and do
not reorder a loop nest. Do not introduce temporary allocations, automatic
arrays, or runtime polymorphism into a hot path.

Do not rename a module, type, or public procedure, and do not move files. Do not
change communicator usage or the halo exchange pattern, and do not add a blocking
collective inside a time stepping loop. Do not move an OpenMP or GPU data region.
Do not add an external dependency, and do not modernize to Fortran 2018.

Do not delete a comment describing an equation, an algorithm, or a numerical
assumption.

## Scope

Make the change that was asked for and stop. Provide a minimal diff. Do not
refactor code you happened to read, do not fix unrelated style in a file you are
editing, and do not expand a bug fix into a redesign. If you find a separate
problem, say so in one sentence and record it under "Out of scope" in the pull
request description rather than fixing it here.

New code and bug fixes both need a test. For a bug fix, the test must fail
against the unfixed code, and you should say why the existing suite did not catch
the bug.

Flag uncertainty explicitly. If a change could affect accuracy, stability,
performance, or parallel behaviour and you are not certain, stop and ask.

## Exemplars

Read these before writing:

- `src/SELF_Lagrange_t.f90` for docstring density and argument level `!!` docs
- `src/SELF_LinearEuler2D_t.f90` for a module level LaTeX docstring and inline
  flux annotation
- `src/SELF_LinearShallowWater2D_t.f90` for generic bindings, `do concurrent`,
  and reasoning comments
- `src/SELF_DGModel2D_t.f90` for type and binding block layout
- `docs/Models/linear-shallow-water-model.md` for a model page
- `docs/Tutorials/LinearShallowWater/KelvinWaves.md` for a tutorial page

Do not imitate `src/SELF_SupportRoutines.f90`. Its doxygen markup is a holdover
from an earlier convention and is not the target style.
