---
applyTo: "docs/**/*.md,*.md"
---

# Documentation conventions

SELF documentation is prose. It explains the mathematics in paragraphs, gives the
equations in LaTeX, and shows the calling sequence in a Fortran block. It is not
a slide deck and not a checklist. `docs/Contributing/StyleGuide.md` describes
each rule and names the check that enforces it.

## Voice

- First person plural for what SELF and its authors do: "we use a local upwind
  (Lax-Friedrich's) Riemann solver".
- Second person for the reader: "Assuming you've created interpolant, mesh, and
  geometry objects, you can set the `fCori%interior` attribute".
- Write full sentences in paragraphs. One sentence per line, no hard wrapping.
- Reserve lists for genuine enumerations, such as a parameter list with units.
  Do not convert an explanation into bullets.
- Emphasis is rare. Carry the emphasis in the sentence rather than in `**bold**`.

## Mechanics

- ASCII punctuation only. No em dashes, no typographic quotes, no emoji.
- One level one heading, on the first line of the page. Headings go no deeper
  than level four; a page that needs more depth should be split.
- Mathematics uses MathJax: `$inline$` and display `$$ ... $$`. The established
  notation is `\vec{s}` for the solution vector, `\overleftrightarrow{f}` for the
  conservative flux, and `\vec{q}` for the source terms.
- A new page under `docs/` must be added to the `nav:` block of `mkdocs.yml`, or
  it is published but unreachable.
- Do not add a checklist of work items. Track work in issues.

## Page structure

A model page opens with `## Definition`: a link to the FORD generated module and
type pages, the generic conservation law
$\vec{s}_t + \nabla \cdot \overleftrightarrow{f} = \vec{q}$, and then the
specialization of $\vec{s}$, $\overleftrightarrow{f}$, and $\vec{q}$ for this
model. It continues with `## Implementation`: which base class is extended and
which type-bound procedures are overridden, each subsection closing with a
`fortran` code block.

A tutorial page opens with a paragraph of physical motivation, then
`## Configuration` with the equations, parameters and units, mesh, polynomial
degree, quadrature, and the CFL arithmetic worked out inline, then the initial
condition with figures, then `## Running this example` with shell blocks.

Read `docs/Models/linear-shallow-water-model.md` and
`docs/Tutorials/LinearShallowWater/KelvinWaves.md` before writing a new page.
