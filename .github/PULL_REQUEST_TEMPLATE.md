<!--
Every section below is required and is checked automatically by the
pr-template-check workflow. Replace the guidance in each section; do not
delete the headings.
-->

## Scope

<!--
What this pull request changes, as a bounded list. A reviewer should be able to
read this and know what they are being asked to look at. If a change is not
described here, it does not belong in this pull request.
-->

## Out of scope

<!--
Problems you found while working on this that you are deliberately not fixing
here. Reference an issue for each one, or say that it is not yet filed. Write
"None" if you found nothing.

Anything listed here is off limits for review comments, so be specific.
-->

## Issues resolved

<!--
Link the issue this pull request closes, using a closing keyword so that GitHub
closes it on merge:

    Fixes #123

Every pull request resolves a filed issue. If no issue exists, open one first
and describe the problem there rather than in this description.
-->

## Tests introduced

<!--
Which tests are new, and what each one establishes.

For a bug fix, name the test that now fails against the old code, and say why
the existing suite did not catch the bug. A bug that reached main is evidence
of a gap in the tests, and closing that gap is part of the fix.

For a new feature, give the validation criteria: the quantity being measured
and the tolerance it must meet.

If this pull request genuinely introduces no testable behaviour, say so and
explain why.
-->

## Checklist

- [ ] `fprettify` has been run over the directories this pull request touches
- [ ] `python3 .github/scripts/style_check.py` reports no new violations
- [ ] `python3 .github/scripts/docs_style_check.py` reports no new violations
- [ ] New tests pass, and MPI tests were exercised on at least two ranks
- [ ] No change to the mathematical formulation, discretization order, quadrature,
      nodal ordering, or time integration scheme without maintainer approval
