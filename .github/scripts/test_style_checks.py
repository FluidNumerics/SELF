#!/usr/bin/env python3
"""Tests for the deterministic style checks.

The checkers gate every pull request, so a broken rule is expensive: it either
blocks work that is fine or lets through what it was meant to catch. These
tests pin the behaviour of each rule against a small input rather than against
the repository, so that they keep working as the tree changes.

Run them with:

    python3 -m unittest discover -s .github/scripts -p 'test_*.py'
"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import docs_style_check
import pr_body_check
import style_check
from self_style import load_rules, mask_strings, split_comment

RULES = load_rules(
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "style-rules.json"
    )
)
LICENSE = open(
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "license-header.txt",
    ),
    encoding="utf-8",
).read()


class FortranCase(unittest.TestCase):
    """Base class providing a temporary Fortran file with a valid banner."""

    def rules_for(self, path, body, license_text=LICENSE):
        """Write a source file and return the rule identifiers it trips."""
        with tempfile.TemporaryDirectory() as directory:
            full = os.path.join(directory, path)
            os.makedirs(os.path.dirname(full), exist_ok=True)
            with open(full, "w", encoding="utf-8") as handle:
                handle.write(license_text + "\n" + body)
            found, _ = style_check.check_file(
                full, RULES, style_check.canonical_license()
            )
        return sorted({violation.rule for violation in found})


MODULE_BODY = """
module self_Example
  !! A minimal module used as a fixture for the rule tests.
  implicit none

contains

  subroutine Init_Example(this)
    !! Set the counter to zero.
    implicit none
    class(Example_t),intent(out) :: this
    this%n = 0
  endsubroutine Init_Example

endmodule self_Example
"""


class TestFortranRules(FortranCase):
    def test_clean_module_passes(self):
        self.assertEqual(self.rules_for("src/a.f90", MODULE_BODY), [])

    def test_f001_rejects_a_truncated_banner(self):
        truncated = "\n".join(LICENSE.splitlines()[:12] + [LICENSE.splitlines()[-1]])
        self.assertIn("F001", self.rules_for("src/a.f90", MODULE_BODY, truncated))

    def test_f001_rejects_a_corrupted_banner(self):
        corrupted = LICENSE.replace("BUSINESS INTERRUPTION", "BUsLESS INTERRUPTION")
        self.assertIn("F001", self.rules_for("src/a.f90", MODULE_BODY, corrupted))

    def test_f001_tolerates_the_year_and_quote_variants(self):
        variant = LICENSE.replace("2024", "2026").replace("\u201cAS IS\u201d", '"AS IS"')
        self.assertEqual(self.rules_for("src/a.f90", MODULE_BODY, variant), [])

    def test_f001_reports_a_missing_banner(self):
        self.assertIn("F001", self.rules_for("src/a.f90", MODULE_BODY, ""))

    def test_f002_rejects_the_predocmark(self):
        body = MODULE_BODY.replace(
            "module self_Example", "!> describes the module\nmodule self_Example"
        )
        self.assertIn("F002", self.rules_for("src/a.f90", body))

    def test_f003_rejects_an_em_dash(self):
        body = MODULE_BODY.replace(
            "    this%n = 0", "    ! set the counter — it starts at zero\n    this%n = 0"
        )
        self.assertIn("F003", self.rules_for("src/a.f90", body))

    def test_f003_allows_scientific_superscripts(self):
        body = MODULE_BODY.replace(
            "    this%n = 0", "    this%n = 0 ! wave speed in m s⁻¹"
        )
        self.assertEqual(self.rules_for("src/a.f90", body), [])

    def test_f004_rejects_markdown_in_a_comment(self):
        for comment in ("!! ## Purpose", "!! see **this**", "!! ```fortran"):
            body = MODULE_BODY.replace("module self_Example", comment + "\nmodule self_Example")
            self.assertIn("F004", self.rules_for("src/a.f90", body), comment)

    def test_f004_allows_fortran_exponentiation_in_a_comment(self):
        body = MODULE_BODY.replace(
            "    this%n = 0", "    ! exts(2) = (nhat(1)**2-nhat(2)**2)*s(2)\n    this%n = 0"
        )
        self.assertEqual(self.rules_for("src/a.f90", body), [])

    def test_f005_rejects_a_work_marker(self):
        body = MODULE_BODY.replace("    this%n = 0", "    ! TODO fix this\n    this%n = 0")
        self.assertIn("F005", self.rules_for("src/a.f90", body))

    def test_f006_rejects_a_split_end_keyword(self):
        body = MODULE_BODY.replace("endsubroutine Init_Example", "end subroutine Init_Example")
        self.assertIn("F006", self.rules_for("src/a.f90", body))

    def test_f006_ignores_a_keyword_inside_a_string(self):
        body = MODULE_BODY.replace(
            "    this%n = 0", '    write(*,*) "expected end if here"'
        )
        self.assertNotIn("F006", self.rules_for("src/a.f90", body))

    def test_f007_reports_a_procedure_without_implicit_none(self):
        body = MODULE_BODY.replace("    implicit none\n    class(", "    class(")
        self.assertIn("F007", self.rules_for("src/a.f90", body))

    def test_f007_reports_a_program_without_implicit_none(self):
        body = "\nprogram example\n  integer :: i\n  i = 0\nendprogram example\n"
        self.assertIn("F007", self.rules_for("test/a.f90", body))

    def test_f007_accepts_a_program_with_implicit_none(self):
        body = "\nprogram example\n  implicit none\n  integer :: i\n  i = 0\nendprogram example\n"
        self.assertNotIn("F007", self.rules_for("test/a.f90", body))

    def test_f007_does_not_borrow_implicit_none_from_an_internal_procedure(self):
        body = """
program example
  integer :: i
  i = 0
contains
  subroutine inner()
    implicit none
  endsubroutine inner
endprogram example
"""
        self.assertIn("F007", self.rules_for("test/a.f90", body))

    def test_f007_reports_a_submodule_without_implicit_none(self):
        body = (
            "\nsubmodule (self_Example) self_Example_impl\n"
            "  integer :: i\n"
            "endsubmodule self_Example_impl\n"
        )
        self.assertIn("F007", self.rules_for("src/a.f90", body))

    def test_f007_accepts_a_submodule_with_implicit_none(self):
        body = (
            "\nsubmodule (self_Example) self_Example_impl\n"
            "  !! Implementation of the example interface.\n"
            "  implicit none\n"
            "  integer :: i\n"
            "endsubmodule self_Example_impl\n"
        )
        self.assertNotIn("F007", self.rules_for("src/a.f90", body))

    def test_f101_reports_a_sparsely_commented_core_file(self):
        body = "\nmodule self_Example\n  implicit none\n" + "".join(
            "  integer :: v%d = 0\n" % n for n in range(40)
        ) + "endmodule self_Example\n"
        self.assertIn("F101", self.rules_for("src/a.f90", body))

    def test_f101_does_not_apply_to_the_backend_directories(self):
        body = "\nmodule self_Example\n  implicit none\n" + "".join(
            "  integer :: v%d = 0\n" % n for n in range(40)
        ) + "endmodule self_Example\n"
        self.assertNotIn("F101", self.rules_for("src/gpu/a.f90", body))

    def test_f007_reports_a_separate_module_procedure_without_implicit_none(self):
        body = """
submodule (self_Example) self_Example_impl
  !! Implementation of the example interface.
  implicit none
contains
  module procedure Init_Example
    this%n = 0
  endprocedure Init_Example
endsubmodule self_Example_impl
"""
        self.assertIn("F007", self.rules_for("src/a.f90", body))

    def test_f007_accepts_a_separate_module_procedure_with_implicit_none(self):
        body = """
submodule (self_Example) self_Example_impl
  !! Implementation of the example interface.
  implicit none
contains
  module procedure Init_Example
    implicit none
    this%n = 0
  endprocedure Init_Example
endsubmodule self_Example_impl
"""
        self.assertNotIn("F007", self.rules_for("src/a.f90", body))

    def test_f102_reports_a_long_comment_line(self):
        body = MODULE_BODY.replace("    this%n = 0", "    ! " + "x" * 140 + "\n    this%n = 0")
        self.assertIn("F102", self.rules_for("src/a.f90", body))


class TestSplitComment(unittest.TestCase):
    def test_a_quoted_exclamation_is_not_a_comment(self):
        code, comment = split_comment('  s = "a!b"')
        self.assertIsNone(comment)

    def test_a_trailing_comment_is_returned_verbatim(self):
        code, comment = split_comment("  x = 1 ! density")
        self.assertEqual(comment.strip(), "density")

    def test_masking_preserves_line_length(self):
        line = '  write(*,*) "hello"'
        self.assertEqual(len(mask_strings(line)), len(line))


PROSE = "\n".join(
    "This sentence explains one aspect of the method in a reasonable number of words."
    for _ in range(30)
)


class TestDocsRules(unittest.TestCase):
    def check(self, text, nav=(), path="docs/Models/page.md"):
        with tempfile.TemporaryDirectory() as directory:
            full = os.path.join(directory, path)
            os.makedirs(os.path.dirname(full), exist_ok=True)
            with open(full, "w", encoding="utf-8") as handle:
                handle.write(text)
            found, _ = docs_style_check.check_file(full, RULES, set(nav), "docs")
        return sorted({violation.rule for violation in found})

    def test_clean_page_passes(self):
        page = "# Title\n\n" + PROSE + "\n"
        self.assertEqual(self.check(page, nav=("Models/page.md",)), [])

    def test_d001_rejects_an_emoji(self):
        page = "# Title\n\nStatus ✅ done.\n\n" + PROSE + "\n"
        self.assertIn("D001", self.check(page, nav=("Models/page.md",)))

    def test_d002_rejects_an_em_dash(self):
        page = "# Title\n\nThe flux — which is conservative — is computed here.\n\n" + PROSE
        self.assertIn("D002", self.check(page, nav=("Models/page.md",)))

    def test_d002_ignores_a_fenced_block(self):
        page = "# Title\n\n```\na — b\n```\n\n" + PROSE + "\n"
        self.assertNotIn("D002", self.check(page, nav=("Models/page.md",)))

    def test_d003_reports_a_missing_title(self):
        self.assertIn("D003", self.check("## Section\n\n" + PROSE, nav=("Models/page.md",)))

    def test_d003_allows_leading_blank_lines(self):
        page = "\n\n# Title\n\n" + PROSE + "\n"
        self.assertNotIn("D003", self.check(page, nav=("Models/page.md",)))

    def test_d003_reports_a_title_below_content(self):
        page = "Some text first.\n\n# Title\n\n" + PROSE + "\n"
        self.assertIn("D003", self.check(page, nav=("Models/page.md",)))

    def test_d004_rejects_a_level_five_heading(self):
        page = "# Title\n\n##### Too deep\n\n" + PROSE + "\n"
        self.assertIn("D004", self.check(page, nav=("Models/page.md",)))

    def test_d005_rejects_a_task_list(self):
        page = "# Title\n\n- [ ] do the thing\n\n" + PROSE + "\n"
        self.assertIn("D005", self.check(page, nav=("Models/page.md",)))

    def test_d006_reports_a_page_missing_from_the_nav(self):
        self.assertIn("D006", self.check("# Title\n\n" + PROSE, nav=()))

    def test_d100_ignores_inline_code(self):
        page = "# Title\n\n" + "\n".join(
            "Use the `__shared__` qualifier for this buffer in the kernel body."
            for _ in range(30)
        )
        self.assertNotIn("D100", self.check(page, nav=("Models/page.md",)))

    def test_d100_reports_heavy_emphasis(self):
        page = "# Title\n\n" + "\n".join(
            "**Purpose** the routine does something for this case." for _ in range(30)
        )
        self.assertIn("D100", self.check(page, nav=("Models/page.md",)))


COMPLETE_BODY = """## Scope
Adds a thing.

## Out of scope
None

## Issues resolved
Fixes #185

## Tests introduced
A unit test, because nothing covered this before.
"""


class TestNavRegressions(unittest.TestCase):
    """The --nav-only baseline filter, which decides what the gate fails on."""

    def violation(self, page):
        from self_style import Violation

        return Violation(
            "D006", "docs/" + page, None,
            "page is not in the mkdocs.yml nav tree, so it is published but "
            "unreachable; add an entry for %r" % page,
        )

    def test_page_orphaned_by_this_change_is_a_regression(self):
        found = [self.violation("Models/one.md")]
        regressions, skipped = docs_style_check.nav_regressions(
            found, {"Models/one.md", "index.md"}
        )
        self.assertEqual(len(regressions), 1)
        self.assertEqual(skipped, 0)

    def test_page_already_orphaned_is_left_to_the_backlog(self):
        found = [self.violation("Models/one.md")]
        regressions, skipped = docs_style_check.nav_regressions(found, {"index.md"})
        self.assertEqual(regressions, [])
        self.assertEqual(skipped, 1)

    def test_the_two_are_separated_in_one_pass(self):
        found = [self.violation("Models/one.md"), self.violation("Models/two.md")]
        regressions, skipped = docs_style_check.nav_regressions(
            found, {"Models/one.md"}
        )
        self.assertEqual([docs_style_check.violation_page(v) for v in regressions],
                         ["Models/one.md"])
        self.assertEqual(skipped, 1)


class TestNavParsing(unittest.TestCase):
    def parse(self, text):
        with tempfile.NamedTemporaryFile("w", suffix=".yml", delete=False) as handle:
            handle.write(text)
            path = handle.name
        try:
            return docs_style_check.nav_pages(path)
        finally:
            os.unlink(path)

    def test_nested_entries_are_found(self):
        nav = self.parse(
            "site_name: x\nnav:\n"
            "  - Home: index.md\n"
            "  - Models:\n"
            "    - One: Models/one.md\n"
            "theme:\n  name: material\n"
        )
        self.assertEqual(nav, {"index.md", "Models/one.md"})

    def test_commented_entries_are_ignored(self):
        nav = self.parse("nav:\n  - Home: index.md\n  # - Old: Models/old.md\n")
        self.assertEqual(nav, {"index.md"})

    def test_keys_after_the_nav_block_are_ignored(self):
        nav = self.parse("nav:\n  - Home: index.md\nextra:\n  - Other: nope.md\n")
        self.assertEqual(nav, {"index.md"})


class TestPullRequestBody(unittest.TestCase):
    def problems(self, body, labels=()):
        return pr_body_check.check(body, set(labels), "skip-pr-checks")

    def test_complete_body_passes(self):
        self.assertEqual(self.problems(COMPLETE_BODY), [])

    def test_missing_section_is_reported(self):
        body = COMPLETE_BODY.split("## Tests introduced")[0]
        self.assertEqual(len(self.problems(body)), 1)

    def test_bare_issue_reference_is_rejected(self):
        body = COMPLETE_BODY.replace("Fixes #185", "#185")
        self.assertTrue(any("closing keyword" in p for p in self.problems(body)))

    def test_absent_issue_reference_is_rejected(self):
        body = COMPLETE_BODY.replace("Fixes #185", "It fixes the flaky test.")
        self.assertTrue(any("does not reference an issue" in p for p in self.problems(body)))

    def test_none_is_accepted_for_out_of_scope_only(self):
        self.assertEqual(self.problems(COMPLETE_BODY.replace("Adds a thing.", "None")) and True, True)
        body = COMPLETE_BODY.replace("Adds a thing.", "None")
        self.assertTrue(any("'## Scope' section is empty" in p for p in self.problems(body)))

    def test_heading_inside_a_fence_is_not_a_section(self):
        body = COMPLETE_BODY.split("## Tests introduced")[0] + (
            "Example:\n\n```\n## Tests introduced\nnot a real answer\n```\n"
        )
        self.assertTrue(any("Tests introduced" in p for p in self.problems(body)))

    def test_html_comment_guidance_does_not_count_as_an_answer(self):
        body = "## Scope\n<!-- describe the scope -->\n" + COMPLETE_BODY.split("## Out of scope")[1]
        self.assertTrue(any("'## Scope'" in p for p in self.problems(body)))

    def test_skip_label_bypasses_every_check(self):
        self.assertEqual(self.problems("", labels=("skip-pr-checks",)), [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
