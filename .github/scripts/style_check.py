#!/usr/bin/env python3
"""Deterministic coding style checks for SELF Fortran source.

fprettify already enforces layout: indentation, keyword case, and whitespace.
These checks cover the conventions fprettify cannot see, namely how comments
and FORD docstrings are written, so that contributions keep the voice of the
original source rather than drifting toward generated prose.

Every rule was calibrated against the tree at f3e1e57c, the last commit before
the first AI assisted contribution. A rule exists only where that tree shows a
consistent convention; conventions the original authors did not follow
uniformly, such as bullet lists in comments, are deliberately not enforced.
The rules are described in docs/Contributing/StyleGuide.md.

Usage:
    style_check.py                     check src/, test/ and examples/
    style_check.py --files a.f90 b.f90 check an explicit list
    style_check.py --stats [paths]     report the metrics without failing
"""

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from self_style import Violation, comment_marker, load_rules, report, split_comment

DEFAULT_ROOTS = ("src", "test", "examples")

RULER = re.compile(r"^!\s*/{40,}\s*!\s*$")
LICENSE_MARKERS = ("Maintainers : support@fluidnumerics.com", "Fluid Numerics LLC")

# A procedure definition, as opposed to a type bound procedure declaration or
# an END statement. The optional prefix covers typed function results such as
# "real(prec) function Foo(...)".
PROCEDURE = re.compile(
    r"^\s*(?:(?:pure|elemental|impure|recursive|module)\s+)*"
    r"(?:[a-z]+\s*(?:\([^)]*\))?\s*(?:,\s*[a-z_]+\s*)*)?"
    r"\b(subroutine|function)\s+([A-Za-z_][A-Za-z0-9_]*)",
    re.IGNORECASE,
)
MODULE = re.compile(r"^\s*module\s+([A-Za-z_][A-Za-z0-9_]*)\s*$", re.IGNORECASE)
INTERFACE_OPEN = re.compile(r"^\s*(abstract\s+)?interface\b", re.IGNORECASE)
INTERFACE_CLOSE = re.compile(r"^\s*end\s*interface\b", re.IGNORECASE)
PROCEDURE_CLOSE = re.compile(r"^\s*end\s*(subroutine|function)\b", re.IGNORECASE)
IMPLICIT_NONE = re.compile(r"^\s*implicit\s+none\b", re.IGNORECASE)

SPLIT_END = re.compile(
    r"\bend\s+(module|submodule|subroutine|function|type|do|if|interface|program|"
    r"select|associate|where|block|forall|enum|procedure)\b",
    re.IGNORECASE,
)

# Characters that never appear in hand written SELF source below the license
# banner. Scientific notation such as the superscripts in "m s^-1" and the dot
# operator is part of the original style and is deliberately permitted; what is
# banned is typographic punctuation and pictographs, which only ever arrive
# through generated prose.
BANNED_CHARACTERS = {
    "—": "em dash",
    "–": "en dash",
    "‒": "figure dash",
    "―": "horizontal bar",
    "…": "horizontal ellipsis",
    "‘": "left single quotation mark",
    "’": "right single quotation mark",
    "“": "left double quotation mark",
    "”": "right double quotation mark",
    " ": "non-breaking space",
    "•": "bullet",
    "→": "rightwards arrow",
    "←": "leftwards arrow",
    "⇒": "rightwards double arrow",
    "✓": "check mark",
    "✔": "heavy check mark",
    "✗": "ballot x",
    "❌": "cross mark",
    "✅": "white heavy check mark",
    "⚠": "warning sign",
}
EMOJI_RANGES = ((0x1F000, 0x1FAFF), (0x2600, 0x27BF), (0xFE0F, 0xFE0F))

MARKDOWN_HEADING = re.compile(r"^\s*#{1,6}\s+\S")
# Markdown emphasis, written so that Fortran exponentiation such as
# "nhat(1)**2-nhat(2)**2" inside a commented out expression is not matched:
# the opening marker must follow whitespace or start of line, and the closing
# marker must be followed by whitespace or sentence punctuation.
MARKDOWN_BOLD = re.compile(r"(?:^|(?<=\s))\*\*(?=\S)[^*]+?\S\*\*(?=$|[\s.,;:)])")
MARKDOWN_FENCE = re.compile(r"```")
BANNED_WORDS = re.compile(
    r"\b(TODO|FIXME|XXX|HACK|Claude|Copilot|ChatGPT|Co-Authored|AI-generated)\b",
    re.IGNORECASE,
)


def emoji(char):
    """True when the character falls in a pictograph or dingbat block."""
    point = ord(char)
    return any(low <= point <= high for low, high in EMOJI_RANGES)


class FortranFile:
    """A parsed Fortran source file: lines, license extent, and comment map."""

    def __init__(self, path):
        self.path = path
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            self.lines = handle.read().splitlines()
        self.license_end = self._license_end()

    def _license_end(self):
        """Index of the last line of the BSD license banner, or -1 if absent."""
        if not self.lines or not RULER.match(self.lines[0]):
            return -1
        for index in range(1, min(len(self.lines), 60)):
            if RULER.match(self.lines[index]):
                return index
        return -1

    def body(self):
        """Yield (line number, text) for every line after the license banner."""
        start = self.license_end + 1
        for offset, text in enumerate(self.lines[start:], start=start + 1):
            yield offset, text

    def comments(self):
        """Yield (line number, marker, text) for every comment after the banner."""
        for number, text in self.body():
            marker = comment_marker(text)
            if marker is None:
                continue
            _, comment = split_comment(text)
            yield number, marker, (comment or "").rstrip()

    def code(self):
        """Yield (line number, code) for every line with code after the banner."""
        for number, text in self.body():
            code, _ = split_comment(text)
            if code.strip():
                yield number, code


def check_license(source, found):
    """F001: every source file opens with the project BSD-3 banner."""
    if source.license_end < 0:
        found.append(
            Violation(
                "F001",
                source.path,
                1,
                "missing the BSD-3 license banner; copy the 25 line header "
                "verbatim from an existing source file",
            )
        )
        return
    header = "\n".join(source.lines[: source.license_end + 1])
    for marker in LICENSE_MARKERS:
        if marker not in header:
            found.append(
                Violation("F001", source.path, 1, "license banner is missing %r" % marker)
            )


def check_docmark(source, found):
    """F002: FORD documentation uses the post mark !!, never the pre mark !>."""
    for number, marker, _ in source.comments():
        if marker == "!>":
            found.append(
                Violation(
                    "F002",
                    source.path,
                    number,
                    "!> is not used in SELF; document with !! placed after the "
                    "statement being described",
                )
            )


def check_characters(source, found):
    """F003: no typographic punctuation or pictographs below the banner."""
    for number, text in source.body():
        for char in sorted(set(text)):
            if char in BANNED_CHARACTERS:
                found.append(
                    Violation(
                        "F003",
                        source.path,
                        number,
                        "%s (U+%04X); use ASCII punctuation"
                        % (BANNED_CHARACTERS[char], ord(char)),
                    )
                )
            elif emoji(char):
                found.append(
                    Violation(
                        "F003",
                        source.path,
                        number,
                        "pictograph U+%04X does not belong in source" % ord(char),
                    )
                )


def check_markdown(source, found):
    """F004: comments are prose, not markdown."""
    for number, _, text in source.comments():
        if MARKDOWN_HEADING.match(text):
            found.append(
                Violation(
                    "F004",
                    source.path,
                    number,
                    "markdown heading in a comment; describe the routine in "
                    "sentences instead",
                )
            )
        if MARKDOWN_BOLD.search(text):
            found.append(
                Violation("F004", source.path, number, "markdown emphasis in a comment")
            )
        if MARKDOWN_FENCE.search(text):
            found.append(
                Violation("F004", source.path, number, "markdown code fence in a comment")
            )


def check_banned_words(source, found):
    """F005: no work in progress markers and no tool attribution in comments."""
    for number, _, text in source.comments():
        match = BANNED_WORDS.search(text)
        if match:
            found.append(
                Violation(
                    "F005",
                    source.path,
                    number,
                    "%r does not belong in source; open an issue instead"
                    % match.group(0),
                )
            )


def check_end_keywords(source, found):
    """F006: closing keywords are fused, as in endsubroutine and enddo."""
    for number, code in source.code():
        match = SPLIT_END.search(code)
        if match:
            found.append(
                Violation(
                    "F006",
                    source.path,
                    number,
                    "write %r as %r" % (match.group(0), "end" + match.group(1).lower()),
                )
            )


def procedures(source):
    """Yield (line index, name) for each procedure defined outside an interface.

    Interface bodies are skipped because they declare a procedure without
    providing one, so they carry neither implicit none nor a docstring.
    """
    depth = 0
    for index, text in enumerate(source.lines):
        code, _ = split_comment(text)
        if INTERFACE_OPEN.match(code):
            depth += 1
            continue
        if INTERFACE_CLOSE.match(code):
            depth = max(0, depth - 1)
            continue
        if depth or PROCEDURE_CLOSE.match(code):
            continue
        match = PROCEDURE.match(code)
        if match:
            yield index, match.group(2)


def check_implicit_none(source, found):
    """F007: implicit none appears at module scope and in every procedure.

    This is a hard requirement of CLAUDE.md section 2. Parts of the legacy
    source predate it, so this rule reports genuine pre-existing gaps as well
    as new ones.
    """
    if any(MODULE.match(split_comment(t)[0]) for t in source.lines):
        module_scope = []
        for text in source.lines:
            code, _ = split_comment(text)
            if PROCEDURE.match(code) and not PROCEDURE_CLOSE.match(code):
                break
            module_scope.append(code)
        if not any(IMPLICIT_NONE.match(c) for c in module_scope):
            found.append(
                Violation("F007", source.path, 1, "module scope is missing implicit none")
            )

    for index, name in procedures(source):
        body = []
        for text in source.lines[index + 1 :]:
            code, _ = split_comment(text)
            if PROCEDURE_CLOSE.match(code):
                break
            body.append(code)
        if not any(IMPLICIT_NONE.match(c) for c in body):
            found.append(
                Violation(
                    "F007",
                    source.path,
                    index + 1,
                    "procedure %r is missing implicit none" % name,
                )
            )


def measure(source):
    """Return the numeric style metrics for one file.

    docstring_coverage is the fraction of procedures whose signature is
    followed by a !! docstring, comment_density is comment lines over code
    lines with the license banner excluded, and max_comment_length is the
    longest comment bearing line.
    """
    total = 0
    documented = 0
    for index, _ in procedures(source):
        total += 1
        for follower in source.lines[index + 1 : index + 3]:
            if comment_marker(follower) == "!!":
                documented += 1
                break

    comment_lines = 0
    longest = 0
    for number, marker, _ in source.comments():
        comment_lines += 1
        text = source.lines[number - 1]
        if text.lstrip().startswith("!"):
            longest = max(longest, len(text))
    code_lines = sum(1 for _ in source.code())

    return {
        "procedures": total,
        "docstring_coverage": round(documented / total, 3) if total else None,
        "comment_density": round(comment_lines / code_lines, 3) if code_lines else None,
        "max_comment_length": longest,
    }


def in_scope(path, config):
    """True when the file matches a metric's applies_to and no excludes entry.

    Paths are matched on the repository relative form so that the checks give
    the same answer whether they are run from the repository root, from a
    pre-commit hook, or against an extracted reference tree.
    """
    normalized = path.replace(os.sep, "/")
    while normalized.startswith("./"):
        normalized = normalized[2:]
    for scope in config["applies_to"]:
        if normalized.startswith(scope) or ("/" + scope) in normalized:
            break
    else:
        return False
    for scope in config.get("excludes", ()):
        if normalized.startswith(scope) or ("/" + scope) in normalized:
            return False
    return True


def check_metrics(source, metrics, rules, found):
    """F101 and F102: the numeric thresholds, applied per configured scope."""
    thresholds = rules["fortran"]["metrics"]

    density = thresholds["F101"]
    if (
        metrics["comment_density"] is not None
        and in_scope(source.path, density)
        and metrics["comment_density"] < density["min"]
    ):
        found.append(
            Violation(
                "F101",
                source.path,
                None,
                "comment density %.3f is below the minimum %.3f; explain what "
                "the routines compute, as the surrounding core source does"
                % (metrics["comment_density"], density["min"]),
            )
        )

    limit = thresholds["F102"]["max"]
    for number, _, _ in source.comments():
        text = source.lines[number - 1]
        if not text.lstrip().startswith("!"):
            # Trailing comments share a line with code, whose length fprettify
            # already governs; only whole line comments are measured here.
            continue
        if len(text) > limit:
            found.append(
                Violation(
                    "F102",
                    source.path,
                    number,
                    "comment line is %d characters, the limit is %d"
                    % (len(text), limit),
                )
            )


BINARY_CHECKS = (
    ("F001", check_license),
    ("F002", check_docmark),
    ("F003", check_characters),
    ("F004", check_markdown),
    ("F005", check_banned_words),
    ("F006", check_end_keywords),
    ("F007", check_implicit_none),
)


def check_file(path, rules):
    """Run every enabled rule against one file, returning violations and metrics."""
    source = FortranFile(path)
    found = []
    disabled = set(rules["fortran"].get("disabled", ()))
    for rule, check in BINARY_CHECKS:
        if rule not in disabled:
            check(source, found)
    metrics = measure(source)
    check_metrics(source, metrics, rules, found)
    return [v for v in found if v.rule not in disabled], metrics


def collect(paths):
    """Expand the given paths into the list of Fortran files to check."""
    files = []
    for path in paths:
        if os.path.isdir(path):
            for root, _, names in os.walk(path):
                files.extend(
                    os.path.join(root, name)
                    for name in sorted(names)
                    if name.endswith(".f90")
                )
        elif path.endswith(".f90") and os.path.isfile(path):
            files.append(path)
    return sorted(files)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", help="files or directories to check")
    parser.add_argument("--files", nargs="*", default=None, help="explicit file list")
    parser.add_argument("--stats", action="store_true", help="report metrics, exit 0")
    parser.add_argument("--rules", default=None, help="path to style-rules.json")
    args = parser.parse_args(argv)

    targets = args.files if args.files is not None else (args.paths or list(DEFAULT_ROOTS))
    files = collect(targets)
    if not files:
        print("No Fortran files to check.")
        return 0

    rules = load_rules(args.rules)
    violations = []
    stats = {}
    for path in files:
        found, metrics = check_file(path, rules)
        violations.extend(found)
        stats[path] = metrics

    if args.stats:
        report([], stats=stats)
        return 0
    return report(violations)


if __name__ == "__main__":
    sys.exit(main())
