"""Shared helpers for the SELF deterministic style checks.

This module holds the pieces that both style_check.py and docs_style_check.py
need: violation reporting, the rule threshold configuration, and the small
amount of Fortran lexing required to tell code apart from comments.

Only the Python standard library is used so that the checks run in the same
bare CI runner that hosts the fprettify lint job, and in a pre-commit hook,
without installing anything.
"""

import json
import os
import re
import sys

DEFAULT_RULES_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "style-rules.json"
)


class Violation:
    """A single rule failure, anchored to a file and (optionally) a line."""

    def __init__(self, rule, path, line, message):
        self.rule = rule
        self.path = path
        self.line = line
        self.message = message

    def __str__(self):
        where = "%s:%d" % (self.path, self.line) if self.line else self.path
        return "%s: %s: %s" % (where, self.rule, self.message)


def load_rules(path=None):
    """Read the rule configuration, falling back to the checked-in defaults."""
    with open(path or DEFAULT_RULES_PATH, "r", encoding="utf-8") as handle:
        return json.load(handle)


def report(violations, stats=None, stream=sys.stdout):
    """Print violations grouped by file, then a one line summary.

    Returns the process exit status: 0 when there are no violations.
    """
    if stats:
        for path in sorted(stats):
            measured = stats[path]
            summary = ", ".join("%s=%s" % (k, measured[k]) for k in sorted(measured))
            stream.write("%s: %s\n" % (path, summary))

    by_file = {}
    for violation in violations:
        by_file.setdefault(violation.path, []).append(violation)

    for path in sorted(by_file):
        for violation in sorted(by_file[path], key=lambda v: (v.line or 0, v.rule)):
            stream.write("%s\n" % violation)

    if violations:
        stream.write(
            "\n%d style violation(s) in %d file(s). "
            "See docs/Contributing/StyleGuide.md for the rule descriptions.\n"
            % (len(violations), len(by_file))
        )
        return 1

    return 0


# Fortran string literals are masked before comment detection so that an
# exclamation point inside a character constant is not mistaken for a comment.
_STRING = re.compile(r"'[^']*'|\"[^\"]*\"")


def mask_strings(line):
    """Return the line with the contents of character constants blanked out.

    Any check that looks for Fortran syntax must run against the masked form,
    otherwise a keyword quoted inside a string, as in
    ``write(*,*) "expected end if"``, is read as code.
    """
    return _STRING.sub(lambda m: m.group(0)[0] + " " * (len(m.group(0)) - 2) +
                       m.group(0)[-1], line)


def split_comment(line):
    """Split a Fortran source line into its code and comment halves.

    Returns (code, comment) where code is the masked source preceding the
    first comment marker and comment is the text following it with the marker
    removed, or None when the line has no comment. Character constants are
    masked so that neither half is confused by a quoted exclamation point or a
    quoted keyword.
    """
    masked = mask_strings(line)
    index = masked.find("!")
    if index < 0:
        return masked, None
    return masked[:index], line[index:].lstrip("!")


def comment_marker(line):
    """Return the comment marker used on a line: '!!', '!>', '!' or None."""
    masked = _STRING.sub(lambda m: " " * len(m.group(0)), line)
    index = masked.find("!")
    if index < 0:
        return None
    rest = line[index:]
    if rest.startswith("!!"):
        return "!!"
    if rest.startswith("!>"):
        return "!>"
    return "!"
