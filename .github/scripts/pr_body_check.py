#!/usr/bin/env python3
"""Check that a pull request description follows the repository template.

A pull request is reviewable only when its scope is stated. This check reads
the description from the GitHub event payload and confirms that the author
answered the four questions the template asks: what is in scope, what is
knowingly left out, which filed issue is being resolved, and which tests were
added and why.

The check is skipped when the pull request carries the label named by
--skip-label, so that a maintainer can land an urgent fix.

Usage:
    pr_body_check.py                    read $GITHUB_EVENT_PATH
    pr_body_check.py --event event.json read an explicit payload
"""

import argparse
import json
import os
import re
import sys

REQUIRED_SECTIONS = ("Scope", "Out of scope", "Issues resolved", "Tests introduced")

# Sections whose only legitimate short answer is that there is nothing to say.
ALLOW_EMPTY_ANSWER = {"Out of scope": ("none", "nothing", "n/a")}

PLACEHOLDERS = ("n/a", "na", "none", "tbd", "todo", "-", ".", "xxx")
CLOSING_KEYWORD = re.compile(
    r"\b(clos(?:e|es|ed)|fix(?:e[sd])?|resolv(?:e|es|ed))\b\s*:?\s*"
    r"(?:#\d+|https?://\S+/issues/\d+)",
    re.IGNORECASE,
)
BARE_ISSUE = re.compile(r"(?:#\d+|https?://\S+/issues/\d+)")
HTML_COMMENT = re.compile(r"<!--.*?-->", re.DOTALL)
HEADING = re.compile(r"^\s{0,3}#{1,6}\s+(.*?)\s*#*\s*$")


def sections(body):
    """Split a pull request body into a mapping of heading to section text.

    HTML comments are removed first so that the guidance carried by the
    template does not count as an answer.
    """
    text = HTML_COMMENT.sub("", body or "")
    found = {}
    current = None
    for line in text.splitlines():
        match = HEADING.match(line)
        if match:
            current = match.group(1).strip()
            found.setdefault(current, [])
            continue
        if current is not None:
            found[current].append(line)
    return {name: "\n".join(lines).strip() for name, lines in found.items()}


def is_answered(name, content):
    """True when a section carries a real answer rather than a placeholder."""
    if not content:
        return False
    stripped = content.strip().strip("*_`").rstrip(".").strip().lower()
    if stripped in ALLOW_EMPTY_ANSWER.get(name, ()):
        return True
    if stripped in PLACEHOLDERS:
        return False
    # A section holding nothing but an unticked checklist has not been answered.
    return bool(re.sub(r"^\s*[-*+]\s*\[\s*\]\s*$", "", content, flags=re.M).strip())


def check(body, labels, skip_label):
    """Return the list of problems with a pull request description."""
    if skip_label and skip_label in labels:
        return []

    problems = []
    present = sections(body)
    normalized = {name.lower(): name for name in present}

    for name in REQUIRED_SECTIONS:
        actual = normalized.get(name.lower())
        if actual is None:
            problems.append(
                "the '## %s' section is missing; copy the headings from "
                ".github/PULL_REQUEST_TEMPLATE.md" % name
            )
            continue
        content = present[actual]
        if not is_answered(name, content):
            allowed = ALLOW_EMPTY_ANSWER.get(name)
            hint = (
                " (write %s if there is nothing to record)"
                % " or ".join(repr(word) for word in allowed)
                if allowed
                else ""
            )
            problems.append("the '## %s' section is empty%s" % (name, hint))
            continue

        if name == "Issues resolved" and not CLOSING_KEYWORD.search(content):
            if BARE_ISSUE.search(content):
                problems.append(
                    "'## Issues resolved' references an issue but without a "
                    "closing keyword; write 'Fixes #123' so the issue closes on "
                    "merge"
                )
            else:
                problems.append(
                    "'## Issues resolved' does not reference an issue; open an "
                    "issue describing the problem and link it as 'Fixes #123'"
                )

    return problems


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--event",
        default=os.environ.get("GITHUB_EVENT_PATH"),
        help="path to the GitHub event payload",
    )
    parser.add_argument(
        "--skip-label",
        default="skip-pr-checks",
        help="label that bypasses this check",
    )
    args = parser.parse_args(argv)

    if not args.event or not os.path.isfile(args.event):
        print("No event payload to read; pass --event or set GITHUB_EVENT_PATH.")
        return 2

    with open(args.event, "r", encoding="utf-8") as handle:
        payload = json.load(handle)

    pull_request = payload.get("pull_request") or {}
    body = pull_request.get("body") or ""
    labels = {label.get("name") for label in pull_request.get("labels") or []}

    if args.skip_label in labels:
        print("Label %r is set; skipping the description check." % args.skip_label)
        return 0

    problems = check(body, labels, args.skip_label)
    if not problems:
        print("Pull request description is complete.")
        return 0

    print("This pull request description is incomplete:\n")
    for problem in problems:
        print("  - %s" % problem)
    print(
        "\nSee .github/PULL_REQUEST_TEMPLATE.md. Edit the description and this "
        "check will run again."
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
