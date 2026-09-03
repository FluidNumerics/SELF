#!/usr/bin/env python3
"""Deterministic documentation style checks for the SELF handbook.

SELF documentation is written as continuous prose: paragraphs that explain the
mathematics, with LaTeX for the equations and code blocks for the calling
sequence. Generated documentation drifts toward a different shape, dense with
bold labels, short fragments, and deep heading trees, and these checks measure
that difference rather than judging the content.

Every limit is calibrated against the docs at f3e1e57c, the last commit before
the first AI assisted contribution, and each rule is described in
docs/Contributing/StyleGuide.md.

Usage:
    docs_style_check.py                    check docs/ and the root pages
    docs_style_check.py --files a.md b.md  check an explicit list
    docs_style_check.py --stats [paths]    report the metrics without failing
"""

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from self_style import FenceTracker, Violation, load_rules, report

# The handbook plus the markdown pages at the repository root. This mirrors the
# path filter the style-check workflow uses to pick changed files.
DEFAULT_ROOTS = ("docs",) + tuple(
    sorted(name for name in os.listdir(".") if name.endswith(".md"))
) if os.path.isdir("docs") else ("docs",)

# Typographic punctuation that the hand written documentation never uses. The
# ASCII replacement is given so that the failure message is actionable.
BANNED_CHARACTERS = {
    "—": ("em dash", "-- or a comma"),
    "–": ("en dash", "-"),
    "‒": ("figure dash", "-"),
    "―": ("horizontal bar", "-"),
    "…": ("horizontal ellipsis", "..."),
    "‘": ("left single quotation mark", "'"),
    "’": ("right single quotation mark", "'"),
    "“": ("left double quotation mark", '"'),
    "”": ("right double quotation mark", '"'),
    " ": ("non-breaking space", "a plain space"),
    "•": ("bullet", "a markdown list"),
}
EMOJI_RANGES = ((0x1F000, 0x1FAFF), (0x2600, 0x27BF), (0xFE0F, 0xFE0F))

FENCE = re.compile(r"^\s*(```|~~~)")
HEADING = re.compile(r"^(#{1,10})\s+(\S.*)$")
TASK_LIST = re.compile(r"^\s*[-*+]\s+\[[ xX]\]\s")
BULLET = re.compile(r"^\s*(?:[-*+]|\d+\.)\s+\S")
TABLE_ROW = re.compile(r"^\s*\|")
HTML_BLOCK = re.compile(r"^\s*<")
BOLD = re.compile(r"\*\*(?=\S)[^*]+?\S\*\*|__(?=\S)[^_]+?\S__")
INLINE_MATH = re.compile(r"\$[^$]*\$")
INLINE_CODE = re.compile(r"`[^`]*`")
LINK = re.compile(r"\[([^\]]*)\]\([^)]*\)")
SENTENCE_END = re.compile(r"[.!?](?:\s|$)")


def emoji(char):
    """True when the character falls in a pictograph or dingbat block."""
    point = ord(char)
    return any(low <= point <= high for low, high in EMOJI_RANGES)


class MarkdownPage:
    """A markdown page split into headings, prose, and fenced code."""

    def __init__(self, path):
        self.path = path
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            self.lines = handle.read().splitlines()
        self.front_matter = self._front_matter()
        self.fenced = self._fenced()

    def _front_matter(self):
        """Return the number of leading lines taken by a YAML front matter block.

        mkdocs and the agent instruction files both allow a page to open with a
        "---" delimited metadata block. It is metadata, not content, so it does
        not count against the title.
        """
        if not self.lines or self.lines[0].strip() != "---":
            return 0
        for index in range(1, len(self.lines)):
            if self.lines[index].strip() in ("---", "..."):
                return index + 1
        return 0

    def _fenced(self):
        """Return the set of line indices that lie inside a fenced code block."""
        inside = set(range(self.front_matter))
        tracker = FenceTracker()
        for index, text in enumerate(self.lines):
            if index < self.front_matter:
                continue
            if tracker.feed(text):
                inside.add(index)
                continue
            if tracker.inside:
                inside.add(index)
        return inside

    def headings(self):
        """Yield (line number, level, text) for every ATX heading."""
        for index, text in enumerate(self.lines):
            if index in self.fenced:
                continue
            match = HEADING.match(text)
            if match:
                yield index + 1, len(match.group(1)), match.group(2)

    def prose(self):
        """Yield (line number, text) for the narrative lines of the page.

        Code blocks, headings, tables, and raw HTML are excluded so that the
        metrics measure the prose a reader actually reads.
        """
        for index, text in enumerate(self.lines):
            if index in self.fenced or not text.strip():
                continue
            if HEADING.match(text) or TABLE_ROW.match(text) or HTML_BLOCK.match(text):
                continue
            yield index + 1, text


def check_pictographs(page, found):
    """D001: documentation carries no emoji or other pictographs."""
    for number, text in enumerate(page.lines, start=1):
        for char in sorted(set(text)):
            if emoji(char):
                found.append(
                    Violation(
                        "D001",
                        page.path,
                        number,
                        "pictograph U+%04X; say in words what the symbol stands for"
                        % ord(char),
                    )
                )


def check_punctuation(page, found):
    """D002: documentation uses ASCII punctuation."""
    for number, text in enumerate(page.lines, start=1):
        if number - 1 in page.fenced:
            continue
        for char in sorted(set(text)):
            if char in BANNED_CHARACTERS:
                name, replacement = BANNED_CHARACTERS[char]
                found.append(
                    Violation(
                        "D002",
                        page.path,
                        number,
                        "%s (U+%04X); write %s instead" % (name, ord(char), replacement),
                    )
                )


def check_title(page, found):
    """D003: a page opens with exactly one level one heading."""
    top = [(number, text) for number, level, text in page.headings() if level == 1]
    if not top:
        found.append(
            Violation("D003", page.path, 1, "page has no level one heading")
        )
        return
    if len(top) > 1:
        for number, text in top[1:]:
            found.append(
                Violation(
                    "D003",
                    page.path,
                    number,
                    "second level one heading %r; a page has one title" % text,
                )
            )
    # Leading blank lines and a front matter block are harmless, so the title
    # must be the first content on the page rather than literally line one.
    first = next(
        (
            index + 1
            for index, text in enumerate(page.lines)
            if index >= page.front_matter and text.strip()
        ),
        1,
    )
    if top[0][0] != first:
        found.append(
            Violation(
                "D003",
                page.path,
                top[0][0],
                "the title is not the first content on the page; %d line(s) of "
                "content precede it" % (top[0][0] - first),
            )
        )


def check_heading_depth(page, found, limit):
    """D004: headings go no deeper than the configured level."""
    for number, level, text in page.headings():
        if level > limit:
            found.append(
                Violation(
                    "D004",
                    page.path,
                    number,
                    "heading level %d (%r) is deeper than the limit of %d; "
                    "a page that needs this depth should be split"
                    % (level, text, limit),
                )
            )


def check_task_lists(page, found):
    """D005: documentation is not a checklist."""
    for number, text in page.prose():
        if TASK_LIST.match(text):
            found.append(
                Violation(
                    "D005",
                    page.path,
                    number,
                    "task list checkbox; track work in issues, not in the handbook",
                )
            )


def nav_pages(mkdocs_path):
    """Return the set of pages referenced by the mkdocs nav tree.

    The nav block is read with a small line scanner rather than a YAML parser
    so that the check has no dependency outside the standard library. Only the
    document paths are needed, and every one of them appears as the value of a
    "- Title: path" entry.
    """
    pages = set()
    if not os.path.isfile(mkdocs_path):
        return pages
    with open(mkdocs_path, "r", encoding="utf-8") as handle:
        inside = False
        for line in handle:
            stripped = line.strip()
            if not inside:
                inside = line.startswith("nav:")
                continue
            if line[:1] not in (" ", "\t", "#", "") and stripped:
                break
            if stripped.startswith("#"):
                continue
            match = re.search(r":\s*(\S+\.md)\s*$", stripped)
            if match:
                pages.add(match.group(1))
            elif re.match(r"^-\s+(\S+\.md)$", stripped):
                pages.add(stripped.split(None, 1)[1])
    return pages


def check_nav(page, found, nav, docs_dir):
    """D006: every handbook page is reachable from the mkdocs nav tree."""
    normalized = page.path.replace(os.sep, "/")
    while normalized.startswith("./"):
        normalized = normalized[2:]
    prefix = docs_dir.rstrip("/") + "/"
    if prefix not in normalized:
        return
    relative = normalized.split(prefix, 1)[1]
    if relative.startswith("ford/"):
        return
    if relative not in nav:
        found.append(
            Violation(
                "D006",
                page.path,
                None,
                "page is not in the mkdocs.yml nav tree, so it is published "
                "but unreachable; add an entry for %r" % relative,
            )
        )


def nav_regressions(violations, baseline):
    """Split D006 violations into the ones this change caused and the backlog.

    A page that the baseline navigation already failed to reach was broken
    before this change and belongs to the backlog. A page the baseline did
    reach and this change does not is a regression, and is what the check
    fails on.
    """
    regressions = [v for v in violations if violation_page(v) in baseline]
    return regressions, len(violations) - len(regressions)


def violation_page(violation):
    """Recover the nav-relative page path a D006 violation names."""
    return violation.message.rsplit("'", 2)[-2]


def measure(page):
    """Return the numeric prose metrics for one page.

    bold_per_100_lines counts emphasis spans against the prose line count,
    bullet_fraction and mean_sentence_words are reported for information
    only. Neither is enforced: measured across the documentation at f3e1e57c
    they overlap completely with the current tree, so no threshold on them
    separates hand written prose from generated prose.
    """
    prose = list(page.prose())
    if not prose:
        return {
            "prose_lines": 0,
            "bold_per_100_lines": None,
            "bullet_fraction": None,
            "mean_sentence_words": None,
        }

    # Inline code is neutralized first: identifiers such as `__shared__` and
    # `**kwargs` are literals, not emphasis, and counting them as emphasis
    # would fail a page that contains none.
    bold = sum(len(BOLD.findall(INLINE_CODE.sub(" ", text))) for _, text in prose)
    bullets = sum(1 for _, text in prose if BULLET.match(text))

    words = 0
    sentences = 0
    for _, text in prose:
        clean = INLINE_MATH.sub(" ", text)
        clean = INLINE_CODE.sub(" ", clean)
        clean = LINK.sub(lambda m: m.group(1), clean)
        clean = BULLET.sub("", clean).strip()
        if not clean:
            continue
        words += len(clean.split())
        # A line that ends without terminal punctuation still closes a thought
        # in this style, where paragraphs are written one sentence per line.
        sentences += max(1, len(SENTENCE_END.findall(clean)))

    return {
        "prose_lines": len(prose),
        "bold_per_100_lines": round(100.0 * bold / len(prose), 2),
        "bullet_fraction": round(float(bullets) / len(prose), 3),
        "mean_sentence_words": round(float(words) / sentences, 2) if sentences else None,
    }


def check_metrics(page, metrics, rules, found):
    """D100: the numeric prose threshold."""
    thresholds = rules["docs"]["metrics"]
    if metrics["prose_lines"] < thresholds.get("min_prose_lines", 20):
        # A page too short to measure is left to review rather than to metrics.
        return

    bold = thresholds["D100"]["max"]
    if metrics["bold_per_100_lines"] > bold:
        found.append(
            Violation(
                "D100",
                page.path,
                None,
                "%.2f bold spans per 100 prose lines exceeds the limit of %.2f; "
                "emphasis is rare in SELF documentation, carry the emphasis in "
                "the sentence" % (metrics["bold_per_100_lines"], bold),
            )
        )


def check_file(path, rules, nav, docs_dir):
    """Run every enabled rule against one page, returning violations and metrics."""
    page = MarkdownPage(path)
    found = []
    disabled = set(rules["docs"].get("disabled", ()))
    if "D001" not in disabled:
        check_pictographs(page, found)
    if "D002" not in disabled:
        check_punctuation(page, found)
    if "D003" not in disabled:
        check_title(page, found)
    if "D004" not in disabled:
        check_heading_depth(page, found, rules["docs"]["max_heading_level"])
    if "D005" not in disabled:
        check_task_lists(page, found)
    if "D006" not in disabled:
        check_nav(page, found, nav, docs_dir)
    metrics = measure(page)
    check_metrics(page, metrics, rules, found)
    return [v for v in found if v.rule not in disabled], metrics


def collect(paths, skip):
    """Expand the given paths into the list of markdown pages to check."""
    files = []
    for path in paths:
        if os.path.isdir(path):
            for root, names, filenames in os.walk(path):
                names[:] = [n for n in names if n not in skip]
                files.extend(
                    os.path.join(root, name)
                    for name in sorted(filenames)
                    if name.endswith(".md")
                )
        elif path.endswith(".md") and os.path.isfile(path):
            files.append(path)
    return sorted(files)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", help="files or directories to check")
    parser.add_argument("--files", nargs="*", default=None, help="explicit file list")
    parser.add_argument("--stats", action="store_true", help="report metrics, exit 0")
    parser.add_argument("--rules", default=None, help="path to style-rules.json")
    parser.add_argument("--mkdocs", default="mkdocs.yml", help="path to mkdocs.yml")
    parser.add_argument(
        "--nav-only",
        action="store_true",
        help="check only D006, over every page; use when mkdocs.yml changes",
    )
    parser.add_argument(
        "--baseline-mkdocs",
        default=None,
        help="the mkdocs.yml this change starts from; with --nav-only, only "
        "pages that this change orphans are reported, leaving pages that were "
        "already unreachable to the backlog",
    )
    args = parser.parse_args(argv)

    rules = load_rules(args.rules)
    skip = set(rules["docs"].get("skip_directories", ()))
    skip_files = set(rules["docs"].get("skip_files", ()))

    targets = args.files if args.files is not None else (args.paths or list(DEFAULT_ROOTS))
    files = []
    for path in collect(targets, skip):
        normalized = path.replace(os.sep, "/")
        while normalized.startswith("./"):
            normalized = normalized[2:]
        if any(("/%s/" % name) in normalized for name in skip):
            continue
        if normalized in skip_files:
            continue
        files.append(path)
    if not files:
        print("No documentation pages to check.")
        return 0

    nav = nav_pages(args.mkdocs)
    docs_dir = rules["docs"].get("docs_dir", "docs")

    if args.nav_only:
        # A change to the navigation tree can orphan a page it does not touch,
        # so reachability is checked across every page rather than across the
        # changed ones.
        violations = []
        for path in files:
            check_nav(MarkdownPage(path), violations, nav, docs_dir)
        if args.baseline_mkdocs:
            violations, skipped = nav_regressions(
                violations, nav_pages(args.baseline_mkdocs)
            )
            if skipped:
                print(
                    "%d page(s) were already unreachable before this change and "
                    "are left to the backlog.\n" % skipped
                )
        return report(violations)

    violations = []
    stats = {}
    for path in files:
        found, metrics = check_file(path, rules, nav, docs_dir)
        violations.extend(found)
        stats[path] = metrics

    if args.stats:
        report([], stats=stats)
        return 0
    return report(violations)


if __name__ == "__main__":
    sys.exit(main())
