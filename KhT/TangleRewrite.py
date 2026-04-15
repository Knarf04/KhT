# Reidemeister slice-rewriting scaffolding for OPEN_QUESTIONS.md item 6.
#
# Some tangles (notably TH-pattern) never pass through (1, 3) mid-
# computation because of how their slice sequence was written down.  The
# same tangle up to ambient isotopy may admit a Reidemeister-equivalent
# slicing that does reach (1, 3), which would let the existing
# intermediate-cleanup infrastructure fire.
#
# This module is a stub.  The rewriter currently returns the input
# unchanged.  The intended future implementation would apply
# Reidemeister I/II moves to the slice string to hunt for width
# reductions — that's a research-grade question in itself (no
# known complete heuristic), so the library just exposes the hook.
#
# Usage:
#     from TangleRewrite import rewrite_slices
#     candidates = rewrite_slices(slice_str, target_checkpoints=[3])
#     # candidates[0] is always the original.


def rewrite_slices(slice_str, target_checkpoints=(3,)):
    """Return a list of slice strings that are Reidemeister-equivalent
    to ``slice_str`` and (ideally) pass through widths listed in
    ``target_checkpoints``.

    Current implementation: identity.  The input is returned as the
    sole candidate.  See OPEN_QUESTIONS.md item 6 — populating this with
    real Reidemeister moves is a math deliverable.
    """
    return [slice_str]


def slice_widths(slice_str, start=1):
    """Return the running tangle-width after each slice of ``slice_str``,
    where width = number of ends at the bottom of the tangle after that
    slice.  Cheap helper for heuristics that want to decide whether a
    rewrite has achieved a (1, 3) checkpoint.

    Convention: ``slice_str`` is the left-to-right slice-from-bottom
    format used by BNbracket; ``start`` is the number of ends at the
    top (default 1, matching the standard (1, *)-tangle setting).
    """
    width = start
    widths = []
    for word in slice_str.split('.'):
        op = word[:3]
        if op == "cup":
            width += 2
        elif op == "cap":
            width -= 2
        # pos / neg preserve width.
        widths.append(width)
    return widths
