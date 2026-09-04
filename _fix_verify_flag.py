"""Replace the invalid `abaqus verify -user_explicit` with `-user_exp`.

    python _fix_verify_flag.py           # show what would change
    python _fix_verify_flag.py --apply

The wrong flag is not an Abaqus verify option and never was. Running it aborts
the launcher before anything is submitted:

    Abaqus Error: Command line option "user_explicit" may not be used with
    "verify". Valid options for use with "verify" are:
    [... 'user_std', 'user_exp', ...]

A VUMAT is an EXPLICIT user subroutine, so the flag is `user_exp`; `user_std`
is the Standard/implicit equivalent and `exp` on its own verifies the solver
rather than the Fortran toolchain.

This was project-wide -- every launcher in RUN_ME, RUN_ME_SIC, PROBE_* and
RUN_SAG, every generator that writes one, both notebooks and three READMEs, 46
files -- so it is fixed at the generators AND in the already-written files.

Note the token below is assembled from two halves rather than written out. This
script sweeps the whole tree, and on its first run it rewrote its OWN docstring
and constant, leaving a file that claimed to replace a string with itself.
"""
from __future__ import annotations

import io
import os
import sys

# Assembled, so a sweep of this tree cannot rewrite the thing being
# searched for -- which is exactly what happened the first time.
BAD = "user_" + "explicit"
GOOD = "user_exp"

SKIP_DIRS = {".git", "__pycache__", "_colabcheck", "node_modules"}
# Binary-ish or generated things a text swap must not touch.
SKIP_EXT = {".glb", ".png", ".jpg", ".jpeg", ".tif", ".pkl", ".odb", ".gz",
            ".tar", ".zip", ".pyc", ".inp"}


def files():
    for root, dirs, names in os.walk("."):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for n in names:
            if os.path.splitext(n)[1].lower() in SKIP_EXT:
                continue
            p = os.path.join(root, n)
            try:
                if os.path.getsize(p) > 40_000_000:
                    continue
                s = io.open(p, encoding="utf-8", errors="strict").read()
            except (UnicodeDecodeError, OSError):
                continue
            # Skip this script: its docstring quotes the Abaqus error
            # message verbatim, which is documentation, not a launcher.
            if os.path.abspath(p) == os.path.abspath(__file__):
                continue
            if BAD in s:
                yield p, s


def main(argv):
    apply = "--apply" in argv
    hits = list(files())
    if not hits:
        print("no occurrences of %r left" % BAD)
        return 0

    print("%d file(s) contain %r:" % (len(hits), BAD))
    by_kind = {}
    for p, s in hits:
        ext = os.path.splitext(p)[1] or "(none)"
        by_kind[ext] = by_kind.get(ext, 0) + 1
        print("   %-58s %d occurrence(s)" % (p, s.count(BAD)))
    print()
    print("by type: %s" % ", ".join("%s %d" % (k, v)
                                    for k, v in sorted(by_kind.items())))

    if not apply:
        print()
        print("dry run. Re-run with --apply to make the change.")
        return 0

    # Preserve each file's own line endings: run.bat is CRLF on purpose, and
    # rewriting it as LF would break `cmd` on some setups.
    changed = 0
    for p, s in hits:
        raw = io.open(p, "rb").read()
        crlf = raw.count(b"\r\n")
        lf = raw.count(b"\n") - crlf
        nl = "\r\n" if crlf > lf else "\n"
        io.open(p, "w", encoding="utf-8", newline=nl).write(s.replace(BAD,
                                                                      GOOD))
        changed += 1
    print()
    print("rewrote %d file(s)" % changed)

    left = list(files())
    if left:
        print("STILL PRESENT in %d file(s) -- the swap did not take:" % len(left))
        for p, _ in left:
            print("   %s" % p)
        return 1
    print("no occurrences of %r remain" % BAD)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
