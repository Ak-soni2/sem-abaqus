"""Prefix every `abaqus ...` line in a run.bat with `call`.

    python _fix_bat_call.py           # show what would change
    python _fix_bat_call.py --apply

On Windows `abaqus` is `abaqus.bat`. Running one batch file from another
WITHOUT `call` transfers control permanently -- the rest of the caller never
executes. So a launcher written as

    abaqus verify -user_exp
    abaqus job=... datacheck
    abaqus job=... interactive

runs the verify, and then simply stops. No error, no message, and the two
commands that actually matter never run. That is exactly what happened:

    D:\\temp\\RUN_SAG\\6um>run.bat
    ... Abaqus/Explicit with user subroutines verification ... result : PASS
    Verification procedure complete

    D:\\temp\\RUN_SAG\\6um>          <- back at the prompt, nothing submitted

`call abaqus ...` runs it as a subroutine and returns, so the script continues
and `errorlevel` is the one the command actually set.

This affects run.bat only. A POSIX shell does not have this behaviour, so
run.sh is left alone.
"""
from __future__ import annotations

import io
import os
import re
import sys

SKIP_DIRS = {".git", "__pycache__", "_colabcheck", "node_modules"}
# `abaqus` at the start of a line (any indent), not already called.
PAT = re.compile(r"^(\s*)(abaqus\s)", re.M | re.I)


def bats():
    for root, dirs, names in os.walk("."):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for n in sorted(names):
            if n.endswith(".bat"):
                p = os.path.join(root, n)
                try:
                    raw = io.open(p, "rb").read()
                    yield p, raw.decode("utf-8", "replace"), raw
                except OSError:
                    continue


def main(argv):
    apply = "--apply" in argv
    hits = []
    for p, s, raw in bats():
        n = len(PAT.findall(s))
        if n:
            hits.append((p, s, raw, n))
    if not hits:
        print("every `abaqus` line in every run.bat is already `call`ed")
        return 0

    total = sum(h[3] for h in hits)
    print("%d .bat file(s), %d uncalled `abaqus` line(s):" % (len(hits), total))
    for p, _, _, n in hits:
        print("   %-58s %d" % (p, n))

    if not apply:
        print()
        print("dry run. Re-run with --apply.")
        return 0

    for p, s, raw, _ in hits:
        # Keep the file's own line endings -- a LF-only .bat misparses.
        crlf = raw.count(b"\r\n")
        lf = raw.count(b"\n") - crlf
        nl = "\r\n" if crlf > lf else "\n"
        io.open(p, "w", encoding="utf-8", newline=nl).write(
            PAT.sub(r"\1call \2", s))
    print()
    print("rewrote %d file(s)" % len(hits))

    left = [q for q, body, _ in bats() if PAT.search(body)]
    if left:
        print("STILL UNCALLED in: %s" % ", ".join(left))
        return 1
    print("no uncalled `abaqus` lines remain")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
