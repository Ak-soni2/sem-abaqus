"""Give the datacheck its own job name, so it cannot lock out the solve.

    python _fix_datacheck_lock.py           # show what would change
    python _fix_datacheck_lock.py --apply

A datacheck is a real Abaqus submission. It writes ``<job>.lck`` and LEAVES it,
because the job never completed -- there is nothing to complete. A solve
submitted under the same job name then finds that lock and refuses:

    Abaqus Error: Detected lock file micro_15um.lck. Please confirm that no
    other applications are attempting to write to the output database
    associated with this job before removing the lock file and resubmitting.

on a job that has never been run. Every launcher in this project used one job
name for both, so every one of them failed on its first use.

The fix is to submit the datacheck as ``<job>_check``. That is better than
deleting the lock between the two, which is the other obvious option: a lock
file is sometimes REAL -- another process genuinely is writing that .odb -- and
a launcher that unconditionally removes it would clobber a live job. Renaming
never creates the conflict at all, and it keeps the datacheck's own artefacts
under a name nobody mistakes for the real run.
"""
from __future__ import annotations

import io
import os
import re
import sys

SKIP_DIRS = {".git", "__pycache__", "_colabcheck", "node_modules"}
# `abaqus job=NAME ... datacheck`, where NAME does not already end in _check
PAT = re.compile(r"(abaqus\s+job=)([\w.\-]+?)(?<!_check)(\s+[^\n]*?datacheck)")


def launchers():
    for root, dirs, names in os.walk("."):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for n in sorted(names):
            if n in ("run.bat", "run.sh"):
                p = os.path.join(root, n)
                raw = io.open(p, "rb").read()
                yield p, raw.decode("utf-8", "replace"), raw


def main(argv):
    apply = "--apply" in argv
    hits = [(p, s, raw) for p, s, raw in launchers() if PAT.search(s)]
    if not hits:
        print("every datacheck already runs under its own job name")
        return 0

    print("%d launcher(s) share a job name between datacheck and solve:"
          % len(hits))
    for p, s, _ in hits:
        for m in PAT.finditer(s):
            print("   %-52s job=%s -> job=%s_check"
                  % (p, m.group(2), m.group(2)))

    if not apply:
        print()
        print("dry run. Re-run with --apply.")
        return 0

    for p, s, raw in hits:
        crlf = raw.count(b"\r\n")
        lf = raw.count(b"\n") - crlf
        nl = "\r\n" if crlf > lf else "\n"
        out = PAT.sub(r"\1\2_check\3", s)
        # point any "read <job>.dat" advice at the file that will exist
        for m in PAT.finditer(s):
            job = m.group(2)
            out = out.replace("%s.dat" % job, "%s_check.dat" % job)
        io.open(p, "w", encoding="utf-8", newline=nl).write(out)
    print()
    print("rewrote %d launcher(s)" % len(hits))

    left = [p for p, s, _ in launchers() if PAT.search(s)]
    if left:
        print("STILL SHARING a job name: %s" % ", ".join(left))
        return 1
    print("no launcher shares a job name between datacheck and solve")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
