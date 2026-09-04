"""Check every run.bat / run.sh before anyone waits in a queue to find out.

    python verify_launchers.py

A launcher is the one file in a deck package that nothing else validates. The
deck itself is checked by four independent gates; the two-line script that
submits it was checked by nobody, and a wrong flag there wastes the whole trip
to the machine.

That is not hypothetical. Every launcher in this project carried

    abaqus verify -user_explicit

which is not an Abaqus option and never was, so `run.bat` aborted before
submitting anything:

    Abaqus Error: Command line option "user_explicit" may not be used with
    "verify". Valid options for use with "verify" are:
    [... 'user_std', 'user_exp', ...]

46 files, every package, both notebooks. One typo, propagated by the
generators, and no gate looked at it.

WHAT IS CHECKED
---------------
* Only real ``abaqus verify`` options are used, against the list Abaqus itself
  prints, and a VUMAT deck must verify ``user_exp`` -- not ``exp``, which
  checks the solver rather than the Fortran toolchain.
* ``double=both`` on every solve. It is not a preference here: h and dc are
  compared at tens of nanometres against a millimetre geometry, and single
  precision does not have the digits. The failure is silent.
* The subroutine named actually exists NEXT TO the launcher, because these
  folders are meant to be copied to a work directory and a relative path out
  of the folder breaks the moment they are.
* The deck named exists, or is documented as regenerable.
* A datacheck runs before the solve, and the script stops if it fails --
  `errorlevel` on Windows, `||` or `set -e` on POSIX.
* run.bat is CRLF. A LF-only .bat misparses on some Windows setups.
* Every ``abaqus`` line in a run.bat is ``call``ed. On Windows ``abaqus`` is
  ``abaqus.bat``, and running one batch file from another WITHOUT ``call``
  transfers control permanently -- the caller never resumes. This bit the
  project immediately after the flag fix above: the verify ran, printed
  ``result : PASS``, and the script simply ended, having submitted nothing
  and printed no error. 57 lines across all 18 .bat files.
"""

from __future__ import annotations

import io
import os
import re
import sys

# What `abaqus verify` actually accepts, from the error message Abaqus prints
# when given something else.
VERIFY_OPTS = {
    "all", "install", "std", "exp", "foundation", "user_std", "user_exp",
    "param", "ams", "design", "moldflow", "tosca", "cae", "viewer", "noGUI",
    "adams", "dcatiav5", "parasolid", "proe", "swi", "docUrl", "cPerf",
    "ioPerf", "make", "parallel", "scripting", "noComp", "retainFiles",
    "noProd", "log", "noLic", "noRestrict", "help", "verbose",
    "remoteRunDirectory",
}

SKIP_DIRS = {".git", "__pycache__", "_colabcheck", "node_modules"}

FAIL = []
PASS = 0


def chk(what, ok, detail=""):
    global PASS
    if ok:
        PASS += 1
    else:
        FAIL.append("%s%s" % (what, ("  -- " + detail) if detail else ""))
    return ok


def launchers():
    for root, dirs, names in os.walk("."):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for n in sorted(names):
            if n in ("run.bat", "run.sh"):
                yield os.path.join(root, n)


def check(path):
    raw = io.open(path, "rb").read()
    s = raw.decode("utf-8", "replace")
    d = os.path.dirname(path)
    rel = path.replace(os.sep, "/").lstrip("./")
    bat = path.endswith(".bat")
    bad0 = len(FAIL)

    # --- verify options -------------------------------------------------
    for m in re.finditer(r"abaqus\s+verify\s+((?:-\w+\s*)+)", s):
        for opt in re.findall(r"-(\w+)", m.group(1)):
            chk("%s: 'abaqus verify -%s' is not a real option" % (rel, opt),
                opt in VERIFY_OPTS,
                "Abaqus accepts %s" % ", ".join(sorted(VERIFY_OPTS))[:90])

    uses_user = "user=" in s
    if uses_user:
        verifies = re.findall(r"abaqus\s+verify\s+-(\w+)", s)
        chk("%s: a user-subroutine deck verifies user_exp" % rel,
            "user_exp" in verifies or not verifies,
            "found %r -- 'exp' checks the solver, 'user_exp' checks that "
            "Abaqus can BUILD a subroutine, which is the thing that fails"
            % (verifies or "no verify step"))

    # --- every solve must be double precision ---------------------------
    jobs = re.findall(r"abaqus\s+(?:job=|-job\s)[^\n]*", s)
    for j in jobs:
        if "datacheck" in j:
            continue
        chk("%s: 'double=both' on every solve" % rel, "double=both" in j,
            "h and dc are compared at nanometres on a millimetre geometry; "
            "single precision has ~7 digits and the failure is SILENT")

    # --- the subroutine must be beside the launcher ---------------------
    for sub in set(re.findall(r"user=([^\s]+)", s)):
        sub = sub.strip('"')
        chk("%s: the subroutine %s is in the folder" % (rel, sub),
            os.path.exists(os.path.join(d, sub)),
            "these folders get copied to a work directory, so a path out of "
            "the folder breaks")

    # --- the deck must exist, or be documented as regenerable -----------
    for inp in set(re.findall(r"input=([^\s]+)", s)):
        inp = inp.strip('"')
        here = os.path.exists(os.path.join(d, inp))
        readme = os.path.join(d, "README.md")
        excused = (os.path.exists(readme)
                   and re.search(r"rebuild|regener|gitignor",
                                 io.open(readme, encoding="utf-8",
                                         errors="replace").read(), re.I))
        chk("%s: the deck %s is present or documented as regenerable"
            % (rel, inp), here or bool(excused),
            "neither on disk nor explained in a README")

    # --- a datacheck, and a stop if it fails ----------------------------
    if jobs:
        has_dc = any("datacheck" in j for j in jobs)
        chk("%s: a datacheck runs before the solve" % rel, has_dc,
            "seconds, and it catches the keyword and material-card errors "
            "that otherwise surface after the queue")
        if has_dc:
            guards = (re.search(r"if\s+errorlevel\s+1", s, re.I) if bat
                      else re.search(r"\|\||set\s+-e", s))
            chk("%s: it stops if the datacheck fails" % rel, bool(guards),
                "otherwise the solve starts anyway and the datacheck was "
                "decoration")

    # --- `call` on every abaqus line, in .bat only ----------------------
    if bat:
        uncalled = re.findall(r"(?m)^\s*(abaqus\s+\S+)", s)
        chk("%s: every abaqus line is `call`ed" % rel, not uncalled,
            "%s -- on Windows abaqus is abaqus.bat, and one .bat running "
            "another without `call` never returns, so everything after the "
            "first abaqus line silently never runs"
            % (uncalled[:2] if uncalled else ""))

    # --- line endings ----------------------------------------------------
    if bat:
        crlf = raw.count(b"\r\n")
        chk("%s: CRLF line endings" % rel,
            crlf > 0 and raw.count(b"\n") == crlf,
            "a LF-only .bat misparses on some Windows setups")

    return len(FAIL) == bad0


def main(argv):
    files = list(launchers())
    if not files:
        print("no run.bat / run.sh found")
        return 1
    print("checking %d launcher(s)" % len(files))
    print()
    ok_files = 0
    for f in files:
        if check(f):
            ok_files += 1
    print("  %d of %d launchers clean" % (ok_files, len(files)))
    print()
    print("=" * 74)
    if FAIL:
        print("  %d check(s) passed, %d FAILED" % (PASS, len(FAIL)))
        for f in FAIL:
            print("    - %s" % f)
        return 1
    print("  ALL %d CHECKS PASSED across %d launchers" % (PASS, len(files)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
