"""Verification gate for vumat_grind.for, the hybrid ductile/brittle VUMAT.

Compiles the subroutine with a single-material-point driver and checks it
against three independent things: closed-form algebra, the published JH-2
benchmarks, and the original ``vumat_jh2.for`` run through the same driver.

    python verify_vumat_grind.py            # run everything
    python verify_vumat_grind.py -v         # print every check

Exits non-zero on the first failed check, so it can gate a build.

Why a separate gate. ``verify_all.py`` and the two deck verifiers check the
geometry half of this project; none of them can see whether the constitutive
law is right. The strongest check here is B1: with PROPS(56) = 3 the hybrid
routine must reproduce ``vumat_jh2.for`` to the last bit on identical input.
That is what lets the JH-2 branch inherit the validation the original already
carries, instead of claiming it.
"""

from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "_hybrid_test")
SRC_GRIND = os.path.join(HERE, "vumat_grind.for")
SRC_JH2 = os.path.join(HERE, "vumat_jh2.for")

# How many constants vumat_grind.for reads. Kept here so the grind2 gate can
# slice a grind2 card back to a grind card without hard-coding the number in
# two files that must agree.
N_PROPS_GRIND = 56

VERBOSE = "-v" in sys.argv or "--verbose" in sys.argv

_PASS = 0
_FAIL: list[str] = []


# --------------------------------------------------------------------------
# harness
# --------------------------------------------------------------------------

def find_gfortran() -> str:
    exe = shutil.which("gfortran")
    if exe:
        return exe
    root = os.path.expandvars(
        r"%LOCALAPPDATA%\Microsoft\WinGet\Packages")
    for dirpath, _dirs, files in os.walk(root) if os.path.isdir(root) else []:
        if "gfortran.exe" in files:
            return os.path.join(dirpath, "gfortran.exe")
    raise SystemExit(
        "gfortran not found. Install one (winget install "
        "BrechtSanders.WinLibs.POSIX.UCRT) or put it on PATH; without a "
        "compiler this file can check the source text but not the physics.")


def check(name: str, ok: bool, detail: str = "") -> None:
    global _PASS
    if ok:
        _PASS += 1
        if VERBOSE:
            print("  PASS  %-58s %s" % (name, detail))
    else:
        _FAIL.append(name)
        print("  FAIL  %-58s %s" % (name, detail))


def close(name: str, got: float, want: float, rtol: float = 1e-9,
          atol: float = 0.0) -> None:
    err = abs(got - want)
    tol = atol + rtol * abs(want)
    check(name, err <= tol, "got %.12g want %.12g (err %.3g, tol %.3g)"
          % (got, want, err, tol))


class Driver:
    """One compiled VUMAT plus the material-point driver."""

    def __init__(self, fc: str, source: str, exe: str):
        self.exe = os.path.join(WORK, exe)
        local = os.path.join(WORK, os.path.basename(exe).replace(".exe", ".f"))
        shutil.copyfile(source, local)
        # Static: the driver has to run from anywhere without the compiler's
        # runtime DLLs being on PATH, or this gate passes on the machine that
        # built it and fails everywhere else.
        cmd = [fc, "-std=legacy", "-ffixed-form", "-ffixed-line-length-72",
               "-Wall", "-Wconversion", "-O2", "-static", "-I", WORK,
               "-o", self.exe, os.path.join(WORK, "driver.f"), local]
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=WORK)
        if r.returncode != 0:
            raise SystemExit("compile failed for %s:\n%s" % (source, r.stderr))
        self.warnings = [
            ln for ln in r.stderr.splitlines()
            if "Warning:" in ln
            and "Unused dummy argument" not in ln
            and "Unused parameter" not in ln
        ]

    def run(self, props, segs, *, nstatev=20, fields=(), coord=(0.0, 0.0, 0.0),
            dt=1e-8, nout=1, charlen=1.0e-3, density=2.35e-9):
        """segs = [(nstep, (de11, de22, de33, de12, de23, de13)), ...]"""
        lines = [str(len(props)),
                 " ".join("%.17g" % v for v in props),
                 "%d %d" % (nstatev, len(fields))]
        if fields:
            lines.append(" ".join("%.17g" % v for v in fields))
        lines.append(" ".join("%.17g" % v for v in coord))
        lines.append("%.17g %.17g" % (charlen, density))
        lines.append("%.17g %d %d" % (dt, len(segs), nout))
        for n, de in segs:
            lines.append("%d %s" % (n, " ".join("%.17g" % v for v in de)))
        r = subprocess.run([self.exe], input="\n".join(lines) + "\n",
                           capture_output=True, text=True, cwd=WORK)
        if r.returncode != 0:
            raise SystemExit("driver failed:\n%s\n%s" % (r.stdout[-2000:],
                                                         r.stderr[-2000:]))
        out = []
        for ln in r.stdout.split("\n"):
            f = ln.split()
            if len(f) < 8:
                continue
            v = [float(x) for x in f[1:]]
            out.append({"t": v[0], "s": v[1:7], "sdv": v[7:7 + nstatev],
                        "eint": v[7 + nstatev], "einel": v[8 + nstatev]})
        return out


# --------------------------------------------------------------------------
# material cards
# --------------------------------------------------------------------------

# Baranowski et al. sandstone, exactly the set vumat_jh2.for falls back to.
JH2 = [3735.6, 2686.0, 1982.0, 1374.0, 8.0, 0.71, 0.30, 0.022,
       0.55, 0.40, 1.0, 0.002, 1.20, 9000.0, 22000.0, 0.25, 912.0]

# Johnson-Cook for Ti-6Al-4V (Yadav et al. 2022, Table 4) in mm-MPa-tonne-s,
# with the damage set from Johnson & Cook / Lesuer. Used here to exercise the
# algebra, not because titanium is what gets ground.
JC = dict(A=1098.0, B=1092.0, n=0.93, C=0.014, m=1.1, edot0=1.0,
          E=113000.0, nu=0.34, rho=4.43e-9, cp=5.26e8, betaq=0.9,
          T0=293.15, Tm=1933.0,
          b=2.95e-7, M=3.0, alpha=0.3, lam=1.0, rprime=2.0, gsge=0.0,
          D1=-0.09, D2=0.25, D3=-0.5, D4=0.014, D5=3.87, dcrit=1.0)


def card(*, jh2=None, jc=None, dcut=0.0, lamc=0.15, hard=0.0, kic=0.0,
         idcf=1, thc=0.0, h0=0.0, hg=0.0, rtip=0.0, ihmode=0,
         edot0=1.0, edmin=1.0, itcut=1.0, fsmax=-1.0):
    j = list(jh2 if jh2 is not None else JH2)
    c = dict(JC if jc is None else jc)
    p = j + [edot0, edmin, itcut, fsmax]
    p += [c["A"], c["B"], c["n"], c["C"], c["m"], c["edot0"], c["E"], c["nu"],
          c["rho"], c["cp"], c["betaq"], c["T0"], c["Tm"], c["b"], c["M"],
          c["alpha"], c["lam"], c["rprime"], c["gsge"]]
    p += [c["D1"], c["D2"], c["D3"], c["D4"], c["D5"], c["dcrit"]]
    p += [dcut, lamc, hard, kic, float(idcf), thc, h0, hg, rtip, float(ihmode)]
    assert len(p) == 56, len(p)
    return p


def sge_flow(ep, dep, dt, hlen, thom, c):
    """The closed form the Fortran must reproduce. Written from the papers,
    not from the Fortran, so agreement means something."""
    epn = max(ep + dep, 0.0)
    hard = c["A"] + c["B"] * epn ** c["n"] if epn > 0 else c["A"]
    edr = max((dep / dt) / c["edot0"], 1.0) if dt > 0 else 1.0
    fthm = max(1.0 - min(max(thom, 0.0), 1.0) ** c["m"], 1e-16)
    sjc = max(hard * (1.0 + c["C"] * math.log(edr)) * fthm, 1e-16)
    g = c["gsge"] if c["gsge"] > 0 else c["E"] / (2.0 * (1.0 + c["nu"]))
    hl = max(hlen, c["b"])
    eta = 4.0 * epn / hl
    sgec = c["rprime"] * c["b"] * (c["M"] * c["alpha"] * g) ** 2
    arg = sgec * eta / (sjc * sjc)
    fsge = math.sqrt(1.0 + arg ** c["lam"]) if arg > 0 else 1.0
    return sjc * fsge, sjc, fsge, eta


def dc_formula(lamc, H, E, Kc, form):
    if form == 2:
        return lamc * (E / H) * (Kc / H) ** 2
    return lamc * math.sqrt(H / E) * (Kc / H) ** 2


# --------------------------------------------------------------------------
# A. source hygiene
# --------------------------------------------------------------------------

def test_source(grind: Driver) -> None:
    print("A. source hygiene")
    text = open(SRC_GRIND, encoding="ascii").read()
    lines = text.split("\n")

    long_code = []
    bad_cont = []
    for i, ln in enumerate(lines, start=1):
        if not ln.strip():
            continue
        if ln[0] in "Cc*!":
            continue
        # Fixed form ignores everything past column 72 SILENTLY. gfortran
        # truncates without complaint and so does ifort, so a 73-column line
        # is a bug that compiles and then computes the wrong number.
        if len(ln.rstrip()) > 72:
            long_code.append(i)
        if len(ln) >= 6 and ln[5] not in (" ", "0"):
            if ln[:5].strip():
                bad_cont.append(i)
    check("no code line exceeds column 72", not long_code, str(long_code[:8]))
    check("continuation lines have columns 1-5 blank", not bad_cont,
          str(bad_cont[:8]))
    check("ascii only, unix line endings", "\r" not in text
          and all(ord(ch) < 128 for ch in text))
    check("compiles with no conversion/typing warnings",
          not grind.warnings, "; ".join(grind.warnings[:3]))
    # The implicit typing trap that already bit once: any variable whose name
    # starts i-n is an INTEGER under vaba_param.inc, so a real constant read
    # into one is silently truncated.
    import re
    bad = set()
    for m in re.finditer(r"^ {6,}([a-z][a-z0-9_]*)\s*=\s*props\(", text,
                         re.MULTILINE | re.IGNORECASE):
        if m.group(1)[0].lower() in "ijklmn":
            bad.add(m.group(1))
    check("no real prop is read into an implicitly-integer name", not bad,
          str(sorted(bad)))


# --------------------------------------------------------------------------
# B. the JH-2 branch
# --------------------------------------------------------------------------

def jh2_histories():
    """(name, props tail, segments) shared by the hybrid and the original."""
    de = 1.0e-5
    return [
        ("uniaxial strain compression 5% then unload",
         [(5000, (-de, 0, 0, 0, 0, 0)), (5000, (de, 0, 0, 0, 0, 0))]),
        ("uniaxial strain tension 1%",
         [(1000, (de, 0, 0, 0, 0, 0))]),
        ("hydrostatic compression 3%",
         [(3000, (-de, -de, -de, 0, 0, 0))]),
        ("pure shear to 2%",
         [(2000, (0, 0, 0, de, 0, 0))]),
        ("triaxial: hydrostatic then deviatoric",
         [(1000, (-de, -de, -de, 0, 0, 0)),
          (4000, (-de, 0.5 * de, 0.5 * de, 0, 0, 0))]),
    ]


def test_jh2_identity(grind: Driver, jh2: Driver) -> None:
    print("B. JH-2 branch is the original, bit for bit")
    for name, segs in jh2_histories():
        # PROPS(56) = 3 forces the brittle branch everywhere.
        p_h = card(ihmode=3)
        p_j = JH2 + [1.0, 1.0, 1.0, -1.0]          # the original reads 1..21
        a = grind.run(p_h, segs, nstatev=20, nout=250)
        b = jh2.run(p_j, segs, nstatev=12, nout=250)
        check("same number of output rows: " + name, len(a) == len(b),
              "%d vs %d" % (len(a), len(b)))
        if len(a) != len(b):
            continue
        worst_s = 0.0
        worst_v = 0.0
        for ra, rb in zip(a, b):
            for x, y in zip(ra["s"], rb["s"]):
                worst_s = max(worst_s, abs(x - y))
            for x, y in zip(ra["sdv"][:12], rb["sdv"][:12]):
                worst_v = max(worst_v, abs(x - y))
        check("stress identical to vumat_jh2.for: " + name, worst_s == 0.0,
              "max |dsigma| = %.3g" % worst_s)
        check("SDV 1-12 identical to vumat_jh2.for: " + name, worst_v == 0.0,
              "max |dSDV| = %.3g" % worst_v)


def test_jh2_benchmarks(grind: Driver) -> None:
    print("B. JH-2 published benchmarks")
    a, b, cc, rn, rm = JH2[5], JH2[6], JH2[7], JH2[8], JH2[9]
    phel, rt, sfmax, sighel = JH2[3], JH2[4], JH2[15], JH2[16]
    tstar = rt / phel

    # The strength surfaces themselves, read out of the running code at a
    # spread of pressures. A hydrostatic path sets P without needing any
    # stress control, and SDV6/SDV7 are the normalised intact and fractured
    # strengths the routine actually used.
    segs = [(4000, (-1e-5, -1e-5, -1e-5, 0, 0, 0))]
    r = grind.run(card(ihmode=3), segs, nout=50, dt=1.0)
    worst_i = worst_f = 0.0
    for row in r:
        p = row["sdv"][2]
        pit = max(p / phel + tstar, 1e-16)
        pf = max(p / phel, 0.0)
        wi = a * pit ** rn                      # rate factor is 1: eds -> EDMIN
        wf = min(b * pf ** rm, sfmax)
        worst_i = max(worst_i, abs(row["sdv"][5] - wi))
        worst_f = max(worst_f, abs(row["sdv"][6] - wf))
    check("intact surface equals A(P*+T*)^N", worst_i < 1e-12,
          "max err %.3g over P up to %.1f MPa" % (worst_i, r[-1]["sdv"][2]))
    check("fractured surface equals min(B P*^M, SFMAX)", worst_f < 1e-12,
          "max err %.3g" % worst_f)

    # The header's quasi-static uniaxial compressive strength, 90.000 MPa, is
    # where the uniaxial-stress elastic path meets the intact surface: at an
    # axial stress s, P = s/3 and q = s, so s = A((s/3 + T)/PHEL)^N * SIGHEL.
    # Solve that here and confirm the routine puts its surface in the same
    # place at the same pressure -- which is the published claim, without
    # needing a stress-controlled driver to walk there.
    def uniaxial_root(sign):
        lo, hi = 1e-6, 1.0e4
        for _ in range(200):
            s = 0.5 * (lo + hi)
            p = sign * s / 3.0
            q = a * max(p / phel + tstar, 1e-30) ** rn * sighel
            if q > s:
                lo = s
            else:
                hi = s
        return 0.5 * (lo + hi)

    close("closed-form uniaxial compressive strength", uniaxial_root(+1.0),
          90.0, rtol=2e-4)
    close("closed-form uniaxial tensile strength", uniaxial_root(-1.0),
          17.932, rtol=2e-4)
    for conf, want in ((10.0, 119.886), (17.0, 138.380), (25.0, 158.049)):
        # Triaxial: P = (s + 2*conf)/3 with q = s - conf.
        lo, hi = 1e-6, 1.0e4
        for _ in range(200):
            s = 0.5 * (lo + hi)
            p = (s + 2.0 * conf) / 3.0
            q = a * max(p / phel + tstar, 1e-30) ** rn * sighel
            if q > s - conf:
                lo = s
            else:
                hi = s
        # The paper quotes the AXIAL stress at failure, not the deviator.
        close("closed-form triaxial strength at %.0f MPa confinement" % conf,
              0.5 * (lo + hi), want, rtol=2e-4)

    # And the surface the code evaluates at the uniaxial-compression pressure
    # must be that same 90 MPa.
    p_at = 90.0 / 3.0
    idx = min(range(len(r)), key=lambda i: abs(r[i]["sdv"][2] - p_at))
    p_got = r[idx]["sdv"][2]
    q_code = r[idx]["sdv"][5] * sighel
    q_want = a * (p_got / phel + tstar) ** rn * sighel
    close("q_limit from the running code at P = s/3", q_code, q_want,
          rtol=1e-12)
    close("that limit is the published 90 MPa", q_code, 90.0, rtol=3e-2)

    # Fully fractured residual strength, triaxial. The reference paper reports
    # ~60% of peak for the single-element triaxial test.
    segs = [(2000, (-1e-5, -1e-5, -1e-5, 0, 0, 0)),
            (60000, (-1e-5, 5e-6, 5e-6, 0, 0, 0))]
    r = grind.run(card(ihmode=3), segs, nout=100, dt=1.0)
    qs = [row["sdv"][3] for row in r]
    dmg = [row["sdv"][0] for row in r]
    peak = max(qs)
    resid = qs[-1]
    check("damage reaches 1 in triaxial compression", max(dmg) >= 0.999,
          "D_max = %.6f" % max(dmg))
    check("fully fractured residual is 40-80%% of peak",
          0.40 <= resid / peak <= 0.80,
          "resid/peak = %.3f (%.4g / %.4g)" % (resid / peak, resid, peak))

    # Bulking: pressure must exceed the pure-EOS pressure once damage grows,
    # and only in compression.
    dp = [row["sdv"][10] for row in r]
    check("bulking pressure is non-negative and non-decreasing",
          all(b >= -1e-12 for b in dp)
          and all(dp[i + 1] >= dp[i] - 1e-12 for i in range(len(dp) - 1)),
          "max deltaP = %.4g MPa" % max(dp))
    r2 = grind.run(card(ihmode=3), [(3000, (1e-5, 0, 0, 0, 0, 0))],
                   nout=100, dt=1.0)
    check("no bulking in tension",
          max(row["sdv"][10] for row in r2) <= 1e-12,
          "max deltaP = %.3g" % max(row["sdv"][10] for row in r2))


# --------------------------------------------------------------------------
# C. the ductile JC + SGE branch
# --------------------------------------------------------------------------

def test_ductile(grind: Driver) -> None:
    print("C. ductile branch: Johnson-Cook + strain-gradient enhancement")
    c = dict(JC)
    g = c["E"] / (2.0 * (1.0 + c["nu"]))
    k = c["E"] / (3.0 * (1.0 - 2.0 * c["nu"]))

    # -- C1 elasticity, before yield -------------------------------------
    de = 1.0e-7
    n = 200
    r = grind.run(card(ihmode=2), [(n, (de, 0, 0, 0, 0, 0))], nout=n, dt=1e-9)
    eps = n * de
    close("uniaxial-strain s11 = (K+4G/3) eps", r[-1]["s"][0],
          (k + 4.0 * g / 3.0) * eps, rtol=1e-12)
    close("uniaxial-strain s22 = (K-2G/3) eps", r[-1]["s"][1],
          (k - 2.0 * g / 3.0) * eps, rtol=1e-12)
    close("no plastic strain while elastic", r[-1]["sdv"][1], 0.0, atol=0.0)

    # -- C2 yield point in pure shear ------------------------------------
    # q = sqrt(3)|s12|, so first yield is at s12 = A/sqrt(3) with eta = 0.
    dg = 1.0e-5
    r = grind.run(card(ihmode=2), [(3000, (0, 0, 0, dg, 0, 0))], nout=1,
                  dt=1e-9)
    first = next(i for i, row in enumerate(r) if row["sdv"][1] > 0)
    s12_at_yield = r[first - 1]["s"][3]
    close("pure-shear yield at s12 = A/sqrt(3)", s12_at_yield * math.sqrt(3.0),
          c["A"], rtol=2e-3)
    close("pressure stays zero in pure shear", r[-1]["sdv"][2], 0.0, atol=1e-9)

    # -- C3 the flow stress the return map lands on ----------------------
    # After each increment q must sit exactly on sigma_e(ep, dep/dt, T), with
    # ep and T taken at the START of that increment -- which is what the
    # routine uses, and what the previous output row holds.
    def start_state(rows, i):
        if i == 0:
            return 0.0, c["T0"]
        return rows[i - 1]["sdv"][1], rows[i - 1]["sdv"][15]

    for hloc, tag in ((1.0e-2, "h = 10 um, SGE negligible"),
                      (1.0e-4, "h = 0.1 um, SGE active"),
                      (2.0e-6, "h = 2 nm, SGE strong")):
        p = card(ihmode=2, h0=hloc)
        r = grind.run(p, [(3000, (0, 0, 0, dg, 0, 0))], nout=1, dt=1e-9)
        worst = 0.0
        worst_f = 0.0
        worst_q = 0.0
        for i, row in enumerate(r):
            dep = row["sdv"][8]
            if dep <= 0:
                continue
            ep0, t0 = start_state(r, i)
            thom = (t0 - c["T0"]) / (c["Tm"] - c["T0"])
            # The cap, rebuilt from this card's own JH-2 constants at the
            # pressure and rate the routine actually saw (SDV3, SDV5). Derived
            # independently of the Fortran, so this now checks the cap as well
            # as the SGE algebra.
            sy, sjc, fsge, eta = sge_flow(ep0, dep, 1e-9, hloc,
                                          thom, c)
            worst = max(worst, abs(row["sdv"][17] - sy) / sy)
            worst_f = max(worst_f, abs(row["sdv"][18] - fsge) / fsge)
            # and the stress actually carried is that surface, degraded
            # by the damage reached on this increment
            qd = (1.0 - row["sdv"][0]) * sy
            worst_q = max(worst_q, abs(row["sdv"][3] - qd)
                          / max(qd, 1e-30))
        check("q lies on the SGE flow surface: " + tag, worst < 1e-9,
              "max rel err %.3g" % worst)
        check("SDV19 equals the closed-form SGE factor: " + tag,
              worst_f < 1e-12, "max rel err %.3g" % worst_f)
        check("q equals (1-D) times the flow surface: " + tag,
              worst_q < 1e-9, "max rel err %.3g" % worst_q)
        check("SGE amplification recorded: " + tag,
              r[-1]["sdv"][18] >= 1.0,
              "fsge = %.6f at ep = %.4f" % (r[-1]["sdv"][18], r[-1]["sdv"][1]))

    # SGE must be monotone in 1/h and must vanish as h grows.
    fs = []
    for hloc in (1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6):
        r = grind.run(card(ihmode=2, h0=hloc),
                      [(3000, (0, 0, 0, dg, 0, 0))], nout=3000, dt=1e-9)
        fs.append(r[-1]["sdv"][18])
    check("SGE factor increases as h falls",
          all(fs[i + 1] > fs[i] for i in range(len(fs) - 1)),
          " -> ".join("%.4f" % v for v in fs))
    check("SGE factor -> 1 for a coarse chip", abs(fs[0] - 1.0) < 5e-3,
          "fsge(h=10um) = %.6f" % fs[0])

    # The three papers' algebra must agree at Lambda = 1: the peening form
    # sqrt(1 + 2 eta b (M a G)^2/s^2) and the blanking form
    # sqrt(s^2 + r' eta b (M a G)^2)/s are the same number when r' = 2.
    ep_t, h_t = 0.05, 1.0e-5
    sjc_t = c["A"] + c["B"] * ep_t ** c["n"]
    eta_t = 4.0 * ep_t / h_t
    f_pee = math.sqrt(1.0 + 2.0 * eta_t * c["b"]
                      * (c["M"] * c["alpha"] * g) ** 2 / sjc_t ** 2)
    f_bla = math.sqrt(sjc_t ** 2 + 2.0 * eta_t * c["b"]
                      * (c["M"] * c["alpha"] * g) ** 2) / sjc_t
    close("peening eq.25 and blanking eq.7 are the same formula",
          f_pee, f_bla, rtol=1e-14)

    # -- C4 strain-rate term ---------------------------------------------
    fast = grind.run(card(ihmode=2, h0=1.0), [(2000, (0, 0, 0, dg, 0, 0))],
                     nout=1, dt=1e-9)
    slow = grind.run(card(ihmode=2, h0=1.0), [(2000, (0, 0, 0, dg, 0, 0))],
                     nout=1, dt=1e-6)
    check("rate hardening raises the flow stress",
          fast[-1]["sdv"][3] > slow[-1]["sdv"][3],
          "q_fast %.4f > q_slow %.4f" % (fast[-1]["sdv"][3],
                                         slow[-1]["sdv"][3]))
    pc = card(ihmode=2, h0=1.0)
    for rows, dtv, tag in ((fast, 1e-9, "1e-9 s"), (slow, 1e-6, "1e-6 s")):
        ep0, t0 = rows[-2]["sdv"][1], rows[-2]["sdv"][15]
        thom = (t0 - c["T0"]) / (c["Tm"] - c["T0"])
        sy, _, _, _ = sge_flow(ep0, rows[-1]["sdv"][8], dtv, 1.0,
                               thom, c)
        close("rate term matches 1 + C ln(edot/edot0), dt = " + tag,
              rows[-1]["sdv"][17], sy, rtol=1e-9)

    # -- C5 adiabatic heating --------------------------------------------
    r = grind.run(card(ihmode=2, h0=1.0), [(4000, (0, 0, 0, dg, 0, 0))],
                  nout=1, dt=1e-9)
    tprev = c["T0"]
    dprev_d = 0.0
    worst = 0.0
    for row in r:
        dep = row["sdv"][8]
        if dep > 0:
            qpl = (1.0 - dprev_d) * row["sdv"][17]
            want = tprev + c["betaq"] * qpl * dep / (c["rho"] * c["cp"])
            worst = max(worst, abs(row["sdv"][15] - want))
        tprev = row["sdv"][15]
        dprev_d = row["sdv"][0]
    check("temperature follows beta*q*dep/(rho cp)", worst < 1e-9,
          "max |dT| error = %.3g K" % worst)
    check("temperature actually rises", r[-1]["sdv"][15] > c["T0"] + 1e-6,
          "T = %.4f K after ep = %.4f" % (r[-1]["sdv"][15], r[-1]["sdv"][1]))
    # Thermal softening must actually reduce the flow stress: sigma_jc at the
    # end has to sit below the isothermal value at the same plastic strain
    # and rate.
    ep_end, dep_end = r[-2]["sdv"][1], r[-1]["sdv"][8]
    iso, _, _, _ = sge_flow(ep_end, dep_end, 1e-9, 1.0, 0.0, c)
    check("thermal softening lowers the flow stress",
          0.0 < r[-1]["sdv"][16] < iso,
          "sigma_jc %.4f < isothermal %.4f MPa" % (r[-1]["sdv"][16], iso))

    # -- C6 Johnson-Cook damage and deletion ------------------------------
    # Pure shear: triaxiality is 0, so epsilon_f = (D1 + D2)(rate)(temp).
    r = grind.run(card(ihmode=2, h0=1.0), [(40000, (0, 0, 0, dg, 0, 0))],
                  nout=1, dt=1e-9)
    worst = 0.0
    for i, row in enumerate(r):
        dep = row["sdv"][8]
        dprev = r[i - 1]["sdv"][0] if i else 0.0
        if dep > 0 and row["sdv"][0] < 1.0:
            _, t0 = start_state(r, i)
            thom = (t0 - c["T0"]) / (c["Tm"] - c["T0"])
            edr = max((dep / 1e-9) / c["edot0"], 1.0)
            trix = -row["sdv"][2] / row["sdv"][3] if row["sdv"][3] > 0 else 0.0
            epf = ((c["D1"] + c["D2"] * math.exp(c["D3"] * trix))
                   * (1.0 + c["D4"] * math.log(edr))
                   * (1.0 + c["D5"] * thom))
            worst = max(worst, abs(row["sdv"][0] - (dprev + dep / epf)))
    check("JC damage accumulates as dep/eps_f", worst < 1e-9,
          "max |dD| error = %.3g" % worst)
    check("damage monotonically increases",
          all(r[i + 1]["sdv"][0] >= r[i]["sdv"][0] - 1e-15
              for i in range(len(r) - 1)))
    check("element deletes when D reaches DCRIT",
          r[-1]["sdv"][0] >= 1.0 and r[-1]["sdv"][11] == 0.0,
          "D = %.6f, STATUS = %.0f" % (r[-1]["sdv"][0], r[-1]["sdv"][11]))
    first_del = next((i for i, row in enumerate(r) if row["sdv"][11] == 0.0),
                     None)
    check("STATUS stays 0 once set", first_del is not None
          and all(row["sdv"][11] == 0.0 for row in r[first_del:]))
    check("STATUS is 1 before damage completes",
          r[0]["sdv"][11] == 1.0)

    # Compression must be far more ductile than tension: that asymmetry is
    # the whole reason a grit can plough without cracking.
    r_t = grind.run(card(ihmode=2, h0=1.0),
                    [(40000, (1e-6, -0.3e-6, -0.3e-6, 0, 0, 0))], nout=200,
                    dt=1e-9)
    r_c = grind.run(card(ihmode=2, h0=1.0),
                    [(40000, (-1e-6, 0.3e-6, 0.3e-6, 0, 0, 0))], nout=200,
                    dt=1e-9)
    check("both paths yield", min(r_t[-1]["sdv"][1], r_c[-1]["sdv"][1]) > 1e-4,
          "ep_t %.4f, ep_c %.4f" % (r_t[-1]["sdv"][1], r_c[-1]["sdv"][1]))
    check("tension damages faster than compression",
          r_t[-1]["sdv"][0] > r_c[-1]["sdv"][0],
          "D_tension %.4f > D_compression %.4f"
          % (r_t[-1]["sdv"][0], r_c[-1]["sdv"][0]))
    # NOTE for whoever revisits the constitutive model: the +/-1.5 triaxiality
    # clamp in vumat_grind.for was ADDED by the derived routine -- the reference
    # 'vumat_jc_damage (1).for':262-265 feeds the triaxiality straight into
    # eps_f with no clamp. Under a grit, P of order H against q of a few hundred
    # MPa gives trix ~ -4, well outside the clamp, where the unclamped law gives
    # eps_f = 0.15*exp(6) ~ 60 (no ductile fracture) against the clamp's 1.42.
    # That is a real question about the model, NOT a bug to be fixed silently:
    # it decides whether the zone this model labels ductile ploughs or fails.
    # Recorded here, deliberately not acted on.

    # -- C7 the damaged point must still carry compression ----------------
    # A grit sits on top of material that has already failed in shear. If the
    # damage factor were applied to the compressive pressure as well, that
    # material would carry nothing and the grit would sink through elements
    # that are still in the mesh.
    # Compression plus shear: the compression keeps the triaxiality negative
    # (so the failure strain stays long, as it must) while the shear supplies
    # the plastic strain to get there in a reasonable number of increments.
    r_cc = grind.run(card(ihmode=2, h0=1.0),
                     [(60000, (-1e-5, 0.3e-5, 0.3e-5, 1e-5, 0, 0))], nout=500,
                     dt=1e-9)
    check("that path stays in compression", r_cc[-1]["sdv"][2] > 0.0,
          "P = %.4g MPa" % r_cc[-1]["sdv"][2])
    check("compression reaches full damage", r_cc[-1]["sdv"][0] >= 1.0,
          "D = %.6f" % r_cc[-1]["sdv"][0])
    check("a fully damaged point still carries compressive pressure",
          r_cc[-1]["sdv"][2] > 0.0,
          "P = %.6g MPa at D = %.4f" % (r_cc[-1]["sdv"][2],
                                        r_cc[-1]["sdv"][0]))


# --------------------------------------------------------------------------
# D. the switch
# --------------------------------------------------------------------------

def test_switch(grind: Driver) -> None:
    print("D. the depth-of-cut switch")
    c = dict(JC)
    seg = [(10, (0, 0, 0, 1e-9, 0, 0))]

    # -- D1 h(u) from the coordinates -------------------------------------
    # h must stay positive across the whole probe, or the clamp at zero
    # masks what is being measured.
    thc, h0, hg, rtip = 0.35, 4.0e-3, -0.05, 25.0
    worst = 0.0
    for u in (-0.024, -0.01, 0.0, 0.01, 0.024):
        x = (-math.sin(thc) * u + 25.0 * math.cos(thc),
             math.cos(thc) * u + 25.0 * math.sin(thc), 0.004)
        # A radial offset must not change h: only the tangential station does.
        r = grind.run(card(dcut=1.0, thc=thc, h0=h0, hg=hg, rtip=rtip),
                      seg, coord=x, nout=10)
        want = h0 + hg * u - u * u / (2.0 * rtip)
        worst = max(worst, abs(r[-1]["sdv"][13] - want))
    check("h(u) matches H0 + HG u - u^2/(2 RTIP)", worst < 1e-15,
          "max |dh| = %.3g mm" % worst)

    u = 0.02
    x = (-math.sin(thc) * u + 25.0 * math.cos(thc),
         math.cos(thc) * u + 25.0 * math.sin(thc), 0.0)
    r = grind.run(card(dcut=1.0, thc=thc, h0=h0, hg=hg, rtip=0.0),
                  seg, coord=x, nout=10)
    check("RTIP <= 0 drops the curvature term",
          abs(r[-1]["sdv"][13] - (h0 + hg * u)) < 1e-15,
          "%.17g vs %.17g" % (r[-1]["sdv"][13], h0 + hg * u))

    r = grind.run(card(dcut=1.0, thc=0.0, h0=-1.0, hg=0.0), seg, nout=10)
    close("h is clamped at zero", r[-1]["sdv"][13], 0.0, atol=0.0)

    # -- D2 dc, both published forms --------------------------------------
    H, Kc = 1000.0, 0.3
    for form in (1, 2):
        want = dc_formula(0.15, H, c["E"], Kc, form)
        r = grind.run(card(dcut=0.0, lamc=0.15, hard=H, kic=Kc, idcf=form,
                           h0=0.0), seg, nout=10)
        close("dc, form %d" % form, r[-1]["sdv"][14], want, rtol=1e-13)
    d1 = dc_formula(0.15, H, c["E"], Kc, 1)
    d2 = dc_formula(0.15, H, c["E"], Kc, 2)
    close("the two dc forms differ by (E/H)^1.5", d2 / d1,
          (c["E"] / H) ** 1.5, rtol=1e-12)
    r = grind.run(card(dcut=7.5e-6, hard=H, kic=Kc), seg, nout=10)
    close("an explicit dc overrides the formula", r[-1]["sdv"][14], 7.5e-6,
          rtol=0.0, atol=0.0)

    # -- D3 the branch actually taken -------------------------------------
    dcv = 1.0e-4
    for h, want in ((0.5 * dcv, 1), (0.999 * dcv, 1), (dcv, 2), (2.0 * dcv, 2)):
        r = grind.run(card(dcut=dcv, h0=h), seg, nout=10)
        check("h = %.4g mm -> mode %d" % (h, want),
              int(round(r[-1]["sdv"][12])) == want,
              "got mode %d" % int(round(r[-1]["sdv"][12])))

    # The two branches must actually behave differently, or the switch is
    # decoration. Same history, one either side of dc.
    hist = [(4000, (-1e-6, 0, 0, 0, 0, 0))]
    duct = grind.run(card(dcut=dcv, h0=0.1 * dcv), hist, nout=4000)
    brit = grind.run(card(dcut=dcv, h0=10.0 * dcv), hist, nout=4000)
    check("the two branches give different stress",
          abs(duct[-1]["s"][0] - brit[-1]["s"][0]) > 1.0,
          "s11 ductile %.4g vs brittle %.4g MPa"
          % (duct[-1]["s"][0], brit[-1]["s"][0]))
    check("ductile branch leaves the JH-2 state variables untouched",
          duct[-1]["sdv"][9] == 0.0 and duct[-1]["sdv"][10] == 0.0)
    check("brittle branch leaves the JC state variables untouched",
          brit[-1]["sdv"][16] == 0.0 and brit[-1]["sdv"][17] == 0.0)

    # -- D4 latching -------------------------------------------------------
    long_hist = [(200, (0, 0, 0, 1e-7, 0, 0)), (200, (0, 0, 0, -1e-7, 0, 0)),
                 (200, (0, 0, 0, 1e-7, 0, 0))]
    for h in (0.5 * dcv, 2.0 * dcv):
        r = grind.run(card(dcut=dcv, h0=h), long_hist, nout=1)
        modes = {int(round(row["sdv"][12])) for row in r}
        hs = {round(row["sdv"][13], 15) for row in r}
        check("mode is latched for h = %.4g" % h, len(modes) == 1, str(modes))
        check("h is latched for h = %.4g" % h, len(hs) == 1, str(hs))
        check("INIT is set for h = %.4g" % h, r[-1]["sdv"][19] == 1.0)

    # -- D5 the overrides --------------------------------------------------
    r = grind.run(card(dcut=1e-9, h0=1.0, ihmode=2), seg, nout=10)
    check("IHMODE 2 forces ductile even when h >> dc",
          int(round(r[-1]["sdv"][12])) == 1)
    r = grind.run(card(dcut=1.0, h0=1e-9, ihmode=3), seg, nout=10)
    check("IHMODE 3 forces brittle even when h << dc",
          int(round(r[-1]["sdv"][12])) == 2)

    # -- D6 field-variable source -----------------------------------------
    for hv, want in ((0.5 * dcv, 1), (2.0 * dcv, 2)):
        r = grind.run(card(dcut=dcv, h0=99.0, ihmode=1), seg, fields=(hv,),
                      nout=10)
        check("field variable 1 supplies h = %.4g" % hv,
              abs(r[-1]["sdv"][13] - hv) < 1e-18
              and int(round(r[-1]["sdv"][12])) == want)


# --------------------------------------------------------------------------
# E. things that must hold whichever branch runs
# --------------------------------------------------------------------------

def test_shared(grind: Driver) -> None:
    print("E. shared invariants")
    dcv = 1.0e-4
    hist = [(3000, (-1e-6, 2e-7, 2e-7, 1e-7, 0, 0))]
    for h, tag in ((0.1 * dcv, "ductile"), (10.0 * dcv, "brittle")):
        r = grind.run(card(dcut=dcv, h0=h), hist, nout=1)
        check("%s: no NaN or Inf in the stress" % tag,
              all(math.isfinite(v) for row in r for v in row["s"]))
        check("%s: no NaN or Inf in the state" % tag,
              all(math.isfinite(v) for row in r for v in row["sdv"]))
        check("%s: damage stays in [0, 1]" % tag,
              all(-1e-15 <= row["sdv"][0] <= 1.0 + 1e-15 for row in r))
        check("%s: plastic strain never decreases" % tag,
              all(r[i + 1]["sdv"][1] >= r[i]["sdv"][1] - 1e-15
                  for i in range(len(r) - 1)))
        check("%s: STATUS is 0 or 1" % tag,
              all(row["sdv"][11] in (0.0, 1.0) for row in r))
        check("%s: inelastic energy never decreases" % tag,
              all(r[i + 1]["einel"] >= r[i]["einel"] - 1e-20
                  for i in range(len(r) - 1)))
        check("%s: dc is recorded on every point" % tag,
              all(row["sdv"][14] == dcv for row in r))

    # A zero increment must be a no-op, which is what Abaqus does on its very
    # first call to size the stable increment.
    r = grind.run(card(dcut=dcv, h0=0.1 * dcv),
                  [(1, (0, 0, 0, 0, 0, 0)), (1, (0, 0, 0, 0, 0, 0))], nout=1)
    check("a zero strain increment leaves the stress at zero",
          all(abs(v) < 1e-30 for row in r for v in row["s"]))

    # nstatev below the layout must not crash: the routine guards every write.
    for nsv in (12, 13, 15, 20):
        r = grind.run(card(dcut=dcv, h0=0.1 * dcv),
                      [(50, (0, 0, 0, 1e-7, 0, 0))], nstatev=nsv, nout=50)
        check("runs with nstatev = %d" % nsv, len(r) == 1 and
              all(math.isfinite(v) for v in r[-1]["s"]))


# --------------------------------------------------------------------------

def main() -> int:
    if not os.path.isdir(WORK):
        os.makedirs(WORK)
    for need in ("driver.f", "vaba_param.inc"):
        if not os.path.exists(os.path.join(WORK, need)):
            raise SystemExit("missing %s in %s" % (need, WORK))
    fc = find_gfortran()
    print("=" * 78)
    print("vumat_grind.for verification")
    print("  compiler: %s" % fc)
    print("=" * 78)
    grind = Driver(fc, SRC_GRIND, "grind.exe")
    jh2 = Driver(fc, SRC_JH2, "jh2.exe")

    test_source(grind)
    test_jh2_identity(grind, jh2)
    test_jh2_benchmarks(grind)
    test_ductile(grind)
    test_switch(grind)
    test_shared(grind)

    print("=" * 78)
    if _FAIL:
        print("%d passed, %d FAILED" % (_PASS, len(_FAIL)))
        for f in _FAIL:
            print("   - " + f)
        return 1
    print("ALL %d CHECKS PASSED" % _PASS)
    return 0


if __name__ == "__main__":
    sys.exit(main())
