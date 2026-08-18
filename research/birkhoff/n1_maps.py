"""Session brief N, part 1 -- which induced maps on the cone are LINEAR,
and do any of them contract the Hilbert (= exchange) metric?

Setting (see FINDINGS.md of this directory for the statements).

    K   = cone of positive functions on the compactified [0, inf], pointwise
          order.  Its Hilbert projective metric is
              d(Phi, Psi) = log sup(Psi/Phi) + log sup(Phi/Psi) .
    C   = { Phi convex, nondecreasing, Phi >= Lam_Phi * beta }  subset of K,
          the projective closure of { F_a = log Z_a }.
    d(a,b) = d(F_a, F_b)      (brief I, OBSTRUCTION.md Sec. 1.1)

Birkhoff--Hopf needs a map that is LINEAR on the cone whose order defines the
metric.  The resource operations are:

    tensor        a (x) c            F -> F + F_c            translation
    Cartesian     a^{(x)k}           F -> k F                scalar
    fibre power   a^[m] = (a_i^m)    F -> F(m beta)          dilation
    disjoint      a + c              F -> log(e^F + e^{F_c}) log-sum-exp

This script checks, numerically and exhaustively where possible:

  (1) the dictionary itself, against genuine integer signatures;
  (2) additivity and homogeneity of each induced map (the linearity audit);
  (3) whether each map contracts, preserves or EXPANDS d;
  (4) the mediant proposition: tensoring is d-nonexpansive exactly on the
      mutually non-dominating pairs, and expands without bound off them;
  (5) the same sweep on the F_3 / F_5 quadratic and F_3 cubic map classes.

    python research/birkhoff/n1_maps.py
"""
from __future__ import annotations

import itertools
import math
import random
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "realizability"))
sys.path.insert(0, str(HERE.parent.parent / "src"))

import i_cone as T                                        # noqa: E402
import common as CM                                       # noqa: E402
from optimizers import differential_evolution, pattern_search  # noqa: E402

LOG2 = math.log(2.0)


# ---------------------------------------------------------------------------
# cone operations
# ---------------------------------------------------------------------------

def tensor(a: T.Trop, b: T.Trop) -> T.Trop:
    """Phi_a + Phi_b.  max_i(c_i+x_i b) + max_j(c'_j+x'_j b) = max_ij(...)."""
    c = (a.c[:, None] + b.c[None, :]).ravel()
    x = (a.x[:, None] + b.x[None, :]).ravel()
    return T.Trop(c, x)


def power(a: T.Trop, k: float) -> T.Trop:
    """Cartesian power: Phi -> k Phi."""
    return T.Trop(k * a.c, k * a.x)


def dilate(a: T.Trop, m: float) -> T.Trop:
    """Fibre power a_i -> a_i^m: Phi -> Phi(m beta)."""
    return T.Trop(a.c, m * a.x)


def mixture(a: T.Trop, ws, ms) -> T.Trop:
    """T_mu Phi = sum_j w_j Phi(m_j beta): a tensor of fibre powers."""
    out = power(dilate(a, ms[0]), ws[0])
    for w, m in zip(ws[1:], ms[1:]):
        out = tensor(out, power(dilate(a, m), w))
    return out


def flatten(a: T.Trop) -> T.Trop:
    """L Phi = Phi(0) + Lam_Phi beta -- the flat signature with the same (r,M).

    Linear:  Phi(0) and Lam_Phi = lim Phi' are both linear functionals on C.
    """
    return T.Trop([a.R], [a.Lam])


def corner(a: T.Trop) -> T.Trop:
    """N Phi = max(Phi(0), Lam_Phi beta) -- the OTHER retraction onto the
    sigma-line.  Not linear (a max of two linear functionals is not additive)."""
    return T.Trop([a.R, 0.0], [0.0, a.Lam])


def disjoint(a: T.Trop, b: T.Trop, grid) -> np.ndarray:
    """log(e^Phi_a + e^Phi_b) on a grid (not a cone operation: it leaves C)."""
    fa, fb = a.val(grid), b.val(grid)
    m = np.maximum(fa, fb)
    return m + np.log(np.exp(fa - m) + np.exp(fb - m))


def disjoint_d(a, b, c, grid):
    """d(a + c, b + c) with the two endpoints beta = 0, infinity included."""
    fa = disjoint(a, c, grid)
    fb = disjoint(b, c, grid)
    v = np.log(fb) - np.log(fa)
    e0 = math.log(a.R + c.R)                    # F(0) = log(r_a + r_c)
    e0 = math.log(math.log(math.exp(b.R) + math.exp(c.R))) - \
        math.log(math.log(math.exp(a.R) + math.exp(c.R)))
    e1 = math.log(max(b.Lam, c.Lam)) - math.log(max(a.Lam, c.Lam))
    hi = max(float(v.max()), e0, e1)
    lo = min(float(v.min()), e0, e1)
    return hi - lo


# ---------------------------------------------------------------------------
# random cone elements
# ---------------------------------------------------------------------------

def rand_trop(rng, k=None, hi=5.0):
    while True:
        kk = k if k is not None else rng.randint(1, 4)
        try:
            t = T.Trop([rng.uniform(0.02, hi) for _ in range(kk)],
                       [rng.uniform(0.02, hi) for _ in range(kk)])
            return t
        except ValueError:
            continue


# ---------------------------------------------------------------------------
# (1) the dictionary, against genuine integer signatures
# ---------------------------------------------------------------------------

def sig_tensor(a, b):
    return tuple(sorted((x * y for x in a for y in b), reverse=True))


def sig_union(s1, s2):
    """a + b: disjoint union of signatures, on compressed Sig objects."""
    return CM.Sig.from_logs(list(s1.xs) + list(s2.xs),
                            list(s1.mults) + list(s2.mults))


def check_dictionary(rng, trials=40):
    print("=== (1) dictionary: resource operation  ->  action on F = log Z ===")
    grid = np.linspace(0.0, 1000.0, 4001)          # beta horizon 10^3
    wt = wp = wd = wu = 0.0
    for _ in range(trials):
        a = tuple(sorted((rng.randint(2, 9) for _ in range(rng.randint(2, 4))),
                         reverse=True))
        b = tuple(sorted((rng.randint(2, 9) for _ in range(rng.randint(2, 4))),
                         reverse=True))
        Sa, Sb = CM.Sig.of(a), CM.Sig.of(b)
        Fa, Fb = Sa.F(grid), Sb.F(grid)
        # tensor
        St = CM.Sig.of(sig_tensor(a, b))
        wt = max(wt, float(np.abs(St.F(grid) - (Fa + Fb)).max()))
        # Cartesian power k = 3
        Sp = CM.Sig.of(sig_tensor(sig_tensor(a, a), a))
        wp = max(wp, float(np.abs(Sp.F(grid) - 3.0 * Fa).max()))
        # fibre power m = 2  (a_i -> a_i^2)
        Sd = CM.Sig.of(tuple(v * v for v in a))
        wd = max(wd, float(np.abs(Sd.F(grid) - Sa.F(2.0 * grid)).max()))
        # disjoint union
        Su = CM.Sig.of(tuple(sorted(a + b, reverse=True)))
        m = np.maximum(Fa, Fb)
        wu = max(wu, float(np.abs(Su.F(grid)
                                  - (m + np.log(np.exp(Fa - m) + np.exp(Fb - m)))).max()))
    print(f"  a (x) b        F_a + F_b            max error {wt:.3e}")
    print(f"  a^(x)3         3 F_a                max error {wp:.3e}")
    print(f"  a^[2]          F_a(2 beta)          max error {wd:.3e}")
    print(f"  a + b          log(e^F_a + e^F_b)   max error {wu:.3e}")
    print(f"  ({trials} random integer signature pairs, beta grid [0, 10^3])")


# ---------------------------------------------------------------------------
# (2) the linearity audit
# ---------------------------------------------------------------------------

def check_linearity(rng, trials=400):
    print("\n=== (2) linearity audit on the cone: T(Phi+Psi) = T Phi + T Psi ?"
          "   T(lam Phi) = lam T Phi ? ===")
    grid = np.linspace(0.0, 1000.0, 4001)
    ops = {}
    fixed = rand_trop(rng, 2)
    ops["tensor with a fixed c"] = lambda p: tensor(p, fixed)
    ops["Cartesian power k=3"] = lambda p: power(p, 3.0)
    ops["fibre power m=2"] = lambda p: dilate(p, 2.0)
    ops["dilation mixture"] = lambda p: mixture(p, [1.0, 2.0], [0.5, 3.0])
    ops["flattening L"] = flatten
    ops["corner N"] = corner
    rows = []
    for name, op in ops.items():
        add = hom = 0.0
        for _ in range(trials):
            p, q = rand_trop(rng), rand_trop(rng)
            lam = rng.uniform(0.3, 3.0)
            lhs = op(tensor(p, q)).val(grid)
            rhs = op(p).val(grid) + op(q).val(grid)
            scale = max(1.0, float(np.abs(rhs).max()))
            add = max(add, float(np.abs(lhs - rhs).max()) / scale)
            lhs = op(power(p, lam)).val(grid)
            rhs = lam * op(p).val(grid)
            scale = max(1.0, float(np.abs(rhs).max()))
            hom = max(hom, float(np.abs(lhs - rhs).max()) / scale)
        rows.append((name, add, hom))
    # the disjoint union, on the same footing
    add = hom = 0.0
    for _ in range(trials):
        p, q = rand_trop(rng), rand_trop(rng)
        lam = rng.uniform(0.3, 3.0)
        lhs = disjoint(tensor(p, q), fixed, grid)
        rhs = disjoint(p, fixed, grid) + disjoint(q, fixed, grid)
        add = max(add, float(np.abs(lhs - rhs).max()) / max(1.0, float(np.abs(rhs).max())))
        lhs = disjoint(power(p, lam), fixed, grid)
        rhs = lam * disjoint(p, fixed, grid)
        hom = max(hom, float(np.abs(lhs - rhs).max()) / max(1.0, float(np.abs(rhs).max())))
    rows.append(("disjoint union with a fixed c", add, hom))
    print(f"  {'operation':<32} {'rel. additivity defect':>22} "
          f"{'rel. homogeneity defect':>24}")
    for name, add, hom in rows:
        print(f"  {name:<32} {add:22.3e} {hom:24.3e}")
    print("  (a defect at machine precision means LINEAR; O(1) means not)")


# ---------------------------------------------------------------------------
# (3)+(4) contraction audit
# ---------------------------------------------------------------------------

def contact(a: T.Trop, b: T.Trop):
    """(inf, sup) of Phi_a/Phi_b over beta in [0, inf]."""
    bs = T.candidates(a, b)
    v = -T.phi(a, b, bs)                      # log(Phi_a/Phi_b)
    e_inf = math.log(a.Lam) - math.log(b.Lam)
    return min(float(v.min()), e_inf), max(float(v.max()), e_inf)


def check_contraction(rng, trials=200000):
    print("\n=== (3) does the map contract d?   e = d(Ta,Tb) - d(a,b),  "
          "kappa = sup d(Ta,Tb)/d(a,b) ===")
    names = ["tensor with fixed c", "fibre power m=2.7", "Cartesian power k=3",
             "dilation mixture", "flattening L", "corner N"]
    hi = {n: -1e9 for n in names}
    lo = {n: +1e9 for n in names}
    kap = {n: 0.0 for n in names}
    worst_nd = -1e9
    n_nd = 0
    sgrid = np.exp(np.linspace(-30.0, 30.0, 3001))
    worst_disj = -1e9
    for t in range(trials):
        a, b = rand_trop(rng), rand_trop(rng)
        d0 = T.hilbert(a, b)
        c = rand_trop(rng)
        img = {
            "tensor with fixed c": (tensor(a, c), tensor(b, c)),
            "fibre power m=2.7": (dilate(a, 2.7), dilate(b, 2.7)),
            "Cartesian power k=3": (power(a, 3.0), power(b, 3.0)),
            "dilation mixture": (mixture(a, [1.0, 2.0], [0.5, 3.0]),
                                 mixture(b, [1.0, 2.0], [0.5, 3.0])),
            "flattening L": (flatten(a), flatten(b)),
            "corner N": (corner(a), corner(b)),
        }
        for n, (ta, tb) in img.items():
            d1 = T.hilbert(ta, tb)
            hi[n] = max(hi[n], d1 - d0)
            lo[n] = min(lo[n], d1 - d0)
            if d0 > 1e-9:
                kap[n] = max(kap[n], d1 / d0)
        m, M = contact(a, b)
        if m <= 0.0 <= M:                      # both exchange rates <= 1
            n_nd += 1
            worst_nd = max(worst_nd, img["tensor with fixed c"] and
                           T.hilbert(*img["tensor with fixed c"]) - d0)
        if t < 3000:
            worst_disj = max(worst_disj, disjoint_d(a, b, c, sgrid) - d0)
    print(f"  {trials} random cone triples (1-4 lines each)")
    print(f"  {'map':<22} {'max e':>12} {'min e':>12} {'sup ratio':>12}")
    for n in names:
        print(f"  {n:<22} {hi[n]:+12.6f} {lo[n]:+12.6f} {kap[n]:12.9f}")
    print(f"  tensoring restricted to the {n_nd} mutually non-dominating pairs "
          f"(both rates <= 1): max e = {worst_nd:+.3e}")
    print(f"  disjoint union with c: max e = {worst_disj:+.6f}   "
          f"(3000 triples, beta in [e^-30, e^30])")


# ---------------------------------------------------------------------------
# the unbounded-expansion family, in closed form and in signatures
# ---------------------------------------------------------------------------

def unbounded_family():
    print("\n=== (4) tensoring expands d without bound: a and a^(x)k are at "
          "d = 0, their tensors with c are not ===")
    from mpmath import mp, mpf, log as mplog
    mp.dps = 45
    print("  Phi_a = max(1, beta),  Phi_{a^k} = k Phi_a  (so d = 0 exactly),")
    print("  Phi_c = N + eps*beta  (a flat resource).  Exact:")
    print("    d = log[ (k+eps)(1+N) / ((1+eps)(k+N)) ]  ->  log k .")
    print(f"  {'k':>4} {'N':>10} {'eps':>10} {'d(a,a^k)':>10} "
          f"{'d(a(x)c, a^k(x)c)':>34} {'log k - d':>12}")
    for k in (2, 3, 8):
        for e10 in (2, 10, 40):
            N = mpf(10) ** e10
            eps = mpf(10) ** (-e10)
            a = T.Trop([1.0, 0.0], [0.0, 1.0])
            ak = power(a, float(k))
            d0 = T.hilbert(a, ak)
            d1 = mplog((k + eps) * (1 + N) / ((1 + eps) * (k + N)))
            print(f"  {k:>4} {'1e%d' % e10:>10} {'1e-%d' % e10:>10} "
                  f"{d0:10.2e} {mp.nstr(d1, 30):>34} "
                  f"{float(mplog(k) - d1):12.3e}")
    # numerical confirmation in the Trop machinery at one point of the family
    a = T.Trop([1.0, 0.0], [0.0, 1.0])
    ak = power(a, 3.0)
    c = T.Trop([1e6], [1e-6])
    # the same phenomenon for the disjoint union, on genuine integer-style
    # signatures: Z_c dominates Z_a at small beta and is dominated at large beta
    print("  the disjoint union does the same, on genuine signatures:")
    grid = np.linspace(-25.0, 25.0, 200001)

    def sig_power(sig, k):
        out = {0.0: 1.0}
        for _ in range(k):
            new = {}
            for x0, m0 in out.items():
                for x, m in zip(sig.xs, sig.mults):
                    new[x0 + x] = new.get(x0 + x, 0.0) + m0 * m
            out = new
        return CM.Sig.from_logs(list(out.keys()), list(out.values()))

    sa = CM.Sig.from_logs([1.0, 0.0], [1.0, math.e ** 3])
    for k in (2, 3):
        for N in (10, 30, 60):
            sak = sig_power(sa, k)
            sc = CM.Sig.from_logs([1e-3], [math.exp(N)])
            d0, _ = CM.d_and_A(sa, sak, grid)
            d1, _ = CM.d_and_A(sig_union(sa, sc), sig_union(sak, sc), grid)
            print(f"    k={k} log r_c={N:>3}: d(a, a^(x)k) = {d0:.2e},  "
                  f"d(a+c, a^(x)k+c) = {d1:.9f},  log k = {math.log(k):.9f}")
    print(f"  Trop cross-check k=3, N=1e6, eps=1e-6:  "
          f"d = {T.hilbert(tensor(a, c), tensor(ak, c)):.15f}   "
          f"exact = {float(mplog((3+mpf('1e-6'))*(1+mpf('1e6'))/((1+mpf('1e-6'))*(3+mpf('1e6'))))):.15f}")


# ---------------------------------------------------------------------------
# hill-climb: how large can the expansion be at a given d(a,b) > 0 ?
# ---------------------------------------------------------------------------

def _mk(z, k, n):
    z = np.asarray(z, float).reshape(n, k, 2)
    out = []
    for row in z:
        c = np.maximum(row[:, 0], 0.0)
        x = np.maximum(row[:, 1], 0.0)
        if c.max() <= 1e-9:
            c = c + 1e-3
        if x.max() <= 1e-9:
            x = x + 1e-3
        out.append(T.Trop(c, x))
    return out


def _neg_expand_nd(z, k):
    """expansion under tensoring, restricted to non-dominating pairs."""
    try:
        a, b, c = _mk(z, k, 3)
    except ValueError:
        return 1e3
    lo, hi = contact(a, b)
    if not (lo <= 1e-12 <= hi + 1e-12):
        return 1e3 + abs(lo) + abs(hi)
    return -(T.hilbert(tensor(a, c), tensor(b, c)) - T.hilbert(a, b))


def climb_nondominating(seed=3, k=3, restarts=5):
    print("\n  hill-climb of the expansion over NON-DOMINATING pairs "
          "(Prop. 2.3 says it must be <= 0):")
    bounds = [(0.0, 8.0)] * (3 * k * 2)
    best = math.inf
    for t in range(restarts):
        z, f = differential_evolution(_neg_expand_nd, bounds, args=(k,),
                                      seed=seed + 977 * t, maxiter=350,
                                      popsize=14, F=(0.3, 1.2), CR=0.9)
        for step in (0.5, 0.1, 0.02, 4e-3, 8e-4, 1.6e-4, 3e-5, 6e-6, 1e-6):
            z, f = pattern_search(_neg_expand_nd, z, args=(k,), step=step,
                                  min_step=1e-13, maxiter=30000, bounds=bounds)
        best = min(best, f)
    print(f"    max expansion found = {-best:+.12f}   (k = {k} lines each)")


# ---------------------------------------------------------------------------
# (5) the finite-field map classes
# ---------------------------------------------------------------------------

def field_classes(q, degree):
    mons = [(i, j) for i in range(degree + 1) for j in range(degree + 1)
            if i + j <= degree]
    seen = set()
    for co in itertools.product(range(q), repeat=len(mons)):
        counts = {}
        for x in range(q):
            for y in range(q):
                v = sum(c * pow(x, i, q) * pow(y, j, q)
                        for c, (i, j) in zip(co, mons)) % q
                counts[v] = counts.get(v, 0) + 1
        sig = tuple(sorted(counts.values(), reverse=True))
        if len(sig) > 1:
            seen.add(sig)
    return sorted(seen, reverse=True)


def sig_to_trop(sig):
    """F_a is not piecewise linear; use the exact signature machinery instead."""
    return CM.Sig.of(sig)


def check_classes(label, pool):
    print(f"\n  --- {label}: {len(pool)} classes {pool} ---")
    S = [CM.Sig.of(s) for s in pool]
    n = len(pool)
    worst = -1e9
    worst_arg = None
    n_dom = n_nd = 0
    for i, j, kk in itertools.product(range(n), repeat=3):
        if i >= j:
            continue
        a, b, c = S[i], S[j], S[kk]
        d0, _ = CM.d_and_A(a, b)
        at = CM.Sig.of(sig_tensor(pool[i], pool[kk]))
        bt = CM.Sig.of(sig_tensor(pool[j], pool[kk]))
        d1, _ = CM.d_and_A(at, bt)
        e0 = math.log(b.R) - math.log(a.R)
        e1 = math.log(b.Lam) - math.log(a.Lam)
        hi, lo, _, _ = CM.extrema(a, b)
        dominating = not (lo <= 0.0 <= hi)
        n_dom += dominating
        n_nd += (not dominating)
        if d1 - d0 > worst:
            worst = d1 - d0
            worst_arg = (pool[i], pool[j], pool[kk], d0, d1, dominating)
    print(f"    triples (a<b, c): {n*(n-1)//2*n}   dominating pairs used "
          f"{n_dom}, non-dominating {n_nd}")
    print(f"    max d(a(x)c, b(x)c) - d(a,b) = {worst:+.9f}")
    if worst_arg is not None:
        a, b, c, d0, d1, dom = worst_arg
        print(f"      at a={a}, b={b}, c={c}:  {d0:.6f} -> {d1:.6f}"
              f"   (one resource dominates: {dom})")


def thompson(a: T.Trop, b: T.Trop):
    """Thompson's part metric  max( log sup(Phi_b/Phi_a), log sup(Phi_a/Phi_b) )."""
    bs = T.candidates(a, b)
    v = T.phi(a, b, bs)
    e_inf = math.log(b.Lam) - math.log(a.Lam)
    hi = max(float(v.max()), e_inf)
    lo = min(float(v.min()), e_inf)
    return max(hi, -lo)


def check_thompson(rng, trials=100000):
    """Nussbaum's hypothesis is order-preserving + homogeneous of degree 1.
    Tensoring is order-preserving but NOT homogeneous, which is exactly why it
    expands d.  It is however Thompson-nonexpansive."""
    print("\n=== (3b) the same maps in Thompson's part metric ===")
    worst_t = worst_u = -1e9
    sgrid = np.exp(np.linspace(-30.0, 30.0, 3001))
    for t in range(trials):
        a, b, c = rand_trop(rng), rand_trop(rng), rand_trop(rng)
        t0 = thompson(a, b)
        worst_t = max(worst_t, thompson(tensor(a, c), tensor(b, c)) - t0)
        if t < 3000:
            fa, fb, = disjoint(a, c, sgrid), disjoint(b, c, sgrid)
            v = np.log(fb) - np.log(fa)
            e0 = (math.log(math.log(math.exp(b.R) + math.exp(c.R)))
                  - math.log(math.log(math.exp(a.R) + math.exp(c.R))))
            e1 = math.log(max(b.Lam, c.Lam)) - math.log(max(a.Lam, c.Lam))
            hi = max(float(v.max()), e0, e1)
            lo = min(float(v.min()), e0, e1)
            worst_u = max(worst_u, max(hi, -lo) - t0)
    print(f"  tensor with a fixed c:  max d_T(Ta,Tb) - d_T(a,b) = {worst_t:+.3e} "
          f"  ({trials} triples)   NONEXPANSIVE")
    print(f"  disjoint union with c:  max d_T(Ta,Tb) - d_T(a,b) = {worst_u:+.6f} "
          f"  (3000 triples)")


def birkhoff_diameters(rng, trials=2000):
    """Delta(T) = projective diameter of the image.  Birkhoff ratio tanh(Delta/4)."""
    print("\n=== (6) the Birkhoff projective diameters ===")
    ws, ms = [1.0, 2.0], [0.5, 3.0]
    shift = math.log(sum(ws)) - math.log(sum(w * m for w, m in zip(ws, ms)))
    worst = 0.0
    for _ in range(trials):
        p = rand_trop(rng)
        worst = max(worst, abs(mixture(p, ws, ms).sigma - p.sigma - shift))
    print(f"  dilation mixture: sigma(T Phi) - sigma(Phi) is constant "
          f"= {shift:+.12f}; max deviation over {trials} draws = {worst:.3e}")
    print("  hence d(T Phi, T Psi) >= |sigma_Phi - sigma_Psi|, unbounded on C:")
    print(f"  {'W':>8} {'d(T flat_0, T flat_W)':>24} {'tanh(Delta/4) bound':>22}")
    for W in (1.0, 10.0, 100.0, 600.0):
        f0 = T.Trop([1.0], [1.0])                      # sigma = 0
        f1 = T.Trop([math.exp(W)], [1.0])              # sigma = W
        print(f"  {W:8.0f} {T.hilbert(mixture(f0, ws, ms), mixture(f1, ws, ms)):24.9f}"
              f" {'1 (Delta = inf)':>22}")


def main():
    rng = random.Random(20260818)
    check_dictionary(rng)
    check_linearity(rng)
    check_contraction(rng)
    check_thompson(rng)
    unbounded_family()
    climb_nondominating()
    birkhoff_diameters(rng)
    print("\n=== (5) the finite-field map classes ===")
    check_classes("F_3 quadratic", field_classes(3, 2))
    check_classes("F_3 cubic", field_classes(3, 3))
    check_classes("F_5 quadratic", field_classes(5, 2))


if __name__ == "__main__":
    main()
