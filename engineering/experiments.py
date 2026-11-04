"""Open-track experiments: noise-precision-query laws for relation-lattice recovery.

    python engineering/experiments.py e1        # cost law: smallest b for >= 95 % 4-sample recovery
    python engineering/experiments.py e2        # noise sweep (full solver) + corrupted-shot rate
    python engineering/experiments.py e3        # mitigation: 3-fold majority vote on stages m <= 19
    python engineering/experiments.py e4        # group structure from the Smith normal form
    python engineering/experiments.py plots     # PNGs from the CSVs in results/
    python engineering/experiments.py all

Development code only: it uses instance answers / true lattices to SCORE runs, never inside solve().
Every statement in REPORT_NOTES.md about these results is labelled theorem / simulator fact /
empirical observation.
"""

from __future__ import annotations

import csv
import importlib.util
import math
import multiprocessing as mp
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "engineering"))
sys.path.insert(0, str(ROOT / "tests"))

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister  # noqa: E402

from silentgroup import MockGroup, Oracle, oracle_gate  # noqa: E402

RES = ROOT / "results"
SOLVER = ROOT / "submission" / "my_solver.py"
WORKERS = max(1, mp.cpu_count() - 2)


def load(path=SOLVER):
    spec = importlib.util.spec_from_file_location("solver_exp", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = load()


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, mid - half), min(1.0, mid + half))


def write_csv(name, rows):
    RES.mkdir(exist_ok=True)
    with open(RES / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def read_csv(name):
    with open(RES / name) as f:
        return list(csv.DictReader(f))


def true_hnf(inst):
    """Exact relation lattice of (A, B, C, g) from the group (test helper), canonical HNF."""
    from test_lattice import relation_lattice_basis
    els = [inst.elements[e] for e in S.ELS]
    return S.hnf(relation_lattice_basis(list(inst.group.moduli), els))


def cyclic_modulus(k):
    return (1 << k) - 3                                # odd, ~2^k; MockGroup picks random g, a, b, c


# ============================================================================ E1: cost law
E1_SIZES = (20, 30, 40, 50, 58)
E1_BITS = tuple(range(14, 45, 2))
E1_SEEDS = 30


def e1_task(args):
    k, b, seed = args
    inst = MockGroup([cyclic_modulus(k)], seed=seed)
    orc = Oracle(inst, seed=10_000 + seed)
    ctx = S._Ctx(orc)
    qc = S.build_pe_circuit(b)
    samples = [S.decode(next(iter(orc.run(qc))), b) for _ in range(4)]
    res = S.relation_lattice(ctx, samples, b)           # 4 samples, no consensus, no tests
    ok = res is not None and res[0] == true_hnf(inst)
    return {"log2G": k, "b": b, "seed": seed, "ok": int(ok), "log2H": round(math.log2(S.det_of(true_hnf(inst))), 2)}


def run_e1():
    tasks = [(k, b, s) for k in E1_SIZES for b in E1_BITS for s in range(E1_SEEDS)
             if b >= 10 + k // 4 - 2 and b <= 24 + k // 4 + 4]   # skip hopeless / trivially-ok b
    with mp.Pool(WORKERS) as pool:
        rows = pool.map(e1_task, tasks, chunksize=4)
    write_csv("e1_cost_law_raw.csv", rows)
    e1_summary(rows)


def e1_summary(rows=None):
    """b50 = interpolated b where 4-sample success first reaches 50 %; plateau = mean success over the
    top b values (last 6 bits of the grid); b95 = smallest b with success >= 95 % (often not reached)."""
    if rows is None:
        rows = [{k: float(v) for k, v in r.items()} for r in read_csv("e1_cost_law_raw.csv")]
    out = []
    for k in E1_SIZES:
        pts = {}
        for r in rows:
            if int(r["log2G"]) == k:
                pts.setdefault(int(r["b"]), []).append(int(r["ok"]))
        bs = sorted(pts)
        rate = [sum(pts[b]) / len(pts[b]) for b in bs]
        b50 = None
        for i in range(1, len(bs)):
            if rate[i - 1] < 0.5 <= rate[i]:
                b50 = bs[i - 1] + (0.5 - rate[i - 1]) / (rate[i] - rate[i - 1]) * (bs[i] - bs[i - 1])
                break
        topb = [b for b in bs if b >= bs[-1] - 6]
        k_top = sum(sum(pts[b]) for b in topb)
        n_top = sum(len(pts[b]) for b in topb)
        lo, hi = wilson(k_top, n_top)
        b95 = next((b for b, x in zip(bs, rate) if x >= 0.95), None)
        log2H = statistics.mean(float(r["log2H"]) for r in rows if int(r["log2G"]) == k)
        pred = 16 + log2H / 4
        out.append({"log2G": k, "log2H": round(log2H, 2), "b50": round(b50, 2) if b50 else "", "b_pred": round(pred, 2),
                    "plateau": round(k_top / n_top, 3), "plateau_lo": round(lo, 3), "plateau_hi": round(hi, 3),
                    "b95": b95 or "not reached", "per_b": " ".join(f"{b}:{x:.2f}" for b, x in zip(bs, rate)),
                    "q4_at_b50": round(8 * b50 * (b50 + 1)) if b50 else "", "q4_pred": round(8 * pred * (pred + 1))})
    write_csv("e1_cost_law.csv", out)
    xs = [r["log2H"] for r in out]
    ys = [r["b50"] for r in out]
    mx, my = statistics.mean(xs), statistics.mean(ys)
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    for r in out:
        print(f"log2|H|={r['log2H']:5}: b50={r['b50']:6} pred={r['b_pred']:6} plateau={r['plateau']} "
              f"[{r['plateau_lo']},{r['plateau_hi']}] b95={r['b95']}  q4(b50)={r['q4_at_b50']} q4(pred)={r['q4_pred']}")
    print(f"fit b50 = {my - slope * mx:.1f} + {slope:.3f} * log2|H|   (heuristic: 16 + 0.25 * log2|H|)")


# ============================================================================ E2: noise sweep
E2_NOISE = (0.002, 0.004, 0.008, 0.02, 0.05)
E2_SHAPES = {"cyclic_2^58": ([(1 << 58) - 27], {}), "public1_shape": ([6145295, 53834531687], {})}
E2_SEEDS = 20
E2_SHOTS = 300


def e2_solve_task(args):
    shape, nu, seed = args
    mod, kw = E2_SHAPES[shape]
    inst = MockGroup(mod, seed=seed, **kw)
    orc = Oracle(inst, budget=80_000, seed=20_000 + seed, noise=nu)
    t0 = time.perf_counter()
    try:
        res, err = S.solve(orc), ""
    except Exception as exc:                            # record, never hide
        res, err = None, repr(exc)[:120]
    return {"shape": shape, "noise": nu, "seed": seed, "ok": int(res == inst.answer),
            "order_ok": int(isinstance(res, dict) and res.get("order") == inst.answer["order"]),
            "false_none": int(isinstance(res, dict) and res.get("secret") == "none" and inst.answer["secret"] != "none"),
            "queries": orc.queries_used, "seconds": round(time.perf_counter() - t0, 2), "error": err}


def corrupted(sample, R, b):
    """A shot is corrupted if some TRUE short relation has a tail beyond the solver's tolerance."""
    tol = 1 << (b - 12)
    for x in R:
        v = sum(a * c for a, c in zip(x, sample)) % (1 << b)
        if min(v, (1 << b) - v) > tol:
            return True
    return False


def e2_shot_task(args):
    nu, chunk, builder, b = args
    inst = MockGroup([cyclic_modulus(40)], seed=7)
    R = S.lll(true_hnf(inst))
    orc = Oracle(inst, seed=30_000 + chunk, noise=nu)
    qc, dec = BUILDERS[builder](b)
    bad = 0
    for _ in range(E2_SHOTS // 10):
        bad += corrupted(dec(next(iter(orc.run(qc)))), R, b)
    return {"noise": nu, "builder": builder, "shots": E2_SHOTS // 10, "corrupted": bad,
            "queries_per_shot": orc.cost(qc, 1)}


def run_e2():
    tasks = [(sh, nu, s) for sh in E2_SHAPES for nu in E2_NOISE for s in range(E2_SEEDS)]
    with mp.Pool(WORKERS) as pool:
        rows = pool.map(e2_solve_task, tasks, chunksize=1)
        write_csv("e2_noise_solve.csv", rows)
        shots = pool.map(e2_shot_task, [(nu, c, "plain", 36) for nu in E2_NOISE for c in range(10)])
    agg = []
    for nu in E2_NOISE:
        rs = [r for r in shots if r["noise"] == nu]
        k, n = sum(r["corrupted"] for r in rs), sum(r["shots"] for r in rs)
        lo, hi = wilson(k, n)
        agg.append({"noise": nu, "shots": n, "corrupted": k, "rate": round(k / n, 4), "lo": round(lo, 4), "hi": round(hi, 4)})
    write_csv("e2_corruption.csv", agg)
    for sh in E2_SHAPES:
        for nu in E2_NOISE:
            rs = [r for r in rows if r["shape"] == sh and r["noise"] == nu]
            k = sum(r["ok"] for r in rs)
            lo, hi = wilson(k, len(rs))
            print(f"{sh:14s} nu={nu:<6} ok={k}/{len(rs)} [{lo:.2f},{hi:.2f}] medQ={statistics.median(r['queries'] for r in rs):.0f} "
                  f"false_none={sum(r['false_none'] for r in rs)} maxT={max(r['seconds'] for r in rs)}")
    for a in agg:
        print(f"nu={a['noise']}: corrupted-shot rate {a['rate']} [{a['lo']},{a['hi']}]")
    fit = [(math.log(a["noise"]), math.log(a["rate"])) for a in agg if 0 < a["corrupted"] and a["rate"] < 0.6]
    if len(fit) >= 2:
        mx = statistics.mean(x for x, _ in fit)
        my = statistics.mean(y for _, y in fit)
        slope = sum((x - mx) * (y - my) for x, y in fit) / sum((x - mx) ** 2 for x, _ in fit)
        print(f"log-log slope of corrupted-shot rate vs noise = {slope:.2f} (hypothesis: 2)")


# ============================================================================ E3: mitigation
VOTE_M = 19                                            # stages m <= 19 (powers <= 2^18) are voted


def build_voted_circuit(b, trunc=12, vote_m=VOTE_M):
    """Like S.build_pe_circuit, but stages m <= vote_m are measured 3 times and the MAJORITY is fed back.

    Each stage has its own classical register (1 or 3 bits); feedback on a voted stage is applied when
    the register holds 3, 5, 6 or 7 (majority 1). if_tests without oracle gates cost no queries."""
    q = QuantumRegister(1, "q")
    regs = {}
    qc = QuantumCircuit(q)
    for ei in range(4):
        for m in range(b, 0, -1):
            regs[ei, m] = ClassicalRegister(3 if m <= vote_m else 1, f"c{ei}_{m}")
            qc.add_register(regs[ei, m])
    for ei, e in enumerate(S.ELS):
        for m in range(b, 0, -1):
            reg = regs[ei, m]
            for copy in range(len(reg)):
                qc.h(0)
                qc.append(oracle_gate(e, 1 << (m - 1)), [0])
                for l in range(m + 1, min(b, m + trunc) + 1):
                    src, ang = regs[ei, l], -2 * math.pi / 2 ** (l - m + 1)
                    if len(src) == 1:
                        with qc.if_test((src[0], 1)):
                            qc.p(ang, 0)
                    else:
                        for v in (3, 5, 6, 7):          # majority of the three copies is 1
                            with qc.if_test((src, v)):
                                qc.p(ang, 0)
                qc.h(0)
                qc.measure(0, reg[copy])
                qc.reset(0)
    index = {bit: i for i, bit in enumerate(qc.clbits)}

    def decode(key):
        bits = key[::-1]
        out = []
        for ei in range(4):
            val = 0
            for m in range(1, b + 1):
                reg = regs[ei, m]
                ones = sum(int(bits[index[bt]]) for bt in reg)
                val |= (1 if 2 * ones > len(reg) else 0) << (b - m)
            out.append(val)
        return out
    return qc, decode


def build_plain(b):
    qc = S.build_pe_circuit(b)
    return qc, lambda key: S.decode(key, b)


BUILDERS = {"plain": build_plain, "voted": build_voted_circuit}
E3_NOISE = (0.002, 0.008, 0.02)


def e3_recovery_task(args):
    nu, builder, seed = args
    inst = MockGroup([cyclic_modulus(58)], seed=seed)
    orc = Oracle(inst, seed=40_000 + seed, noise=nu)
    qc, dec = BUILDERS[builder](36)
    ctx = S._Ctx(orc)
    samples = [dec(next(iter(orc.run(qc)))) for _ in range(4)]
    res = S.relation_lattice(ctx, samples, 36)
    return {"noise": nu, "builder": builder, "seed": seed, "ok": int(res is not None and res[0] == true_hnf(inst)),
            "queries": orc.queries_used}


def run_e3():
    with mp.Pool(WORKERS) as pool:
        shots = pool.map(e2_shot_task, [(nu, c, bl, 36) for nu in E3_NOISE for bl in BUILDERS for c in range(10)])
        rec = pool.map(e3_recovery_task, [(nu, bl, s) for nu in E3_NOISE for bl in BUILDERS for s in range(40)])
    rows = []
    for nu in E3_NOISE:
        for bl in BUILDERS:
            rs = [r for r in shots if r["noise"] == nu and r["builder"] == bl]
            k, n = sum(r["corrupted"] for r in rs), sum(r["shots"] for r in rs)
            rr = [r for r in rec if r["noise"] == nu and r["builder"] == bl]
            ok = sum(r["ok"] for r in rr)
            q4 = rr[0]["queries"]
            lo, hi = wilson(k, n)
            slo, shi = wilson(ok, len(rr))
            rows.append({"noise": nu, "builder": bl, "queries_per_shot": rs[0]["queries_per_shot"],
                         "corrupt_rate": round(k / n, 4), "corrupt_lo": round(lo, 4), "corrupt_hi": round(hi, 4),
                         "recovery4": round(ok / len(rr), 3), "rec_lo": round(slo, 3), "rec_hi": round(shi, 3),
                         "recovery_per_10k_queries": round(1e4 * ok / len(rr) / q4, 4)})
    write_csv("e3_mitigation.csv", rows)
    for r in rows:
        print(r)


# ============================================================================ E4: group structure
def minors_gcd(M, k):
    """gcd of all k x k minors of the 4x4 integer matrix M (determinantal divisor d_k)."""
    from itertools import combinations
    g = 0
    for rows in combinations(range(4), k):
        for cols in combinations(range(4), k):
            sub = [[M[r][c] for c in cols] for r in rows]
            g = math.gcd(g, int_det(sub))
    return g


def int_det(A):
    """Exact integer determinant (Bareiss fraction-free elimination)."""
    A = [row[:] for row in A]
    n, sign, prev = len(A), 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k]), None)
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def group_structure(H):
    """Invariant factors d1 | d2 | ... of Z^4 / L (= the subgroup <A, B, C, g>), from the SNF of L.

    s_k = D_k / D_{k-1}, D_k = gcd of k x k minors (theorem: determinantal divisors). No queries."""
    D = [1] + [minors_gcd(H, k) for k in range(1, 5)]
    return [D[k] // D[k - 1] for k in range(1, 5) if D[k] // D[k - 1] > 1]


E4_SHAPES = {"cyclic_2^58": [(1 << 58) - 27], "2^20x6x3^8": [2**20, 6, 3**8], "2^29_sq": [2**29, 2**29],
             "forged2_AB": [10**8, 3 * 10**8], "forged3_ABC": [2**28, 3 * 2**28], "public2_shape": [10628367750, 12297285]}
E4_FORGED = {"forged2_AB": "AB", "forged3_ABC": "ABC"}


def e4_task(args):
    shape, seed = args
    inst = MockGroup(E4_SHAPES[shape], seed=seed, forged=E4_FORGED.get(shape, ""))
    captured = {}
    orig = S.answer

    def spy(H, trusted):                               # capture the verified lattice (dev only)
        captured["H"] = H
        return orig(H, trusted)
    S.answer = spy
    try:
        res = S.solve(Oracle(inst, budget=80_000, seed=50_000 + seed))
    finally:
        S.answer = orig
    truth = group_structure(true_hnf(inst))
    got = group_structure(captured["H"]) if "H" in captured else None
    G = inst.group
    exponent = math.lcm(*(G.element_order(inst.elements[e]) for e in S.ELS))   # independent check
    return {"shape": shape, "seed": seed, "solve_ok": int(res == inst.answer), "structure_true": truth,
            "structure_got": got, "match": int(got == truth), "exponent_check": int(bool(truth) and truth[-1] == exponent)}


def run_e4():
    tasks = [(sh, s) for sh in E4_SHAPES for s in range(10)]
    with mp.Pool(WORKERS) as pool:
        rows = pool.map(e4_task, tasks)
    write_csv("e4_structure.csv", rows)
    for sh in E4_SHAPES:
        rs = [r for r in rows if r["shape"] == sh]
        print(f"{sh:14s} structure match {sum(r['match'] for r in rs)}/{len(rs)}  exponent check "
              f"{sum(r['exponent_check'] for r in rs)}/{len(rs)}  e.g. Z_{' x Z_'.join(map(str, rs[0]['structure_true']))}")


# ============================================================================ plots
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"


def style(ax, title, xlabel, ylabel):
    ax.set_facecolor(SURF)
    ax.set_title(title, color=INK, fontsize=11, loc="left")
    ax.set_xlabel(xlabel, color=INK2)
    ax.set_ylabel(ylabel, color=INK2)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(INK2)
    ax.tick_params(colors=INK2)


def noise_ticks(ax, nus):
    """Log x-axis with one plain label per measured noise level (no overlapping minor labels)."""
    from matplotlib.ticker import NullFormatter, NullLocator
    ax.set_xscale("log")
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks(nus)
    ax.set_xticklabels([f"{nu:g}" for nu in nus])


def plots():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if (RES / "e1_cost_law.csv").exists():
        rows = read_csv("e1_cost_law.csv")
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=SURF)
        shades = ["#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#104281"]   # one hue, light -> dark = size
        for r, c in zip(rows, shades):
            pts = [p.split(":") for p in r["per_b"].split()]
            xs, ys = [int(b) for b, _ in pts], [float(x) for _, x in pts]
            a1.plot(xs, ys, color=c, linewidth=2, marker="o", markersize=5, label=f"|G| ≈ 2^{r['log2G']}")
        a1.legend(frameon=False, labelcolor=INK2, fontsize=8, loc="lower right")
        a1.axhline(0.95, color=INK2, linewidth=1, linestyle="--")
        style(a1, "4-sample recovery vs precision b (30 seeds/pt)", "bits per element b", "P(lattice exact)")
        xs = [float(r["log2H"]) for r in rows]
        a2.plot(xs, [float(r["b_pred"]) for r in rows], color=ORANGE, linewidth=2, label="heuristic 16 + log2|H|/4")
        a2.plot(xs, [float(r["b50"]) for r in rows], color=BLUE, linewidth=2, marker="o", markersize=8, label="measured b50")
        a2.legend(frameon=False, labelcolor=INK2)
        style(a2, "Precision where 4-sample recovery reaches 50 %", "log2 |H|", "b")
        fig.tight_layout()
        fig.savefig(RES / "e1_cost_law.png", dpi=150)
        plt.close(fig)

    if (RES / "e2_noise_solve.csv").exists():
        rows = read_csv("e2_noise_solve.csv")
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=SURF)
        for sh, c in zip(E2_SHAPES, (BLUE, ORANGE)):
            nus = sorted({float(r["noise"]) for r in rows if r["shape"] == sh})
            rate, lo, hi, medq = [], [], [], []
            for nu in nus:
                rs = [r for r in rows if r["shape"] == sh and float(r["noise"]) == nu]
                k = sum(int(r["ok"]) for r in rs)
                w = wilson(k, len(rs))
                rate.append(k / len(rs))
                lo.append(k / len(rs) - w[0])
                hi.append(w[1] - k / len(rs))
                medq.append(statistics.median(int(r["queries"]) for r in rs))
            a1.errorbar(nus, rate, yerr=[lo, hi], color=c, linewidth=2, marker="o", markersize=7, capsize=3, label=sh)
            a2.plot(nus, medq, color=c, linewidth=2, marker="o", markersize=7, label=sh)
        for a in (a1, a2):
            noise_ticks(a, list(E2_NOISE))
            a.legend(frameon=False, labelcolor=INK2)
        style(a1, "Full solver: success vs noise (20 seeds, Wilson 95 %)", "noise ν (turns/√query)", "P(order & secret right)")
        style(a2, "Full solver: median queries vs noise", "noise ν (turns/√query)", "queries")
        fig.tight_layout()
        fig.savefig(RES / "e2_noise.png", dpi=150)
        plt.close(fig)

    if (RES / "e2_corruption.csv").exists():
        rows = read_csv("e2_corruption.csv")
        fig, ax = plt.subplots(figsize=(5.8, 4.2), facecolor=SURF)
        nus = [float(r["noise"]) for r in rows]
        rate = [max(float(r["rate"]), 1e-4) for r in rows]
        ax.errorbar(nus, rate, yerr=[[max(0, float(r["rate"]) - float(r["lo"])) for r in rows],
                                     [float(r["hi"]) - float(r["rate"]) for r in rows]],
                    color=BLUE, linewidth=2, marker="o", markersize=7, capsize=3, label="measured (300 shots/pt)")
        ref = rate[0] if rate[0] > 1e-4 else rate[1] / (nus[1] / nus[0]) ** 2
        ax.plot(nus, [min(1.0, ref * (nu / nus[0]) ** 2) for nu in nus], color=ORANGE, linewidth=2, linestyle="--",
                label="∝ ν² (hypothesis, capped at 1)")
        noise_ticks(ax, nus)
        ax.set_yscale("log")
        ax.legend(frameon=False, labelcolor=INK2)
        style(ax, "Corrupted-shot rate at b=36 (|G|≈2^40)", "noise ν", "fraction of shots corrupted")
        fig.tight_layout()
        fig.savefig(RES / "e2_corruption.png", dpi=150)
        plt.close(fig)

    if (RES / "e3_mitigation.csv").exists():
        rows = read_csv("e3_mitigation.csv")
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=SURF)
        for bl, c in (("plain", BLUE), ("voted", ORANGE)):
            rs = [r for r in rows if r["builder"] == bl]
            nus = [float(r["noise"]) for r in rs]
            a1.plot(nus, [float(r["corrupt_rate"]) for r in rs], color=c, linewidth=2, marker="o", markersize=7,
                    label=f"{bl} ({rs[0]['queries_per_shot']} q/shot)")
            a2.plot(nus, [float(r["recovery_per_10k_queries"]) for r in rs], color=c, linewidth=2, marker="o", markersize=7, label=bl)
        for a in (a1, a2):
            noise_ticks(a, list(E3_NOISE))
            a.legend(frameon=False, labelcolor=INK2)
        style(a1, "Corrupted-shot rate: plain vs 3-vote on m ≤ 19", "noise ν", "fraction of shots corrupted")
        style(a2, "4-sample recovery per 10k queries (40 seeds)", "noise ν", "recoveries / 10k queries")
        fig.tight_layout()
        fig.savefig(RES / "e3_mitigation.png", dpi=150)
        plt.close(fig)
    print("plots written to", RES)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    steps = {"e1": run_e1, "e1sum": e1_summary, "e2": run_e2, "e3": run_e3, "e4": run_e4, "plots": plots}
    for name in (steps if which == "all" else [which]):
        if name in steps:
            t0 = time.time()
            steps[name]()
            print(f"[{name} done in {time.time() - t0:.0f}s]")
