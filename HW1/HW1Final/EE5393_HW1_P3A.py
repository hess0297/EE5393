"""
EE5393 HW1 P3A — Gillespie SSA for the original CRN.
Runs NUM_RUNS trials and reports mean/std of final (w,z).
"""

import math, random, statistics

# ---- SETTINGS ----
NUM_RUNS = 100
BASE_SEED = 100
PRINT_EVERY = 10  # 0 disables
T_END = 200000.0
MAX_STEPS = 10_000_000

# initial counts (no scaling)
INIT = {"a": 0, "b": 1, "c": 0, "y": 8192, "yP": 0, "w": 0, "wP": 0, "x": 200, "d": 0, "z": 0}

k = dict(
    rb=0.02, r1=8.0, r2=12.0, r3=60.0, r4=1.5, r5=0.034,
    r6=0.002, r7=18000.0, r8=0.091, r9=0.3
)

RXNS = [
    ("rb", {"b": 1},         {"a": 1, "b": 1},              "rb"),
    ("r1", {"a": 1, "y": 2}, {"c": 1, "yP": 1, "a": 1},     "r1"),
    ("r2", {"c": 2},         {"c": 1},                      "r2"),
    ("r3", {"a": 1},         {},                            "r3"),
    ("r4", {"yP": 1},        {"y": 1},                      "r4"),
    ("r5", {"c": 1},         {"w": 1},                      "r5"),
    ("r6", {"x": 1},         {"d": 1},                      "r6"),
    ("r7", {"d": 1, "w": 1}, {"d": 1, "wP": 1, "z": 1},     "r7"),
    ("r8", {"d": 1},         {},                            "r8"),
    ("r9", {"wP": 1},        {"w": 1},                      "r9"),
]

def nCk(n, r):
    if r < 0 or r > n: return 0
    r = min(r, n - r)
    num = den = 1
    for i in range(1, r + 1):
        num *= n - (r - i)
        den *= i
    return num // den

def propensity(c, R, rate):
    a = rate
    for sp, m in R.items():
        n = c.get(sp, 0)
        if n < m: return 0.0
        a *= nCk(n, m)
    return float(a)

def fire(c, R, P):
    for sp, m in R.items(): c[sp] -= m
    for sp, m in P.items(): c[sp] = c.get(sp, 0) + m

def done(c):
    return (
        c.get("y", 0) <= 1 and c.get("x", 0) == 0 and c.get("d", 0) == 0 and
        c.get("wP", 0) == 0 and c.get("c", 0) == 0 and c.get("yP", 0) == 0 and
        c.get("a", 0) == 0
    )

def ssa(seed, init):
    random.seed(seed)
    c = dict(init)
    t = 0.0

    for _ in range(MAX_STEPS):
        if done(c): return c, "done"

        props, a0 = [], 0.0
        for _, R, _, rk in RXNS:
            ai = propensity(c, R, k[rk])
            props.append(ai)
            a0 += ai
        if a0 <= 0.0: return c, "no reactions possible"

        dt = -math.log(random.random()) / a0
        if t + dt > T_END: return c, "reached T_END"
        t += dt

        r = random.random() * a0
        s = 0.0
        for i, ai in enumerate(props):
            s += ai
            if s >= r:
                _, R, P, _ = RXNS[i]
                fire(c, R, P)
                break

    return c, "reached MAX_STEPS"

def mean_std(xs):
    return statistics.mean(xs), (statistics.stdev(xs) if len(xs) > 1 else 0.0)

def main():
    target_w = math.log2(INIT["y"])
    target_z = INIT["x"] * target_w
    print(f"Target z = {target_z}\n")

    finals, reasons = [], {}
    for i in range(NUM_RUNS):
        final, reason = ssa(BASE_SEED + i, INIT)
        finals.append(final)
        reasons[reason] = reasons.get(reason, 0) + 1
        if PRINT_EVERY and (i + 1) % PRINT_EVERY == 0:
            print(f"Run {i+1:3d}/{NUM_RUNS}: z={final.get('z',0)} w={final.get('w',0)} reason={reason}")

    w = [f.get("w", 0) for f in finals]
    z = [f.get("z", 0) for f in finals]

    mw, sw = mean_std(w)
    mz, sz = mean_std(z)

    print("\n--- Summary ---")
    print(f"Stop reasons: {reasons}")

    print(f"\n--- Statistics over {NUM_RUNS} runs ---")

    print("\nW results:")
    print(f"  target(w)  = {target_w:.6f}")
    print(f"  mean(w)    = {mw:.6f}")
    print(f"  std(w)     = {sw:.6f}")

    print("\nZ results:")
    print(f"  target(z)  = {target_z:.6f}")
    print(f"  mean(z)    = {mz:.6f}")
    print(f"  std(z)     = {sz:.6f}")

    species = sorted({sp for f in finals for sp in f})
    print("\n--- Mean final counts (all species) ---")
    for sp in species:
        print(f"{sp:>3} : {statistics.mean([f.get(sp, 0) for f in finals]):.3f}")

if __name__ == "__main__":
    main()