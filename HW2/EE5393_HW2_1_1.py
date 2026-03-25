import math
import random

# -------------------- User settings --------------------
SEED = 1
MAX_TIME = 1e6
MAX_STEPS = 100000

# bootstrap + 10 Fibonacci stages = 11 reactions total
RATES = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05]
# ------------------------------------------------------


def nCk(n, k):
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    k = min(k, n - k)
    out = 1
    for i in range(1, k + 1):
        out = out * (n - (k - i)) // i
    return out


def propensity(counts, reactants, rate):
    a = rate
    for sp, m in reactants.items():
        n = counts.get(sp, 0)
        c = nCk(n, m)
        if c == 0:
            return 0.0
        a *= c
    return float(a)


def run_fibonacci_ssa():
    if SEED is not None:
        random.seed(SEED)

    counts = {}
    counts["S"] = 1

    # recorded Fibonacci species
    for i in range(1, 13):
        counts[f"X{i}"] = 0

    # active pair carriers
    for i in range(1, 11):
        counts[f"A{i}"] = 0
        counts[f"B{i}"] = 0

    reactions = []

    # Bootstrap from a single input
    reactions.append((
        {"S": 1},
        {"A1": 1, "B1": 1, "X1": 1, "X2": 1},
        RATES[0]
    ))

    # Stage 1 through stage 9
    for n in range(1, 10):
        reactions.append((
            {f"A{n}": 1},
            {f"B{n+1}": 1, f"X{n+2}": 1},
            RATES[n]
        ))
        reactions.append((
            {f"B{n}": 1},
            {f"A{n+1}": 1, f"B{n+1}": 1, f"X{n+2}": 1},
            RATES[n]
        ))

    # Final stage n = 10 creates X12 only
    reactions.append((
        {"A10": 1},
        {"X12": 1},
        RATES[10]
    ))
    reactions.append((
        {"B10": 1},
        {"X11": 0, "X12": 1},   # X11 already built earlier; keep syntax uniform
        RATES[10]
    ))

    t = 0.0

    print("=" * 72)
    print("FIBONACCI SSA (single input, X1 through X12)")
    print("=" * 72)
    print(f"Initial: S={counts['S']}")

    for step in range(1, MAX_STEPS + 1):
        if t >= MAX_TIME:
            print("\nStopped: MAX_TIME reached.")
            break

        props = []
        a0 = 0.0

        for reactants, products, rate in reactions:
            a = propensity(counts, reactants, rate)
            props.append(a)
            a0 += a

        if a0 <= 0.0:
            print(f"\nStopped at step {step-1}: no reactions can fire.")
            break

        r1 = random.random()
        dt = -math.log(max(r1, 1e-300)) / a0
        t += dt

        r2 = random.random() * a0
        s = 0.0
        idx = 0
        for i, a in enumerate(props):
            s += a
            if r2 <= s:
                idx = i
                break

        reactants, products, _ = reactions[idx]

        for sp, m in reactants.items():
            counts[sp] -= m
        for sp, m in products.items():
            counts[sp] = counts.get(sp, 0) + m

        x_state = ", ".join(
            f"X{i}={counts[f'X{i}']}" for i in range(1, 13) if counts[f"X{i}"] > 0
        )

    print("\nFinal:")
    for i in range(1, 13):
        print(f"X{i} = {counts[f'X{i}']}")

    print("=" * 72)


if __name__ == "__main__":
    run_fibonacci_ssa()