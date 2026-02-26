import random
import math

# -------------------- User settings --------------------
TRIALS = 5000       # increase for better accuracy (runtime increases)
N_STEPS = 30000     # steps per trial; increase if events seem rare
SEED = 1            # set None for random seed
# -------------------------------------------------------

# Reaction rate constants
k1, k2, k3 = 1.0, 2.0, 3.0

# Initial state
x1_0, x2_0, x3_0 = 110, 26, 55


def choose2(n):
    """n choose 2 (how many pairs can be chosen from n items)."""
    if n < 2:
        return 0
    return n * (n - 1) / 2.0


def one_trial():
    """Run one stochastic trajectory for N_STEPS reaction firings."""
    x1, x2, x3 = x1_0, x2_0, x3_0

    hit_c1 = False
    hit_c2 = False
    hit_c3 = False

    for _ in range(N_STEPS):
        # Propensities (reaction 'weights')
        a1 = k1 * choose2(x1) * x2           # R1: 2X1 + X2 -> 4X3
        a2 = k2 * x1 * choose2(x3)           # R2: X1 + 2X3 -> 3X2
        a3 = k3 * x2 * x3                    # R3: X2 + X3 -> 2X1
        a0 = a1 + a2 + a3

        # If nothing can fire, stop
        if a0 == 0:
            break

        # Pick which reaction fires (Gillespie: choose by relative propensity)
        r = random.random() * a0
        if r < a1:
            # R1: (-2, -1, +4)
            x1 -= 2
            x2 -= 1
            x3 += 4
        elif r < a1 + a2:
            # R2: (-1, +3, -2)
            x1 -= 1
            x2 += 3
            x3 -= 2
        else:
            # R3: (+2, -1, -1)
            x1 += 2
            x2 -= 1
            x3 -= 1

        # Check outcomes after each firing
        if x1 >= 150:
            hit_c1 = True
        if x2 < 10:
            hit_c2 = True
        if x3 > 100:
            hit_c3 = True

        # If we already hit all outcomes at least once, we can stop early
        if hit_c1 and hit_c2 and hit_c3:
            break

    return hit_c1, hit_c2, hit_c3


def main():
    if SEED is not None:
        random.seed(SEED)

    c1_hits = 0
    c2_hits = 0
    c3_hits = 0

    for _ in range(TRIALS):
        h1, h2, h3 = one_trial()
        c1_hits += 1 if h1 else 0
        c2_hits += 1 if h2 else 0
        c3_hits += 1 if h3 else 0

    print("EE 5393 HW1 — Problem 1(a) :contentReference[oaicite:1]{index=1}")
    print(f"TRIALS={TRIALS}, N_STEPS={N_STEPS}, SEED={SEED}")
    print(f"Start state S0 = [{x1_0}, {x2_0}, {x3_0}]")
    print()
    print("Estimated probabilities (event hit at least once within N_STEPS):")
    print(f"Pr(C1: x1 >= 150) = {c1_hits / TRIALS:.5f}")
    print(f"Pr(C2: x2 < 10)   = {c2_hits / TRIALS:.5f}")
    print(f"Pr(C3: x3 > 100)  = {c3_hits / TRIALS:.5f}")


if __name__ == "__main__":
    main()
