import random
import numpy as np

# -------------------- Settings --------------------
TRIALS = 200000   # Number of trials
STEPS = 7         # Steps per trial
SEED = 1          # Seed
# --------------------------------------------------

k1, k2, k3 = 1.0, 2.0, 3.0
x0 = (9, 8, 7)


def choose2(n: int) -> float:
    return 0.0 if n < 2 else n * (n - 1) / 2.0


def step(x1: int, x2: int, x3: int):
    # propensities (weights)
    a1 = k1 * choose2(x1) * x2
    a2 = k2 * x1 * choose2(x3)
    a3 = k3 * x2 * x3
    a0 = a1 + a2 + a3

    if a0 == 0.0:
        return None  # no reaction can fire

    r = random.random() * a0
    if r < a1:
        return (x1 - 2, x2 - 1, x3 + 4)   # R1
    elif r < a1 + a2:
        return (x1 - 1, x2 + 3, x3 - 2)   # R2
    else:
        return (x1 + 2, x2 - 1, x3 - 1)   # R3


def run_one():
    x1, x2, x3 = x0
    for _ in range(STEPS):
        nxt = step(x1, x2, x3)
        if nxt is None:
            break
        x1, x2, x3 = nxt
    return x1, x2, x3


def main():
    if SEED is not None:
        random.seed(SEED)

    finals = np.zeros((TRIALS, 3), dtype=np.int64)

    for i in range(TRIALS):
        finals[i, :] = run_one()

    means = finals.mean(axis=0)
    vars_ = finals.var(axis=0, ddof=0)  # population variance estimate

    print(f"Start: {list(x0)}, Steps: {STEPS}, Trials: {TRIALS}, Seed: {SEED}\n")

    print(f"Mean(X1) = {means[0]:.6f}    Var(X1) = {vars_[0]:.6f}")
    print(f"Mean(X2) = {means[1]:.6f}    Var(X2) = {vars_[1]:.6f}")
    print(f"Mean(X3) = {means[2]:.6f}    Var(X3) = {vars_[2]:.6f}")


if __name__ == "__main__":
    main()