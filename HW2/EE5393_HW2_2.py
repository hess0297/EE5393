import math
import random
from typing import Dict, List, Tuple

INPUT_SEQUENCE = [100, 5, 500, 20, 250]

MAX_TIME_PER_PHASE = 10000.0
MAX_STEPS_PER_PHASE = 1_000_000
SEED = 1

K_SLOW = 0.01
K_FAST = 100.0

Stoich = Dict[str, int]
Reaction = Tuple[Stoich, Stoich, float]

RXNS_BLUE_RED: List[Reaction] = [
    ({"X": 1}, {"A": 1, "R_1": 1}, K_SLOW),
    ({"B_1": 1}, {"F": 1, "C": 1, "R_2": 1}, K_SLOW),
    ({"B_2": 1}, {"H": 1, "E": 1}, K_SLOW),

    ({"A": 2}, {"A_1": 1}, K_FAST),
    ({"A_1": 2}, {"A_2": 1}, K_FAST),
    ({"A_2": 2}, {"Y": 1}, K_FAST),

    ({"F": 2}, {"F_1": 1}, K_FAST),
    ({"F_1": 2}, {"F_2": 1}, K_FAST),
    ({"F_2": 2}, {"X": 1}, K_FAST),

    ({"C": 2}, {"C_1": 1}, K_FAST),
    ({"C_1": 2}, {"C_2": 1}, K_FAST),
    ({"C_2": 2}, {"Y": 1}, K_FAST),

    ({"H": 2}, {"H_1": 1}, K_FAST),
    ({"H_1": 2}, {"H_2": 1}, K_FAST),
    ({"H_2": 2}, {"X": 1}, K_FAST),

    ({"E": 2}, {"E_1": 1}, K_FAST),
    ({"E_1": 2}, {"E_2": 1}, K_FAST),
    ({"E_2": 2}, {"Y": 1}, K_FAST),
]

RXNS_RED_GREEN: List[Reaction] = [
    ({"R_1": 1}, {"G_1": 1}, K_SLOW),
    ({"R_2": 1}, {"G_2": 1}, K_SLOW),
]

RXNS_GREEN_BLUE: List[Reaction] = [
    ({"G_1": 1}, {"B_1": 1}, K_SLOW),
    ({"G_2": 1}, {"B_2": 1}, K_SLOW),
]

INIT_COUNTS: Dict[str, int] = {
    "X": 0,
    "Y": 0,
    "A": 0,   "A_1": 0, "A_2": 0,
    "F": 0,   "F_1": 0, "F_2": 0,
    "C": 0,   "C_1": 0, "C_2": 0,
    "H": 0,   "H_1": 0, "H_2": 0,
    "E": 0,   "E_1": 0, "E_2": 0,
    "R_1": 0, "G_1": 0, "B_1": 0,
    "R_2": 0, "G_2": 0, "B_2": 0,
}


def nCk(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    k = min(k, n - k)
    out = 1
    for i in range(1, k + 1):
        out = out * (n - (k - i)) // i
    return out


def propensity(counts: Dict[str, int], reactants: Stoich, rate: float) -> float:
    a = rate
    for sp, m in reactants.items():
        c = nCk(counts.get(sp, 0), m)
        if c == 0:
            return 0.0
        a *= c
    return float(a)


def fire_reaction(counts: Dict[str, int], reactants: Stoich, products: Stoich) -> None:
    for sp, m in reactants.items():
        counts[sp] = counts.get(sp, 0) - m
    for sp, m in products.items():
        counts[sp] = counts.get(sp, 0) + m


def run_phase(rxns: List[Reaction], counts_in: Dict[str, int]) -> Dict[str, int]:
    counts = dict(counts_in)
    t = 0.0

    for _ in range(MAX_STEPS_PER_PHASE):
        if t >= MAX_TIME_PER_PHASE:
            break

        props: List[float] = []
        a0 = 0.0

        for reactants, products, rate in rxns:
            a = propensity(counts, reactants, rate)
            props.append(a)
            a0 += a

        if a0 <= 0.0:
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

        reactants, products, _ = rxns[idx]
        fire_reaction(counts, reactants, products)

    return counts


def run_one_cycle(counts_in: Dict[str, int], x_value: int) -> Dict[str, int]:
    counts = dict(counts_in)

    counts["Y"] = 0
    counts["X"] = x_value

    counts = run_phase(RXNS_BLUE_RED, counts)
    counts["Y_recorded"] = counts.get("Y", 0)

    counts = run_phase(RXNS_RED_GREEN, counts)

    counts["Y"] = 0

    counts = run_phase(RXNS_GREEN_BLUE, counts)

    return counts


def main() -> None:
    if SEED is not None:
        random.seed(SEED)

    counts = dict(INIT_COUNTS)

    print(" i |   X |   Y |  B1 |  B2")
    print("---------------------------")

    for i, xval in enumerate(INPUT_SEQUENCE, start=1):
        counts = run_one_cycle(counts, xval)
        print(f"{i:>2} | {xval:>4} | {counts['Y_recorded']:>4} | "
              f"{counts['B_1']:>4} | {counts['B_2']:>4}")

    print("\nFinal counts:")
    for sp in sorted(counts):
        print(f"{sp:>10s} : {counts[sp]}")


if __name__ == "__main__":
    main()