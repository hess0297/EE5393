# EE5393 HW1 P2 - Lambda stochastic simulation (Gillespie SSA)
# Assumes the following 3 files are in the SAME folder as this .py file:
#   1) lambda_r.txt   (reactions: "reactants : products : rate")
#   2) lambda_in.txt  (initial counts: "Species  Value  ...")
#   3) this script

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


# -------------------- User settings --------------------
TRIALS_PER_MOI = 50
MOI_VALUES = range(1, 11)

MAX_TIME  = 5000.0
MAX_STEPS = 5_000_000
SEED = 1

STEALTH_THRESHOLD = 145  # stealth when cI2 > 145
HIJACK_THRESHOLD = 55    # hijack when Cro2 > 55

REACTIONS_FILENAME = "lambda_r.txt"
INIT_FILENAME = "lambda_in.txt"
# ------------------------------------------------------


Stoich = Dict[str, int]
Reaction = Tuple[Stoich, Stoich, float]  # (reactants, products, rate)


def nCk(n: int, k: int) -> int:
    """Compute binomial coefficient C(n,k) for integers n>=0."""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    k = min(k, n - k)
    out = 1
    for i in range(1, k + 1):
        out = out * (n - (k - i)) // i
    return out


def parse_stoich(side: str) -> Stoich:
    side = side.strip() # Remove leading and trailing spaces
    if not side: # If line is empty return blank
        return {}
    
    toks = side.split() # Split each part of the line at spaces (create tokens)
    if len(toks) % 2 != 0: # If spaces does not come in a pair of specie and amount then fail program
        raise ValueError(f"Bad stoichiometry side (odd token count): {side!r}")

    d: Stoich = {}
    for i in range(0, len(toks), 2): # Loop over tokens two at a time (gather the specie and number values togethersdfasd)
        sp = toks[i]
        m = int(toks[i + 1])
        d[sp] = d.get(sp, 0) + m
    return d


def load_reactions(path: Path) -> List[Reaction]:
    rxns: List[Reaction] = []
    with path.open("r", encoding="utf-8", errors="replace") as f: # Open path, UTF-8 encoding, replace characters with placeholder
        for raw in f: # Each line in the reaction file is saved under raw, this goes until there is not another line to save (f has been through)
            line = raw.strip() # Remove leading and trailing spaces and newline characters
            
            if not line or line.startswith("#"): # Ignore blank lines
                continue
            
            parts = [p.strip() for p in line.split(":")] # Isolate the parts of each line (reactants:products:rate)
            if len(parts) != 3: # If there are not 3 parts in the line indicate that something is wrong
                raise ValueError(f"Bad reaction line (need 2 colons): {line}")
            
            reactants = parse_stoich(parts[0])
            products = parse_stoich(parts[1])
            rate = float(parts[2])
            rxns.append((reactants, products, rate))

    if not rxns:
        raise ValueError(f"No reactions loaded from {path}")
    return rxns


def load_initial_counts(path: Path) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            toks = line.split()
            if len(toks) < 2:
                continue
            sp = toks[0]
            val = int(toks[1])
            counts[sp] = val

    if not counts:
        raise ValueError(f"No initial values loaded from {path}")
    return counts


def propensity(counts: Dict[str, int], reactants: Stoich, rate: float) -> float:
    a = rate
    for sp, m in reactants.items():
        n = counts.get(sp, 0)
        c = nCk(n, m)
        if c == 0:
            return 0.0
        a *= c
    return float(a)


def classify(counts: Dict[str, int]) -> str | None:
    stealth = counts.get("cI2", 0) > STEALTH_THRESHOLD # If dictionary value of cI2 is already greater than STEALTH_THRESHOLD then stealth = 1, otherwise 0
    hijack = counts.get("Cro2", 0) > HIJACK_THRESHOLD # If dictionary value of Cro2 is already greater than HIJACK_THRESHOLD then stealth = 1, otherwise 0
    
    # Return the condition of the state
    if stealth and hijack:
        return "tie"
    if stealth:
        return "stealth"
    if hijack:
        return "hijack"
    return None


def run_one(rxns: List[Reaction], init_counts: Dict[str, int], moi_value: int) -> str:
    """
    Run one SSA trajectory until:
      - stealth/hijack/tie reached, or
      - MAX_TIME reached, or
      - MAX_STEPS reached, or
      - no reactions can fire
    """
    counts = dict(init_counts) # Copy the initial counts dictionary (number of molecules of each specie)
    counts["MOI"] = moi_value  # Override MOI from init file
    t = 0.0 # Reset time

    out0 = classify(counts) # Run clasify on counts
    if out0 is not None: # If lambda is in a terminal fate state then give that state
        return out0

    for _step in range(MAX_STEPS): 
        if t >= MAX_TIME: # End if max time has been reached
            break

        props: List[float] = [] # Create props as a list of float variables
        a0 = 0.0 # Intialize a0
        for reactants, products, rate in rxns: # Loop through each reactiong
            a = propensity(counts, reactants, rate) # Get the propensity for each reaction
            props.append(a) # Assign the propensity for each reaction
            a0 += a # Create a running total of reaction rates (TOTAL REACTION RATE)

        if a0 <= 0.0: # Break if no reaction can fire
            break

        # sample time increment using Gillespie algorithm
        r1 = random.random() # Make r1 a random number between 0 and 1
        dt = -math.log(max(r1, 1e-300)) / a0 # Timestep according to Gillespie algorithm
        t += dt # Increment time by time step

        # Choose which reaction fires
        r2 = random.random() * a0 # Make r2 a random 0 to 1 a0
        s = 0.0 # Accumulates propensities
        idx = 0 # Stores chosen reaction index
        for i, a in enumerate(props): # Loop through reaction propensities; i - reaction index, a - reactions propensity
            s += a # Add propensities to s
            if r2 <= s: # Go until the random r2 is less than the accumulated propensities
                idx = i # Assign index of breaking reaction to be the reaction that fires
                break

        reactants, products, _rate = rxns[idx]

        # Apply stoichiometry (reactants consumed, products produced) for each reaction
        for sp, m in reactants.items(): 
            counts[sp] = counts.get(sp, 0) - m 
        for sp, m in products.items():
            counts[sp] = counts.get(sp, 0) + m

        out = classify(counts) # Return the current state of lambda
        if out is not None:
            return out

    return "neither"


def main() -> None:
    if SEED is not None:
        random.seed(SEED)

    here = Path(__file__).resolve().parent
    reactions_path = here / REACTIONS_FILENAME
    init_path = here / INIT_FILENAME

    rxns = load_reactions(reactions_path) # Read reactions file
    init_counts = load_initial_counts(init_path) # Read input file

    print(f"Trials/MOI={TRIALS_PER_MOI}, MAX_TIME={MAX_TIME}, MAX_STEPS={MAX_STEPS}, SEED={SEED}")
    print("MOI   P(stealth_first)   P(hijack_first)   P(tie)   P(neither)")
    print("----  -----------------  ---------------   ------   ---------")

    for moi in MOI_VALUES: # Run through MOI values 1 to 10
        nS = nH = nT = nN = 0 # 
        for _ in range(TRIALS_PER_MOI):
            out = run_one(rxns, init_counts, moi) 
            if out == "stealth":
                nS += 1
            elif out == "hijack":
                nH += 1
            elif out == "tie":
                nT += 1
            else:
                nN += 1

        pS = nS / TRIALS_PER_MOI
        pH = nH / TRIALS_PER_MOI
        pT = nT / TRIALS_PER_MOI
        pN = nN / TRIALS_PER_MOI

        print(f"{moi:>3d}   {pS:>16.4f}     {pH:>13.4f}   {pT:>6.4f}   {pN:>8.4f}")


if __name__ == "__main__":
    main()