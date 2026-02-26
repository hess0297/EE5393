"""
Deterministic (sequenced) CRN simulator.
Simulates the chemical reaction network until:
- target value (y) is reached.
- No more reactions can occur.
"""

# -----------------------------
# Initial conditions
# -----------------------------
species = {
    "a": 0,
    "b": 1,
    "c": 0,
    "x": 1234567,  # Initial x value
    "xP": 0,
    "w": 0,
    "y": 1,
    "d": 0,
    "yP": 0,
}

# -----------------------------
# Reaction functions
# -----------------------------
def rb_b_to_a_and_b(s):
    # b -> a + b
    s["a"] += 1
    s["b"] += 1

def r1_a_and_2x_to_c_and_xp_and_a(s):
    # a + 2x -> c + xP + a   (a catalytic)
    s["c"] += 1
    s["xP"] += 1
    s["x"] -= 2

def r2_2c_to_c(s):
    # 2c -> c
    s["c"] -= 1

def r3_a_to_null(s):
    # a -> ∅ (remove one)
    if s["a"] > 0:
        s["a"] -= 1

def r4_xp_to_x(s):
    # xP -> x
    s["xP"] -= 1
    s["x"] += 1

def r5_c_to_w(s):
    # c -> w
    s["w"] += 1
    s["c"] -= 1

def r6_w_to_d(s):
    # w -> d
    s["w"] -= 1
    s["d"] += 1

def r7_d_and_y_to_d_and_2yP(s):
    # d + y -> d + 2yP  (d catalytic)
    s["y"] -= 1
    s["yP"] += 2

def r8_d_to_null(s):
    # d -> ∅ (remove one)
    if s["d"] > 0:
        s["d"] -= 1

def r9_yP_to_y(s):
    # yP -> y
    s["yP"] -= 1
    s["y"] += 1

# -----------------------------
# Simulation
# -----------------------------
def simulate_crn(s, max_steps=2_000_000_000):
    step = 0
    target_y = s["x"]  # Set target y equal to initial x value

    while s["y"] < target_y and step < max_steps:
        step += 1
        reactions_happened = False  # Track if any reaction occurred during this step

        # Apply reactions in a fixed sequence each step
        if s["b"] > 0:  # b -> a + b
            rb_b_to_a_and_b(s)
            reactions_happened = True

        if s["a"] > 0 and s["x"] >= 2:  # a + 2x -> c + xP + a
            r1_a_and_2x_to_c_and_xp_and_a(s)
            reactions_happened = True

        if s["c"] >= 2:  # 2c -> c
            r2_2c_to_c(s)
            reactions_happened = True

        if s["a"] > 0:  # a -> ∅
            r3_a_to_null(s)
            reactions_happened = True

        if s["xP"] > 0:  # xP -> x
            r4_xp_to_x(s)
            reactions_happened = True

        if s["c"] > 0:  # c -> w
            r5_c_to_w(s)
            reactions_happened = True

        if s["w"] > 0:  # w -> d
            r6_w_to_d(s)
            reactions_happened = True

        if s["d"] > 0 and s["y"] > 0:  # d + y -> d + 2yP
            r7_d_and_y_to_d_and_2yP(s)
            reactions_happened = True

        if s["d"] > 0:  # d -> ∅
            r8_d_to_null(s)
            reactions_happened = True

        if s["yP"] > 0:  # yP -> y
            r9_yP_to_y(s)
            reactions_happened = True

        # If no reactions happened and we haven't reached the target, it means we are stuck
        if not reactions_happened:
            print(f"System is stuck at step {step}: y={s['y']} x={s['x']} w={s['w']} a={s['a']} d={s['d']} c={s['c']} yP={s['yP']}")
            
            # Allow one last artificial reaction to kick-start the system
            if s["y"] == target_y and s["d"] == 0:  # Trigger if y is at target, and d is missing
                print("Artificial reaction triggered to re-enable yP growth.")
                # Allow an artificial reaction: Add d (for example, artificially create d to restart reactions)
                s["d"] = 1  # Now `d` is available
                s["y"] -= 1  # Subtract from y to trigger d + y -> d + 2yP
                s["yP"] += 2  # Create yP
                reactions_happened = True
            else:
                break

        # Stop if target is reached and no reactions occurred
        if s["y"] >= target_y and not reactions_happened:
            print(f"Target reached at step {step}, no more reactions.")
            break

    if s["y"] < target_y:
        raise RuntimeError(f"Hit max_steps={max_steps} with y={s['y']} < target_y={target_y}")

    return s, step

# -----------------------------
# Run
# -----------------------------
if __name__ == "__main__":
    final_species, steps = simulate_crn(species)

    print(f"Reached target y={target_y} in {steps} steps\n")
    for name in ["a", "b", "c", "x", "xP", "w", "y", "d", "yP"]:
        print(f"{name}: {final_species[name]}")