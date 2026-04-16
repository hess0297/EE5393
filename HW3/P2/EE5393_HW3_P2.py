def p_and(a, b):
    return a * b

def p_not(a):
    return 1 - a

def p_or(a, b):
    return 1 - (1 - a) * (1 - b)


# ----------------------------
# Base probabilities
# ----------------------------
A = 0.4
B = 0.5

p04 = A
p05 = B
p06 = p_not(A)          # 0.6
p02 = p_and(A, B)       # 0.2
p03 = p_and(p_not(A), p_not(B))  # 0.3
p07 = p_not(p03)        # 0.7
p08 = p_not(p02)        # 0.8
p01 = p_and(p02, p05)   # 0.1
p09 = p_not(p01)        # 0.9


# ----------------------------
# Problem 2(a)(i): 0.8881188
# ----------------------------
p1 = p_or(p04, p03)     # 0.58
p2 = p_or(p08, p1)      # 0.916
p3 = p_or(p03, p2)      # 0.9412
p4 = p_and(p05, p3)     # 0.4706
p5 = p_or(p03, p4)      # 0.62942
p6 = p_and(p07, p5)     # 0.440594
z1 = p_or(p08, p6)      # 0.8881188

print("(a)(i) P(z) =", round(z1, 7))
print("target      =", 0.8881188)
print()


# ----------------------------
# Problem 2(a)(ii): 0.2119209
# ----------------------------
q1 = p_and(p01, p05)    # 0.05
q2 = p_or(p09, q1)      # 0.905
q3 = p_and(p04, q2)     # 0.362
q4 = p_and(p03, q3)     # 0.1086
q5 = p_or(p05, q4)      # 0.5543
q6 = p_or(p03, q5)      # 0.68801
q7 = p_and(p06, q6)     # 0.412806
q8 = p_or(p05, q7)      # 0.706403
z2 = p_and(p03, q8)     # 0.2119209

print("(a)(ii) P(z) =", round(z2, 7))
print("target       =", 0.2119209)
print()


# ----------------------------
# Problem 2(a)(iii): 0.5555555
# ----------------------------
r1 = p_or(p09, p07)     # 0.97
r2 = p_and(p02, r1)     # 0.194
r3 = p_and(p03, r2)     # 0.0582
r4 = p_or(p05, r3)      # 0.5291
r5 = p_and(p07, r4)     # 0.37037
r6 = p_and(p03, r5)     # 0.111111
z3 = p_or(p05, r6)      # 0.5555555

print("(a)(iii) P(z) =", round(z3, 7))
print("target        =", 0.5555555)
print()

# ----------------------------
# Problem 2(b)
# ----------------------------

# ----------------------------
# Problem 2(b): binary probabilities using S = {0.5}
# ----------------------------

p05 = 0.5

# ----------------------------
# Problem 2(b)(i): 0.1011111_2
# ----------------------------
z7 = p05
z6 = p_or(p05, z7)
z5 = p_or(p05, z6)
z4 = p_or(p05, z5)
z3 = p_or(p05, z4)
z2 = p_and(p05, z3)
z1 = p_or(p05, z2)

print("(b)(i) P(z) =", round(z1, 7))
print("target      =", int("1011111", 2) / (2**7))
print()


# ----------------------------
# Problem 2(b)(ii): 0.1101111_2
# ----------------------------
z7 = p05
z6 = p_or(p05, z7)
z5 = p_or(p05, z6)
z4 = p_or(p05, z5)
z3 = p_and(p05, z4)
z2 = p_or(p05, z3)
z1 = p_or(p05, z2)

print("(b)(ii) P(z) =", round(z1, 7))
print("target       =", int("1101111", 2) / (2**7))
print()


# ----------------------------
# Problem 2(b)(iii): 0.1010111_2
# ----------------------------
z7 = p05
z6 = p_or(p05, z7)
z5 = p_or(p05, z6)
z4 = p_and(p05, z5)
z3 = p_or(p05, z4)
z2 = p_and(p05, z3)
z1 = p_or(p05, z2)

print("(b)(iii) P(z) =", round(z1, 7))
print("target        =", int("1010111", 2) / (2**7))
print()