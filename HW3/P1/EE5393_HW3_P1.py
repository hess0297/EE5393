import numpy as np
from math import comb

def poly_to_bernstein_exact(coeffs):
    n = len(coeffs) - 1
    b = np.zeros(n + 1, dtype=float)

    for k in range(n + 1):
        for i in range(k + 1):
            b[k] += coeffs[i] * comb(k, i) / comb(n, i)

    return b

def bernstein_output(b, x):
    """
    b : list of Bernstein coefficients [b0, b1, ..., bn]
    x : input value in [0,1]

    returns: circuit output P(out=1)
    """
    n = len(b) - 1

    # compute Bernstein sum
    y = 0.0
    for k in range(n + 1):
        # binomial coefficient
        c = 1
        for i in range(1, k + 1):
            c = c * (n - i + 1) / i

        y += b[k] * c * (x**k) * ((1 - x)**(n - k))

    return y

# ----------------------------
# Problem 1(a)
# f(x) = x - x^2/4
# ----------------------------

poly_coeffs_1a = [0, 1, -1/4]   # a0 + a1*x + a2*x^2
bern_coeffs_1a = poly_to_bernstein_exact(poly_coeffs_1a)

print("Problem 1a:")
print(bern_coeffs_1a)
print()

# ----------------------------
# Problem 1(b)
# cos(x) approximated by p(x) = 1 - x^2/2
# ----------------------------

poly_coeffs_1b = [1, 0, -1/2]
bern_coeffs_1b = poly_to_bernstein_exact(poly_coeffs_1b)

print("Problem 1b:")
print(bern_coeffs_1b)
print()

# ----------------------------
# Problem 1(c)
# p(t) = 31/32 t^5 + 5/32 t^4 - 5/8 t^3 + 5/4 t^2 - 5/4 t + 1/2
# ----------------------------

poly_coeffs_1c = [1/2, -5/4, 5/4, -5/8, 5/32, 31/32]
bern_coeffs_1c = poly_to_bernstein_exact(poly_coeffs_1c)

print("Problem 1c:")
print(bern_coeffs_1c)
print()

print("X = 0:")
print(bernstein_output(bern_coeffs_1c, 0))
print()

print("X = 0.25")
print(bernstein_output(bern_coeffs_1c, 0.25))
print()

print("X = 0.5")
print(bernstein_output(bern_coeffs_1c, 0.5))
print()

print("X = 0.75")
print(bernstein_output(bern_coeffs_1c, 0.75))
print()

print("X = 1")
print(bernstein_output(bern_coeffs_1c, 1))
print()


