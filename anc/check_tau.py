#!/usr/bin/env python3
"""Exact rational enclosures for the two arithmetic-rate endpoint bounds."""
from fractions import Fraction as F


def log_interval(x, terms=32):
    """Return rational lower/upper bounds for log(x), for 1 <= x <= 2."""
    x = F(x)
    assert 1 <= x <= 2
    z = (x - 1) / (x + 1)
    lower = 2 * sum((z ** (2 * k + 1)) / (2 * k + 1)
                    for k in range(terms))
    remainder = 2 * z ** (2 * terms + 1) / ((2 * terms + 1) * (1 - z * z))
    return lower, lower + remainder


def log_integer(n, log2):
    """Enclose log(n) for 64 <= n <= 127 by reducing to [1, 2]."""
    assert 64 <= n <= 127
    lo, hi = log_interval(F(n, 64))
    return 6 * log2[0] + lo, 6 * log2[1] + hi


def main():
    log2 = log_interval(F(2))
    log54 = log_interval(F(5, 4))
    s_lo = F(643, 96) - F(9, 4) * log2[1] - F(25, 8) * log54[1]
    s_hi = F(643, 96) - F(9, 4) * log2[0] - F(25, 8) * log54[0]

    log66, log67 = log_integer(66, log2), log_integer(67, log2)
    log70, log71 = log_integer(70, log2), log_integer(71, log2)
    tau67_hi = s_hi + F(9, 16) * log66[1] - F(1, 8) * log67[0]
    tau_neg70_hi = s_hi + F(9, 16) * log71[1] - F(1, 8) * log70[0]

    assert tau67_hi < F(6272107, 10**6)
    assert tau_neg70_hi < F(6307708, 10**6)
    print("tau(67) < 6272107/1000000: PASS")
    print("tau(-70) < 6307708/1000000: PASS")
    print("all comparisons used exact Fraction arithmetic")


if __name__ == "__main__":
    main()
