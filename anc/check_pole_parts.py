#!/usr/bin/env python3
"""Pure-Python exact replay of the removable-pole P'' formulas."""
from fractions import Fraction as Q
from math import factorial

def add(a, b):
    n = max(len(a), len(b)); out = [Q(0)] * n
    for i, x in enumerate(a): out[i] += x
    for i, x in enumerate(b): out[i] += x
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out

def scale(a, c): return [x*c for x in a]

def mul(a, b):
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] += x*y
    return out

def deriv(a, k=1):
    for _ in range(k):
        a = [Q(i+1)*a[i+1] for i in range(len(a)-1)] or [Q(0)]
    return a

def evalp(a, x):
    z = Q(0)
    for c in reversed(a): z = z*x + c
    return z

def div_monic(a, q):
    a = a[:]; d = len(q)-1; out = [Q(0)] * (len(a)-d)
    for k in range(len(a)-1, d-1, -1):
        c = a[k]; out[k-d] = c
        if c:
            for j in range(d+1): a[k-d+j] -= c*q[j]
    if any(a[k] for k in range(d)): raise AssertionError(("remainder", a[:d]))
    return out

def run(m):
    g = m // 2
    base = [Q(1)]
    for s in range(1, g+1): base = mul(base, [Q(-s), Q(1)])
    w = scale(mul(mul(base, base), base), Q(1, factorial(g)**3))
    wd = [deriv(w, k) for k in range(4)]
    checks = 0
    for i in range(m):
        for j in range(m, 2*m):
            d = Q(j-i)
            wi, wj = evalp(w, -i), evalp(w, -j)
            ci, cj = wi/d, -wj/d
            # P = (w - ci*(r+j) - cj*(r+i))/((r+i)(r+j)).
            num = w[:]
            num = add(num, [-(ci*j + cj*i), -(ci+cj)])
            den = [Q(i*j), Q(i+j), Q(1)]
            p = div_monic(num, den)
            p2 = deriv(p, 2)
            lhs_i, lhs_j = evalp(p2, -i), evalp(p2, -j)
            rhs_i = (evalp(wd[3], -i)/(3*d) - evalp(wd[2], -i)/d**2
                     + 2*evalp(wd[1], -i)/d**3
                     + 2*(wj-wi)/d**4)
            rhs_j = (-evalp(wd[3], -j)/(3*d) - evalp(wd[2], -j)/d**2
                     - 2*evalp(wd[1], -j)/d**3
                     + 2*(wi-wj)/d**4)
            if lhs_i != rhs_i or lhs_j != rhs_j:
                raise AssertionError((m, i, j, lhs_i, rhs_i, lhs_j, rhs_j))
            checks += 2
    return checks

total = 0
for m in range(4, 33, 4):
    c = run(m); total += c
    print(f"m={m} checks={c} PASS", flush=True)
print(f"ALL_PASS total={total}")
