#!/usr/bin/env python3
"""Recompute the exact seed determinant and check Table B.1."""

from fractions import Fraction as Q
from itertools import permutations
from pathlib import Path
import argparse
import json
import re


HERE = Path(__file__).resolve().parent


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(a, b):
    out = [Q(0)] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return trim(out)


def scale(a, c):
    return trim([x * c for x in a])


def mul(a, b):
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def deriv(a, order=1):
    for _ in range(order):
        a = [Q(i + 1) * a[i + 1] for i in range(len(a) - 1)] or [Q(0)]
    return trim(a)


def value(a, x):
    out = Q(0)
    for c in reversed(a):
        out = out * x + c
    return out


def power(a, n):
    out = [Q(1)]
    for _ in range(n):
        out = mul(out, a)
    return out


def div_monic(a, divisor):
    a = list(a)
    degree = len(divisor) - 1
    quotient = [Q(0)] * max(1, len(a) - degree)
    for k in range(len(a) - 1, degree - 1, -1):
        c = a[k]
        quotient[k - degree] = c
        if c:
            for j, d in enumerate(divisor):
                a[k - degree + j] -= c * d
    if any(a[i] for i in range(min(degree, len(a)))):
        raise AssertionError("polynomial division remainder")
    return trim(quotient)


def quotient_monic(a, divisor):
    a = list(a)
    degree = len(divisor) - 1
    quotient = [Q(0)] * max(1, len(a) - degree)
    for k in range(len(a) - 1, degree - 1, -1):
        c = a[k]
        quotient[k - degree] = c
        if c:
            for j, d in enumerate(divisor):
                a[k - degree + j] -= c * d
    return trim(quotient)


def determinant(matrix):
    total = [Q(0)]
    for permutation in permutations(range(len(matrix))):
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(len(matrix))
                         for j in range(i + 1, len(matrix)))
        term = [Q((-1) ** inversions)]
        for i, j in enumerate(permutation):
            term = mul(term, matrix[i][j])
        total = add(total, term)
    return total


def seed_polynomial():
    variable = [Q(0), Q(1)]
    weight = scale(power([Q(2), Q(-3), Q(1)], 3), Q(1, 8))
    base = add(variable, [Q(-1)])
    matrix = []
    for i in range(4):
        row = []
        for j in range(4, 8):
            denominator = [Q(i * j), Q(i + j), Q(1)]
            polynomial = quotient_monic(weight, denominator)
            second = deriv(polynomial, 2)
            values = [value(second, Q(r)) for r in range(1, len(second) + 1)]
            differences = []
            while values:
                differences.append(values[0])
                values = [values[k + 1] - values[k] for k in range(len(values) - 1)]

            entry = [Q(0)]
            for order, coefficient in enumerate(differences):
                entry = add(entry, scale(power(base, 2 - order), coefficient))

            delta = j - i
            ci = value(weight, Q(-i)) / delta
            cj = -value(weight, Q(-j)) / delta
            for pole, residue in ((i, ci), (j, cj)):
                finite = [Q(0)]
                for s in range(1, pole + 1):
                    monomial = [Q(0)] * (pole - s + 1)
                    monomial[pole - s] = residue * 2 / (s ** 3)
                    finite = add(finite, monomial)
                entry = add(entry, scale(mul(finite, power(base, 3)), Q(-1)))
            row.append(entry)
        matrix.append(row)

    denominator = 2 ** 18 * 3 ** 9 * 5 ** 12 * 7 ** 4
    numerator = scale(determinant(matrix), Q(denominator))
    divisor = power(base, 8)
    result = div_monic(numerator, divisor)
    if any(c.denominator != 1 for c in result):
        raise AssertionError("seed coefficients are not integral")
    return [int(c) for c in result]


def tex_source(path, seen=None):
    seen = seen or set()
    path = path.resolve()
    if path in seen:
        raise RuntimeError("repeated TeX input")
    seen.add(path)
    text = path.read_text()
    pattern = re.compile(r"\\input\{([^}]+)\}")
    def expand(match):
        name = match.group(1)
        if not name.endswith(".tex"):
            name += ".tex"
        return tex_source(path.parent / name, seen)
    return pattern.sub(expand, text)


def integer_field(value):
    value = value.strip().replace(r"\,", "").replace(r"\thinspace", "")
    value = value.replace("$", "").replace("{", "").replace("}", "")
    value = value.replace(" ", "").replace("~", "")
    value = value.replace(r"\textminus", "-").replace("−", "-")
    if not re.fullmatch(r"-?\d+", value):
        return None
    return int(value)


def parse_table(tex):
    label = tex.index(r"\label{tab:seed-coefficients}")
    candidates = []
    for env in ("tabular", "longtable"):
        start = tex.rfind(r"\begin{" + env + "}", 0, label)
        if start >= 0:
            end = tex.find(r"\end{" + env + "}", label)
            if end >= 0:
                candidates.append((start, end))
    if not candidates:
        raise AssertionError("seed table not found")
    start, end = max(candidates)
    rows = []
    for line in tex[start:end].splitlines():
        line = line.split("%", 1)[0].strip()
        if not line.endswith(r"\\"):
            continue
        fields = line[:-2].split("&")
        if len(fields) != 2:
            continue
        j, coefficient = (integer_field(x) for x in fields)
        if j is not None and coefficient is not None:
            rows.append((j, coefficient))
    if [j for j, _ in rows] != list(range(23)):
        raise AssertionError("seed table rows are not 0 through 22")
    return [coefficient for _, coefficient in rows]


def main():
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument("--tex", help="main.tex for an alternate source layout")
    parser.add_argument("--seed-json", help="optional reference seed record")
    args = parser.parse_args()
    source = Path(args.tex) if args.tex else HERE.parent / "main.tex"
    if not source.exists():
        source = HERE / "main.tex"
    printed = parse_table(tex_source(source))
    exact = seed_polynomial()
    if printed != exact:
        raise AssertionError("printed seed coefficients differ from exact determinant")
    residues = [sum(c * pow(a, j, 13) for j, c in enumerate(exact)) % 13
                for a in range(13)]
    if any(value == 0 for value in residues):
        raise AssertionError("a seed residue vanishes modulo 13")
    if args.seed_json:
        reference = json.loads(Path(args.seed_json).read_text())
        if (reference["coefficients_ascending"] != exact or
                reference["mod13_values"] != residues):
            raise AssertionError("exact seed differs from reference record")
        print("reference seed record=PASS")
    print("seed degree=22 coefficients=PASS residues_mod_13=PASS")


if __name__ == "__main__":
    main()
