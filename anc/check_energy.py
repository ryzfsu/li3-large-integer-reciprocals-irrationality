#!/usr/bin/env python3
"""Check the printed density tables and recompute the directed certificate."""

from fractions import Fraction
from pathlib import Path
import argparse
import json
import re
import time

from cert import DEN, TAU, certify


HERE = Path(__file__).resolve().parent
EXPECTED = {
    67: ("[6.2732861583045446758 +/- 2.16e-20]",
         "[6.2733351903561357608 +/- 3.57e-20]"),
    70: ("[6.3100053850648586710 +/- 2.56e-20]",
         "[6.3100585701803409215 +/- 6.02e-21]"),
}


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


def source_text(arg):
    if arg:
        return tex_source(Path(arg))
    default = HERE.parent / "main.tex"
    if default.exists():
        return tex_source(default)
    local = HERE / "main.tex"
    if local.exists():
        return tex_source(local)
    raise SystemExit("main.tex not found next to anc")


def block(tex, label, environments):
    at = tex.index(r"\label{" + label + "}")
    candidates = []
    for env in environments:
        marker = r"\begin{" + env + "}"
        finish = r"\end{" + env + "}"
        starts = []
        before = tex.rfind(marker, 0, at)
        after = tex.find(marker, at)
        if before >= 0:
            starts.append(before)
        if after >= 0:
            starts.append(after)
        for start in starts:
            end = tex.find(finish, max(start, at))
            if end >= 0:
                candidates.append((abs(start - at), start, end + len(finish)))
    if not candidates:
        raise AssertionError("table block not found")
    _, start, end = min(candidates)
    return tex[start:end]


def integer_field(value):
    value = value.strip()
    value = value.replace(r"\,", "").replace(r"\thinspace", "")
    value = value.replace("$", "").replace("{", "").replace("}", "")
    value = value.replace(" ", "").replace("~", "")
    value = value.replace(r"\textminus", "-").replace("−", "-")
    if not re.fullmatch(r"-?\d+", value):
        return None
    return int(value)


def row_fields(line):
    line = line.split("%", 1)[0].strip()
    if not line.endswith(r"\\"):
        return None
    return [integer_field(v) for v in line[:-2].split("&")]


def parse_density(tex):
    values = {67: [0] * 256, 70: [0] * 256}
    seen = set()

    first = block(tex, "tab:energy-density", ("longtable", "tabular"))
    for line in first.splitlines():
        fields = row_fields(line)
        if fields is None:
            continue
        triples = []
        if len(fields) == 3:
            triples = [fields]
        elif len(fields) == 6:
            triples = [fields[:3], fields[3:]]
        for triple in triples:
            if any(x is None for x in triple):
                continue
            i, n67, n70 = triple
            if i < 0 or i >= 256 or i in seen:
                raise AssertionError("bad or duplicate density index")
            seen.add(i)
            values[67][i] = n67
            values[70][i] = n70

    if "tab:energy-tail" in tex:
        tail = block(tex, "tab:energy-tail", ("tabular", "longtable"))
        for line in tail.splitlines():
            fields = row_fields(line)
            if fields is None or len(fields) != 12:
                continue
            for j in range(0, 12, 3):
                triple = fields[j:j + 3]
                if any(x is None for x in triple):
                    continue
                i, n67, n70 = triple
                if i < 146 or i >= 256 or i in seen:
                    raise AssertionError("bad or duplicate tail index")
                seen.add(i)
                values[67][i] = n67
                values[70][i] = n70

    missing = set(range(256)) - seen
    if any(i < 146 for i in missing):
        raise AssertionError("density indices do not cover the stated table range")
    if any(values[q][i] != 0 for q in (67, 70) for i in missing):
        raise AssertionError("omitted tail entries must be zero")
    return values


def parse_level(tex, a):
    match = re.search(rf"\\(?:ell|Lambda)_\{{{a}\}}\s*=\s*(\d+)\s*/\s*(\d+)", tex)
    if not match:
        raise AssertionError("missing certificate level")
    return Fraction(int(match.group(1)), int(match.group(2)))


def main():
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument("--tex", help="main.tex for an alternate source layout")
    parser.add_argument("--density-json", help="optional reference integer table")
    args = parser.parse_args()
    started = time.perf_counter()
    tex = source_text(args.tex)
    values = parse_density(tex)
    if args.density_json:
        reference = json.loads(Path(args.density_json).read_text())
        rows = {record["a"]: record["n"] for record in reference["endpoints"]}
        if reference["sigma_den"] != DEN or rows != values:
            raise AssertionError("printed density tables differ from reference integers")
        print("reference density integers=PASS")
    levels = {a: parse_level(tex, a) for a in (67, 70)}
    for a in (67, 70):
        n = values[a]
        if len(n) != 256 or sum(n) != 128 * DEN:
            raise AssertionError("density mass failed")
        if any(x < 0 or x > 3 * DEN for x in n):
            raise AssertionError("density cap failed")
        result = certify(a, 128, 256, n, levels[a], 256, TAU[a])
        if (result["lower_certified_str"], result["upper_certified_str"]) != EXPECTED[a]:
            raise AssertionError("directed certificate changed")
        if not result["tail_ok"] or not result["margin_positive"]:
            raise AssertionError("certificate margin or tail failed")
        print(f"a={a} rows=256 mass=PASS cap=PASS tail=PASS intervals=PASS")
        print(f"  lower {result['lower_certified_str']}")
        print(f"  upper {result['upper_certified_str']}")
    print(f"energy VERIFY_OK elapsed_seconds={time.perf_counter() - started:.3f}")


if __name__ == "__main__":
    main()
