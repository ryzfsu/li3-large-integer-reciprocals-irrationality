# Irrationality of the trilogarithm at reciprocals of large integers

Ruiyi Zhang

The manuscript proves that `Li_3(1/a)` is irrational for every integer `a >= 67` and every
integer `a <= -70`.

For a multiple `m` of 4 an `m x m` determinant is formed whose entries are linear forms in 1
and `Li_3(1/a)`. A common denominator of the coefficients is bounded, nonvanishing is proved
along a prime-power subsequence, and the size of the determinant is bounded by a logarithmic
energy. Exact density tables with directed rounding give the energy bounds at `a = 67` and
`a = -70`; a slope comparison extends them to the full ranges.

## Contents

- `main.pdf` (29 pages)
- `main.tex`, `body.tex`, `appendices.tex` - LaTeX source
- `anc/` - the ancillary checks of Appendix C, their expected outputs and `results.txt`.
  See `anc/README.md`.
- `SHA256SUMS`

## Build

    pdflatex main.tex
    pdflatex main.tex

## Checks

Python 3.12 with python-flint 0.9.0 (needed by `check_energy.py` and `cert.py` only).
Run from `anc/`; the scripts read the tables from `../main.tex` and its input files:

    cd anc
    python3 check_energy.py      # energy VERIFY_OK (about 5 minutes)
    python3 check_seed.py
    python3 check_tau.py
    python3 check_pole_parts.py  # diagnostic
    python3 check_tower.py       # diagnostic
