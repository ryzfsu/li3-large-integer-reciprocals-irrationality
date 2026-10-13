# Ancillary checks

These scripts verify the finite computations listed in the paper. Run them
from this directory in Python 3.12.3. `check_energy.py` and `cert.py` require
python-flint 0.9.0 for Arb interval arithmetic. The other checks use only the
Python standard library. In the source archive, `anc/` is beside `main.tex`;
the table parsers read `../main.tex` and its input files automatically.

| Check | Command | Expected result | Approximate runtime |
| --- | --- | --- | --- |
| Energy tables and directed certificate | `python3 check_energy.py` | `energy VERIFY_OK`; 256 rows, mass, cap, tail, intervals pass at both endpoints | 5 min |
| Seed identity; seed table and residues | `python3 check_seed.py` | `seed degree=22 coefficients=PASS residues_mod_13=PASS` | under 1 s |
| Arithmetic rates | `python3 check_tau.py` | both exact τ bounds pass | under 1 s |
| Pole finite parts (diagnostic) | `python3 check_pole_parts.py` | `ALL_PASS total=6528` | about 7 s |
| Prime-power tower (diagnostic) | `python3 check_tower.py` | `all_pass: true` for both endpoints | about 5 s |

The seed command combines the exact determinant and the printed coefficient
table checks. The energy command checks both printed density tables, the
exact mass and cap, and the directed K=256 certificate. The table layout has
two density blocks: omitted indices in the tail block are interpreted as
zero, as stated in the paper. The pole and tower checks are finite diagnostics;
neither replaces the corresponding uniform proof.

The v8 source hashes and consolidated outputs are in `results.txt`.
