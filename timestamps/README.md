# Timestamps

[OpenTimestamps](https://opentimestamps.org) proofs for the assets of release `v1`.
Each proof commits to the SHA-256 of the file below and is anchored in Bitcoin block
971169 (2026-10-13 01:55:28 UTC).

| file | SHA-256 |
|---|---|
| `li3-large-integer-reciprocals-irrationality-v1.pdf` | `5dfc5f9677eb41e036e61397978cc0083897d09950a32acc88e0bdc284a513bb` |
| `li3-large-integer-reciprocals-irrationality-v1.tar.gz` | `feb8b11d8046b4942377447f01395e68b04f66554920b852bf92159255129956` |

To verify, download the release asset into this directory and run

    pip install opentimestamps-client
    ots verify timestamps/<file>.ots
