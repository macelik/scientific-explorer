# Archived vs live flank scores (pilot chromosomes)

Command (server on port 8767):

    BASE_URL=http://127.0.0.1:8767 python3 tests/compare_archived_live.py

Pilot chromosomes: chr1, chr2, chr5, chr7, chr8, chr9, chr13. Full output: `archived_vs_live.json`.

| Archived column | Live column | Matched / archived rows | Compared values | Max abs diff |
|---|---|---|---|---|
| flank_scores.delta | delta_median | 162/162 | 324 | 5e-05 |
| flank_scores.delta_mean | delta_mean | 162/162 | 324 | 4.97e-05 |
| flank_scores.delta_trim | delta_iqr_mean | 162/162 | 324 | 4.99e-05 |
| flank_scores.srd_phi | srd_phi | 162/162 | 324 | 4.99e-05 |
| screening.log2fc | delta_mean | 196/196 | 392 | 4.99e-05 |
| screening.srd | srd | 196/196 | 392 | 4.94e-05 |

## Conclusion

Every archived segment was matched by a live segment (no unmatched segments), and every compared value agrees.
Each pair differs by at most about 5e-05, consistent with the archived tables being rounded to four decimals
(for example archived 1.7282 vs live 1.72815003 for cluster3 chr2 81-125, left). No pair agrees to the strict
1e-6 floating-point-noise threshold, so strictly each differs by up to ~5e-05; the differences are rounding-sized
and no larger. Track values on pilot chromosomes therefore change by at most ~5e-05 when the tracks move to
live scores.
