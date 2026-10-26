# Rule 0 Genomic-Gap Fix Validation

## Real-Data Counts

Command:
```bash
python3 - <<'EOF'
import h5py, numpy as np
h = h5py.File('pyepi_results/count_matrix.h5ad', 'r')
s = h['var/start'][...].astype(np.int64); e = h['var/end'][...].astype(np.int64)
seq = h['var/seqnames'][...] if 'seqnames' in h['var'] else None
d = s[1:] - e[:-1]
print('boundaries', len(d), 'old rule (start != end)', int((s[1:] != e[:-1]).sum()), 'new rule (start > end+1)', int((d > 1).sum()))
EOF
```

Output:
```
boundaries 25336 old rule (start != end) 25336 new rule (start > end+1) 528
```

## Analysis

- **Total boundaries**: 25336
- **Old rule flagged as gap**: 25336 (every boundary)
- **New rule flagged as gap**: 528

Before this fix every internal boundary was flagged as a genomic gap; revision 1's 'blocked by gap: 13' on chr8 came from this bug.
