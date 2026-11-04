// tests/cn_recall_helpers.ts
// Run with: web/node_modules/.bin/esbuild tests/cn_recall_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/cn_recall_helpers.mjs
//           node tests/_bundle/cn_recall_helpers.mjs
import { UNCHANGED, karyoCode, clusterStateRow, stateLabel, transitionText, recallProgressText, tsv, recallHeader, mergedBreakpointsTsv, recalledClusterSegmentsTsv, recalledPercellTsv, changeSummaryTsv, type RecallStatus } from '../web/src/cnRecall'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}

eq(karyoCode(2), 2, 'no before: after code')
eq(karyoCode(2, 2), UNCHANGED, 'same state: unchanged')
eq(karyoCode(3, 2), 3, 'changed: after code')
eq(karyoCode(255, 255), 255, 'missing in both stays missing')
eq(karyoCode(255, 2), 255, 'called -> missing shows missing')

const segs = [
  { group: 'cluster0', n_cells: 5, best_s: 0.5, segment: 0, chromosome: 'chr1', bin_start: 0, bin_end: 2, n_bins: 3, continuous_cn: 1, integer_cn: 2, watson: 1, cn_call: 1 },
  { group: 'cluster0', n_cells: 5, best_s: 0.5, segment: 1, chromosome: 'chr1', bin_start: 3, bin_end: 4, n_bins: 2, continuous_cn: 1.5, integer_cn: 3, watson: 1, cn_call: 1.5 },
  { group: 'all', n_cells: 9, best_s: 0.6, segment: 0, chromosome: 'chr1', bin_start: 0, bin_end: 5, n_bins: 6, continuous_cn: 0, integer_cn: 1, watson: 0, cn_call: 0 },
]
eq([...clusterStateRow(segs, 'cluster0', 6)], [2, 2, 2, 3, 3, 255], 'cluster state row')
eq(stateLabel(3), '1.5 (weak gain)', 'state label')
eq(stateLabel(255), 'not called', 'missing label')
eq(transitionText({ before: 1, after: 1.5, bins: 12 }), '1 (baseline) → 1.5 (weak gain): 12 bins', 'transition text')
eq(transitionText({ before: null, after: 1, bins: 1 }), 'not called → 1 (baseline): 1 bins', 'transition from missing')

const st: RecallStatus = {
  id: 'j1', hash: 'h', status: 'running', stage: 'percell', progress: { done: 1240, total: 3460 }, error: null,
  settings: { rules: [{ id: 'rule0' }], allow_gap_crossing: false }, provenance: { scope: 's' },
  merge: { pairs: [], merged_breakpoints: [{ source: 'cluster0', chromosome: 'chr1', absolute_bin: 3, recursion_depth: 1 }] },
  cluster_segments: segs,
  diff: {
    sources: [{ source: 'cluster0', segments_before: 3, segments_after: 2, merges: 1, phase2: 0, best_s_before: 0.51, best_s_after: 0.5, bins_changed: 2, bins_called: 6, pct_changed: 33.333 }],
    cells: { n_cells: 2, n_changed_cells: 1, median_pct_changed: 5, best_s_change: { median: 0, min: -0.01, max: 0.01 },
      per_cell: [{ cell: 'a', wnn_cluster: 0, best_s_before: 0.5, best_s_after: 0.49, bins_changed: 3 }, { cell: 'b', wnn_cluster: 1, best_s_before: null, best_s_after: 0.6, bins_changed: 0 }] },
    transitions: { clusters: [{ before: 1, after: 1.5, bins: 2 }], cells: [{ before: 0, after: 1, bins: 3 }] },
  },
  percell_ready: false,
}
eq(recallProgressText(st), 'merging ✓ → cluster CN ✓ → per-cell 1,240 / 3,460', 'progress percell')
eq(recallProgressText({ ...st, stage: 'merge', progress: { done: 3, total: 144 } }), 'merging… 3 / 144', 'progress merge')
eq(recallProgressText({ ...st, stage: 'cluster', progress: { done: 250, total: 3460 } }), 'merging ✓ → cluster CN… 250 / 3,460', 'progress cluster')
eq(recallProgressText({ ...st, status: 'done', stage: 'done' }), 'merging ✓ → cluster CN ✓ → per-cell ✓ · done', 'progress done')
eq(recallProgressText({ ...st, status: 'error', error: 'percell: boom' }), 'error: percell: boom', 'progress error')
eq(recallProgressText(null), '', 'progress null')

eq(tsv(['a: 1'], ['x', 'y'], [{ x: 1, y: null }, { x: 'q', y: 2 }]), '# a: 1\nx\ty\n1\t\nq\t2\n', 'tsv writer')
eq(recallHeader(st), ['rules: [{"id":"rule0"}]', 'allow_gap_crossing: false', 'provenance: {"scope":"s"}'], 'header lines')
eq(mergedBreakpointsTsv(st).split('\n').slice(3, 5), ['source\tchromosome\tabsolute_bin\trecursion_depth', 'cluster0\tchr1\t3\t1'], 'merged breakpoints tsv')
eq(recalledClusterSegmentsTsv(st).split('\n')[3], 'group\tn_cells\tbest_s\tsegment\tchromosome\tbin_start\tbin_end\tn_bins\tcontinuous_cn\tinteger_cn\twatson\tcn_call', 'cluster segments columns')
eq(recalledPercellTsv(st).split('\n').slice(3, 6), ['cell\twnn_cluster\tbest_s_before\tbest_s_after\tbins_changed', 'a\t0\t0.5\t0.49\t3', 'b\t1\t\t0.6\t0'], 'percell tsv')
const summary = changeSummaryTsv(st).split('\n')
eq(summary[3], 'source\tsegments_before\tsegments_after\tmerges\tphase2\tbest_s_before\tbest_s_after\tbins_changed\tbins_called\tpct_changed', 'summary columns')
eq(summary.slice(summary.indexOf('scope\tbefore\tafter\tbins')), ['scope\tbefore\tafter\tbins', 'clusters\t1\t1.5\t2', 'cells\t0\t1\t3', ''], 'summary transitions')

if (failures) { console.error(`${failures} failure(s)`); process.exit(1) } else console.log('all cn_recall helper checks passed')
