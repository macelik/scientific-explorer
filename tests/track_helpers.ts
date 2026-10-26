// tests/track_helpers.ts
// Run with: web/node_modules/.bin/esbuild tests/track_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/track_helpers.mjs
//           node tests/_bundle/track_helpers.mjs
import { archivedColourToTrack, colourValue, derivedRejectionReason, diagnosticsEffectMetric, migrateSegSettings } from '../web/src/trackMetrics'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}
const row: any = { start: 10, end: 15, phi: 1, gap_left: false, gap_right: false,
  srd_left: -4, srd_right: 1.5, srd_phi_left: null, srd_phi_right: 2.0,
  delta_mean_left: -0.3, delta_mean_right: 0.1, delta_median_left: 0.2, delta_median_right: -0.2,
  delta_iqr_mean_left: null, delta_iqr_mean_right: null }

eq(colourValue(row, 'delta_mean'), { value: 0.1, side: 'right' }, 'effect colour: signed value on smaller flank')
eq(colourValue(row, 'delta_median'), { value: 0.2, side: 'left' }, 'effect colour tie goes left')
eq(colourValue(row, 'srd'), { value: 1.5, side: 'right' }, 'stat colour: absolute value on smaller flank')
eq(colourValue(row, 'srd_phi'), { value: 2.0, side: 'right' }, 'missing side loses')
eq(colourValue(row, 'delta_iqr_mean'), { value: null, side: null }, 'both missing -> N/A')

eq(['delta_median', 'delta_mean', 'log2fc', 'delta_trim', 'srd', 'srd_phi', 'bogus', undefined].map(archivedColourToTrack),
  ['delta_median', 'delta_mean', 'delta_mean', 'delta_iqr_mean', 'srd', 'srd_phi', 'delta_mean', 'delta_mean'], 'archived colour ids map to live')
eq(['delta_median', 'delta_mean', 'delta_iqr_mean', 'srd', 'srd_phi'].map((m) => diagnosticsEffectMetric(m as any)),
  ['delta_median', 'delta_mean', 'delta_trim', 'delta_mean', 'delta_mean'], 'diagnostics effect metric mapping')

const old = { chrom: 'chr9', metric: 'log2fc', effectMetric: 'delta_trim', flank: 'min_srd', tau: 0.3, labels: { consistency: false, proposal: true, cn: true } }
eq(migrateSegSettings(old), { chrom: 'chr9', tau: 0.3, labels: { cn: true, signs: false }, colorMetric: 'delta_mean' }, 'old session migrates')
eq(migrateSegSettings({ colorMetric: 'srd', labels: { signs: true } }), { colorMetric: 'srd', labels: { signs: true } }, 'new session passes through')
eq(migrateSegSettings(undefined), {}, 'missing seg -> empty')

const base = { selected_for_significance: true, accepted: false, rejection_reasons: null, p_local: 0.002997, p_global: null }
eq(derivedRejectionReason(base), 'local permutation p = 0.003 ≥ 0.001 (cutoff inferred: all accepted candidates have p < 0.001)', 'local p reason')
eq(derivedRejectionReason({ ...base, p_local: 0.000999, p_global: 0.0123 }), 'global permutation p = 0.0123 ≥ 0.001 (cutoff inferred: all accepted candidates have p < 0.001)', 'global p reason')
eq(derivedRejectionReason({ ...base, p_local: 0.0005, p_global: 0.0005 }), 'not accepted (no recorded reason) (cutoff inferred: all accepted candidates have p < 0.001)', 'no reason found')
eq(derivedRejectionReason({ ...base, rejection_reasons: 'left_child_too_short' }), 'left_child_too_short', 'recorded reason kept')
eq(derivedRejectionReason({ ...base, accepted: true }), null, 'accepted -> no reason')
eq(derivedRejectionReason({ ...base, selected_for_significance: false }), null, 'not tested and no recorded reason -> null')

if (failures) { console.error(`${failures} failure(s)`); process.exit(1) }
console.log('all track_helpers checks passed')
