// tests/allelic_helpers.ts
// Run with: web/node_modules/.bin/esbuild tests/allelic_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/allelic_helpers.mjs
//           node tests/_bundle/allelic_helpers.mjs
import { DEFAULT_SETTINGS, DEFAULT_VIEW, agreementRows, allelicFromSession, allelicSegmentRows, bpDistance, callShares, comparisonTsv, depthSegmentRows, matchBreakpoints, sameRun, settingsProblems, sharesText, type AllelicBreakpoint, type AllelicSegment } from '../web/src/allelic'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}
const abp = (source: string, chromosome: string, bin: number, gap: [number, number] | null = null): AllelicBreakpoint => ({ source, chromosome, absolute_bin: bin, bp_genomic_start: 0, bp_genomic_end: 0, gap_start: gap ? gap[0] : null, gap_end: gap ? gap[1] : null, recursion_depth: 0, ad_score: 1, p_local: 0.001, p_global: 0.001 })
const dbp = (source: string, chromosome: string, bin: number) => ({ source, chromosome, absolute_bin: bin, bp_genomic_start: 0, bp_genomic_end: 0, recursion_depth: 0, ad_score: 1, p_local: 0.001, p_global: 0.001, n_cells: 1 })

eq(settingsProblems(DEFAULT_SETTINGS), [], 'defaults are valid')
eq(settingsProblems({ ...DEFAULT_SETTINGS, alpha: 0, k: 2.5 }), ['k: must be a whole number', 'alpha: must be > 0 and ≤ 0.1'], 'settings problems name the field')
eq(sameRun({ signal: 'abs_dbaf', settings: DEFAULT_SETTINGS }, 'abs_dbaf', { ...DEFAULT_SETTINGS }), true, 'same run')
eq(sameRun({ signal: 'abs_dbaf', settings: DEFAULT_SETTINGS }, 'dbaf', DEFAULT_SETTINGS), false, 'different signal')
eq(sameRun(null, 'dbaf', DEFAULT_SETTINGS), false, 'no run')
eq(sameRun({ signal: 'dbaf', settings: DEFAULT_SETTINGS }, 'dbaf', { ...DEFAULT_SETTINGS, global_pool: 'chromosome' }), false, 'pool is part of the run')
eq(settingsProblems({ ...DEFAULT_SETTINGS, global_pool: 'cell' as any }), ['global_pool: must be genome or chromosome'], 'pool checked')

eq(bpDistance(50, { bin: 53, gap_start: null, gap_end: null }), 3, 'plain distance')
eq(bpDistance(47, { bin: 50, gap_start: 45, gap_end: 50 }), 0, 'inside a gap is 0')
eq(bpDistance(41, { bin: 50, gap_start: 45, gap_end: 50 }), 4, 'left of a gap: distance to its start')
eq(bpDistance(53, { bin: 50, gap_start: 45, gap_end: 50 }), 3, 'right of a gap: distance to its end')

const m = matchBreakpoints([10, 12, 40], [{ bin: 11, gap_start: null, gap_end: null }, { bin: 44, gap_start: null, gap_end: null }, { bin: 90, gap_start: null, gap_end: null }], 3)
eq(m.pairs.map((p) => [p.d, p.a, p.dist]), [[0, 0, 1]], 'one-to-one: tie on distance goes to the lower depth bin')
eq([m.depthOnly, m.allelicOnly], [[1, 2], [1, 2]], 'unmatched on both sides')
const g = matchBreakpoints([47, 60], [{ bin: 50, gap_start: 45, gap_end: 50 }, { bin: 58, gap_start: null, gap_end: null }], 3)
eq(g.pairs.map((p) => [p.d, p.a, p.dist]), [[0, 0, 0], [1, 1, 2]], 'gap match at distance 0, then nearest')

const depthBps = [dbp('cluster1', 'chr1', 10), dbp('cluster1', 'chr1', 40), dbp('cluster1', 'chr2', 5), dbp('cluster2', 'chr1', 7)]
const allBps = [abp('cluster1', 'chr1', 11), abp('cluster1', 'chr2', 30), abp('cluster1', 'chr2', 6)]
eq(agreementRows(['cluster1', 'cluster2'], ['chr1', 'chr2'], depthBps, allBps, 3), [
  { source: 'cluster1', depth: 3, allelic: 3, matched: 2, depthOnly: 1, allelicOnly: 1, medianOffset: 1 },
  { source: 'cluster2', depth: 1, allelic: 0, matched: 0, depthOnly: 1, allelicOnly: 0, medianOffset: null }], 'agreement rows over chromosomes')
eq(agreementRows(['cluster1'], ['chr2'], depthBps, allBps, 3)[0].matched, 1, 'agreement rows for one chromosome')

const segs = [{ group: 'cluster1', s: 0, e: 20, cn_call: 1 }, { group: 'cluster1', s: 20, e: 60, cn_call: 1.5 }]
const values: Record<string, (number | null)[]> = { cluster1: Array.from({ length: 60 }, (_, i) => (i < 20 ? (i % 2 ? 0.1 : null) : 0.3)) }
eq(depthSegmentRows(segs, values, [abp('cluster1', 'chr1', 40), abp('cluster1', 'chr1', 22)], 'chr1', 3), [
  { source: 'cluster1', bin_start: 0, bin_end: 19, cn_call: 1, kept: 10, mean: 0.1, sd: 0, allelicInside: 0 },
  { source: 'cluster1', bin_start: 20, bin_end: 59, cn_call: 1.5, kept: 40, mean: 0.3, sd: 0, allelicInside: 1 }], 'depth segments: kept, mean, sd, allelic bps more than N from both edges')

eq(callShares({ bin_start: 10, bin_end: 29 }, segs), [{ call: 1, frac: 0.5 }, { call: 1.5, frac: 0.5 }], 'shares: equal split, smaller call first')
eq(sharesText([{ call: 1, frac: 0.8 }, { call: 1.5, frac: 0.2 }]), '1.0: 80% · 1.5: 20%', 'shares text')
eq(sharesText([]), 'no depth call', 'no overlap')
const aseg: AllelicSegment = { source: 'cluster1', chromosome: 'chr1', bin_start: 0, bin_end: 29, n_bins: 30, kept_bins: 25, mean: 0.12, median: 0.1, at_sum: 1, ha_sum: 0, hb_sum: 1 }
eq(allelicSegmentRows([aseg, { ...aseg, chromosome: 'chr2' }], segs, ['cluster1'], 'chr1'), [{ source: 'cluster1', bin_start: 0, bin_end: 29, kept: 25, mean: 0.12, spanned: '1.0: 67% · 1.5: 33%' }], 'allelic segment rows for the chosen chromosome')

const tsv = comparisonTsv(['signal: abs_dbaf'], agreementRows(['cluster2'], ['chr1'], depthBps, allBps, 3), [], [])
eq(tsv.split('\n').slice(0, 4), ['# signal: abs_dbaf', '# table 1: breakpoint agreement', 'source\tdepth_bps\tallelic_bps\tmatched\tdepth_only\tallelic_only\tmedian_offset', 'cluster2\t1\t0\t0\t1\t0\tNA'], 'TSV: header, table 1')
eq(tsv.includes('# table 2: depth segments seen by the allelic signal') && tsv.includes('# table 3: allelic segments seen by depth'), true, 'TSV: all three tables')

eq(allelicFromSession(undefined), { view: DEFAULT_VIEW, notice: null }, 'no saved state: defaults')
const saved = { ...DEFAULT_VIEW, signal: 'dbaf', chrom: 'chr3', sources: ['cluster4'], depthSource: 'original', matchN: 5, scope: 'genome', settings: { ...DEFAULT_SETTINGS, alpha: 0.01 } }
eq(allelicFromSession(JSON.parse(JSON.stringify(saved))), { view: saved, notice: null }, 'round trip')
const junk = allelicFromSession({ signal: 'baf', settings: { alpha: 7 }, chrom: 3, matchN: -1 })
eq(junk.view, DEFAULT_VIEW, 'junk: defaults')
eq(junk.notice?.startsWith('Saved allelic settings could not be applied: '), true, 'junk: notice')

if (failures) { console.error(`${failures} failure(s)`); process.exit(1) } else console.log('all allelic helper checks passed')
