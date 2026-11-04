// tests/merge_rounds_helpers.ts
// Run with: web/node_modules/.bin/esbuild tests/merge_rounds_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rounds_helpers.mjs
//           node tests/_bundle/merge_rounds_helpers.mjs
import { wireRounds, sameRounds, roundLabel, roundStripText, baseHeading, karyoModes, noiseRows, noiseCellsN, noiseTsv, roundsHeader, roundsFromSession, NOISE_DEFINITIONS, type Noise, type Round } from '../web/src/mergeRounds'
import { recallHeader, type RecallStatus } from '../web/src/cnRecall'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}

const r0 = { id: 'rule0', type: 'flank_threshold', tau: 0.5, srdRef: 8, combine: 'any' }
const noise: Noise = { cluster2: {
  chr1: { cluster: { segments: 10, zigzags: 2, zigzag_median_bins: 12.5, gain_loss_adjacent: 1 }, cells: null },
  genome: { cluster: { segments: 144, zigzags: 33, zigzag_median_bins: 40, gain_loss_adjacent: 6 }, cells: { n_cells: 612, segments_mean: 160.5, zigzags_mean: 31, zigzag_median_bins: null, gain_loss_adjacent_mean: 13.42 } } } }
const round: Round = { rules: [r0, { id: 'u1', enabled: true }, { id: 'u2', enabled: false }], allow_gap_crossing: false, label: 'x', merges: 358, noise }

eq(wireRounds([round]), [{ rules: round.rules, allow_gap_crossing: false }], 'wire drops display fields')
eq(sameRounds([round], [{ allow_gap_crossing: false, rules: round.rules } as any]), true, 'same recipe regardless of extra fields and key order')
eq(sameRounds([round], []), false, 'different recipe')
eq(roundLabel([r0]), 'τ 0.5 · ref 8 · any', 'label: Rule 0 only')
eq(roundLabel(round.rules as any[]), 'τ 0.5 · ref 8 · any + 1 user rule', 'label: enabled user rules counted')
eq(roundStripText([]), 'Original', 'strip: original')
eq(roundStripText([{ ...round, label: 'τ 0.5 · ref 8 · any' }]), 'Original → Round 1 (τ 0.5 · ref 8 · any, 358 merges)', 'strip: one round')
eq(baseHeading(0), 'Base: original', 'heading 0')
eq(baseHeading(1), 'Base: round 1 (1 round of merges on the original)', 'heading 1')
eq(baseHeading(2), 'Base: round 2 (2 rounds of merges on the original)', 'heading 2')
eq(karyoModes(0, false), [{ mode: 'before', label: 'Before', enabled: true }, { mode: 'after', label: 'After', enabled: false }, { mode: 'changed', label: 'Changed', enabled: false }], 'modes: no rounds, no re-call')
eq(karyoModes(1, true).map((m) => [m.mode, m.enabled]), [['original', true], ['before', true], ['after', true], ['changed', true]], 'modes: rounds and a re-call')

const cols = [{ name: 'Original', noise }, { name: 'After', noise: { cluster2: { genome: { cluster: { segments: 97, zigzags: 7, zigzag_median_bins: 107, gain_loss_adjacent: 1 }, cells: null } } }, cellsPending: true }]
const rows = noiseRows(cols, 'cluster2', 'genome')
eq(rows.map((r) => [r.level, r.measure, ...r.values]), [
  ['cluster', 'segments', '144', '97'], ['cluster', 'zig-zags', '33', '7'], ['cluster', 'zig-zag median length (bins)', '40', '107'], ['cluster', 'gain next to loss', '6', '1'],
  ['cells', 'segments per cell', '160.5', 're-calling…'], ['cells', 'zig-zags per cell', '31.0', 're-calling…'], ['cells', 'zig-zag median length (bins)', '—', 're-calling…'], ['cells', 'gain next to loss per cell', '13.4', 're-calling…'],
], 'noise rows: formatting, missing median, pending cells')
eq(noiseRows(cols, 'cluster2', 'chr1')[2].values, ['12.5', '—'], 'noise rows: half-bin median and absent key')
eq(noiseRows([{ name: 'Original', noise: null }], 'cluster2', 'genome')[0].values, ['—'], 'noise rows: no data')
eq(noiseCellsN(cols, 'cluster2', 'genome'), 612, 'cells n from first column with cells')

const hdr = roundsHeader([round])
eq(hdr, [`round 1 rules: ${JSON.stringify(round.rules)} · allow_gap_crossing: false`], 'rounds header')
const t = noiseTsv([{ name: 'original', noise }, { name: 'round1', noise }], ['cluster2'], hdr).split('\n')
eq(t[0], hdr.map((h) => `# ${h}`)[0], 'tsv header line')
eq(t[1], 'source\tchromosome\tlevel\tmeasure\toriginal\tround1', 'tsv columns')
eq(t[2], 'cluster2\tchr1\tcluster\tsegments\t10\t10', 'tsv first row')
eq(t.filter((l) => l.startsWith('cluster2\tgenome\tcells\tzig-zag median')).length, 1, 'tsv cells rows present')
eq(t.find((l) => l.startsWith('cluster2\tchr1\tcells\tsegments')), 'cluster2\tchr1\tcells\tsegments per cell\t\t', 'tsv empty cells values')

eq(roundsFromSession({}), [], 'session without rounds')
eq(roundsFromSession({ rounds: [round] }), [round], 'session with rounds')
eq(roundsFromSession({ rounds: 'x' }), [], 'session with garbage rounds')
const rs = roundsFromSession({ rounds: [null, { allow_gap_crossing: false }, { rules: [r0], allow_gap_crossing: 'no' }, { rules: [r0], allow_gap_crossing: true }, round] })
eq(rs.length, 2, 'session rounds: invalid elements dropped')
eq([rs[0].label, rs[0].merges, rs[0].noise], ['τ 0.5 · ref 8 · any', 0, null], 'session rounds: missing display fields filled')
eq(rs[1], round, 'session rounds: complete round kept')
eq(NOISE_DEFINITIONS.length, 4, 'four definitions')
eq(NOISE_DEFINITIONS[1].startsWith('Zig-zag: an internal segment whose state differs from both neighbours'), true, 'zig-zag definition text')

const st = { settings: { rules: [r0], allow_gap_crossing: false, rounds: [{ rules: [r0], allow_gap_crossing: false }] }, provenance: {} } as unknown as RecallStatus
eq(recallHeader(st)[0], `base rounds: ${JSON.stringify([{ rules: [r0], allow_gap_crossing: false }])}`, 'recall header lists base rounds first')
eq(recallHeader({ ...st, settings: { rules: [r0], allow_gap_crossing: false } } as RecallStatus)[0].startsWith('rules: '), true, 'no base rounds: header unchanged')

if (failures) { console.error(`${failures} failure(s)`); process.exit(1) } else console.log('all merge_rounds helper checks passed')
