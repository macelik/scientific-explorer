import { useRef, useState } from 'react'
import DeferredPlot from './DeferredPlot'
import MergeRulesEditor from './MergeRulesEditor'
import { useIntegration } from '../integrationStore'
import { useStore } from '../store'
import { runHatchMerge, hatchShapes, type MergeResult } from '../hatchMerge'
import { canRunMerge, describeCondition } from '../mergeRules'

function download(text: string, name: string, type: string) {
  const url = URL.createObjectURL(new Blob([text], { type })); const a = document.createElement('a')
  a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000)
}

export default function HatchMergeView() {
  const ig = useIntegration(), meta = ig.meta!
  const ex = useStore((s) => s.meta)!, presentation = useStore((s) => s.presentation)
  const [source, setSource] = useState('cluster4')
  const [result, setResult] = useState<MergeResult | null>(null)
  const [stepIdx, setStepIdx] = useState(0)
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState('')
  const token = useRef(0)
  const chrom = ig.seg.chrom
  const runnable = canRunMerge(ig.mergeRules)
  const request = { chrom, source, rules: ig.mergeRules, allow_gap_crossing: ig.allowGapCrossing }
  const signature = JSON.stringify(request)
  const [frozen, setFrozen] = useState('')
  const stale = !!result && frozen !== signature

  const run = async () => {
    if (!runnable) return
    const id = ++token.current; setBusy(true); setError('')
    try {
      const data = await runHatchMerge(request)
      if (token.current !== id) return
      setResult(data); setFrozen(signature); setStepIdx(data.steps.length - 1)
    } catch (e: any) { if (token.current === id) setError(String(e.message)) }
    finally { if (token.current === id) setBusy(false) }
  }
  const step = result?.steps[stepIdx]
  // Describe the displayed run from result.settings, never from the live editor state.
  const ruleSummary = result ? result.settings.rules.filter((r) => r.enabled).map((r, i) => `${i + 1}. ${r.name}: ${r.conditions.map(describeCondition).join(' AND ')}`).join(' | ') : ''

  const layout = { height: 260, margin: { l: 60, r: 20, t: 30, b: 45 }, paper_bgcolor: '#fff', plot_bgcolor: '#fff', font: { family: 'system-ui', size: 12 }, hovermode: presentation ? false : 'closest' as const }
  const segmentTrace = (segs: { start: number; end: number }[], startBp: number[], endBp: number[], color: string, name: string) => ({
    type: 'bar', orientation: 'h', name, y: segs.map(() => name), base: segs.map((s) => startBp[s.start] / 1e6),
    x: segs.map((s) => (endBp[s.end - 1] - startBp[s.start]) / 1e6), marker: { color }, hovertemplate: '%{x:.3f} Mb<extra></extra>',
    customdata: segs.map((s) => `${s.start},${s.end}`),
  })
  const shapesForSegments = (segs: { start: number; end: number }[], startBp: number[], endBp: number[]) =>
    segs.map((s, i) => ({ type: 'rect', xref: 'x', yref: 'paper', x0: startBp[s.start] / 1e6, x1: endBp[s.end - 1] / 1e6, y0: 0, y1: 1,
      fillcolor: i % 2 === 0 ? 'rgba(150,73,29,0.18)' : 'rgba(40,120,160,0.18)', line: { width: 1, color: '#334155' } }))

  return <section className="panel hatch-merge">
    <h3>2. Merge rules — exploratory adjacent-segment merge</h3>
    <p className="lead">Build an ordered list of rules. A rule is a set of conditions on a segment's left/right flank scores that must ALL hold; a segment is flagged when ANY enabled rule matches, and each rule's hatch on the Chromosome tracks above shows exactly the segments it flags. A merge run repeatedly joins one flagged segment to its more similar neighbour (the flank with the smaller |score| under that rule's "merge into" metric): topmost rule first, then smallest score, then leftmost, recomputing every score, including the pooled dispersion φ, after each merge. These are user-defined exploratory diagnostics, not validated decision rules. This is not the original bootstrap/chain-guard merge prototype, and CN calls are not recomputed on merged boundaries.</p>
    <MergeRulesEditor />
    <div className="controls-row">
      <label>Source<select aria-label="Merge source" value={source} onChange={(e) => setSource(e.target.value)}>{meta.sources.map((s) => <option key={s}>{s}</option>)}</select></label>
      <label>Chromosome<select aria-label="Merge chromosome" value={chrom} onChange={(e) => ig.setSeg({ chrom: e.target.value })}>{ex.chromosomes.map((c) => <option key={c.name}>{c.name}</option>)}</select></label>
      <label className="chk-inline"><input type="checkbox" aria-label="Allow merging across a genomic gap" checked={ig.allowGapCrossing} onChange={(e) => ig.setAllowGapCrossing(e.target.checked)} />allow merging across a genomic gap (off by default)</label>
      <button className="btn primary" disabled={busy || !runnable} onClick={run}>{busy ? 'Merging…' : 'Run merge'}</button>
    </div>
    {!runnable && <p className="muted small">Enable at least one rule, and complete every enabled rule (each needs at least one condition with valid values), to run a merge.</p>}
    {error && <p role="alert">{error}</p>}
    {result && step && <>
      <p className={stale ? 'notice' : 'muted'} role="status">{stale ? 'Settings changed; showing the previous run. ' : ''}{result.source} {result.chrom}: {result.original.length} original segments → {result.final.length} after {result.steps.length - 1} merge(s). Rules: {ruleSummary}{result.settings.allow_gap_crossing ? ' · gap crossing allowed' : ''}.</p>
      <div className="controls-row">
        <label>Step<input aria-label="Merge step" type="range" min={0} max={result.steps.length - 1} value={stepIdx} onChange={(e) => setStepIdx(+e.target.value)} /></label>
        <span className="muted small">step {stepIdx} of {result.steps.length - 1}{step.removed_boundary !== null ? ` — removed boundary at bin ${step.removed_boundary}: rule "${step.rule_name}" merged segment [${step.segment?.[0]}, ${step.segment?.[1]}) ${step.side}${step.score !== null ? `, |score| ${step.score.toFixed(3)}` : ''}` : ' — original segmentation'}</span>
        <button className="btn-sm" onClick={() => setStepIdx(0)} disabled={stepIdx === 0}>Reset to original</button>
        <button className="btn-sm" onClick={() => setStepIdx((i) => Math.max(0, i - 1))} disabled={stepIdx === 0}>Undo one step</button>
        <button className="btn-sm" onClick={() => setStepIdx(result.steps.length - 1)} disabled={stepIdx === result.steps.length - 1}>Jump to final</button>
        <button className="btn-sm" onClick={() => download(JSON.stringify(result, null, 2), `hatch-merge-${result.chrom}-${result.source}.json`, 'application/json')}>Download JSON</button>
        <button className="btn-sm" onClick={() => {
          const rows = [['step', 'source', 'chromosome', 'merged_start', 'merged_end', 'removed_boundary', 'segment_start', 'segment_end', 'rule_id', 'rule_name', 'side', 'score', 'phi'],
            ...result.steps.filter((s) => s.removed_boundary !== null).map((s) => [s.step, result.source, result.chrom, s.merged_interval?.[0], s.merged_interval?.[1], s.removed_boundary, s.segment?.[0], s.segment?.[1], s.rule_id, JSON.stringify(s.rule_name ?? ''), s.side, s.score ?? '', s.phi ?? ''])]
          download(rows.map((r) => r.join(',')).join('\n') + '\n', `hatch-merge-history-${result.chrom}-${result.source}.csv`, 'text/csv')
        }}>Download merge-history CSV</button>
      </div>
      <DeferredPlot data={[segmentTrace(result.original, result.start_bp, result.end_bp, '#94a3b8', 'Original')]}
        layout={{ ...layout, title: { text: 'Original segmentation', font: { size: 12 } }, xaxis: { title: `${result.chrom} position (Mb)` }, yaxis: { visible: false }, shapes: shapesForSegments(result.original, result.start_bp, result.end_bp) }} />
      <DeferredPlot data={[segmentTrace(step.segments, result.start_bp, result.end_bp, '#96491d', `Step ${stepIdx}`)]}
        layout={{ ...layout, title: { text: `Segmentation at step ${stepIdx} (${step.segments.length} segments)`, font: { size: 12 } }, xaxis: { title: `${result.chrom} position (Mb)` }, yaxis: { visible: false }, shapes: [...shapesForSegments(step.segments, result.start_bp, result.end_bp), ...(step.merged_interval ? hatchShapes(result.start_bp[step.merged_interval[0]] / 1e6, result.end_bp[step.merged_interval[1] - 1] / 1e6, 'y', '/', '#111827') : [])] }} />
      <div className="table-scroll"><table className="tbl"><thead><tr><th>step</th><th>source</th><th>chromosome</th><th>merged interval (bins)</th><th>removed boundary</th><th>segment</th><th>rule</th><th>side</th><th>score</th><th>φ</th></tr></thead><tbody>
        {result.steps.filter((s) => s.removed_boundary !== null).map((s) => <tr key={s.step} className={s.step === stepIdx ? 'sel' : 'clickable'} onClick={() => setStepIdx(s.step)}>
          <td>{s.step}</td><td>{result.source}</td><td>{result.chrom}</td><td>[{s.merged_interval?.[0]}, {s.merged_interval?.[1]})</td><td>{s.removed_boundary}</td>
          <td>[{s.segment?.[0]}, {s.segment?.[1]})</td><td>{s.rule_name}</td><td>{s.side}</td><td>{s.score !== null ? s.score.toFixed(3) : '–'}</td><td>{s.phi != null ? s.phi.toFixed(3) : '–'}</td>
        </tr>)}
      </tbody></table></div>
      <p className="small muted">The rule hatches on the Chromosome tracks are computed from the original accepted-breakpoint segmentation, not from this panel's selected step. This panel shows the merged segmentation's boundaries (the two plots and the history table); it does not display recomputed per-segment scores for the merged partition. CN calls are not recomputed on merged boundaries.</p>
      <details><summary>Merge run provenance and settings</summary><pre style={{ whiteSpace: 'pre-wrap', overflowWrap: 'anywhere' }}>{JSON.stringify({ settings: result.settings, provenance: result.provenance }, null, 2)}</pre></details>
    </>}
  </section>
}
