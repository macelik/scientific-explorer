import { useRef, useState } from 'react'
import DeferredPlot from './DeferredPlot'
import { useIntegration } from '../integrationStore'
import { useStore } from '../store'
import { runHatchMerge, hatchShapes, type MergeResult, type MergeSettings } from '../hatchMerge'

const METRIC_NAMES: Record<string, string> = { srd: 'SRD', srd_phi: 'SRD/√φ', delta_log2fc: 'delta log2FC' }
const FIRED_NAMES: Record<string, string> = { transition: 'transition zone', proposal: 'merge proposal', ambiguous: 'ambiguous' }
function download(text: string, name: string, type: string) {
  const url = URL.createObjectURL(new Blob([text], { type })); const a = document.createElement('a')
  a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000)
}
function describeLayer(name: string, settings: MergeSettings): string | null {
  if (name === 'transition' && settings.transition) return `transition zone (${METRIC_NAMES[settings.transition.metric]}${settings.transition.metric === 'delta_log2fc' ? `, ${settings.transition.estimator}` : ''})`
  if (name === 'proposal' && settings.proposal) return `merge proposal (${METRIC_NAMES[settings.proposal.metric]}${settings.proposal.metric === 'delta_log2fc' ? `, ${settings.proposal.estimator}` : ''}, threshold ${settings.proposal.threshold})`
  if (name === 'ambiguous' && settings.ambiguous) return `ambiguous (${settings.ambiguous.estimator})`
  return null
}

export default function HatchMergeView() {
  const ig = useIntegration(), meta = ig.meta!
  const ex = useStore((s) => s.meta)!, presentation = useStore((s) => s.presentation)
  const [source, setSource] = useState('cluster4')
  const [smallMaxBins, setSmallMaxBins] = useState<number | ''>('')
  const [allowGapCrossing, setAllowGapCrossing] = useState(false)
  const [result, setResult] = useState<MergeResult | null>(null)
  const [stepIdx, setStepIdx] = useState(0)
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState('')
  const token = useRef(0)
  const chrom = ig.seg.chrom
  const hatchLayers = ig.hatchLayers
  // Eligibility reuses the SAME hatch-layer toggles/metrics/thresholds already
  // configured in the Chromosome tracks controls above, so there is no separate
  // copy of this settings here to go stale relative to those controls.
  const transition = hatchLayers.transition.on ? { metric: hatchLayers.transition.metric, estimator: hatchLayers.transition.estimator } : undefined
  const proposal = hatchLayers.proposal.on ? { metric: hatchLayers.proposal.metric, estimator: hatchLayers.proposal.estimator, threshold: hatchLayers.proposal.threshold } : undefined
  const ambiguous = hatchLayers.ambiguous.on ? { estimator: hatchLayers.ambiguous.estimator } : undefined
  const anyLayerOn = !!(transition || proposal || ambiguous)
  const request = { chrom, source, transition, proposal, ambiguous,
                    small_max_bins: smallMaxBins === '' ? undefined : smallMaxBins, allow_gap_crossing: allowGapCrossing }
  const signature = JSON.stringify(request)
  const [frozen, setFrozen] = useState('')
  const stale = !!result && frozen !== signature

  const run = async () => {
    if (!anyLayerOn) return
    const id = ++token.current; setBusy(true); setError('')
    try {
      const data = await runHatchMerge(request)
      if (token.current !== id) return
      setResult(data); setFrozen(signature); setStepIdx(data.steps.length - 1)
    } catch (e: any) { if (token.current === id) setError(String(e.message)) }
    finally { if (token.current === id) setBusy(false) }
  }
  const reset = () => setStepIdx(0)
  const undo = () => setStepIdx((i) => Math.max(0, i - 1))
  const step = result?.steps[stepIdx]
  // Read the *displayed* run's settings from result.settings, not the live hatch-layer
  // state above: once a layer's config changes after a run, the live state no longer
  // describes the scores actually shown (see the `stale` notice) and must not be used
  // for the summary/history/export of that still-displayed result.
  const rSettings = result?.settings
  const activeLayerNames = rSettings ? (['transition', 'proposal', 'ambiguous'] as const).filter((n) => rSettings[n]) : []

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
    <h3>Exploratory adjacent-segment merge</h3>
    <p className="lead">A new, independent exploratory rule: a boundary is eligible whenever ANY currently-enabled hatch layer flags it (transition zone, merge proposal, or ambiguous — tick them in the Chromosome tracks controls above; each uses its own metric/threshold/estimator there). Repeatedly merges one eligible boundary — the smallest merge-proposal |score| if any eligible boundary has one, else the leftmost — recomputing every score after each merge, until none remain eligible. This is not the original bootstrap/chain-guard merge prototype (see its results and controls further below) and it does not recompute CN calls on the merged boundaries.</p>
    <div className="controls-row">
      <label>Source<select aria-label="Merge source" value={source} onChange={(e) => setSource(e.target.value)}>{meta.sources.map((s) => <option key={s}>{s}</option>)}</select></label>
      <label>Chromosome<select aria-label="Merge chromosome" value={chrom} onChange={(e) => ig.setSeg({ chrom: e.target.value })}>{ex.chromosomes.map((c) => <option key={c.name}>{c.name}</option>)}</select></label>
      <label>Restrict eligibility to boundaries touching a segment ≤ <input aria-label="Small-segment restriction" type="number" min={1} placeholder="off" value={smallMaxBins} onChange={(e) => setSmallMaxBins(e.target.value === '' ? '' : +e.target.value)} style={{ width: 64 }} /> bins</label>
      <label className="chk-inline"><input type="checkbox" checked={allowGapCrossing} onChange={(e) => setAllowGapCrossing(e.target.checked)} />allow merging across a genomic gap (off by default)</label>
      <button className="btn primary" disabled={busy || !anyLayerOn} onClick={run}>{busy ? 'Merging…' : 'Run merge'}</button>
    </div>
    {!anyLayerOn && <p className="muted small">Tick at least one hatch layer (transition zone / merge proposal / ambiguous) in the Chromosome tracks controls above to enable a merge run.</p>}
    {anyLayerOn && <p className="muted small">Will use: {(['transition', 'proposal', 'ambiguous'] as const).filter((n) => request[n]).map((n) => FIRED_NAMES[n]).join(', ')}.</p>}
    {error && <p role="alert">{error}</p>}
    {result && step && <>
      <p className={stale ? 'notice' : 'muted'} role="status">{stale ? 'Settings changed; showing the previous run. ' : ''}{result.source} {result.chrom}: {result.original.length} original segments → {result.final.length} after {result.steps.length - 1} merge(s). Active layers: {activeLayerNames.map((n) => describeLayer(n, rSettings!)).filter(Boolean).join('; ')}.</p>
      <div className="controls-row">
        <label>Step<input aria-label="Merge step" type="range" min={0} max={result.steps.length - 1} value={stepIdx} onChange={(e) => setStepIdx(+e.target.value)} /></label>
        <span className="muted small">step {stepIdx} of {result.steps.length - 1}{step.removed_boundary !== null ? ` — removed boundary at bin ${step.removed_boundary}, fired ${step.fired.join('+')}${step.score !== null ? `, score ${step.score.toFixed(3)}` : ''}` : ' — original segmentation'}</span>
        <button className="btn-sm" onClick={reset} disabled={stepIdx === 0}>Reset to original</button>
        <button className="btn-sm" onClick={undo} disabled={stepIdx === 0}>Undo one step</button>
        <button className="btn-sm" onClick={() => setStepIdx(result.steps.length - 1)} disabled={stepIdx === result.steps.length - 1}>Jump to final</button>
        <button className="btn-sm" onClick={() => download(JSON.stringify(result, null, 2), `hatch-merge-${result.chrom}-${result.source}.json`, 'application/json')}>Download JSON</button>
        <button className="btn-sm" onClick={() => {
          const rows = [['step', 'source', 'chromosome', 'merged_start', 'merged_end', 'removed_boundary', 'fired', 'score', 'phi'],
            ...result.steps.filter((s) => s.removed_boundary !== null).map((s) => [s.step, result.source, result.chrom, s.merged_interval?.[0], s.merged_interval?.[1], s.removed_boundary, s.fired.join('+'), s.score ?? '', s.phi ?? ''])]
          download(rows.map((r) => r.join(',')).join('\n') + '\n', `hatch-merge-history-${result.chrom}-${result.source}.csv`, 'text/csv')
        }}>Download merge-history CSV</button>
      </div>
      <DeferredPlot data={[segmentTrace(result.original, result.start_bp, result.end_bp, '#94a3b8', 'Original')]}
        layout={{ ...layout, title: { text: 'Original segmentation', font: { size: 12 } }, xaxis: { title: `${result.chrom} position (Mb)` }, yaxis: { visible: false }, shapes: shapesForSegments(result.original, result.start_bp, result.end_bp) }} />
      <DeferredPlot data={[segmentTrace(step.segments, result.start_bp, result.end_bp, '#96491d', `Step ${stepIdx}`)]}
        layout={{ ...layout, title: { text: `Segmentation at step ${stepIdx} (${step.segments.length} segments)`, font: { size: 12 } }, xaxis: { title: `${result.chrom} position (Mb)` }, yaxis: { visible: false }, shapes: [...shapesForSegments(step.segments, result.start_bp, result.end_bp), ...(step.merged_interval ? hatchShapes(result.start_bp[step.merged_interval[0]] / 1e6, result.end_bp[step.merged_interval[1] - 1] / 1e6, 'y', '/', '#111827') : [])] }} />
      <div className="table-scroll"><table className="tbl"><thead><tr><th>step</th><th>source</th><th>chromosome</th><th>merged interval (bins)</th><th>removed boundary</th><th>fired</th><th>score</th><th>φ</th></tr></thead><tbody>
        {result.steps.filter((s) => s.removed_boundary !== null).map((s) => <tr key={s.step} className={s.step === stepIdx ? 'sel' : 'clickable'} onClick={() => setStepIdx(s.step)}>
          <td>{s.step}</td><td>{result.source}</td><td>{result.chrom}</td><td>[{s.merged_interval?.[0]}, {s.merged_interval?.[1]})</td><td>{s.removed_boundary}</td>
          <td>{s.fired.map((f) => FIRED_NAMES[f]).join(', ')}</td><td>{s.score !== null ? s.score.toFixed(3) : '–'}</td><td>{s.phi != null ? s.phi.toFixed(3) : '–'}</td>
        </tr>)}
      </tbody></table></div>
      <p className="small muted">Segment summaries and hatch classifications for this step's segmentation use the "new exploratory hatch layers" controls above once you tick a layer — those recompute from live boundaries and are not tied to this panel's step selection. This panel itself shows only the merged segmentation's boundaries (original vs. current-step spans in the two plots above and the history table); it does not display recomputed per-segment summaries or scores for the merged partition — that combination is not yet implemented here. CN calls are not recomputed on merged boundaries; any CN shown elsewhere in this tab still reflects the original accepted-breakpoint segmentation.</p>
      <details><summary>Merge run provenance and settings</summary><pre style={{ whiteSpace: 'pre-wrap', overflowWrap: 'anywhere' }}>{JSON.stringify({ settings: result.settings, provenance: result.provenance }, null, 2)}</pre></details>
    </>}
  </section>
}
