import { useMemo } from 'react'
import { nodeKey, useStore } from '../store'
import { selectionCoords } from '../coords'
import { computeWindow } from './WindowWorkspace'
import { fmt } from '../stats'

export default function Inspector() {
  const meta = useStore((s) => s.meta)!
  const timings = useStore((s) => s.timings)
  const activeId = useStore((s) => s.activeSelectionId)
  const sel = useStore((s) => s.selections.find((x) => x.id === s.activeSelectionId))
  const node = useStore((s) => (sel ? s.nodes[nodeKey(sel.cell, sel.chrom, sel.nodeStart, sel.nodeEnd)] : undefined))
  const events = useStore((s) => s.events)
  const w = useMemo(() => (sel && node ? computeWindow(sel, node) : null), [sel, node])
  const c = w ? selectionCoords(node!, sel!.s, sel!.e) : null
  return (
    <aside className="inspector">
      <h3>Inspector</h3>
      {w && c ? (
        <div className="insp-block">
          <div className="tbl-title" style={{ color: sel!.color }}>{sel!.name}</div>
          <table className="kv"><tbody>
            <tr><td>node</td><td>{sel!.nodeSide} [{sel!.nodeStart}, {sel!.nodeEnd}) v{sel!.nodeVersion}</td></tr>
            <tr><td>node-local</td><td>[{c.s}, {c.e}) · {c.nBins} bins</td></tr>
            <tr><td>chromosome-local</td><td>[{c.chromStart}, {c.chromEnd})</td></tr>
            <tr><td>genome-wide</td><td>[{c.genomeStart}, {c.genomeEnd})</td></tr>
            <tr><td>genomic (1-based incl.)</td><td>{c.startBp.toLocaleString()}–{c.endBp.toLocaleString()}</td></tr>
            <tr><td>span</td><td>{c.spanMb.toFixed(3)} Mb ({c.spanBp.toLocaleString()} bp)</td></tr>
            <tr><td>retained</td><td>{c.retainedMb.toFixed(3)} Mb · {c.gaps} coordinate gap(s)</td></tr>
            <tr><td>AD positions</td><td>{c.adCount}{c.adFirst !== null ? ` (j ${c.adFirst}…${c.adLast})` : ''}</td></tr>
            <tr><td>anchor</td><td>{w.anchorKind} b={w.anchor}</td></tr>
            <tr><td>AD</td><td>mean {fmt(w.adStats.mean)} · CV {fmt(w.adStats.cv)} · lag1 {fmt(w.adStats.lag1)}</td></tr>
            <tr><td>X</td><td>mean {fmt(w.xStats.mean)} · CV {fmt(w.xStats.cv)} · lag1 {fmt(w.xStats.lag1)}</td></tr>
          </tbody></table>
        </div>
      ) : <div className="muted small">Select a window (click its row or drag on a track) to see exact coordinates in all three systems.</div>}
      <details open><summary>Definitions</summary>
        <ul className="small">
          <li>AD curve: pyEpi <code>seq_dist_ad</code> on the node's own X; boundary <i>j</i> splits [0,j) | [j,n), value d[j−1], drawn at the start of bin j.</li>
          <li>Selection: node-local half-open [s,e); AD sample included when max(s,1) ≤ j &lt; e. A 40-bin window has 40 AD positions unless it touches node bin 0.</li>
          <li>Stats: mean, median, population std (ddof 0), CV = std/mean (N/A if mean 0), lag-k = Pearson(v[:-k], v[k:]) needing ≥3 pairs and non-constant arrays. Descriptive only.</li>
          <li>log2FC: log2(rolling mean of X / reference); smoothing before log2; zero windows are −∞ shown clipped; zero reference = undefined (no substitution).</li>
          <li>Baselines: figure-compatible <code>trimmed_mean_iqr(lb=False)</code> or literal zero-inclusive median, genome-wide and for the current node.</li>
          <li>Tests: <code>permutation_test_ad</code> + <code>global_permutation_test_ad</code> (full-genome X row), 1000 permutations, seed 42, strict p &lt; 0.001 on both.</li>
          <li>CN: leaves spliced into the cell's production segmentation, <code>assign_gainloss_new</code> on the full X row; 1 = base/diploid.</li>
        </ul>
      </details>
      <details><summary>Dataset</summary>
        <div className="small">
          <div>identity <code>{meta.dataset.identity.hash}</code></div>
          <div>{meta.dataset.identity.n_cells} cohort cells of {meta.dataset.identity.n_cells_total_h5ad} in h5ad · {meta.dataset.identity.n_bins} bins · {meta.dataset.bin_width_bp} bp</div>
          <div className="muted">{meta.dataset.h5ad}</div>
          <div className="muted">{meta.dataset.cohort_dir}</div>
          <div>signal: {meta.dataset.signal}</div>
          <div>event index: {events.status} {events.total ? `${events.done}/${events.total}` : ''}</div>
        </div>
      </details>
      <details open><summary>Measured timings (this session)</summary>
        <table className="kv small"><tbody>{timings.slice(0, 14).map((t, i) => <tr key={i}><td>{t.label}</td><td>{t.ms} ms</td></tr>)}</tbody></table>
      </details>
    </aside>
  )
}
