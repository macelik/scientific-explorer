import { useIntegration } from '../integrationStore'
import { changeSummaryTsv, downloadText, mergedBreakpointsTsv, recallProgressText, recalledClusterSegmentsTsv, recalledPercellTsv, transitionText } from '../cnRecall'
import { karyoModes, sameRounds } from '../mergeRounds'

const f4 = (v: number | null | undefined) => (v === null || v === undefined ? '—' : v.toFixed(4))

/** Before · After · Changed switch for both karyograms, plus the re-call impact summary and exports. */
export default function RecallImpact() {
  const ig = useIntegration(); const st = ig.recall
  const hasAfter = !!st?.cluster_segments
  const d = st?.diff
  const modes = karyoModes(ig.rounds.length, hasAfter)
  const onBase = !!st && sameRounds(st.settings.rounds ?? [], ig.rounds)
  return <div className="recall-impact" role="region" aria-label="CN re-call impact">
    <div className="controls-row">
      <span className="muted small">CN calls</span>
      <div className="segmented" role="group" aria-label="Karyogram mode">
        {modes.map((m) => <button key={m.mode} className={`btn-sm${ig.karyoMode === m.mode ? ' on' : ''}`} aria-pressed={ig.karyoMode === m.mode} disabled={!m.enabled} onClick={() => ig.set({ karyoMode: m.mode })}>{m.label}</button>)}
      </div>
      {ig.rounds.length > 0 && <span className="muted small">{`Before = round ${ig.rounds.length}; Original = the on-disk calls.`}</span>}
      {st ? <span className="muted small">{recallProgressText(st)}</span> : <span className="muted small">Run "Merge genome &amp; re-call CN" in section 2 (Merge rules) to compare.</span>}
      {ig.karyoMode === 'changed' && <span className="seg-legend"><span><i style={{ background: '#f1f1f1' }} />unchanged</span><span>other colours: the after state of changed bins</span></span>}
    </div>
    {d && <>
      <div className="table-scroll"><table className="tbl small">
        <thead><tr><th>source</th><th>segments before → after</th><th>merges (phase 2)</th><th>best_s before → after</th><th>bins changed</th></tr></thead>
        <tbody>{d.sources.map((s) => <tr key={s.source}><td>{s.source}</td><td>{s.segments_before} → {s.segments_after}</td><td>{s.merges} ({s.phase2})</td><td>{f4(s.best_s_before)} → {f4(s.best_s_after)}</td><td>{s.bins_changed.toLocaleString('en-US')} ({s.pct_changed.toFixed(2)}%)</td></tr>)}</tbody>
      </table></div>
      {d.cells && <p className="small">cells with any change: {d.cells.n_changed_cells.toLocaleString('en-US')} of {d.cells.n_cells.toLocaleString('en-US')} · median {d.cells.median_pct_changed.toFixed(2)}% bins changed · best_s change median {f4(d.cells.best_s_change.median)} [{f4(d.cells.best_s_change.min)}, {f4(d.cells.best_s_change.max)}]</p>}
      <p className="small">Top transitions, clusters: {d.transitions.clusters.slice(0, 5).map(transitionText).join(' · ') || 'none'}{d.transitions.cells ? <> · cells: {d.transitions.cells.slice(0, 5).map(transitionText).join(' · ') || 'none'}</> : null}</p>
      <p className="muted small">A merge on one chromosome can change calls elsewhere: best_s is fitted over the whole genome.</p>
      <div className="controls-row">
        <button className="btn-sm" onClick={() => downloadText(mergedBreakpointsTsv(st!), 'merged_breakpoints.tsv')}>merged_breakpoints.tsv</button>
        <button className="btn-sm" onClick={() => downloadText(recalledClusterSegmentsTsv(st!), 'recalled_cluster_segments.tsv')}>recalled_cluster_segments.tsv</button>
        <button className="btn-sm" disabled={!d.cells} onClick={() => downloadText(recalledPercellTsv(st!), 'recalled_percell_best_s.tsv')}>recalled_percell_best_s.tsv</button>
        <button className="btn-sm" onClick={() => downloadText(changeSummaryTsv(st!), 'cn_change_summary.tsv')}>cn_change_summary.tsv</button>
        {st?.status === 'done' && (onBase
          ? <button className="btn-sm primary" onClick={() => ig.adoptRecall()}>{`Adopt as round ${ig.rounds.length + 1}`}</button>
          : <span className="muted small">computed on a different base</span>)}
      </div>
    </>}
  </div>
}
