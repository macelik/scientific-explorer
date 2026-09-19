import { CN_COLORS, CN_LABELS } from '../colors'
import type { AssignResult, SegmentRow } from '../types'

function Strip({ cn5, label, n }: { cn5: number[]; label: string; n: number }) {
  // contiguous runs of the same categorical state (never averaged)
  const runs: { s: number; e: number; v: number }[] = []
  for (let i = 0; i < cn5.length; i++) { const last = runs[runs.length - 1]; if (last && last.v === cn5[i]) last.e = i + 1; else runs.push({ s: i, e: i + 1, v: cn5[i] }) }
  return (
    <div className="strip-row"><span className="strip-label">{label}</span>
      <div className="strip">{runs.map((r, i) => <span key={i} style={{ left: `${(100 * r.s) / n}%`, width: `${(100 * (r.e - r.s)) / n}%`, background: CN_COLORS[r.v] }} title={`bins [${r.s},${r.e}) · CN ${r.v} = ${CN_LABELS[r.v]}`} />)}</div>
    </div>
  )
}

export default function CNResult({ res, outdated }: { res: AssignResult; outdated: boolean }) {
  const n = res.production.cn5.length
  const prod = res.production, man = res.manual
  const method = res.segment_estimator ?? 'arithmetic'
  const trimmed = method !== 'arithmetic'
  const arithmetic = res.arithmetic_manual ?? man
  const methodLabel = method === 'iqr_upper' ? 'IQR upper tail' : method === 'iqr_two_sided' ? 'IQR two-sided' : 'arithmetic mean'
  const estimatorChanged = (s: SegmentRow) => arithmetic.cn5.slice(s.start_bin,s.end_bin).some((v,i)=>v!==man.cn5[s.start_bin+i])
  const changedSeg = (s: SegmentRow) => { for (let i = s.start_bin; i < s.end_bin; i++) if (prod.cn5[i] !== man.cn5[i]) return true; return false }
  const exportCsv = () => {
    const rows: (string | number | null | undefined)[][] = [
      [`# method=${method}; scope=all full-genome segments; original genome baseline and bin-count weights; five-state 1=base`],
      ['segment_id','start_bin','end_bin_excl','n_bins','start_bp','end_bp','arithmetic_mean_x','selected_summary_x','selected_summary_norm','retained_bins','arithmetic_cn5','selected_cn5','cn_label','cn_int','cn_cont','changed_vs_arithmetic','changed_vs_production'],
      ...man.segments.map(s=>[s.segment_id,s.start_bin,s.end_bin,s.n_bins,s.start_bp,s.end_bp,s.mean_x,s.summary_x??s.mean_x,s.summary_norm??s.mean_norm,s.retained_bins??s.n_bins,arithmetic.cn5[s.start_bin],s.cn5,CN_LABELS[s.cn5],s.cn_int,s.cn_cont,estimatorChanged(s)?'yes':'no',changedSeg(s)?'yes':'no'])]
    const blob = new Blob([rows.map(r=>r.join(',')).join('\n')], {type:'text/csv'})
    const a = document.createElement('a'), url=URL.createObjectURL(blob)
    a.href=url; a.download=`${res.cell}_${res.chrom}_${method}_manual_cn.csv`; a.click()
    setTimeout(()=>URL.revokeObjectURL(url),1000)
  }
  return (
    <div className={`card cn ${outdated ? 'outdated' : ''}`}>
      <div className="row wrap">
        <strong>CN assignment · {res.cell} / {res.chrom}</strong>
        {outdated && <span className="badge fail">OUTDATED — the tree or segment summary changed after this result; press “Finish & assign CN” again</span>}
        <span className="muted">full-genome fit: {res.chrom} leaves spliced into production segmentation ({man.n_segments_genome} genome segments), {methodLabel} · original full-genome baseline {res.genome_baseline_trimmed.toFixed(4)} · {res.compute_ms} ms</span>
        <span className="spacer" />
        <button className="btn-sm" onClick={exportCsv}>CSV</button>
      </div>
      <div className="row wrap">
        <span>production scale factor <b>{prod.best_s.toFixed(6)}</b> {prod.matches_export ? <span className="badge ok">reproduces exported CN row and best_s {prod.best_s_export.toFixed(6)}</span> : <span className="badge fail">differs from export ({prod.best_s_export.toFixed(6)})</span>}</span>
        <span>selected scale factor <b>{man.best_s.toFixed(6)}</b> {man.best_s !== prod.best_s ? <span className="badge warn">shared scale changed</span> : <span className="muted">unchanged</span>}</span>
        <span>changed calls on {res.chrom}: <b>{res.changed_bins_this_chrom}</b> / {n} bins</span>
        <span>elsewhere: {res.changed_elsewhere.length ? res.changed_elsewhere.map((c) => `${c.chrom} ${c.changed_bins}/${c.n}`).join(', ') : 'no other chromosome changed'} <span className="muted small">(other chromosomes keep their boundaries; {trimmed ? 'their summaries and the shared scale can change' : 'only the shared scale can move their calls'})</span></span>
      </div>
      {trimmed && <div className="notice estimator-comparison">Same manual boundaries: arithmetic scale <b>{arithmetic.best_s.toFixed(6)}</b> → {methodLabel} scale <b>{man.best_s.toFixed(6)}</b>. Changed five-state calls from the estimator: <b>{res.estimator_changed_bins ?? 0}</b> / {n} bins on {res.chrom}; {res.estimator_changed_elsewhere?.reduce((sum,c)=>sum+c.changed_bins,0) ?? 0} bins elsewhere. Original segment lengths weight both fits.</div>}
      <Strip cn5={prod.cn5} label="production" n={n} />
      {trimmed && <Strip cn5={arithmetic.cn5} label="manual arithmetic" n={n} />}
      <Strip cn5={man.cn5} label={`manual ${methodLabel}`} n={n} />
      <div className="karyo-legend">{[0, 0.5, 1, 1.5, 2].map((v) => <span key={v}><i style={{ background: CN_COLORS[v] }} /> {v} = {CN_LABELS[v]}</span>)}</div>
      <div style={{overflowX:'auto'}}><table className="tbl stats">
        <thead><tr><th>segment</th><th>bins [start, end)</th><th>Mb</th><th>arithmetic mean X</th><th>selected estimate X</th><th>retained / total</th><th>estimate / baseline</th>{trimmed && <th>arithmetic five-state</th>}<th>selected five-state</th><th>integer CN</th><th>continuous</th><th>vs production</th></tr></thead>
        <tbody>{man.segments.map(s=>{const ch=changedSeg(s);return <tr key={s.segment_id} className={ch?'changed':''}>
          <td>{s.segment_id}</td><td>[{s.start_bin}, {s.end_bin})</td><td>{(s.start_bp/1e6).toFixed(2)}–{(s.end_bp/1e6).toFixed(2)}</td>
          <td>{s.mean_x.toFixed(3)}</td><td>{(s.summary_x??s.mean_x).toFixed(3)}</td><td>{s.retained_bins??s.n_bins} / {s.n_bins}</td><td>{(s.summary_norm??s.mean_norm)?.toFixed(3)}</td>
          {trimmed && <td>{arithmetic.cn5[s.start_bin]} {estimatorChanged(s)&&<span className="badge warn">changed</span>}</td>}
          <td><i className="sw" style={{background:CN_COLORS[s.cn5]}} />{s.cn5} = {CN_LABELS[s.cn5]}</td><td>{s.cn_int}</td><td>{s.cn_cont.toFixed(2)}</td><td>{ch?'changed':'same'}</td>
        </tr>})}</tbody>
      </table></div>
      <details><summary className="muted small">production segments on {res.chrom} ({prod.segments.length})</summary>
        <table className="tbl stats"><thead><tr><th>segment</th><th>bins</th><th>n</th><th>mean X</th><th>five-state</th><th>int</th><th>cont</th></tr></thead>
          <tbody>{prod.segments.map((s) => <tr key={s.segment_id}><td>{s.segment_id}</td><td>[{s.start_bin}, {s.end_bin})</td><td>{s.n_bins}</td><td>{s.mean_x.toFixed(3)}</td><td>{s.cn5}</td><td>{s.cn_int}</td><td>{s.cn_cont.toFixed(2)}</td></tr>)}</tbody></table></details>
    </div>
  )
}
