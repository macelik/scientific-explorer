import { useEffect, useRef, useState } from 'react'
import { useStore } from '../store'
import { experimentKey, sameDraft, useExperimentStore } from '../experimentStore'
import type { ExperimentDraft, ExperimentResult, ExperimentCN } from '../experimentStore'
import { mean, summarize, fmtP } from '../stats'
import { ChildExperimentPlots, ExperimentTracks, EDITED, ORIGINAL } from './ExperimentPlots'

const number = (v: number | null | undefined, digits = 3) => v == null || !Number.isFinite(v) ? 'N/A' : v.toLocaleString(undefined, { maximumFractionDigits: digits })
function download(text: string, name: string, type: string) {
  const url = URL.createObjectURL(new Blob([text], { type })), a = document.createElement('a')
  a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000)
}

function CNStrip({ result }: { result: ExperimentCN }) {
  if (!result.cn5) return <span className="muted">{result.reason}</span>
  const colors: Record<string, string> = { '0': '#336f9a', '0.5': '#a8c4dc', '1': '#f0eee8', '1.5': '#efb4aa', '2': '#bd4a40' }
  const runs: { s: number; e: number; value: number }[] = []
  result.cn5.forEach((value, i) => { const prev = runs[runs.length-1]; if (prev?.value === value) prev.e = i+1; else runs.push({ s: i, e: i+1, value }) })
  return <div className="experiment-cn-strip" aria-label="Five-state CN strip in retained-bin order">{runs.map(r => <span key={r.s} style={{ flex: r.e-r.s, background: colors[r.value] }} title={`bins [${r.s}, ${r.e}) · CN code ${r.value}`} />)}</div>
}

function ExperimentResults({ result, outdated, onRegion }: { result: ExperimentResult; outdated: boolean; onRegion?: (start:number,end:number)=>void }) {
  const { request: r } = result
  const exportSamples = () => {
    const rows = ['chrom,chrom_bin,start_bp,end_bp,X_original,X_edited,delta_X,AD_original_root,AD_edited_root']
    for (let i=result.intervention.event_start; i<result.intervention.event_end; i++) rows.push([r.chrom, i, result.start_bp[i], result.end_bp[i], result.x_original[i], result.x_edited[i], result.x_edited[i]-result.x_original[i],
      i ? result.original.nodes[0].ad[i-1] : '', i ? result.edited.nodes[0].ad[i-1] : ''].join(','))
    download(rows.join('\n'), `${r.cell}-${r.chrom}-experiment-samples.csv`, 'text/csv')
  }
  const eventStart=result.intervention.event_start, eventEnd=result.intervention.event_end
  const x0 = result.x_original.slice(eventStart, eventEnd), x1 = result.x_edited.slice(eventStart, eventEnd)
  const b0 = result.original.nodes[0].argmax_b, b1 = result.edited.nodes[0].argmax_b
  const bp = (b: number | null) => b === null ? 'N/A' : `${number(result.start_bp[b]/1e6)} Mb`
  return <section className="experiment-results" aria-label="Experiment results">
    <div className="section-heading"><div><h2>Original → experiment</h2><span className="muted">{r.cell} / {r.chrom} · event [{eventStart}, {eventEnd}) · {r.operation} {r.operation==='simulate' ? r.cn_state : ['extend','duplicate'].includes(r.operation) ? `${r.direction} ${result.intervention.actual_extension} bins` : r.value} · k={result.parameters.k} · {number(result.compute_ms/1000, 2)} s</span></div>
      <div className="row"><button className="btn-sm" onClick={() => download(JSON.stringify(result, null, 2), `${r.cell}-${r.chrom}-experiment.json`, 'application/json')}>Export experiment JSON</button><button className="btn-sm" onClick={exportSamples}>Export regional CSV</button></div></div>
    {outdated && <div className="notice" role="status">Settings changed. These are the previous run’s results; run the experiment again to update them.</div>}
    <div className="experiment-measures">
      <div><span>Root candidate</span><strong>{bp(b0)} → {bp(b1)}</strong></div>
      <div><span>Mean X in region</span><strong>{number(mean(x0))} → {number(mean(x1))}</strong></div>
      <div><span>Zero fraction in region</span><strong>{number(100*x0.filter(v => v===0).length/x0.length, 1)}% → {number(100*x1.filter(v => v===0).length/x1.length, 1)}%</strong></div>
      <div><span>Genome normalization factor</span><strong>× {number(result.intervention.normalization_factor, 6)}</strong></div>
    </div>
    <ExperimentTracks result={result} starts={result.start_bp} ends={result.end_bp} original={result.x_original} edited={result.x_edited} start={eventStart} end={eventEnd} onRegion={onRegion} />
    <p className="muted">Drag horizontally to select a new source region. Orange marks the edit target; solid lines mark accepted breakpoints and dashed lines mark unaccepted root candidates. Blue = original; brown = edited. All curves use unsmoothed X; coordinate gaps break the lines.</p>
    <details className="scientific-details"><summary>Both children and exact split tests</summary>
      <p>Each node is recomputed on its own interval. Bounds can differ between scenarios. This run uses k={result.parameters.k}: tested depths 0–{result.parameters.k-1}. Immediate hypothetical children may be inspected when the root fails; they do not enter automatic segmentation or expand further. Deeper branches are visited only after accepted parents.</p>
      <ChildExperimentPlots result={result} />
      <div className="table-scroll"><table className="tbl"><thead><tr><th>Scenario / node</th><th>Bins [start, end)</th><th>Candidate</th><th>AD</th><th>Local p</th><th>Global p</th><th>Automatic decision</th></tr></thead><tbody>
        {([['Original', result.original], ['Edited', result.edited]] as const).flatMap(([label, tree]) => tree.nodes.map(node => <tr key={label+node.side}>
          <td>{label} / {node.side} · depth {node.depth ?? (node.side==='root'?0:node.side.split('.').length)}</td><td>[{node.start}, {node.end})</td><td>{bp(node.argmax_b === null ? null : node.start+node.argmax_b)}</td><td>{number(node.argmax_ad)}</td>
          <td>{node.tests ? fmtP(node.tests.p_local, node.tests.p_local_exact) : 'N/A'}</td><td>{node.tests ? fmtP(node.tests.p_global, node.tests.p_global_exact) : 'N/A'}</td>
          <td>{node.argmax_b === null ? 'Cannot split' : node.hypothetical ? `Hypothetical · tests ${node.passes ? 'pass' : 'fail'}` : node.accepted ? 'Accepted' : 'Rejected'}</td>
        </tr>))}
      </tbody></table></div>
      {result.edited_at_original_boundary && <p className="same-boundary">Holding the original root boundary at <b>{bp(b0)}</b>: edited AD {number(result.edited_at_original_boundary.observed_ad)}, local p {fmtP(result.edited_at_original_boundary.p_local, result.edited_at_original_boundary.p_local_exact)}, global p {fmtP(result.edited_at_original_boundary.p_global, result.edited_at_original_boundary.p_global_exact)}.</p>}
      <p className="muted">1,000 permutations · seed 42 · both p &lt; 0.001 · each global test uses that scenario’s full-genome X.</p>
    </details>
    <section className="edge-diagnostics"><h3>Did pyEpi find the event edges?</h3><p className="muted">Exact retained-bin boundaries; chromosome ends require no split. Tests below evaluate the event edge even when the algorithm chose another candidate.</p>
      <div className="table-scroll"><table className="tbl"><thead><tr><th>Event edge</th><th>Accepted by recursion?</th><th>Tested node</th><th>AD rank in node</th><th>Local p</th><th>Global p</th></tr></thead><tbody>
        {result.event_edges.flatMap(edge => edge.nodes.length ? edge.nodes.map(node => <tr key={edge.edge+node.side}><td>{edge.edge} · bin {edge.boundary} · {bp(edge.boundary)}</td><td>{edge.accepted?'Yes':'No'}</td><td>{node.side}</td><td>{node.rank}{node.is_winner?' · argmax':''}</td><td>{fmtP(node.p_local,node.p_local_exact)}</td><td>{fmtP(node.p_global,node.p_global_exact)}</td></tr>) : [<tr key={edge.edge}><td>{edge.edge} · bin {edge.boundary}</td><td>{edge.chromosome_edge?'Chromosome edge':'Not accepted'}</td><td colSpan={4}>No within-node split to test</td></tr>])}
      </tbody></table></div>
    </section>
    <section className="cn-comparison"><h3>Separate the signal effect from the segmentation effect</h3><p className="muted">CN is fitted across the full genome in every row. Only the selected chromosome is resegmented; other chromosomes keep their production boundaries.</p>
      <div className="table-scroll"><table className="tbl"><thead><tr><th>Input and boundaries</th><th>Five-state CN · retained-bin order</th><th>Fitted scale</th><th>Changed bins here</th><th>Changed elsewhere</th></tr></thead><tbody>
        {([['original', 'Original X · production boundaries'], ['fixed', 'Edited X · production boundaries'], ['resegmented', 'Edited X · recomputed boundaries'], ['event_isolated', 'Edited X · event edges imposed (diagnostic)']] as const).map(([key, label]) => {
          const cn = result.cn[key]
          return <tr key={key}><td>{label}</td><td style={{ minWidth: 240 }}><CNStrip result={cn} /></td><td>{number(cn.best_s, 6)}</td><td>{cn.changed_bins ?? 'N/A'}</td><td>{cn.status === 'undefined' ? 'N/A' : cn.changed_elsewhere?.length ? cn.changed_elsewhere.map(v => `${v.chrom}: ${v.changed_bins}`).join(' · ') : '0'}</td></tr>
        })}
      </tbody></table></div>
      <p className="muted">CN code: 0 loss · 0.5 putative loss · 1 base/diploid · 1.5 putative gain · 2 gain. Changes are relative to the original row.</p>
      <details className="scientific-details"><summary>Holmes, Watson and combined calls inside the event</summary><p>The imposed-edge row adds the event edges to production boundaries, isolating CN calling from automatic breakpoint detection. It is not a detected event.</p><table className="tbl"><thead><tr><th>Boundary treatment</th><th>Holmes bins by state</th><th>Watson bins by state</th><th>Combined bins by state</th></tr></thead><tbody>{Object.entries(result.cn).map(([key, cn])=><tr key={key}><td>{key.replace(/_/g,' ')}</td>{['holmes','watson','cn5'].map(kind=><td key={kind}>{cn.event_calls ? Object.entries(cn.event_calls[kind]).map(([state,count])=>`${state}: ${count}`).join(' · ') : 'Undefined'}</td>)}</tr>)}</tbody></table><p className="muted">Holmes/Watson states: 0 loss, 1 base, 2 gain. These are state codes, not integer copy numbers.</p></details>
    </section>
    <details className="scientific-details"><summary>Transformation and provenance</summary><p>{result.intervention.clipped_bins} bins clipped to zero before normalization. Genome total: {number(result.intervention.original_total)} → {number(result.intervention.edited_total)}. {result.intervention.changed_bins_genome} changed bins genome-wide.</p><p>{result.intervention.sampling} {result.intervention.donor_count != null && `Donor pool: ${result.intervention.donor_count} bins; mean ${number(result.intervention.donor_mean)}; zero fraction ${number(100*result.intervention.donor_zero_fraction!,1)}%.`}</p>
      <dl>{Object.entries(result.provenance).map(([k,v]) => <div key={k}><dt>{k.replace(/_/g,' ')}</dt><dd>{v}</dd></div>)}</dl>
    </details>
  </section>
}

export default function ExperimentView() {
  const cell = useStore(s => s.cell), chrom = useStore(s => s.chrom), meta = useStore(s => s.meta), bins = useStore(s => s.bins)
  const cellData = useStore(s => cell ? s.cells[cell] : undefined)
  const selections = useStore(s => s.selections), active = useStore(s => s.activeSelectionId)
  const key = cell && chrom ? experimentKey(cell, chrom) : ''
  const entry = useExperimentStore(s => s.entries[key]), update = useExperimentStore(s => s.update), run = useExperimentStore(s => s.run)
  const initialize = useExperimentStore(s => s.initialize)
  const runSeries = useExperimentStore(s => s.runSeries), cancelSeries=useExperimentStore(s=>s.cancelSeries), selectResult=useExperimentStore(s=>s.selectResult)
  const [previewState, setPreviewState] = useState<{ signature: string; data?: { original: number[]; edited: number[]; intervention: ExperimentResult['intervention'] }; error?: string }>({ signature: '' })
  const cm = meta?.chromosomes.find(c => c.name === chrom)
  const regionInputs = useRef<{ start: HTMLInputElement | null; end: HTMLInputElement | null }>({ start: null, end: null })
  // The packed bin buffer has [start..., end..., chromosome-local index...];
  // coordinates below come from the authoritative node response instead.
  const nodes = useStore(s => s.nodes)
  const root = cell && chrom && cm ? nodes[`${cell}|${chrom}|0|${cm.n}`] : undefined
  useEffect(() => {
    if (!key || !cm || !root) return
    const region = meta?.regions.find(r => r.chrom === chrom && r.id === '9p21')
    const idx = region ? root.start_bp.map((_,i)=>i).filter(i => root.end_bp[i]>=region.start_mb*1e6 && root.start_bp[i]<=region.end_mb*1e6) : []
    initialize(key, { start: idx.length ? idx[0] : 0, end: idx.length ? idx[idx.length-1]+1 : Math.min(40, cm.n), operation: 'multiply', value: 0.5, preserve_total: false, cn_state:'loss',direction:'right',length:40,seed:42,k:2 })
  }, [key, cm, root])
  const signature=JSON.stringify({cell,chrom,...entry?.draft})
  useEffect(()=>{
    if(!entry || !cell || !chrom) return
    const controller=new AbortController()
    const timer=setTimeout(async()=>{
      try {
        const response=await fetch('/api/experiment/preview',{method:'POST',headers:{'Content-Type':'application/json'},body:signature,signal:controller.signal})
        const data=await response.json()
        if(!response.ok) throw new Error(typeof data.detail==='string'?data.detail:JSON.stringify(data.detail))
        setPreviewState({signature,data})
      }catch(e:any){ if(e.name!=='AbortError') setPreviewState({signature,error:e.message}) }
    },150)
    return()=>{clearTimeout(timer);controller.abort()}
  },[signature])
  const preview=previewState.signature===signature ? previewState.data : null
  if (!cell || !chrom) return <div className="empty"><h2>X-value experiments</h2><p>Open a cell from the cohort to edit a region and compare its algorithmic consequences.</p><button className="btn" onClick={() => useStore.getState().setTab('cohort')}>Choose a cell</button></div>
  if (!entry || !cm || !root || !cellData) return <div className="empty">Loading original X…</div>
  const d = entry.draft, set = (patch: Partial<ExperimentDraft>) => update(key, patch)
  const validBounds=Number.isInteger(d.start)&&Number.isInteger(d.end)&&d.start>=0&&d.end<=cm.n&&d.start<d.end
  const windows = selections.filter(w => w.cell===cell && w.chrom===chrom)
  const applyWindow = (id: string) => { const w=windows.find(v=>v.id===id); if(w) set({ start:w.nodeStart+w.s, end:w.nodeStart+w.e }) }
  const chooseMb = (lo: number, hi: number) => {
    const idx=root.start_bp.map((_,i)=>i).filter(i=>root.end_bp[i]>=lo*1e6 && root.start_bp[i]<=hi*1e6)
    if(Number.isFinite(lo)&&Number.isFinite(hi)&&lo<=hi&&idx.length) set({ start:idx[0], end:idx[idx.length-1]+1 })
    else useStore.getState().setToast('This interval contains no retained bins. Enter increasing bounds within this chromosome.')
  }
  const outdated=!!entry.result&&!sameDraft(d,entry.result.request)
  const span = preview ? (root.end_bp[d.end-1]-root.start_bp[d.start]+1)/1e6 : null
  const retained = preview ? root.end_bp.slice(d.start,d.end).reduce((a,v,i)=>a+v-root.start_bp[d.start+i]+1,0)/1e6 : null
  return <div className="experiment-view">
    <div className="section-heading"><div><span className="eyebrow">CONTROLLED REGIONAL EDIT</span><h1>X-value experiments</h1><p>Change normalized X in a region and trace the effect through AD, split tests, and CN assignment.</p></div><span className="source-label">Original data preserved</span></div>
    <fieldset className="experiment-editor" aria-label="Experiment settings" disabled={entry.series?.running}>
      <div className="experiment-control-row">
        <label>Use a region<select aria-label="Experiment region" defaultValue="" onChange={e => { const r=meta!.regions.find(v=>v.id===e.target.value); if(r) chooseMb(r.start_mb,r.end_mb); e.target.value='' }}><option value="">Named regions…</option>{meta!.regions.filter(r=>r.chrom===chrom).map(r=><option key={r.id} value={r.id}>{r.label}</option>)}</select></label>
        <label>Use an explorer window<select aria-label="Experiment window" value="" onChange={e=>applyWindow(e.target.value)}><option value="">{windows.length ? `${windows.length} available windows…` : 'Draw a window in Explore first'}</option>{windows.map(w=><option key={w.id} value={w.id}>{w.name} · {w.e-w.s} bins</option>)}</select></label>
        {active && windows.some(w=>w.id===active) && <button className="btn-sm" onClick={()=>applyWindow(active)}>Use active window</button>}
        <label>Start bin<input aria-label="Experiment start bin" type="number" min={0} max={cm.n-1} value={d.start} onChange={e=>set({start:e.target.valueAsNumber})} /></label>
        <label>End bin (exclusive)<input aria-label="Experiment end bin" type="number" min={1} max={cm.n} value={d.end} onChange={e=>set({end:e.target.valueAsNumber})} /></label>
      </div>
      <div className="experiment-control-row">
        <label>Start Mb<input aria-label="Experiment start Mb" ref={el=>{regionInputs.current.start=el}} key={`s${d.start}`} type="number" step="0.1" defaultValue={validBounds ? root.start_bp[d.start]/1e6 : 0} /></label>
        <label>End Mb<input aria-label="Experiment end Mb" ref={el=>{regionInputs.current.end=el}} key={`e${d.end}`} type="number" step="0.1" defaultValue={validBounds ? root.end_bp[d.end-1]/1e6 : 0} /></label>
        <button className="btn-sm" onClick={()=>chooseMb(regionInputs.current.start?.valueAsNumber ?? NaN,regionInputs.current.end?.valueAsNumber ?? NaN)}>Apply Mb bounds</button>
        <span className="muted">{d.end-d.start} retained bins · {number(span)} Mb span · {number(retained)} Mb covered</span>
      </div>
      <div className="experiment-control-row intervention-row">
        <label>Transformation<select aria-label="X transformation" value={d.operation} onChange={e=>set({operation:e.target.value as ExperimentDraft['operation'],value:e.target.value==='multiply'?0.5:0})}><option value="multiply">Multiply X by</option><option value="add">Add to X (clip at zero)</option><option value="set">Replace X with</option><option value="simulate">Simulate loss / base / gain</option><option value="extend">Extend with similar values</option><option value="duplicate">Duplicate ordered pattern</option></select></label>
        {['multiply','add','set'].includes(d.operation) && <label>{d.operation==='multiply'?'Factor':'Value'}<input aria-label="X edit value" type="number" step="0.1" value={d.value} onChange={e=>set({value:e.target.valueAsNumber})} /></label>}
        {d.operation==='simulate' && <label>Donor class<select aria-label="Simulation class" value={d.cn_state||'loss'} onChange={e=>set({cn_state:e.target.value as ExperimentDraft['cn_state']})}><option value="loss">Loss-called bins</option><option value="base">Base/diploid-called bins</option><option value="gain">Gain-called bins</option></select></label>}
        {['extend','duplicate'].includes(d.operation) && <><label>Direction<select aria-label="Extension direction" value={d.direction||'right'} onChange={e=>set({direction:e.target.value as 'left'|'right'})}><option value="left">Left</option><option value="right">Right</option></select></label><label>Additional bins<input aria-label="Extension bins" type="number" min={1} max={cm.n} value={d.length??40} onChange={e=>set({length:e.target.valueAsNumber})}/></label><button className="btn-xs" onClick={()=>set({length:d.end-d.start})}>One source-width</button></>}
        {['simulate','extend'].includes(d.operation) && <label>Sampling seed<input aria-label="Sampling seed" type="number" min={0} max={4294967295} value={d.seed??42} onChange={e=>set({seed:e.target.valueAsNumber})}/></label>}
        {d.operation==='multiply' && <div className="preset-factors">{[0,0.5,1,2].map(v=><button key={v} className={`btn-xs ${d.value===v?'on':''}`} onClick={()=>set({value:v})}>×{v}</button>)}</div>}
        <label className="chk-inline"><input type="checkbox" checked={d.preserve_total} onChange={e=>set({preserve_total:e.target.checked})} />Preserve original genome total</label>
        <label>Recursion limit k<select aria-label="Experiment recursion limit" value={d.k??2} onChange={e=>set({k:+e.target.value})}>{Array.from({length:10},(_,i)=>i+1).map(k=><option key={k} value={k}>{k}{k===2?' · production default':''}</option>)}</select></label>
        <button className="btn primary" disabled={!preview||entry.status==='running'} onClick={()=>run(cell,chrom)}>{entry.status==='running'?'Computing exact tests…':'Run experiment'}</button>
        <button className="btn" disabled={!preview||entry.status==='running'} onClick={()=>runSeries(cell,chrom,cm.n)}>Run width series</button>
      </div>
      <p className="editor-note">{d.preserve_total?'After the edit, every genome bin is rescaled to restore the original total.':'Genome total is allowed to change.'} Every run starts from original X. {['extend','duplicate'].includes(d.operation)?'Neighbouring bins are overwritten; genomic coordinates stay fixed.':'These edits do not rerun GC correction.'}</p>
      {d.operation==='simulate' && <p className="editor-note">Seeded sampling from this cell’s genome-wide bins with exact production CN code {d.cn_state==='gain'?'2':d.cn_state==='base'?'1':'0'}. Call-conditioned simulation is not biological ground truth and does not preserve spatial correlation.</p>}
      {preview?.intervention.donor_count!=null && <p className="editor-note">Donor pool: {preview.intervention.donor_count} bins · mean {number(preview.intervention.donor_mean)} · {number(100*preview.intervention.donor_zero_fraction!,1)}% zeros. {preview.intervention.actual_extension!=null&&`Extension: ${preview.intervention.actual_extension}/${preview.intervention.requested_extension} requested bins; total event ${preview.intervention.event_end-preview.intervention.event_start} bins.`}</p>}
      <p className="editor-note">Width series: 10, 25, 50, 100, 150 {['extend','duplicate'].includes(d.operation)?'additional bins':'bins centred on the selected region'}. Chromosome-edge clipping is reported by actual event width. Fixed sampling seed; one realization per width.</p>
      <p className="editor-note">k={d.k??2} tests node depths 0–{(d.k??2)-1}; depth 0 is the chromosome root. Deeper recursion requires accepted parents and can take longer. This setting applies to both scenarios and every width-series run.</p>
      {!preview && <p className="notice">{previewState.signature===signature&&previewState.error ? previewState.error : 'Updating X preview…'}</p>}
      {entry.status==='running' && <p className="computation-status" role="status">Computing original and edited nodes, 1,000-permutation tests, and full-genome CN fits. Navigation remains available.</p>}
      {entry.error && <p className="notice" role="alert">{entry.error}</p>}
    </fieldset>
    {entry.series?.running && <div className="computation-status" role="status">Width series: {entry.series.done}/{entry.series.total} complete. <button className="btn-sm" onClick={()=>cancelSeries(key)}>Stop after current run</button></div>}
    {preview && (!entry.result || outdated) && <section className="experiment-preview"><div className="section-heading"><h3>Live X preview</h3><span className="muted">Drag horizontally to select the source region · shading marks the edit target · AD and CN update on Run</span></div><ExperimentTracks starts={root.start_bp} ends={root.end_bp} original={preview.original} edited={preview.edited} start={preview.intervention.target_start} end={preview.intervention.target_end} onRegion={entry.series?.running?undefined:(start,end)=>set({start,end})} /></section>}
    {!!entry.history?.length && <details className="scientific-details" open><summary>Run history and size sensitivity · {entry.history.length} retained runs</summary><p className="muted">Each row is a frozen experiment, not a detection probability. Edge recovery requires the exact bin boundary. Compare rows with the same transformation, seed, and genomic location.</p><div className="table-scroll"><table className="tbl"><thead><tr><th>Experiment</th><th>Event bins / span</th><th>Root candidate</th><th>Left / right edge found</th><th>Combined event calls</th><th></th></tr></thead><tbody>{entry.history.map((v,i)=><tr key={i}><td>{v.request.operation} {v.request.operation==='simulate'?v.request.cn_state:['extend','duplicate'].includes(v.request.operation)?v.request.direction:v.request.value} · seed {v.request.seed??42} · k={v.request.k??2}</td><td>{v.intervention.event_end-v.intervention.event_start} / {number((v.end_bp[v.intervention.event_end-1]-v.start_bp[v.intervention.event_start]+1)/1e6)} Mb</td><td>{v.edited.nodes[0].argmax_b===null?'N/A':`${number(v.start_bp[v.edited.nodes[0].argmax_b!]/1e6)} Mb`}</td><td>{v.event_edges.map(e=>e.chromosome_edge?'chr end':e.accepted?'yes':'no').join(' / ')}</td><td>{v.cn.resegmented.event_calls ? Object.entries(v.cn.resegmented.event_calls.cn5).map(([state,count])=>`${state}: ${count}`).join(' · ') : 'Undefined'}</td><td><button className="btn-xs" disabled={entry.status==='running'||entry.series?.running} onClick={()=>selectResult(key,v)}>Inspect run {i+1}</button></td></tr>)}</tbody></table></div><button className="btn-sm" onClick={()=>download(JSON.stringify(entry.history,null,2),`${cell}-${chrom}-width-experiments.json`,'application/json')}>Export run history</button></details>}
    {entry.result && <ExperimentResults result={entry.result} outdated={outdated} onRegion={entry.series?.running?undefined:(start,end)=>set({start,end})} />}
  </div>
}
