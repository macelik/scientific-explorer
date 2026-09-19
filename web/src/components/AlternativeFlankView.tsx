import PairedFlankBootstrap from './PairedFlankBootstrap'
import { useRef, useState } from 'react'
import Plot from './Plot'
import { useIntegration } from '../integrationStore'
import { useStore } from '../store'

type Metric = 'mean' | 'median' | 'iqr_mean'
type Summary = Record<Metric,number> & { n:number; retained:number; std:number|null; cv:number|null; zero_fraction:number }
type Region = { start:number; end:number; summary:Summary }
type Row = Region & { key:string; source:string; start_bp:number; end_bp:number; left:Region & {effects:Record<Metric,number|null>}; right:Region & {effects:Record<Metric,number|null>} }
type Result = {chrom:string;sources:string[];max_bins:number;rows:Row[];profiles:Record<string,number[]>;start_bp:number[];end_bp:number[];provenance:unknown}
const names:Record<Metric,string>={mean:'Arithmetic mean',median:'Median',iqr_mean:'Two-sided IQR mean'}
const fmt=(x:number|null|undefined)=>x==null?'undefined':x.toFixed(3)
function download(text:string,name:string,type:string){const url=URL.createObjectURL(new Blob([text],{type}));const a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)}

export default function AlternativeFlankView(){
 const ig=useIntegration(), meta=ig.meta!
 const ex=useStore(s=>s.meta)!, presentation=useStore(s=>s.presentation)
 const [source,setSource]=useState('cluster4'),[max,setMax]=useState(100),[metric,setMetric]=useState<Metric>('mean')
 const [result,setResult]=useState<Result|null>(null),[selected,setSelected]=useState(''),[busy,setBusy]=useState(false),[error,setError]=useState('')
 const token=useRef(0)
 const sources=source==='selected'?ig.seg.sources:[source]
 const request={chrom:ig.seg.chrom,sources,max_bins:max}, signature=JSON.stringify(request)
 const [frozen,setFrozen]=useState('')
 const stale=!!result&&frozen!==signature
 const generate=async()=>{
  const id=++token.current;setBusy(true);setError('')
  try{const response=await fetch('/api/integration/flank-view',{method:'POST',headers:{'Content-Type':'application/json'},body:signature});if(!response.ok)throw new Error((await response.json()).detail??response.statusText)
   const data:Result=await response.json();if(token.current!==id)return
   setResult(data);setFrozen(signature);setSelected(data.rows[0]?.key??'')
  }catch(e:any){if(token.current===id)setError(String(e.message))}finally{if(token.current===id)setBusy(false)}
 }
 const row=result?.rows.find(r=>r.key===selected)
 const points=result?.rows.filter(r=>r.left.effects[metric]!==null&&r.right.effects[metric]!==null)??[]
 const layout={height:350,margin:{l:65,r:20,t:25,b:55},paper_bgcolor:'#fff',plot_bgcolor:'#fff',font:{family:'system-ui',size:12},hovermode:presentation?false:'closest'}
 const csv=()=>{
  if(!result)return
  const rows=[['source','chrom','start_bin','end_bin_exclusive','side','flank_start','flank_end_exclusive','estimator','segment_level','flank_level','log2_effect','segment_bins','segment_iqr_retained','flank_bins','flank_iqr_retained'],
   ...result.rows.flatMap(r=>(['left','right'] as const).flatMap(side=>(Object.keys(names) as Metric[]).map(m=>[r.source,result.chrom,r.start,r.end,side,r[side].start,r[side].end,m,r.summary[m],r[side].summary[m],r[side].effects[m]??'',r.summary.n,r.summary.retained,r[side].summary.n,r[side].summary.retained])))]
  download(rows.map(r=>r.join(',')).join('\n')+'\n',`alternative-flanks-${result.chrom}.csv`,'text/csv')
 }
 const selectedRegions=row?[{label:'Left flank',color:'#2878a0',...row.left},{label:'Segment',color:'#96491d',...row},{label:'Right flank',color:'#27806d',...row.right}]:[]
 const tracks=selectedRegions.map(r=>{
  const x:(number|null)[]=[],y:(number|null)[]=[]
  for(let i=r.start;i<r.end;i++){if(i>r.start&&result!.start_bp[i]>result!.end_bp[i-1]+1){x.push(null);y.push(null)}x.push(result!.start_bp[i]/1e6);y.push(result!.profiles[row!.source][i])}
  return {type:'scatter',mode:'lines+markers',name:r.label,x,y,marker:{size:3},line:{color:r.color,width:1},connectgaps:false}
 })
 return <section className="panel alternative-flanks">
  <h3>Alternative flank view</h3>
  <p className="lead">Recalculate descriptive effects on fixed accepted segments. Compare both flanks and inspect the actual bin values; no weaker-flank selection is hidden in the plot.</p>
  <div className="controls-row">
   <label>Chromosome<select aria-label="Alternative chromosome" value={ig.seg.chrom} onChange={e=>ig.setSeg({chrom:e.target.value})}>{ex.chromosomes.map(c=><option key={c.name}>{c.name}</option>)}</select></label>
   <label>Source<select aria-label="Alternative source" value={source} onChange={e=>setSource(e.target.value)}>{meta.sources.map(s=><option key={s}>{s}</option>)}<option value="selected">Selected track rows</option></select></label>
   <label>Maximum internal segment bins<input aria-label="Alternative maximum bins" type="number" min={1} max={10000} value={max} onChange={e=>setMax(+e.target.value)} /></label>
   <button className="btn primary" disabled={busy||!sources.length||!Number.isInteger(max)||max<1||max>10000} onClick={generate}>{busy?'Generating…':'Generate alternative flank view'}</button>
  </div>
  {error&&<p role="alert">{error}</p>}
  {result&&<>
   <p className={stale?'notice':'muted'} role="status">{stale?'Settings changed; showing the previous generated result. ':''}{result.rows.length} internal segments on {result.chrom} · {result.sources.join(', ')} · ≤ {result.max_bins} retained bins. Fixed WNN res-0.3 membership.</p>
   <div className="controls-row"><label>Level estimate<select aria-label="Alternative level estimate" value={metric} onChange={e=>setMetric(e.target.value as Metric)}>{(Object.keys(names) as Metric[]).map(m=><option key={m} value={m}>{names[m]}</option>)}</select></label>
    <button className="btn-sm" onClick={csv}>Download alternative CSV</button><button className="btn-sm" onClick={()=>download(JSON.stringify(result,null,2),`alternative-flanks-${result.chrom}.json`,'application/json')}>Download alternative JSON</button>
   </div>
   <p className="small muted">IQR means retain values within Q1 − 1.5 × IQR and Q3 + 1.5 × IQR. Zero-level ratios are undefined without adding a pseudocount. Per-bin spread is not uncertainty of a cluster estimate. SRD and merge decisions are not recomputed. Paired cell bootstrap is available separately for the selected segment.</p>
   {result.rows.length>0?<>
    <Plot data={[{type:'scatter',mode:'markers',x:points.map(r=>r.left.effects[metric]),y:points.map(r=>r.right.effects[metric]),customdata:points.map(r=>r.key),text:points.map(r=>`${r.source} bins [${r.start},${r.end}) · ${r.summary.n} bins`),marker:{color:points.map(r=>r.key===selected?'#96491d':'#2878a0'),size:points.map(r=>r.key===selected?11:7)},hovertemplate:'%{text}<br>left Δ %{x:.3f}<br>right Δ %{y:.3f}<extra></extra>'}]} layout={{...layout,xaxis:{title:'log2(segment / LEFT flank)',zeroline:true},yaxis:{title:'log2(segment / RIGHT flank)',zeroline:true}}} onClick={e=>{const key=e?.points?.[0]?.customdata;if(key)setSelected(key)}} />
    <p className="small muted">Each dot is one segment. Both positive: elevated versus both neighbors; both negative: depressed versus both. Opposite signs: intermediate levels. {result.rows.length-points.length} segments with an undefined axis omitted from this plot; all remain selectable below.</p>
    <label>Inspect segment <select aria-label="Alternative segment" value={selected} onChange={e=>setSelected(e.target.value)}>{result.rows.map(r=><option key={r.key} value={r.key}>{r.source} · [{r.start}, {r.end}) · {(r.start_bp/1e6).toFixed(2)}–{(r.end_bp/1e6).toFixed(2)} Mb</option>)}</select></label>
    {row&&<>
     <Plot data={tracks} layout={{...layout,height:300,xaxis:{title:`${result.chrom} position (Mb)`},yaxis:{title:'Mean normalized X over member cells'},legend:{orientation:'h',y:1.1}}} />
     <Plot data={selectedRegions.map(r=>({type:'box',name:r.label,y:result.profiles[row.source].slice(r.start,r.end),boxpoints:'all',jitter:0,pointpos:0,marker:{size:3,color:r.color},line:{color:r.color}}))} layout={{...layout,height:280,yaxis:{title:'Per-bin pseudobulk X'},showlegend:false}} />
     <div className="table-scroll"><table className="tbl"><thead><tr><th>Region</th><th>Bins [start,end)</th><th>Mean</th><th>Median</th><th>IQR mean</th><th>IQR retained / total</th><th>Sample std</th><th>CV</th><th>Zero fraction</th><th>Δ ({names[metric]})</th></tr></thead><tbody>{selectedRegions.map((r,i)=><tr key={r.label}><td>{r.label}</td><td>[{r.start},{r.end})</td><td>{fmt(r.summary.mean)}</td><td>{fmt(r.summary.median)}</td><td>{fmt(r.summary.iqr_mean)}</td><td>{r.summary.retained} / {r.summary.n}</td><td>{fmt(r.summary.std)}</td><td>{fmt(r.summary.cv)}</td><td>{fmt(r.summary.zero_fraction)}</td><td>{i===1?'reference segment':fmt(row[i===0?'left':'right'].effects[metric])}</td></tr>)}</tbody></table></div>
    <PairedFlankBootstrap key={`${row.key}|${frozen}`} source={row.source} chrom={result.chrom} start={row.start} end={row.end} />
    </>}
   </>:<p>No internal segments satisfy this width limit. This is a selection result, not evidence of no biological events.</p>}
   <details><summary>Generation provenance</summary><pre style={{whiteSpace:'pre-wrap',overflowWrap:'anywhere'}}>{JSON.stringify(result.provenance,null,2)}</pre></details>
  </>}
 </section>
}
