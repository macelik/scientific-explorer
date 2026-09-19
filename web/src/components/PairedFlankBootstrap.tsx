import { useEffect, useRef, useState } from 'react'
import Plot from './Plot'
import { useStore } from '../store'
type Side={observed:number|null;valid:number;invalid:number;ci95:[number,number]|null;se:number|null;sign_agreement:number|null;effects:number[]}
type Result={left:Side;right:Side;n_cells:number;replicates:number;seed:number;interpretation:string}
const fmt=(x:number|null)=>x===null?'undefined':x.toFixed(4)
export default function PairedFlankBootstrap({source,chrom,start,end}:{source:string;chrom:string;start:number;end:number}){
 const [result,setResult]=useState<Result|null>(null),[busy,setBusy]=useState(false),[error,setError]=useState('')
 const controller=useRef<AbortController|null>(null),presentation=useStore(s=>s.presentation)
 useEffect(()=>()=>controller.current?.abort(),[])
 const run=async()=>{setBusy(true);setError('');const c=new AbortController();controller.current=c
  try{const r=await fetch('/api/integration/flank-bootstrap',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({source,chrom,start,end}),signal:c.signal});const body=await r.json();if(!r.ok)throw new Error(body.detail??r.statusText);setResult(body)}catch(e:any){if(e.name!=='AbortError')setError(String(e.message))}finally{setBusy(false)}
 }
 return <div className="paired-bootstrap">
  <h4>Stability across member cells</h4>
  <p className="small muted">Arithmetic-mean effects only. Resample the same member cells jointly for left flank, segment and right flank (200 replicates, seed 42). Boundaries and membership stay fixed. Percentile intervals are exploratory and do not account for selecting these segments.</p>
  <button className="btn" disabled={busy} onClick={run}>{busy?'Resampling member cells…':'Compute paired cell bootstrap'}</button>
  {error&&<p role="alert">{error}</p>}
  {result&&<>
   <p>{result.n_cells} member cells · {result.replicates} paired replicates.</p>
   <table className="tbl"><thead><tr><th>Flank</th><th>Observed mean effect</th><th>95% percentile interval</th><th>Bootstrap SE</th><th>Valid / attempted</th><th>Sign agreement</th></tr></thead><tbody>{(['left','right'] as const).map(side=>{const r=result[side];return <tr key={side}><td>{side}</td><td>{fmt(r.observed)}</td><td>{r.ci95?r.ci95.map(fmt).join(' to '):'undefined'}</td><td>{fmt(r.se)}</td><td>{r.valid} / {result.replicates}</td><td>{fmt(r.sign_agreement)}</td></tr>})}</tbody></table>
   <Plot data={(['left','right'] as const).map(side=>({type:'histogram',name:`${side} flank`,x:result[side].effects,opacity:.65,marker:{color:side==='left'?'#2878a0':'#27806d'}}))} layout={{height:260,margin:{l:50,r:15,t:30,b:50},barmode:'overlay',xaxis:{title:'Bootstrap log2(segment mean / flank mean)'},yaxis:{title:'Replicates'},legend:{orientation:'h'},hovermode:presentation?false:'closest'}} />
   <button className="btn-sm" onClick={()=>{const url=URL.createObjectURL(new Blob([JSON.stringify(result,null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download=`paired-bootstrap-${source}-${chrom}-${start}.json`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)}}>Download bootstrap JSON</button>
  </>}
 </div>
}
