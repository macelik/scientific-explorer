import assert from 'node:assert/strict'
import { experimentOverlays, nodeDepth } from '../web/src/experimentOverlays'
import { sameDraft } from '../web/src/experimentStore'
const node=(start:number,end:number,side:string,b:number,accepted:boolean)=>({start,end,side,argmax_b:b,accepted,hypothetical:!accepted})
const r:any={start_bp:[100,200,500,600,700,800,900,1000],end_bp:[199,299,599,699,799,899,999,1099],
 intervention:{target_start:2,target_end:6,event_start:1,event_end:7},
 original:{nodes:[node(0,4,'left',2,true)],boundaries:[2,4,6]},
 edited:{nodes:[node(0,5,'left',3,false)],boundaries:[]}}
const shapes=experimentOverlays(r,'left')
const rect=shapes.filter(s=>s.type==='rect')
assert.equal(rect.length,1)
assert.equal(rect[0].x0,500/1e6) // target, not larger source+extension event
assert.equal(rect[0].x1,799/1e6) // clipped to union of displayed child bounds
assert.equal(rect[0].yref,'paper')
const accepted=shapes.filter(s=>s.name?.includes('accepted'))
assert.equal(accepted.length,1) // node edges and outside boundaries excluded
assert.equal(accepted[0].x0,500/1e6) // chromosome-local bin at right-bin start
assert.equal(accepted[0].line.dash,'solid')
assert.equal(shapes.filter(s=>s.name?.includes('candidate')).length,1)
r.intervention.target_start=6;r.intervention.target_end=8
assert.equal(experimentOverlays(r,'left').filter(s=>s.type==='rect').length,0)
assert.equal(nodeDepth({side:'left.right'} as any),2)
const draft:any={start:1,end:3,operation:'multiply',value:0.5,preserve_total:false}
assert(sameDraft(draft,{...draft,k:2}))
assert(!sameDraft(draft,{...draft,k:4}))
console.log('PASS clipped edit targets, boundary positions, accepted/candidate distinction, legacy depth and stale settings')
