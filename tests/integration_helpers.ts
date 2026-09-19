import assert from 'node:assert/strict'
import {pickFlank,withGaps,median} from '../web/src/integration'
const row:any={left:{delta_mean:1,srd_phi:null},right:{delta_mean:null,srd_phi:3}}
assert.equal(pickFlank(row,'max_delta','delta_mean'),'left')
assert.equal(median([1,4]),2.5)
const g=withGaps([.000001,.100001,.300001],[1,2,3],[.1,.2,.4])
assert.equal(g.y.filter(v=>v===null).length,1,'A single omitted bin must break the line')
console.log('PASS: missing-flank selection, even median, exact coordinate gaps')
