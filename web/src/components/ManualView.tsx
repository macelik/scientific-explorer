import { useEffect, useMemo } from 'react'
import { BP0, LEFT, NODE_COLOR, RIGHT } from '../colors'
import { boundaryMb } from '../coords'
import { fmtP } from '../stats'
import { leafBoundaries, leaves, manualKey, nodeKey, testKey, useStore, type NodeRefLike } from '../store'
import type { Marker } from '../tracks'
import type { ManualNode, NodeData, SegmentEstimator } from '../types'
import CNResult from './CNResult'
import DisplayControls from './DisplayControls'
import NodeCard from './NodeCard'
import WindowWorkspace from './WindowWorkspace'

function SplitStatus({ node, data }: { node: ManualNode; data: NodeData }) {
  const tests = useStore((s) => s.tests)
  const meta = useStore((s) => s.meta)!
  const runTest = useStore((s) => s.runTest)
  if (node.splitB === null) return null
  const k = testKey(node.cell, node.chrom, node.start, node.end, node.splitB)
  const t = tests[k]
  const b = node.splitB
  const ad = data.ad[b - 1]
  const alpha = meta.defaults.alpha
  return (
    <div className="split-status">
      <div><b style={{ color: BP0 }}>manual split</b> at boundary {b} (chrom bin {node.start + b}) · {boundaryMb(data, b).toFixed(2)} Mb · AD = {ad.toFixed(3)} · children ({b}, {data.n - b}){data.argmax_b === b ? ' · this is the node argmax' : data.argmax_b !== null ? ` · node argmax is boundary ${data.argmax_b} (AD ${data.argmax_ad?.toFixed(2)})` : ''}</div>
      {!t || t.status === 'pending' ? <span className="badge pending">exact permutation tests pending ({meta.defaults.n_permutations} permutations, seed {meta.defaults.seed})…</span>
        : t.status === 'error' ? <span className="badge err">test error: {t.error} <button className="btn-xs" onClick={() => runTest(node.cell, node.chrom, node.start, node.end, b)}>retry</button></span>
          : (() => { const r = t.result!; const pass = r.p_local < alpha && r.p_global < alpha; return <span className={`badge ${pass ? 'ok' : 'fail'}`}>{pass ? 'PASSES' : 'FAILS'} strict p &lt; {alpha} on both · p_local {fmtP(r.p_local, r.p_local_exact)} · p_global {fmtP(r.p_global, r.p_global_exact)} · observed AD {r.observed_ad.toFixed(3)}{t.elapsed !== undefined ? ` · ${(t.elapsed * 1000).toFixed(0)} ms` : ''} · manual label: {pass ? 'would be accepted by production' : 'kept as an intentional manual boundary (tests are descriptive here)'}</span> })()}
    </div>
  )
}

function TreeNav({ tree, activeId, onPick }: { tree: any; activeId: string; onPick: (id: string) => void }) {
  const tests = useStore((s) => s.tests)
  const meta = useStore((s) => s.meta)!
  const render = (id: string): JSX.Element => {
    const n: ManualNode = tree.nodes[id]
    let st = ''
    if (n.splitB !== null) { const t = tests[testKey(n.cell, n.chrom, n.start, n.end, n.splitB)]; st = !t || t.status === 'pending' ? '⏳' : t.status === 'error' ? '⚠' : (t.result!.p_local < meta.defaults.alpha && t.result!.p_global < meta.defaults.alpha) ? '✓' : '✗' }
    return (
      <li key={id}>
        <button className={`tree-node ${id === activeId ? 'active' : ''}`} style={{ borderLeftColor: NODE_COLOR[n.side] }} onClick={() => onPick(id)}>
          <b>{n.label}</b> [{n.start}, {n.end}) · {n.end - n.start} bins {n.splitB !== null ? <span className="muted">· split @{n.start + n.splitB} {st}</span> : <span className="muted">· leaf</span>}
        </button>
        {n.splitB !== null && <ul>{render(n.leftId!)}{render(n.rightId!)}</ul>}
      </li>
    )
  }
  return <ul className="tree">{render(tree.rootId)}</ul>
}

export default function ManualView() {
  const cell = useStore((s) => s.cell); const chrom = useStore((s) => s.chrom)
  const meta = useStore((s) => s.meta)!
  const manual = useStore((s) => (s.cell && s.chrom ? s.manual[manualKey(s.cell, s.chrom)] : undefined))
  const nodes = useStore((s) => s.nodes)
  const cellData = useStore((s) => (s.cell ? s.cells[s.cell] : undefined))
  const display = useStore((s) => s.display)
  const initManual = useStore((s) => s.initManual)
  const place = useStore((s) => s.manualPlace); const remove = useStore((s) => s.manualRemove)
  const undo = useStore((s) => s.manualUndo); const redo = useStore((s) => s.manualRedo); const reset = useStore((s) => s.manualReset)
  const setActive = useStore((s) => s.manualSetActive); const clearNotice = useStore((s) => s.manualClearNotice)
  const finalize = useStore((s) => s.manualFinalize)
  const setEstimator = useStore(s=>s.manualSetEstimator)
  const setMode = useStore((s) => s.setMode); const mode = useStore((s) => s.mode)
  const ensureNode = useStore((s) => s.ensureNode)
  useEffect(() => { if (cell && chrom) initManual(cell, chrom) }, [cell, chrom])
  const tests = useStore((s) => s.tests); const runTest = useStore((s) => s.runTest)
  useEffect(() => {
    if (!manual) return
    for (const id in manual.tree.nodes) {
      const n = manual.tree.nodes[id]
      if (!nodes[nodeKey(n.cell, n.chrom, n.start, n.end)]) ensureNode(n.cell, n.chrom, n.start, n.end)
      // restored sessions: make sure every placed split has its exact tests requested (content-addressed, cached server-side)
      if (n.splitB !== null && !tests[testKey(n.cell, n.chrom, n.start, n.end, n.splitB)]) runTest(n.cell, n.chrom, n.start, n.end, n.splitB)
    }
  }, [manual?.tree.version, manual?.cell, manual?.chrom])
  const active = manual?.tree.nodes[manual.activeNodeId]
  const dataOf = (n?: ManualNode | null) => (n ? nodes[nodeKey(n.cell, n.chrom, n.start, n.end)] : undefined)
  const activeData = dataOf(active)
  const leftN = active?.leftId ? manual!.tree.nodes[active.leftId] : null, rightN = active?.rightId ? manual!.tree.nodes[active.rightId] : null
  const leftD = dataOf(leftN), rightD = dataOf(rightN)
  const crumbs = useMemo(() => { const out: ManualNode[] = []; let n = active; while (n) { out.unshift(n); n = n.parentId ? manual!.tree.nodes[n.parentId] : undefined } return out }, [active, manual?.tree.version])
  if (!cell || !chrom) return <div className="empty"><p>Pick a cell and chromosome first (cohort tab or a reference cell).</p></div>
  if (!manual || !active || !activeData || !cellData) return <div className="empty"><p>Preparing manual segmentation for {cell} / {chrom}…</p></div>
  const prod = cellData.production_breakpoints[chrom] || []
  const ref = (n: ManualNode): NodeRefLike => ({ id: n.id, version: n.version, side: n.side, depth: n.depth })
  const treeBoundaries = leafBoundaries(manual.tree)
  const markersFor = (n: ManualNode, d: NodeData): Marker[] => {
    const m: Marker[] = []
    for (const b of treeBoundaries) if (b > n.start && b < n.end) m.push({ chromB: b, color: BP0, label: b === n.start + (n.splitB ?? -1) ? 'manual split' : 'manual (deeper)', style: 'manual', width: b === n.start + (n.splitB ?? -1) ? 2 : 1.2 })
    if (d.argmax_b !== null && n.splitB !== d.argmax_b) m.push({ chromB: n.start + d.argmax_b, color: NODE_COLOR[n.side], label: 'argmax (not placed)', style: 'candidate' })
    if (display.showProduction) for (const p of prod) if (p.bp > n.start && p.bp < n.end && !treeBoundaries.includes(p.bp)) m.push({ chromB: p.bp, color: '#7c3aed', label: 'production (reference)', style: 'production', width: 1 })
    return m
  }
  const segs = (n: ManualNode, d: NodeData) => [0, ...treeBoundaries.filter((b) => b > n.start && b < n.end).map((b) => b - n.start), d.n]
  const canSplit = activeData.n >= 2
  const childCard = (n: ManualNode | null, d: NodeData | undefined, side: 'left' | 'right') => {
    if (!n) return null
    if (!d) return <div className="card empty">loading {side} child…</div>
    return <NodeCard key={n.id} node={d} nodeRef={ref(n)} title={`${n.label} · ${side} child · [${n.start}, ${n.end}) = ${d.n} bins · ${(d.start_bp[0] / 1e6).toFixed(1)}–${(d.end_bp[d.n - 1] / 1e6).toFixed(1)} Mb`} subtitle={d.argmax_b !== null ? `own argmax · ${boundaryMb(d, d.argmax_b).toFixed(2)} Mb (${d.argmax_b}, ${d.n - d.argmax_b}) · AD ${d.argmax_ad?.toFixed(2)} — recomputed on this child's own X` : 'single bin: terminal segment'} markers={markersFor(n, d)} segmentBoundaries={segs(n, d)} height={440} compact
      extraControls={<>{n.splitB === null && d.argmax_b !== null && <button className="btn-sm" onClick={() => place(cell, chrom, n.id, d.argmax_b!)}>split at argmax</button>}<button className="btn-sm" onClick={() => setActive(cell, chrom, n.id)}>continue here</button></>}
      badge={n.splitB !== null ? <span className="badge">already split @{n.start + n.splitB}</span> : <span className="badge muted">terminal segment</span>} />
  }
  return (
    <div className="manual">
      <DisplayControls manual />
      <div className="card row wrap">
        <strong>Manual segmentation</strong>
        <span className="muted">start: one unsplit chromosome; production boundaries (thin violet) are a reference overlay only</span>
        <span className="spacer" />
        <button className="btn-sm" onClick={() => undo(cell, chrom)} disabled={!manual.history.length}>↶ undo ({manual.history.length})</button>
        <button className="btn-sm" onClick={() => redo(cell, chrom)} disabled={!manual.future.length}>↷ redo ({manual.future.length})</button>
        <button className="btn-sm" onClick={() => reset(cell, chrom)}>reset tree</button>
        <label>Segment summary <select aria-label="CN segment summary" value={manual.cnEstimator??'arithmetic'} onChange={e=>setEstimator(cell,chrom,e.target.value as SegmentEstimator)}><option value="arithmetic">Arithmetic mean (production)</option><option value="iqr_upper">IQR mean · upper tail only, keep zeros</option><option value="iqr_two_sided">IQR mean · two-sided rule</option></select></label>
        <button className="btn" onClick={() => finalize(cell, chrom)} disabled={!leaves(manual.tree).length||manual.cnBusy}>{manual.cnBusy?'Computing CN…':`Finish & assign CN (${leaves(manual.tree).length} segments)`}</button>
      </div>
      {(manual.cnEstimator??'arithmetic')!=='arithmetic' && <div className="notice small">Experimental segment estimator: {manual.cnEstimator==='iqr_upper'?'keep 0 ≤ X ≤ Q3 + 1.5 × IQR':'keep Q1 − 1.5 × IQR ≤ X ≤ Q3 + 1.5 × IQR'} within each full-genome segment. Retained values determine the mean; original bin counts still weight the fit. AD, split tests, and the genome normalization baseline are unchanged. Arithmetic CN on the same manual boundaries is shown alongside.</div>}
      {manual.notice && <div className="notice">{manual.notice} <button className="btn-xs" onClick={() => clearNotice(cell, chrom)}>dismiss</button></div>}
      <div className="manual-grid">
        <div className="card tree-card">
          <div className="tbl-title">Tree navigator · v{manual.tree.version}</div>
          <TreeNav tree={manual.tree} activeId={manual.activeNodeId} onPick={(id) => setActive(cell, chrom, id)} />
          <div className="muted small">leaf boundaries (chrom-local): {treeBoundaries.length ? treeBoundaries.join(', ') : 'none'}</div>
        </div>
        <div className="manual-main">
          <div className="crumbs">{crumbs.map((c, i) => <span key={c.id}>{i > 0 && ' › '}<button className={`btn-xs ${c.id === active.id ? 'on' : ''}`} onClick={() => setActive(cell, chrom, c.id)}>{c.label} [{c.start},{c.end})</button></span>)}</div>
          {mode !== 'place' && canSplit && <div className="notice small">Switch to <b>Place breakpoint</b> mode and click any boundary on the active node's AD curve (not only the argmax) to insert a split.{activeData.argmax_b !== null && <> Or <button className="btn-xs" onClick={() => place(cell, chrom, active.id, activeData.argmax_b!)}>split at argmax ({activeData.argmax_b})</button></>}</div>}
          {mode === 'place' && canSplit && <div className="notice small">Place mode: click a point on the AD curve of <b>{active.label}</b>. {active.splitB !== null ? 'This will replace the current split and rebuild its descendants (undo restores them).' : ''}</div>}
          <NodeCard node={activeData} nodeRef={ref(active)} title={`Active node · ${active.label} · depth ${active.depth} · [${active.start}, ${active.end}) = ${activeData.n} bins`}
            subtitle={activeData.argmax_b !== null ? `argmax (reference only, not placed unless you choose it) · ${boundaryMb(activeData, activeData.argmax_b).toFixed(2)} Mb (${activeData.argmax_b}, ${activeData.n - activeData.argmax_b}) · AD ${activeData.argmax_ad?.toFixed(2)}` : 'single bin — cannot split'}
            markers={markersFor(active, activeData)} segmentBoundaries={segs(active, activeData)} height={520} onPlace={(b) => { place(cell, chrom, active.id, b); setMode('navigate') }}
            extraControls={active.splitB !== null ? <button className="btn-sm" onClick={() => remove(cell, chrom, active.id)}>remove split / subtree</button> : undefined}
            badge={active.splitB !== null ? <span className="badge">split placed</span> : <span className="badge muted">unsplit</span>} />
          <SplitStatus node={active} data={activeData} />
          {active.splitB !== null && <div className="children">{childCard(leftN, leftD, 'left')}{childCard(rightN, rightD, 'right')}</div>}
        </div>
      </div>
      {manual.cnResult && <CNResult res={manual.cnResult} outdated={manual.cnOutdated} />}
      <WindowWorkspace nodes={[{ ref: ref(active), data: activeData }, ...(leftN && leftD ? [{ ref: ref(leftN), data: leftD }] : []), ...(rightN && rightD ? [{ ref: ref(rightN), data: rightD }] : [])]} nodeLabel={(r) => { const n = manual.tree.nodes[r.id]; return n ? `${n.label} [${n.start},${n.end})` : r.side }} />
    </div>
  )
}
export { LEFT as _l, RIGHT as _r }
