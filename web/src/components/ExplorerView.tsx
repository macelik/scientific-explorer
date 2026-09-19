import { useEffect, useMemo, useState } from 'react'
import { BP0, LEFT, NODE_COLOR, RIGHT } from '../colors'
import { boundaryMb } from '../coords'
import { fmtP } from '../stats'
import { explorerNodes, nodeKey, testKey, useStore } from '../store'
import type { Marker } from '../tracks'
import type { NodeData } from '../types'
import NodeCard from './NodeCard'
import WindowWorkspace from './WindowWorkspace'
import DisplayControls from './DisplayControls'

export function TestBadge({ cell, chrom, node, b, inProduction, prodP, parentAccepted, alpha }: { cell: string; chrom: string; node: NodeData; b: number | null; inProduction: boolean; prodP?: { p_local: number; p_global: number }; parentAccepted: boolean; alpha: number }) {
  const tests = useStore((s) => s.tests)
  const runTest = useStore((s) => s.runTest)
  const k = b !== null ? testKey(cell, chrom, node.start, node.end, b) : ''
  const t = k ? tests[k] : undefined
  useEffect(() => { if (b !== null && !inProduction && !t) runTest(cell, chrom, node.start, node.end, b) }, [b, inProduction, k])
  if (b === null) return <span className="badge muted">single bin: cannot split</span>
  if (inProduction && parentAccepted) return <span className="badge ok" title={prodP ? `production p_local=${prodP.p_local.toPrecision(4)}, p_global=${prodP.p_global.toPrecision(4)}` : ''}>accepted by production{prodP ? ` · p_local ${prodP.p_local.toPrecision(3)} · p_global ${prodP.p_global.toPrecision(3)}` : ''}</span>
  const hypo = !parentAccepted ? 'hypothetical (parent split not accepted by production) · ' : ''
  if (!t || t.status === 'pending') return <span className="badge pending">{hypo}tests pending…</span>
  if (t.status === 'error') return <span className="badge err" title={t.error}>test error: {t.error}</span>
  const r = t.result!
  const pass = r.p_local < alpha && r.p_global < alpha
  return <span className={`badge ${pass ? 'ok' : 'fail'}`} title={`observed AD ${r.observed_ad.toFixed(4)} · ${r.n_permutations} permutations, seed ${r.seed}`}>{hypo}{pass ? 'PASSES' : 'FAILS'} both tests · p_local {fmtP(r.p_local, r.p_local_exact)} · p_global {fmtP(r.p_global, r.p_global_exact)}{inProduction ? '' : ' · not a production breakpoint'}</span>
}

export default function ExplorerView() {
  const [compactTracks, setCompactTracks] = useState(true)
  const meta = useStore((s) => s.meta)!
  const cell = useStore((s) => s.cell)
  const chrom = useStore((s) => s.chrom)
  const nodes = useStore((s) => s.nodes)
  const focus = useStore((s) => s.focus)
  const setFocus = useStore((s) => s.setFocus)
  const display = useStore((s) => s.display)
  const cellData = useStore((s) => (s.cell ? s.cells[s.cell] : undefined))
  const ex = useMemo(() => explorerNodes({ cell, chrom, meta, nodes } as any), [cell, chrom, meta, nodes])
  if (!cell || !chrom) return <div className="empty"><p>Pick a cell from the cohort karyogram (or a reference cell) to open the explorer.</p></div>
  if (!ex.root || !cellData) return <div className="empty"><p>Loading {cell} / {chrom}…</p></div>
  const root = ex.root.data, left = ex.left?.data, right = ex.right?.data
  const b0 = root.argmax_b
  const prod = cellData.production_breakpoints[chrom] || []
  const prodOf = (b: number) => prod.find((p) => p.bp === b)
  const rootAccepted = root.argmax_in_production
  const alpha = meta.defaults.alpha
  const fmtB = (n: NodeData, b: number | null) => b === null ? '—' : `${boundaryMb(n, b).toFixed(2)} Mb (${b}, ${n.n - b}) · AD = ${n.argmax_ad?.toFixed(2)}`

  const childMarkers = (side: 'left' | 'right', n: NodeData): Marker[] => n.argmax_b === null ? [] : [{ chromB: n.start + n.argmax_b, color: NODE_COLOR[side], label: `${side} argmax`, style: n.argmax_in_production && rootAccepted ? 'accepted' : 'candidate' }]
  const rootMarkers: Marker[] = []
  if (b0 !== null) rootMarkers.push({ chromB: b0, color: BP0, label: `depth 0 ${rootAccepted ? 'accepted' : 'candidate'}`, style: rootAccepted ? 'accepted' : 'candidate' })
  if (left?.argmax_b != null) rootMarkers.push({ chromB: left.start + left.argmax_b, color: LEFT, label: 'left child argmax', style: 'child' })
  if (right?.argmax_b != null) rootMarkers.push({ chromB: right.start + right.argmax_b, color: RIGHT, label: 'right child argmax', style: 'child' })
  if (display.showProduction) for (const p of prod) if (!rootMarkers.some((m) => m.chromB === p.bp)) rootMarkers.push({ chromB: p.bp, color: BP0, label: 'production', style: 'production', width: 1 })
  const rootSegs = [0, ...(left?.argmax_b != null ? [left.argmax_b] : []), ...(b0 !== null ? [b0] : []), ...(right?.argmax_b != null ? [right.start + right.argmax_b] : []), root.n]
  const prodInChild = (n: NodeData) => display.showProduction ? prod.filter((p) => p.bp > n.start && p.bp < n.end && p.bp !== n.start + (n.argmax_b ?? -1)).map((p) => ({ chromB: p.bp, color: BP0, label: 'production', style: 'production' as const, width: 1 })) : []
  const adRange: [number, number] | null = display.sharedAdScale ? [0, Math.max(root.argmax_ad ?? 0, left?.argmax_ad ?? 0, right?.argmax_ad ?? 0) * 1.12] : null
  const childAdRange: [number, number] | null = display.sharedAdScale ? [0, Math.max(left?.argmax_ad ?? 0, right?.argmax_ad ?? 0) * 1.15] : null

  const childCard = (side: 'left' | 'right', n: NodeData | undefined, ref: any) => {
    if (!n || !ref) return <div className="card empty">no {side} child (depth-0 split unavailable)</div>
    const focused = focus === side
    return <NodeCard node={n} nodeRef={ref} title={`${side === 'left' ? 'Left' : 'Right'} child · depth 1 · [${n.start}, ${n.end}) = ${n.n} bins · ${(n.start_bp[0] / 1e6).toFixed(1)}–${(n.end_bp[n.n - 1] / 1e6).toFixed(1)} Mb`}
      subtitle={`own argmax · ${fmtB(n, n.argmax_b)} · AD recomputed on this child's own X (not a crop of the parent curve)`}
      markers={[...childMarkers(side, n), ...prodInChild(n)]} segmentBoundaries={n.argmax_b !== null ? [0, n.argmax_b, n.n] : [0, n.n]} adRange={childAdRange} height={focused ? 560 : 470} compact
      badge={<TestBadge cell={cell} chrom={chrom} node={n} b={n.argmax_b} inProduction={n.argmax_in_production} prodP={n.argmax_b !== null ? prodOf(n.start + n.argmax_b) : undefined} parentAccepted={rootAccepted} alpha={alpha} />}
      extraControls={<button className="btn-sm" onClick={() => setFocus(focused ? null : side)}>{focused ? 'unfocus' : 'focus'}</button>} />
  }

  return (
    <div className="explorer">
      <div className="workspace-toolbar"><strong>Signal → candidates → recursion</strong><span className="spacer"/><label><input type="checkbox" checked={compactTracks} onChange={e=>setCompactTracks(e.target.checked)}/>Compact chromosome overview</label><button className="btn-sm" onClick={()=>document.getElementById('window-comparisons')?.scrollIntoView({behavior:'smooth'})}>Window comparisons ↓</button><button className="btn-sm" onClick={()=>useStore.getState().setTab('experiment')}>Experiment with X →</button></div>
      <DisplayControls />
      {focus ? (
        <div className="context-strip">
          <span>Whole-chromosome context collapsed: <b style={{ color: BP0 }}>depth 0</b> {fmtB(root, b0)} · {rootAccepted ? 'accepted' : 'candidate'} · <button className="btn-xs" onClick={() => setFocus(null)}>show all</button></span>
          <span className="muted">focus: {focus} child</span>
        </div>
      ) : (
        <NodeCard node={root} nodeRef={ex.root.ref} title={`Whole chromosome · depth 0 · ${chrom} · ${root.n} bins`}
          subtitle={`depth-0 argmax · ${fmtB(root, b0)} · dashed lines = where each child places its own argmax`}
          markers={rootMarkers} segmentBoundaries={rootSegs} adRange={adRange} height={compactTracks ? 380 : 540}
          badge={<TestBadge cell={cell} chrom={chrom} node={root} b={b0} inProduction={rootAccepted} prodP={b0 !== null ? prodOf(b0) : undefined} parentAccepted={true} alpha={alpha} />} />
      )}
      {!rootAccepted && b0 !== null && <div className="notice">Production would stop at this chromosome: the depth-0 candidate is not an accepted production breakpoint. Child inspection below is hypothetical / exploratory.</div>}
      <div className={`children ${focus ? 'focused' : ''}`}>
        {(!focus || focus === 'left') && childCard('left', left, ex.left?.ref)}
        {(!focus || focus === 'right') && childCard('right', right, ex.right?.ref)}
      </div>
      <section id="window-comparisons"><WindowWorkspace nodes={[ex.root, ex.left, ex.right].filter(Boolean) as any} /></section>
    </div>
  )
}
