import { useEffect, useMemo, useRef, useState } from 'react'
import { NODE_COLOR } from '../colors'
import { boundaryMb, mbToBoundary, selectionCoords } from '../coords'
import { useStore, type NodeRefLike } from '../store'
import { buildTracks, type Marker } from '../tracks'
import type { NodeData } from '../types'
import Plot, { downloadPlot } from './Plot'

export interface NodeCardProps {
  node: NodeData; nodeRef: NodeRefLike; title: string; subtitle?: string
  markers: Marker[]; segmentBoundaries: number[]
  adRange?: [number, number] | null; height?: number; compact?: boolean
  onPlace?: (b: number) => void
  extraControls?: React.ReactNode
  badge?: React.ReactNode
}

export default function NodeCard(p: NodeCardProps) {
  const { node, nodeRef } = p
  const display = useStore((s) => s.display)
  const mode = useStore((s) => s.mode)
  const presentation = useStore((s) => s.presentation)
  const meta = useStore((s) => s.meta)!
  const cellData = useStore((s) => s.cells[node.cell])
  const allSelections = useStore((s) => s.selections)
  const activeSelectionId = useStore((s) => s.activeSelectionId)
  const addSelection = useStore((s) => s.addSelection)
  const updateSelection = useStore((s) => s.updateSelection)
  const setActiveSelection = useStore((s) => s.setActiveSelection)
  const setToast = useStore((s) => s.setToast)
  const plotEl = useRef<any>(null)

  // selections shown here: originating on this node or on a descendant node (this node is an ancestor)
  const selections = useMemo(() => allSelections.filter((s) => s.cell === node.cell && s.chrom === node.chrom && s.nodeStart >= node.start && s.nodeEnd <= node.end), [allSelections, node])
  const editable = useMemo(() => new Set(selections.filter((s) => s.nodeId === nodeRef.id).map((s) => s.id)), [selections, nodeRef.id])
  const overlaps = useMemo(() => {
    const own = selections.filter((s) => s.nodeId === nodeRef.id && s.visible)
    const out: { s: number; e: number }[] = []
    for (let i = 0; i < own.length; i++) for (let j = i + 1; j < own.length; j++) {
      const a = own[i], b = own[j]; const s = Math.max(a.s, b.s), e = Math.min(a.e, b.e)
      if (e > s) out.push({ s, e })
    }
    return out
  }, [selections, nodeRef.id])

  const fig = useMemo(() => buildTracks({
    node, side: nodeRef.side, title: p.title, display,
    genomeBaseline: { trimmed: cellData?.genome.baseline_trimmed ?? 0, median: cellData?.genome.baseline_median ?? 0 },
    selections, editableSelectionIds: editable, activeSelectionId,
    markers: p.markers, regions: meta.regions, segmentBoundaries: p.segmentBoundaries, mode, presentation,
    adRange: p.adRange, height: p.height, compact: p.compact, overlapRanges: overlaps,
  }), [node, nodeRef.side, p.title, display, cellData, selections, editable, activeSelectionId, p.markers, meta.regions, p.segmentBoundaries, mode, presentation, p.adRange, p.height, p.compact, overlaps])
  const shapeIds = useRef(fig.shapeSelectionIds)
  shapeIds.current = fig.shapeSelectionIds

  // modern Plotly shape editing: per-shape `editable`, click activates, then drag edges/corners or the interior
  const config = useMemo(() => ({ modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d'], displayModeBar: presentation ? false : 'hover' }), [presentation])

  const makeSel = (s: number, e: number, anchor: number | null, source: 'drawn' | 'peak') => addSelection({
    cell: node.cell, chrom: node.chrom, nodeId: nodeRef.id, nodeVersion: nodeRef.version, nodeStart: node.start, nodeEnd: node.end, nodeSide: nodeRef.side, s, e, anchor, source,
  })

  const onSelected = (e: any) => {
    if (mode !== 'select' || !e || !e.range || !e.range.x) return
    const [m0, m1] = e.range.x as [number, number]
    const s = mbToBoundary(node, Math.min(m0, m1)), en = mbToBoundary(node, Math.max(m0, m1))
    if (en > s) { const sel = makeSel(s, en, null, 'drawn'); if (sel) { const c = selectionCoords(node, sel.s, sel.e); setToast(`${sel.name}: bins [${c.chromStart}, ${c.chromEnd}) · ${c.nBins} bins · ${c.startMb.toFixed(2)}–${c.endMb.toFixed(2)} Mb`) } }
    // clear plotly's own selection box
    try { const el = plotEl.current; if (el) (window as any).Plotly?.restyle?.(el, { selectedpoints: [null] }) } catch { /* ignore */ }
  }
  const onPlainClick = (xMb: number, row: number) => {
    // clicking inside one of this node's own windows activates it (for dragging/inspection) in any mode
    const bin = mbToBoundary(node, xMb)
    const hit = selections.filter((s) => s.nodeId === nodeRef.id && s.visible && bin >= s.s && bin <= s.e)
    if (hit.length) { const h = hit.sort((a, b) => (a.e - a.s) - (b.e - b.s))[0]; setActiveSelection(h.id); return }
    if (mode === 'navigate') return
    const j = mbToBoundary(node, xMb)
    if (j <= 0 || j >= node.n) { setToast(`boundary ${j} would leave an empty child; pick a boundary strictly inside the node`); return }
    if (mode === 'place') { if (row === 0) p.onPlace?.(j); else setToast('Click on the AD row to place a breakpoint'); return }
    if (mode === 'select') {
      const h = meta.defaults.half_window_bins
      const sel = makeSel(j - h, j + h, j, 'peak')
      if (sel) setToast(`${sel.name}: peak window around boundary ${j} (chrom bin ${node.start + j}) · ${sel.e - sel.s} bins${sel.clipped ? ` · clipped at node edge (${sel.clipped.left} left, ${sel.clipped.right} right)` : ''}`)
    }
  }
  // ---- drag handles overlay (own implementation; independent of Plotly's shape editing) ----
  const [tick, setTick] = useState(0)
  const dragging = useRef<{ id: string; kind: 'l' | 'r' | 'm'; startPx: number; s0: number; e0: number; anchor0: number | null } | null>(null)
  const geom = () => {
    const el = plotEl.current; const fl = el?._fullLayout
    if (!fl || !fl.xaxis || !fl.yaxis || !fl.yaxis3) return null
    const xa = fl.xaxis
    return { xa, left: xa._offset, top: fl.yaxis._offset, height: fl.yaxis3._offset + fl.yaxis3._length - fl.yaxis._offset, width: xa._length }
  }
  const own = selections.filter((s) => s.nodeId === nodeRef.id && s.visible)
  const handleGeoms = (() => {
    const g = geom(); if (!g) return []
    return own.map((s) => {
      const x0 = g.left + g.xa.d2p(boundaryMb(node, s.s)), x1 = g.left + g.xa.d2p(boundaryMb(node, s.e))
      return { s, x0: Math.max(g.left, x0), x1: Math.min(g.left + g.width, x1), top: g.top, height: g.height, offscreen: x1 < g.left || x0 > g.left + g.width }
    })
  })()
  const onHandleDown = (e: React.PointerEvent, id: string, kind: 'l' | 'r' | 'm') => {
    const sel = own.find((s) => s.id === id); if (!sel) return
    e.preventDefault(); e.stopPropagation()
    ;(e.target as HTMLElement).setPointerCapture(e.pointerId)
    dragging.current = { id, kind, startPx: e.clientX, s0: sel.s, e0: sel.e, anchor0: sel.anchor }
    setActiveSelection(id)
  }
  const onHandleMove = (e: React.PointerEvent) => {
    const d = dragging.current; if (!d) return
    const g = geom(); if (!g) return
    const dx = e.clientX - d.startPx
    const shiftBoundary = (b: number) => mbToBoundary(node, g.xa.p2d(g.xa.d2p(boundaryMb(node, b)) + dx))
    let s = d.s0, en = d.e0, anchor = d.anchor0
    if (d.kind === 'l') s = Math.max(0, Math.min(d.e0 - 1, shiftBoundary(d.s0)))
    else if (d.kind === 'r') en = Math.max(d.s0 + 1, Math.min(node.n, shiftBoundary(d.e0)))
    else {
      const width = d.e0 - d.s0
      let ns = shiftBoundary(d.s0)
      ns = Math.max(0, Math.min(node.n - width, ns))
      s = ns; en = ns + width
      if (anchor !== null) { anchor = anchor + (ns - d.s0); if (anchor <= 0 || anchor >= node.n) anchor = null }
    }
    if (anchor !== null && (anchor < s || anchor > en)) { /* anchor may leave the window when resizing; keep it as the alignment reference */ }
    const cur = selections.find((x) => x.id === d.id)
    if (cur && (cur.s !== s || cur.e !== en || cur.anchor !== anchor)) updateSelection(d.id, { s, e: en, anchor, clipped: undefined })
  }
  const onHandleUp = (e: React.PointerEvent) => {
    const d = dragging.current; dragging.current = null
    if (!d) return
    const sel = selections.find((x) => x.id === d.id); if (!sel) return
    const c = selectionCoords(node, sel.s, sel.e)
    setToast(`${sel.name} → bins [${c.chromStart}, ${c.chromEnd}) · ${c.nBins} bins · ${c.startMb.toFixed(2)}–${c.endMb.toFixed(2)} Mb (snapped to retained-bin boundaries)`)
  }
  useEffect(() => { const el = plotEl.current; if (!el) return; const h = () => setTick((t) => t + 1); el.on('plotly_afterplot', h); return () => { try { el.removeListener('plotly_afterplot', h) } catch { /* ignore */ } } }, [plotEl.current])

  const color = NODE_COLOR[nodeRef.side]
  return (
    <div className="card node-card" style={{ borderTopColor: color }}>
      <div className="card-head">
        <div>
          <div className="card-title" style={{ color }}>{p.title}</div>
          {p.subtitle && <div className="card-sub">{p.subtitle}</div>}
        </div>
        <div className="card-actions">
          {p.badge}
          {p.extraControls}
          <button className="btn-sm" title="Download this figure as PNG" onClick={() => plotEl.current && downloadPlot(plotEl.current, `${node.cell}_${node.chrom}_${nodeRef.side}_${node.start}-${node.end}`)}>PNG</button>
        </div>
      </div>
      <div className="plot-wrap" data-tick={tick}>
        <Plot data={fig.data} layout={fig.layout} config={config} onSelected={onSelected} onPlainClick={onPlainClick} onReady={(el) => { plotEl.current = el; setTick((t) => t + 1) }} style={{ width: '100%' }} />
        <div className="sel-overlay" onPointerMove={onHandleMove} onPointerUp={onHandleUp} onPointerCancel={onHandleUp}>
          {handleGeoms.filter((h) => !h.offscreen).map((h) => (
            <div key={h.s.id} className={`sel-handles ${h.s.id === activeSelectionId ? 'active' : ''}`} style={{ left: h.x0, width: Math.max(2, h.x1 - h.x0), top: h.top, height: h.height, ['--c' as any]: h.s.color }}>
              <div className="grip l" title={`${h.s.name}: drag to resize the left edge`} onPointerDown={(e) => onHandleDown(e, h.s.id, 'l')} />
              <div className="grip r" title={`${h.s.name}: drag to resize the right edge`} onPointerDown={(e) => onHandleDown(e, h.s.id, 'r')} />
              <div className="grip m" title={`${h.s.name}: drag to move`} onPointerDown={(e) => onHandleDown(e, h.s.id, 'm')}>{h.x1 - h.x0 > 40 ? h.s.name : ''}</div>
            </div>
          ))}
        </div>
      </div>
    </div>
  )
}
