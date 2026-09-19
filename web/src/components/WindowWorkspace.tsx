import { useMemo, useRef, useState } from 'react'
import { LEFT, NODE_COLOR } from '../colors'
import { boundaryMb, selectionCoords } from '../coords'
import { acf, adSamples, ecdf, findCompetingPeaks, fmt, summarize, xSamples } from '../stats'
import { newId, nodeKey, useStore, type NodeRefLike } from '../store'
import { hexA } from '../tracks'
import type { NodeData, Selection, Snapshot, SnapshotWindow, SummaryStats } from '../types'
import Plot, { downloadPlot } from './Plot'

interface NodeEntry { ref: NodeRefLike; data: NodeData }
interface Props { nodes: NodeEntry[]; nodeLabel?: (ref: NodeRefLike) => string }

export interface WindowComputed {
  sel: Selection; node: NodeData; anchor: number; anchorKind: 'peak' | 'midpoint'
  ad: { j: number[]; values: number[] }; x: { i: number[]; values: number[] }
  adOff: number[]; xOff: number[]; adMb: number[]; xMb: number[]
  adStats: SummaryStats; xStats: SummaryStats; adAcf: ReturnType<typeof acf>; xAcf: ReturnType<typeof acf>
  coords: ReturnType<typeof selectionCoords>
  lr: { adL: SummaryStats; adR: SummaryStats; xL: SummaryStats; xR: SummaryStats; adAtAnchor: number | null }
}

export function computeWindow(sel: Selection, node: NodeData): WindowComputed {
  const ad = adSamples(node, sel.s, sel.e), x = xSamples(node, sel.s, sel.e)
  const anchorKind = sel.anchor !== null ? 'peak' : 'midpoint'
  const anchor = sel.anchor !== null ? sel.anchor : Math.floor((sel.s + sel.e) / 2)
  const adL = ad.values.filter((_, k) => ad.j[k] < anchor), adR = ad.values.filter((_, k) => ad.j[k] > anchor)
  const xL = x.values.filter((_, k) => x.i[k] < anchor), xR = x.values.filter((_, k) => x.i[k] >= anchor)
  return {
    sel, node, anchor, anchorKind, ad, x,
    adOff: ad.j.map((j) => j - anchor), xOff: x.i.map((i) => i - anchor + 0.5),
    adMb: ad.j.map((j) => boundaryMb(node, j)), xMb: x.i.map((i) => (node.start_bp[i] + node.end_bp[i]) / 2e6),
    adStats: summarize(ad.values), xStats: summarize(x.values), adAcf: acf(ad.values, 10), xAcf: acf(x.values, 10),
    coords: selectionCoords(node, sel.s, sel.e),
    lr: { adL: summarize(adL), adR: summarize(adR), xL: summarize(xL), xR: summarize(xR), adAtAnchor: anchor > 0 && anchor < node.n ? node.ad[anchor - 1] : null },
  }
}

function csvDownload(name: string, rows: (string | number | null)[][]) {
  const esc = (v: any) => (v === null || v === undefined ? 'N/A' : typeof v === 'string' && /[",\n]/.test(v) ? `"${v.replace(/"/g, '""')}"` : String(v))
  const blob = new Blob([rows.map((r) => r.map(esc).join(',')).join('\n')], { type: 'text/csv' })
  const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000)
}

const statRows = (kind: 'AD' | 'X', w: WindowComputed) => {
  const s = kind === 'AD' ? w.adStats : w.xStats
  return [w.sel.name, kind === 'AD' ? `${w.coords.adCount} AD positions${w.coords.adFirst !== null ? ` (j ${w.coords.adFirst}–${w.coords.adLast})` : ''}` : `${w.coords.nBins} data bins`, fmt(s.mean), fmt(s.median), fmt(s.std), s.cv === null ? `N/A (${s.cvReason})` : fmt(s.cv), s.lag1 === null ? `N/A (${s.lag1Reason})` : fmt(s.lag1)]
}

export default function WindowWorkspace({ nodes, nodeLabel }: Props) {
  const meta = useStore((s) => s.meta)!
  const selections = useStore((s) => s.selections)
  const allNodes = useStore((s) => s.nodes)
  const display = useStore((s) => s.display)
  const activeSelectionId = useStore((s) => s.activeSelectionId)
  const setActiveSelection = useStore((s) => s.setActiveSelection)
  const addSelection = useStore((s) => s.addSelection)
  const updateSelection = useStore((s) => s.updateSelection)
  const removeSelection = useStore((s) => s.removeSelection)
  const clearSelections = useStore((s) => s.clearSelections)
  const compareView = useStore((s) => s.compareView); const setCompareView = useStore((s) => s.setCompareView)
  const compareCoords = useStore((s) => s.compareCoords); const setCompareCoords = useStore((s) => s.setCompareCoords)
  const showLR = useStore((s) => s.showLeftRight); const setShowLR = useStore((s) => s.setShowLeftRight)
  const snapshots = useStore((s) => s.snapshots); const addSnapshot = useStore((s) => s.addSnapshot); const removeSnapshot = useStore((s) => s.removeSnapshot)
  const savedWindows = useStore((s) => s.savedWindows); const saveWindow = useStore((s) => s.saveWindow); const removeSavedWindow = useStore((s) => s.removeSavedWindow); const reapplyWindow = useStore((s) => s.reapplyWindow)
  const setToast = useStore((s) => s.setToast)
  const pushTiming = useStore((s) => s.pushTiming)
  const [peakNode, setPeakNode] = useState<string>('')
  const [nPeaks, setNPeaks] = useState(meta.defaults.n_competing_peaks)
  const [minSep, setMinSep] = useState(meta.defaults.min_peak_separation_mb)
  const [halfWin, setHalfWin] = useState(meta.defaults.half_window_bins)
  const [reapplyTarget, setReapplyTarget] = useState<string>('')
  const [onlyCurrent, setOnlyCurrent] = useState(true)
  const presentation = useStore((s) => s.presentation)
  const curCell = useStore((s) => s.cell); const curChrom = useStore((s) => s.chrom)
  const plotRefs = useRef<Record<string, any>>({})
  const label = nodeLabel ?? ((r: NodeRefLike) => r.side === 'root' ? 'depth 0' : `${r.side} child`)

  // windows relevant to the displayed nodes (plus any whose node is loaded)
  const computed = useMemo(() => {
    const t0 = performance.now()
    const out: WindowComputed[] = []
    for (const sel of selections) {
      if (onlyCurrent && (sel.cell !== curCell || sel.chrom !== curChrom)) continue
      const node = allNodes[nodeKey(sel.cell, sel.chrom, sel.nodeStart, sel.nodeEnd)]
      if (!node) continue
      out.push(computeWindow(sel, node))
    }
    const ms = performance.now() - t0
    if (out.length) setTimeout(() => pushTiming(`selection statistics (${out.length} windows)`, ms), 0)
    return out
  }, [selections, allNodes, onlyCurrent, curCell, curChrom])
  const visible = computed.filter((w) => w.sel.visible)
  const overlaps = useMemo(() => {
    const out: { a: WindowComputed; b: WindowComputed; s: number; e: number; n: number }[] = []
    for (let i = 0; i < computed.length; i++) for (let j = i + 1; j < computed.length; j++) {
      const a = computed[i], b = computed[j]
      if (a.sel.cell !== b.sel.cell || a.sel.chrom !== b.sel.chrom) continue
      const s = Math.max(a.sel.nodeStart + a.sel.s, b.sel.nodeStart + b.sel.s), e = Math.min(a.sel.nodeStart + a.sel.e, b.sel.nodeStart + b.sel.e)
      if (e > s) out.push({ a, b, s, e, n: e - s })
    }
    return out
  }, [computed])

  const currentNodeIds = new Set(nodes.map((n) => n.ref.id))
  const peakTarget = nodes.find((n) => n.ref.id === peakNode) ?? nodes.find((n) => n.ref.side === 'left') ?? nodes[0]

  const addPeakBatch = () => {
    if (!peakTarget) return
    const node = peakTarget.data
    if (node.argmax_b === null) return
    const region = display.activeRegionId ? meta.regions.find((r) => r.id === display.activeRegionId && r.chrom === node.chrom) : null
    const peaks = findCompetingPeaks(node, node.argmax_b, nPeaks, minSep, region ? { startMb: region.start_mb, endMb: region.end_mb } : null, (b) => boundaryMb(node, b))
    let added = 0
    for (const pk of peaks) {
      const w = { s: pk.b - halfWin, e: pk.b + halfWin }  // addSelection clips to the node and records the clipping
      const name = pk.isWinner ? `winner ${pk.mb.toFixed(1)} Mb` : pk.isEvent ? `${region?.label.split(' ')[0]} best ${pk.mb.toFixed(1)} Mb` : `peak ${pk.mb.toFixed(1)} Mb`
      const before = useStore.getState().selections.length
      addSelection({ cell: node.cell, chrom: node.chrom, nodeId: peakTarget.ref.id, nodeVersion: peakTarget.ref.version, nodeStart: node.start, nodeEnd: node.end, nodeSide: peakTarget.ref.side, s: w.s, e: w.e, anchor: pk.b, source: pk.isWinner ? 'winner' : pk.isEvent ? 'event' : 'peak', name, color: pk.isWinner ? NODE_COLOR[peakTarget.ref.side] : undefined, note: pk.relaxed ? `event-region candidate included although closer than ${minSep} Mb to another peak (spacing relaxed)` : pk.isEvent ? 'strongest candidate inside the active event region' : undefined })
      if (useStore.getState().selections.length > before) added++
    }
    setToast(`${added} window(s) added on ${label(peakTarget.ref)} (${peaks.length} candidates, ${peaks.length - added} already present)`)
  }
  const takeSnapshot = () => {
    if (!visible.length) return
    const wins: SnapshotWindow[] = visible.map((w) => ({ selection: { ...w.sel }, adOffsets: w.adOff, adValues: w.ad.values, xOffsets: w.xOff, xValues: w.x.values, adStats: w.adStats, xStats: w.xStats, coords: w.coords, adAcf: w.adAcf, xAcf: w.xAcf }))
    const first = visible[0]
    const snap: Snapshot = { id: newId('snap'), name: `snapshot ${snapshots.length + 1} · ${first.sel.cell} ${first.sel.chrom}`, createdAt: Date.now(), cell: first.sel.cell, chrom: first.sel.chrom, nodeLabel: [...new Set(visible.map((w) => w.sel.nodeSide))].join('+'), windows: wins, provenance: `frozen from ${visible.map((w) => `${w.sel.name}@${w.sel.nodeSide}[${w.sel.nodeStart},${w.sel.nodeEnd}) v${w.sel.nodeVersion}`).join('; ')}` }
    addSnapshot(snap); setToast(`Snapshot frozen with ${wins.length} window(s); it will not change when another cell is opened.`)
  }
  const exportStats = () => {
    const head = ['window', 'cell', 'chrom', 'node', 'node_bounds', 'node_local_s', 'node_local_e', 'chrom_bin_start', 'chrom_bin_end_excl', 'genome_bin_start', 'genome_bin_end_excl', 'start_bp', 'end_bp', 'span_bp', 'retained_bp', 'span_mb', 'retained_mb', 'coordinate_gaps', 'anchor_boundary', 'anchor_kind', 'source', 'AD_count', 'AD_mean', 'AD_median', 'AD_std', 'AD_cv', 'AD_lag1', 'X_count', 'X_mean', 'X_median', 'X_std', 'X_cv', 'X_lag1']
    const rows = computed.map((w) => [w.sel.name, w.sel.cell, w.sel.chrom, w.sel.nodeSide, `[${w.sel.nodeStart},${w.sel.nodeEnd})`, w.sel.s, w.sel.e, w.coords.chromStart, w.coords.chromEnd, w.coords.genomeStart, w.coords.genomeEnd, w.coords.startBp, w.coords.endBp, w.coords.spanBp, w.coords.retainedBp, w.coords.spanMb, w.coords.retainedMb, w.coords.gaps, w.anchor, w.anchorKind, w.sel.source, w.adStats.n, w.adStats.mean, w.adStats.median, w.adStats.std, w.adStats.cv, w.adStats.lag1, w.xStats.n, w.xStats.mean, w.xStats.median, w.xStats.std, w.xStats.cv, w.xStats.lag1])
    csvDownload('window_statistics.csv', [['# selection = node-local half-open [s,e) of retained bins; AD boundary j included when max(s,1)<=j<e (value d[j-1]); std ddof=0; CV=std/mean; lag1 = Pearson(v[:-1], v[1:]); descriptive only'], head, ...rows])
  }
  const exportSamples = () => {
    const rows: any[][] = [['window', 'kind', 'index', 'offset_from_anchor', 'mb', 'value']]
    for (const w of computed) { w.ad.j.forEach((j, k) => rows.push([w.sel.name, 'AD', j, w.adOff[k], w.adMb[k], w.ad.values[k]])); w.x.i.forEach((i, k) => rows.push([w.sel.name, 'X', i, w.xOff[k], w.xMb[k], w.x.values[k]])) }
    csvDownload('window_samples.csv', rows)
  }

  // ---- comparison figures ----------------------------------------------------
  const xAxisTitle = compareCoords === 'offset' ? 'offset from anchor boundary (retained-bin steps; 0 = split boundary)' : 'position (Mb)'
  const overlayFig = (kind: 'AD' | 'X') => {
    const data: any[] = []
    const layoutBase: any = { height: compareView === 'overlay' ? 300 : 120 * Math.max(1, visible.length) + 60, margin: { l: 50, r: 10, t: 30, b: 40 }, paper_bgcolor: 'white', plot_bgcolor: 'white', font: { size: 11 }, showlegend: compareView === 'overlay', legend: { orientation: 'h', y: -0.25, font: { size: 9 } }, hovermode: presentation ? false : 'closest', title: { text: kind === 'AD' ? 'AD samples (from the originating node curve)' : 'X samples (unsmoothed .X under the same bins)', font: { size: 12 }, x: 0 } }
    const shapes: any[] = []
    visible.forEach((w, k) => {
      const xs = kind === 'AD' ? (compareCoords === 'offset' ? w.adOff : w.adMb) : (compareCoords === 'offset' ? w.xOff : w.xMb)
      const ys = kind === 'AD' ? w.ad.values : w.x.values
      const ax = compareView === 'overlay' ? 'y' : `y${k + 1 === 1 ? '' : k + 1}`
      data.push({ type: 'scatter', mode: 'lines+markers', x: xs, y: ys, name: `${w.sel.name} (${w.sel.nodeSide})`, line: { color: w.sel.color, width: w.sel.source === 'winner' ? 2.6 : 1.4 }, marker: { size: w.sel.source === 'winner' ? 5 : 3.5 }, opacity: w.sel.source === 'winner' ? 1 : 0.85, yaxis: ax, hovertemplate: `${w.sel.name} · %{x} · ${kind} %{y:.3f}<extra></extra>` })
      if (compareView === 'small') {
        const dom0 = 1 - (k + 1) / visible.length, dom1 = 1 - k / visible.length - 0.02
        layoutBase[`yaxis${k + 1 === 1 ? '' : k + 1}`] = { domain: [dom0, dom1], title: { text: w.sel.name, font: { size: 9 } }, gridcolor: '#f1f5f9', rangemode: 'tozero' }
        if (compareCoords === 'absolute') shapes.push({ type: 'line', xref: 'x', yref: ax, x0: boundaryMb(w.node, w.anchor), x1: boundaryMb(w.node, w.anchor), y0: 0, y1: 1, line: { color: w.sel.color, dash: 'dot', width: 1 } })
      } else if (compareCoords === 'absolute') shapes.push({ type: 'line', xref: 'x', yref: 'paper', x0: boundaryMb(w.node, w.anchor), x1: boundaryMb(w.node, w.anchor), y0: 0, y1: 1, line: { color: w.sel.color, dash: 'dot', width: 1 } })
    })
    if (compareCoords === 'offset') shapes.push({ type: 'line', xref: 'x', yref: 'paper', x0: 0, x1: 0, y0: 0, y1: 1, line: { color: '#94a3b8', dash: 'dot', width: 1 } })
    const layout = { ...layoutBase, shapes, xaxis: { title: { text: xAxisTitle, font: { size: 10 } }, gridcolor: '#f1f5f9', zeroline: false } }
    if (compareView === 'overlay') layout.yaxis = { title: { text: kind === 'AD' ? 'AD statistic' : 'X (GC-corr.)' }, gridcolor: '#f1f5f9', rangemode: 'tozero' }
    return { data, layout }
  }
  const ecdfFig = (kind: 'AD' | 'X') => {
    const data: any[] = [], shapes: any[] = []
    visible.forEach((w) => {
      const v = kind === 'AD' ? w.ad.values : w.x.values, st = kind === 'AD' ? w.adStats : w.xStats
      const e = ecdf(v)
      data.push({ type: 'scatter', mode: 'lines', line: { shape: 'hv', color: w.sel.color, width: 1.6 }, x: e.x, y: e.y, name: w.sel.name, hovertemplate: `${w.sel.name} · value %{x:.3f} · F %{y:.2f}<extra></extra>` })
      if (st.mean !== null) shapes.push({ type: 'line', xref: 'x', yref: 'paper', x0: st.mean, x1: st.mean, y0: 0, y1: 1, line: { color: w.sel.color, width: 1.2 } })
      if (st.median !== null) shapes.push({ type: 'line', xref: 'x', yref: 'paper', x0: st.median, x1: st.median, y0: 0, y1: 1, line: { color: w.sel.color, width: 1.2, dash: 'dash' } })
    })
    return { data, layout: { hovermode: presentation ? false : 'closest', height: 240, margin: { l: 45, r: 10, t: 30, b: 36 }, paper_bgcolor: 'white', plot_bgcolor: 'white', font: { size: 11 }, showlegend: false, title: { text: `${kind} empirical CDF · solid = mean, dashed = median`, font: { size: 12 }, x: 0 }, shapes, xaxis: { title: { text: kind === 'AD' ? 'AD statistic' : 'X value', font: { size: 10 } }, gridcolor: '#f1f5f9' }, yaxis: { range: [0, 1.02], gridcolor: '#f1f5f9' } } }
  }
  const acfFig = (kind: 'AD' | 'X') => {
    const data: any[] = visible.map((w) => { const a = kind === 'AD' ? w.adAcf : w.xAcf; return { type: 'scatter', mode: 'lines+markers', x: a.lags, y: a.values, text: a.reasons.map((r) => r ? `N/A: ${r}` : ''), name: w.sel.name, line: { color: w.sel.color }, marker: { size: 5 }, connectgaps: false, hovertemplate: `${w.sel.name} · lag %{x} · r %{y:.3f} %{text}<extra></extra>` } })
    return { data, layout: { hovermode: presentation ? false : 'closest', height: 240, margin: { l: 45, r: 10, t: 30, b: 36 }, paper_bgcolor: 'white', plot_bgcolor: 'white', font: { size: 11 }, showlegend: false, title: { text: `${kind} autocorrelation vs lag (Pearson of v[:-k], v[k:]) · lag unit = retained-bin step`, font: { size: 12 }, x: 0 }, xaxis: { title: { text: 'lag (retained bins; gaps flagged in the table)', font: { size: 10 } }, dtick: 1, gridcolor: '#f1f5f9' }, yaxis: { range: [-1.05, 1.05], gridcolor: '#f1f5f9', zeroline: true, zerolinecolor: '#cbd5e1' } } }
  }
  const figs = useMemo(() => ({ ad: overlayFig('AD'), x: overlayFig('X'), eAd: ecdfFig('AD'), eX: ecdfFig('X'), aAd: acfFig('AD'), aX: acfFig('X') }), [visible, compareView, compareCoords, presentation])

  const boundsInput = (w: WindowComputed) => {
    const n = w.node.n
    const set = (s: number, e: number) => { s = Math.max(0, Math.min(n - 1, s)); e = Math.max(s + 1, Math.min(n, e)); updateSelection(w.sel.id, { s, e, clipped: undefined, anchor: w.sel.anchor !== null && (w.sel.anchor <= 0 || w.sel.anchor >= n) ? null : w.sel.anchor }) }
    return (
      <span className="bounds">
        <input type="number" value={w.sel.s} min={0} max={n - 1} onChange={(e) => set(+e.target.value, w.sel.e)} style={{ width: 62 }} title="node-local start bin (inclusive)" />
        <input type="number" value={w.sel.e} min={1} max={n} onChange={(e) => set(w.sel.s, +e.target.value)} style={{ width: 62 }} title="node-local end bin (exclusive)" />
        <button className="btn-xs" title="grow 1 bin left" onClick={() => set(w.sel.s - 1, w.sel.e)}>◀+</button>
        <button className="btn-xs" title="shrink 1 bin left" onClick={() => set(w.sel.s + 1, w.sel.e)}>◀−</button>
        <button className="btn-xs" title="shrink 1 bin right" onClick={() => set(w.sel.s, w.sel.e - 1)}>−▶</button>
        <button className="btn-xs" title="grow 1 bin right" onClick={() => set(w.sel.s, w.sel.e + 1)}>+▶</button>
      </span>
    )
  }

  return (
    <div className="card workspace">
      <div className="row wrap">
        <strong>Selected-window comparison</strong>
        <span className="muted">{computed.length} of {selections.length} window(s) · draw on any track in “Select window” mode, click an AD point for a {2 * halfWin}-bin peak window, or batch-add competing peaks</span>
        <span className="spacer" />
        <label>node <select value={peakTarget?.ref.id ?? ''} onChange={(e) => setPeakNode(e.target.value)}>{nodes.map((n) => <option key={n.ref.id} value={n.ref.id}>{label(n.ref)}</option>)}</select></label>
        <label>peaks <input type="number" min={1} max={12} value={nPeaks} onChange={(e) => setNPeaks(+e.target.value)} style={{ width: 46 }} /></label>
        <label>min sep <input type="number" min={0} step={0.5} value={minSep} onChange={(e) => setMinSep(+e.target.value)} style={{ width: 52 }} /> Mb</label>
        <label>half-width <input type="number" min={1} value={halfWin} onChange={(e) => setHalfWin(+e.target.value)} style={{ width: 46 }} /> bins</label>
        <button className="btn" onClick={addPeakBatch} disabled={!peakTarget}>Compare competing peaks</button>
        <button className="btn-sm" onClick={takeSnapshot} disabled={!visible.length}>Freeze snapshot</button>
        <button className="btn-sm" onClick={exportStats} disabled={!computed.length}>CSV stats</button>
        <button className="btn-sm" onClick={exportSamples} disabled={!computed.length}>CSV samples</button>
        <label className="chk-inline" title="windows from other cells/chromosomes are kept but hidden"><input type="checkbox" checked={onlyCurrent} onChange={(e) => setOnlyCurrent(e.target.checked)} />only this cell/chromosome</label>
        <button className="btn-sm" onClick={() => clearSelections()} disabled={!selections.length}>clear all</button>
      </div>
      {display.activeRegionId && <div className="muted small">Active event region {meta.regions.find((r) => r.id === display.activeRegionId)?.label}: the batch includes its strongest candidate (never replaced by a nearby different peak), relaxing spacing if needed.</div>}

      {computed.length > 0 && (
        <table className="tbl windows">
          <thead><tr><th></th><th>name</th><th>node</th><th>node-local [s,e)</th><th>chrom bins</th><th>Mb</th><th>width</th><th>span / retained</th><th>anchor</th><th>AD / X counts</th><th>edit bounds</th><th></th></tr></thead>
          <tbody>
            {computed.map((w) => (
              <tr key={w.sel.id} className={w.sel.id === activeSelectionId ? 'active' : ''} onClick={() => setActiveSelection(w.sel.id)}>
                <td><input type="color" value={w.sel.color} onChange={(e) => updateSelection(w.sel.id, { color: e.target.value })} title="colour" /> <input type="checkbox" checked={w.sel.visible} onChange={(e) => updateSelection(w.sel.id, { visible: e.target.checked })} title="visible" /></td>
                <td><input className="name" value={w.sel.name} onChange={(e) => updateSelection(w.sel.id, { name: e.target.value })} />{w.sel.note && <div className="muted small">{w.sel.note}</div>}{!currentNodeIds.has(w.sel.nodeId) && <div className="muted small">from another node/cell view</div>}</td>
                <td>{w.sel.nodeSide} [{w.sel.nodeStart},{w.sel.nodeEnd})<br /><span className="muted small">{w.sel.cell}</span></td>
                <td>[{w.sel.s}, {w.sel.e})</td>
                <td>[{w.coords.chromStart}, {w.coords.chromEnd})</td>
                <td>{w.coords.startMb.toFixed(2)}–{w.coords.endMb.toFixed(2)}</td>
                <td>{w.coords.nBins} bins{w.sel.clipped && <span className="warn small"> clipped {w.sel.clipped.left}/{w.sel.clipped.right}</span>}</td>
                <td>{w.coords.spanMb.toFixed(2)} / {w.coords.retainedMb.toFixed(2)} Mb{w.coords.gaps ? <span className="warn small"> · {w.coords.gaps} gap(s)</span> : ''}</td>
                <td>{w.anchorKind === 'peak' ? `b=${w.anchor} (${boundaryMb(w.node, w.anchor).toFixed(2)} Mb)` : <span className="muted">midpoint b={w.anchor}</span>}</td>
                <td>{w.coords.adCount} / {w.coords.nBins}</td>
                <td>{boundsInput(w)}</td>
                <td><button className="btn-xs" onClick={() => saveWindow(w.sel)} title="save as named window">save</button> <button className="btn-xs" onClick={() => removeSelection(w.sel.id)}>✕</button></td>
              </tr>
            ))}
          </tbody>
        </table>
      )}
      {overlaps.length > 0 && <div className="small">Overlaps: {overlaps.map((o, i) => <span key={i} className="pill">{o.a.sel.name} ∩ {o.b.sel.name}: chrom bins [{o.s}, {o.e}) = {o.n} retained bins</span>)}</div>}

      {visible.length > 0 && (
        <>
          <div className="row wrap">
            <span className="seg"><button className={compareView === 'overlay' ? 'on' : ''} onClick={() => setCompareView('overlay')}>overlay</button><button className={compareView === 'small' ? 'on' : ''} onClick={() => setCompareView('small')}>small multiples</button></span>
            <span className="seg"><button className={compareCoords === 'offset' ? 'on' : ''} onClick={() => setCompareCoords('offset')}>offset from anchor</button><button className={compareCoords === 'absolute' ? 'on' : ''} onClick={() => setCompareCoords('absolute')}>absolute Mb</button></span>
            <label className="chk-inline"><input type="checkbox" checked={showLR} onChange={(e) => setShowLR(e.target.checked)} />left/right-of-anchor summaries</label>
            <span className="muted small">AD offsets are integer boundary steps; X bins sit at half-integer offsets so 0 marks the split boundary. Windows without a peak anchor use their midpoint boundary (labelled). Unequal or clipped windows keep their real lengths.</span>
          </div>
          <div className="two-col">
            <div><Plot data={figs.ad.data} layout={figs.ad.layout} onReady={(el) => (plotRefs.current.ad = el)} /><button className="btn-xs" onClick={() => downloadPlot(plotRefs.current.ad, 'ad_comparison')}>PNG</button></div>
            <div><Plot data={figs.x.data} layout={figs.x.layout} onReady={(el) => (plotRefs.current.x = el)} /><button className="btn-xs" onClick={() => downloadPlot(plotRefs.current.x, 'x_comparison')}>PNG</button></div>
          </div>
          <div className="two-col">
            {(['AD', 'X'] as const).map((kind) => (
              <div key={kind}>
                <div className="tbl-title">{kind === 'AD' ? 'AD statistics — selected samples of the originating node’s AD curve (each value still compares all bins on either side of its boundary within the node)' : 'X statistics — corresponding unsmoothed .X values (the data underneath)'}</div>
                <table className="tbl stats">
                  <thead><tr><th>window</th><th>count</th><th>mean</th><th>median</th><th>std (ddof 0)</th><th>CV</th><th>lag-1 autocorr.</th></tr></thead>
                  <tbody>{visible.map((w) => { const r = statRows(kind, w); return <tr key={w.sel.id}><td><i className="sw" style={{ background: w.sel.color }} />{r[0]}</td>{r.slice(1).map((c, i) => <td key={i}>{c}</td>)}</tr> })}</tbody>
                </table>
              </div>
            ))}
          </div>
          {showLR && (
            <table className="tbl stats">
              <thead><tr><th>window</th><th>anchor</th><th>AD at anchor</th><th>AD left (j&lt;b): n / mean / CV</th><th>AD right (j&gt;b): n / mean / CV</th><th>X left [s,b): n / mean / CV</th><th>X right [b,e): n / mean / CV</th></tr></thead>
              <tbody>{visible.map((w) => <tr key={w.sel.id}><td>{w.sel.name}</td><td>{w.anchorKind === 'peak' ? `b=${w.anchor}` : `midpoint b=${w.anchor}`}</td><td>{fmt(w.lr.adAtAnchor)}</td><td>{w.lr.adL.n} / {fmt(w.lr.adL.mean)} / {fmt(w.lr.adL.cv)}</td><td>{w.lr.adR.n} / {fmt(w.lr.adR.mean)} / {fmt(w.lr.adR.cv)}</td><td>{w.lr.xL.n} / {fmt(w.lr.xL.mean)} / {fmt(w.lr.xL.cv)}</td><td>{w.lr.xR.n} / {fmt(w.lr.xR.mean)} / {fmt(w.lr.xR.cv)}</td></tr>)}</tbody>
            </table>
          )}
          <div className="two-col"><Plot data={figs.eAd.data} layout={figs.eAd.layout} /><Plot data={figs.eX.data} layout={figs.eX.layout} /></div>
          <div className="two-col"><Plot data={figs.aAd.data} layout={figs.aAd.layout} /><Plot data={figs.aX.data} layout={figs.aX.layout} /></div>
          <div className="muted small">Statistics are descriptive, not additional significance tests. Lag units are retained-bin steps; a window spanning coordinate gaps (flagged above) does not have a uniform 100 kb lag. Selecting or resizing a window never reruns node AD or permutation tests.</div>
        </>
      )}

      {(snapshots.length > 0 || savedWindows.length > 0) && (
        <div className="two-col">
          <div>
            <div className="tbl-title">Frozen snapshots (keep their original source; never change when a different cell is opened)</div>
            {snapshots.map((s) => (
              <details key={s.id} className="snap">
                <summary><b>{s.name}</b> · {new Date(s.createdAt).toLocaleTimeString()} · {s.windows.length} window(s) · {s.nodeLabel} <button className="btn-xs" onClick={(e) => { e.preventDefault(); removeSnapshot(s.id) }}>✕</button></summary>
                <div className="muted small">{s.provenance}</div>
                <table className="tbl stats"><thead><tr><th>window</th><th>node</th><th>bins</th><th>AD n/mean/median/std/CV/lag1</th><th>X n/mean/median/std/CV/lag1</th></tr></thead>
                  <tbody>{s.windows.map((w) => <tr key={w.selection.id}><td><i className="sw" style={{ background: w.selection.color }} />{w.selection.name}</td><td>{w.selection.cell} {w.selection.chrom} {w.selection.nodeSide} [{w.selection.nodeStart},{w.selection.nodeEnd}) v{w.selection.nodeVersion}</td><td>[{w.coords.chromStart},{w.coords.chromEnd})</td><td>{w.adStats.n} / {fmt(w.adStats.mean)} / {fmt(w.adStats.median)} / {fmt(w.adStats.std)} / {fmt(w.adStats.cv)} / {fmt(w.adStats.lag1)}</td><td>{w.xStats.n} / {fmt(w.xStats.mean)} / {fmt(w.xStats.median)} / {fmt(w.xStats.std)} / {fmt(w.xStats.cv)} / {fmt(w.xStats.lag1)}</td></tr>)}</tbody></table>
                <button className="btn-xs" onClick={() => csvDownload(`${s.name.replace(/[^a-z0-9]+/gi, '_')}.csv`, [['# ' + s.provenance], ['window', 'cell', 'chrom', 'node', 'chrom_bin_start', 'chrom_bin_end_excl', 'AD_n', 'AD_mean', 'AD_median', 'AD_std', 'AD_cv', 'AD_lag1', 'X_n', 'X_mean', 'X_median', 'X_std', 'X_cv', 'X_lag1'], ...s.windows.map((w) => [w.selection.name, w.selection.cell, w.selection.chrom, w.selection.nodeSide, w.coords.chromStart, w.coords.chromEnd, w.adStats.n, w.adStats.mean, w.adStats.median, w.adStats.std, w.adStats.cv, w.adStats.lag1, w.xStats.n, w.xStats.mean, w.xStats.median, w.xStats.std, w.xStats.cv, w.xStats.lag1])])}>CSV</button>
              </details>
            ))}
          </div>
          <div>
            <div className="tbl-title">Saved windows (reapplying a genomic window to another node is explicit; clipping is reported)</div>
            <div className="row wrap"><label>target node <select value={reapplyTarget} onChange={(e) => setReapplyTarget(e.target.value)}><option value="">choose…</option>{nodes.map((n) => <option key={n.ref.id} value={n.ref.id}>{label(n.ref)} [{n.data.start},{n.data.end})</option>)}</select></label></div>
            {savedWindows.map((s) => (
              <div key={s.id} className="row saved">
                <i className="sw" style={{ background: s.color }} /><b>{s.name}</b><span className="muted small">{s.cell} {s.chrom} {s.nodeSide} · chrom bins [{s.nodeStart + s.s},{s.nodeStart + s.e})</span>
                <span className="spacer" />
                <button className="btn-xs" disabled={!reapplyTarget} onClick={() => { const t = nodes.find((n) => n.ref.id === reapplyTarget); if (!t) return; const r = reapplyWindow(s, t.data, t.ref); setToast(r.added ? `Reapplied ${s.name} to ${label(t.ref)}: ${r.report}` : `Could not reapply: ${r.report}`) }}>reapply</button>
                <button className="btn-xs" onClick={() => removeSavedWindow(s.id)}>✕</button>
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  )
}
export { LEFT as _unusedLeft, hexA as _unusedHexA }
