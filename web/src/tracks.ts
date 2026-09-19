/** Builds the three-track Plotly figure (AD / log2FC / X) for one node. */
import { AD0, BP0, EVENT, LOCALBASE, NODE_COLOR, OVERLAP, SEGSTAT } from './colors'
import { boundaryMb, breakGaps, gapMask, MB } from './coords'
import { log2Track, mean, median, rollingMean } from './stats'
import type { Display } from './store'
import type { NodeData, RegionSeed, Selection } from './types'

export interface Marker { chromB: number; color: string; label: string; style: 'accepted' | 'candidate' | 'rejected' | 'manual' | 'child' | 'production'; width?: number }
export interface TrackOptions {
  node: NodeData; side: 'root' | 'left' | 'right'; title: string
  display: Display; genomeBaseline: { trimmed: number; median: number }
  selections: Selection[]; editableSelectionIds: Set<string>; activeSelectionId: string | null
  markers: Marker[]; regions: RegionSeed[]
  segmentBoundaries: number[]  // node-local, sorted, for mean/median overlays
  mode: 'navigate' | 'select' | 'place'
  presentation?: boolean
  adRange?: [number, number] | null
  height?: number
  compact?: boolean
  overlapRanges?: { s: number; e: number }[]  // node-local, hatched overlap of windows
}

export interface TrackFigure { data: any[]; layout: any; shapeSelectionIds: (string | null)[] }

export function referenceBaselines(node: NodeData, genome: { trimmed: number; median: number }, method: 'trimmed' | 'median') {
  const g = method === 'trimmed' ? genome.trimmed : genome.median
  const l = method === 'trimmed' ? (node.baseline_trimmed ?? 0) : (node.baseline_median ?? 0)
  const offset = g > 0 && l > 0 ? Math.log2(l / g) : null
  return { genome: g, node: l, offset }
}

export function buildTracks(o: TrackOptions): TrackFigure {
  const { node, display } = o
  const color = NODE_COLOR[o.side]
  const n = node.n
  const mbBins = node.start_bp.map((s) => s / MB)
  const gaps = gapMask(node.start_bp)
  const x0 = node.start_bp[0] / MB, x1 = node.end_bp[n - 1] / MB
  const data: any[] = []

  // ---- AD row (y) ---------------------------------------------------------
  const adX: number[] = [], adY: number[] = [], adCd: any[] = []
  for (let j = 1; j < n; j++) { adX.push(mbBins[j]); adY.push(node.ad[j - 1]); adCd.push([j, j, n - j, node.start + j]) }
  const adGaps = gaps.slice(1)
  const adB = breakGaps(adX, adY, adGaps)
  const adCdB: any[] = []
  { let k = 0; for (let i = 0; i < adX.length; i++) { if (adGaps[i] && i > 0) adCdB.push(null); adCdB.push(adCd[k++]) } }
  data.push({
    type: 'scatter', mode: 'lines', x: adB.x, y: adB.y, customdata: adCdB, name: 'AD', line: { color, width: 1.4 }, fill: 'tozeroy', fillcolor: hexA(color, 0.12),
    connectgaps: false, xaxis: 'x', yaxis: 'y', hovertemplate: '%{x:.2f} Mb · boundary %{customdata[1]} (chrom bin %{customdata[3]}) · AD %{y:.3f} · (%{customdata[1]}, %{customdata[2]})<extra>AD</extra>',
  })
  if (node.argmax_b !== null) {
    data.push({ type: 'scatter', mode: 'markers', x: [boundaryMb(node, node.argmax_b)], y: [node.argmax_ad], marker: { size: 9, color: o.side === 'root' ? BP0 : color, line: { color: 'white', width: 1.3 } }, name: 'argmax', hoverinfo: 'skip', xaxis: 'x', yaxis: 'y', showlegend: false })
  }

  // ---- log2FC row (y2) -----------------------------------------------------
  const refs = referenceBaselines(node, o.genomeBaseline, display.baseline)
  const roll = rollingMean(node.x, display.rolling, display.minPeriods)
  const tG = log2Track(roll, refs.genome, display.clip)
  const tL = log2Track(roll, refs.node, display.clip)
  const statusText = (t: ReturnType<typeof log2Track>, i: number) => t.status[i] === 'zero' ? ` (zero signal: −∞, drawn at −${display.clip})` : t.status[i] === 'nan' ? ' (insufficient bins)' : (t.values[i] !== null && Math.abs(t.values[i]!) >= display.clip ? ` (clipped at ±${display.clip})` : '')
  if (display.showGenome) {
    const b = breakGaps(mbBins, tG.values, gaps)
    const txt: (string | null)[] = []; { let k = 0; for (let i = 0; i < n; i++) { if (gaps[i] && i > 0) txt.push(null); txt.push(`bin ${k} · ${statusText(tG, k)}`); k++ } }
    data.push({ type: 'scatter', mode: 'lines', x: b.x, y: b.y, text: txt, name: 'log2FC vs genome-wide', line: { color, width: 1.2 }, connectgaps: false, xaxis: 'x', yaxis: 'y2', hovertemplate: '%{x:.2f} Mb · %{y:.3f} · %{text}<extra>vs genome</extra>' })
  }
  if (display.showLocal) {
    const b = breakGaps(mbBins, tL.values, gaps)
    const txt: (string | null)[] = []; { let k = 0; for (let i = 0; i < n; i++) { if (gaps[i] && i > 0) txt.push(null); txt.push(`bin ${k} · ${statusText(tL, k)}`); k++ } }
    data.push({ type: 'scatter', mode: 'lines', x: b.x, y: b.y, text: txt, name: 'log2FC vs this node', line: { color: LOCALBASE, width: 1.0, dash: 'dot' }, connectgaps: false, xaxis: 'x', yaxis: 'y2', hovertemplate: '%{x:.2f} Mb · %{y:.3f} · %{text}<extra>vs node</extra>' })
  }
  // segment mean / median overlays (computed from X, transformed with displayed reference)
  const segRef = refs.genome
  if ((display.showMean || display.showMedian) && segRef > 0) {
    const bs = o.segmentBoundaries
    for (let k = 0; k + 1 < bs.length; k++) {
      const a = bs[k], b = bs[k + 1]
      if (b <= a) continue
      const seg = node.x.slice(a, b)
      const m = mean(seg)!, md = median(seg)!
      const lm = m > 0 ? Math.log2(m / segRef) : null, lmd = md > 0 ? Math.log2(md / segRef) : null
      const overlap = lm !== null && lmd !== null && Math.abs(lm - lmd) < 0.04
      const c = overlap ? OVERLAP : SEGSTAT
      const xa = mbBins[a], xb = b < n ? mbBins[b] : x1
      if (display.showMean) data.push({ type: 'scatter', mode: 'lines', x: [xa, xb], y: [lm ?? -display.clip, lm ?? -display.clip], line: { color: c, width: 1.8 }, xaxis: 'x', yaxis: 'y2', showlegend: false, hovertemplate: `segment [${a},${b}) mean X ${m.toFixed(3)} → log2 ${lm === null ? '−∞ (mean 0)' : lm.toFixed(3)}${overlap ? ' · mean≈median (violet)' : ''}<extra>mean</extra>` })
      if (display.showMedian) data.push({ type: 'scatter', mode: 'lines', x: [xa, xb], y: [lmd ?? -display.clip, lmd ?? -display.clip], line: { color: c, width: 1.4, dash: 'dash' }, xaxis: 'x', yaxis: 'y2', showlegend: false, hovertemplate: `segment [${a},${b}) median X ${md.toFixed(3)} → log2 ${lmd === null ? '−∞ (median 0)' : lmd.toFixed(3)}${overlap ? ' · mean≈median (violet)' : ''}<extra>median</extra>` })
    }
  }

  // ---- X row (y3) ------------------------------------------------------------
  const rawB = breakGaps(mbBins, node.x, gaps)
  const rawCd: any[] = []; { let k = 0; for (let i = 0; i < n; i++) { if (gaps[i] && i > 0) rawCd.push(null); rawCd.push([k, node.start + k, node.offset_genome + node.start + k]); k++ } }
  data.push({ type: 'scatter', mode: 'lines', x: rawB.x, y: rawB.y, customdata: rawCd, name: 'X (per bin)', line: { color, width: 0.7 }, opacity: 0.45, connectgaps: false, xaxis: 'x', yaxis: 'y3', hovertemplate: '%{x:.2f} Mb · node bin %{customdata[0]} · chrom bin %{customdata[1]} · genome %{customdata[2]} · X %{y:.3f}<extra>X</extra>' })
  const rollB = breakGaps(mbBins, roll, gaps)
  data.push({ type: 'scatter', mode: 'lines', x: rollB.x, y: rollB.y, name: `rolling mean (${display.rolling})`, line: { color, width: 1.6 }, connectgaps: false, xaxis: 'x', yaxis: 'y3', hovertemplate: '%{x:.2f} Mb · rolling %{y:.3f}<extra>roll</extra>' })

  // ---- shapes: regions, markers, selections -----------------------------------
  const shapes: any[] = []
  const shapeSelectionIds: (string | null)[] = []
  const annotations: any[] = []
  const push = (s: any, id: string | null = null) => { shapes.push(s); shapeSelectionIds.push(id) }
  let regionCount = 0
  if (display.showRegions) {
    for (const r of o.regions) {
      if (r.chrom !== node.chrom) continue
      const active = r.id === display.activeRegionId
      push({ type: 'rect', xref: 'x', yref: 'paper', x0: r.start_mb, x1: r.end_mb, y0: 0, y1: 1, fillcolor: hexA(EVENT, active ? 0.12 : 0.06), line: { width: active ? 1 : 0, color: EVENT }, layer: 'below' })
      annotations.push({ xref: 'x', yref: 'paper', x: (r.start_mb + r.end_mb) / 2, y: 0.655 + 0.03 * (regionCount % 2), yanchor: 'bottom', text: r.label.replace(/ \(.*\)/, ''), showarrow: false, font: { size: 9, color: EVENT }, opacity: active ? 1 : 0.6 })
      regionCount++
    }
  }
  let markerCount = 0
  for (const m of o.markers) {
    const b = m.chromB - node.start
    if (b <= 0 || b >= n) continue
    const xm = boundaryMb(node, b)
    const dash = m.style === 'child' ? 'dash' : m.style === 'candidate' || m.style === 'rejected' ? 'dot' : 'solid'
    const width = m.width ?? (m.style === 'accepted' || m.style === 'manual' ? 1.8 : 1.3)
    push({ type: 'line', xref: 'x', yref: 'paper', x0: xm, x1: xm, y0: 0, y1: 1, line: { color: m.color, width, dash }, opacity: m.style === 'rejected' ? 0.55 : 0.9, layer: 'above' })
    annotations.push({ xref: 'x', yref: 'paper', x: xm, y: 0.995 - 0.045 * (markerCount % 4), yanchor: 'top', xanchor: 'left', text: ' ' + m.label, showarrow: false, font: { size: 9, color: m.color }, bgcolor: 'rgba(255,255,255,0.6)' })
    markerCount++
  }
  for (const s of o.selections) {
    if (!s.visible) continue
    const sb = s.nodeStart + s.s - node.start, eb = s.nodeStart + s.e - node.start
    if (eb <= 0 || sb >= n) continue
    const xa = boundaryMb(node, Math.max(0, sb)), xb = boundaryMb(node, Math.min(n, eb))
    const editable = o.editableSelectionIds.has(s.id)
    const active = s.id === o.activeSelectionId
    // rects are non-interactive in Plotly; editing is done with the NodeCard's own drag handles overlay
    push({ type: 'rect', xref: 'x', yref: 'paper', x0: xa, x1: xb, y0: 0, y1: 1, fillcolor: hexA(s.color, editable ? 0.2 : 0.14), line: { width: active ? 2 : 1, color: s.color, dash: editable ? 'solid' : 'dot' }, layer: 'below' }, s.id)
    if (s.anchor !== null) {
      const ab = s.nodeStart + s.anchor - node.start
      if (ab > 0 && ab < n) push({ type: 'line', xref: 'x', yref: 'paper', x0: boundaryMb(node, ab), x1: boundaryMb(node, ab), y0: 0, y1: 1, line: { color: s.color, width: 1, dash: 'dot' }, layer: 'above' })
    }
    annotations.push({ xref: 'x', yref: 'paper', x: xa, y: 0.02, yanchor: 'bottom', xanchor: 'left', text: s.name, showarrow: false, font: { size: 9, color: s.color }, bgcolor: 'rgba(255,255,255,0.7)' })
  }
  for (const ov of o.overlapRanges || []) {
    push({ type: 'rect', xref: 'x', yref: 'paper', x0: boundaryMb(node, ov.s), x1: boundaryMb(node, ov.e), y0: 0, y1: 1, fillcolor: 'rgba(0,0,0,0)', line: { width: 1.2, color: '#333', dash: 'dot' }, layer: 'above' })
  }
  // baselines on the X row and zero line on log2FC
  shapes.push({ type: 'line', xref: 'paper', yref: 'y2', x0: 0, x1: 1, y0: 0, y1: 0, line: { color: 'black', width: 0.5, dash: 'dot' }, opacity: 0.35 }); shapeSelectionIds.push(null)
  if (refs.genome > 0) { shapes.push({ type: 'line', xref: 'paper', yref: 'y3', x0: 0, x1: 1, y0: refs.genome, y1: refs.genome, line: { color: 'black', width: 0.8, dash: 'dot' }, opacity: 0.6 }); shapeSelectionIds.push(null) }
  if (refs.node > 0) { shapes.push({ type: 'line', xref: 'paper', yref: 'y3', x0: 0, x1: 1, y0: refs.node, y1: refs.node, line: { color: LOCALBASE, width: 1.1, dash: 'dash' }, opacity: 0.9 }); shapeSelectionIds.push(null) }

  const method = display.baseline === 'trimmed' ? 'trimmed-mean(IQR, lb=0)' : 'median (zeros kept)'
  annotations.push({ xref: 'paper', yref: 'paper', x: 1, y: 0.365, xanchor: 'right', yanchor: 'bottom', showarrow: false, font: { size: 9, color: LOCALBASE }, bgcolor: 'rgba(255,255,255,0.75)',
    text: `${method} · genome-wide ${refs.genome.toFixed(3)} · this node ${refs.node.toFixed(3)} · log2(node/genome) = ${refs.offset === null ? 'undefined' : (refs.offset >= 0 ? '+' : '') + refs.offset.toFixed(3)}` })
  if (tG.undefinedReason || tL.undefinedReason) annotations.push({ xref: 'paper', yref: 'paper', x: 0.5, y: 0.5, showarrow: false, font: { size: 10, color: EVENT }, text: `log2FC undefined: ${tG.undefinedReason || tL.undefinedReason}` })
  const xnf = node.x_nonfinite ? `⚠ ${node.x_nonfinite} non-finite X values` : ''
  if (xnf) annotations.push({ xref: 'paper', yref: 'paper', x: 0, y: 0.3, showarrow: false, font: { size: 9, color: EVENT }, text: xnf })

  const xp99 = percentile(node.x, 0.99)
  const layout: any = {
    height: o.height ?? 520,
    margin: { l: 58, r: 14, t: o.compact ? 30 : 36, b: 36 },
    title: { text: o.title, font: { size: 12 }, x: 0, xanchor: 'left' },
    paper_bgcolor: 'white', plot_bgcolor: 'white',
    font: { family: 'Inter, system-ui, sans-serif', size: 11, color: '#1f2937' },
    hovermode: o.presentation ? false : 'x', hoverdistance: 30, clickmode: 'event',
    dragmode: o.mode === 'navigate' ? 'zoom' : o.mode === 'select' ? 'select' : false,
    selectdirection: 'h',
    showlegend: false,
    xaxis: { range: [x0, x1], domain: [0, 1], anchor: 'y3', title: { text: `${node.chrom} position (Mb)`, standoff: 4 }, showspikes: !o.presentation, spikemode: 'across', spikethickness: 1, spikecolor: '#94a3b8', spikedash: 'dot', gridcolor: '#f1f5f9', zeroline: false },
    yaxis: { domain: [0.64, 1], title: { text: 'AD statistic', standoff: 4 }, gridcolor: '#f1f5f9', zeroline: false, rangemode: 'tozero', ...(o.adRange ? { range: o.adRange } : {}), fixedrange: false },
    yaxis2: { domain: [0.37, 0.61], title: { text: 'log2FC', standoff: 4 }, range: [-display.clip * 0.55, display.clip * 0.55], gridcolor: '#f1f5f9', zeroline: false },
    yaxis3: { domain: [0, 0.34], title: { text: 'X (GC-corr.)', standoff: 4 }, range: [0, Math.max(20, xp99 * 1.4)], gridcolor: '#f1f5f9', zeroline: false },
    shapes, annotations,
  }
  return { data, layout, shapeSelectionIds }
}

export function hexA(hex: string, a: number): string {
  const h = hex.replace('#', '')
  const r = parseInt(h.slice(0, 2), 16), g = parseInt(h.slice(2, 4), 16), b = parseInt(h.slice(4, 6), 16)
  return `rgba(${r},${g},${b},${a})`
}
export function percentile(v: number[], q: number): number {
  const s = [...v].filter(Number.isFinite).sort((a, b) => a - b)
  if (!s.length) return 0
  const i = Math.min(s.length - 1, Math.max(0, Math.floor(q * (s.length - 1))))
  return s[i]
}
export { AD0 }
