import AlternativeFlankView from './AlternativeFlankView'
import { useEffect, useMemo, useState } from 'react'
import Plot, { downloadPlot } from './Plot'
import ClusterKaryogram, { type KaryoRow } from './ClusterKaryogram'
import { CellCard } from './IntegrationClustersView'
import { useIntegration } from '../integrationStore'
import { useStore } from '../store'
import { DEPTH_COLORS, METRICS, NOISE_CV_BINS, SOURCE_COLORS, STORY_CN, STORY_CN_LABELS, VIRIDIS, binColor, binIndex, clusterColor, clusterLevels, clusterOf, clusterSegmentsLocal, consistencyClass, hatchShapes, levelTrack, metricValue, pickFlank, rollingMean, segmentsFromBreakpoints, viridisAt, withGaps, type Metric, type SmallSeg } from '../integration'
import { CN_COLORS } from '../colors'
import { MB } from '../coords'

const fmt = (v: number | null | undefined, d = 3) => (v === null || v === undefined || !Number.isFinite(v) ? 'N/A' : v.toFixed(d))
const PROPOSAL_LABEL: Record<string, string> = { merge_left: 'merge ←', merge_right: 'merge →', ambiguous_keep: 'ambiguous', keep_focal: 'keep focal' }
const PROPOSAL_SHORT: Record<string, string> = { merge_left: 'merge←', merge_right: 'merge→', ambiguous_keep: 'ambig.', keep_focal: 'keep' }
const REASON_LABEL: Record<string, string> = { effect_ge_tau: 'effect ≥ τ', terminal_unchanged: 'terminal', long_long: 'long–long' }
const REASON_SHORT: Record<string, string> = { effect_ge_tau: '≥τ', terminal_unchanged: 'term', long_long: 'LL' }

export default function IntegrationSegmentationView() {
  const ig = useIntegration()
  const meta = ig.meta
  const exMeta = useStore((s) => s.meta)!
  const bins = useStore((s) => s.bins)!
  const presentation = useStore((s) => s.presentation)
  const setToast = useStore((s) => s.setToast)
  const [plotEl, setPlotEl] = useState<any>(null)
  const [chromTable, setChromTable] = useState(true)
  useEffect(() => { ig.load() }, [])
  useEffect(() => { if (meta) { ig.ensurePseudobulk(); ig.ensureKaryo() } }, [meta])
  const seg = ig.seg
  useEffect(() => { if (seg.labels.candidates && meta) for (const s of seg.sources) ig.ensureCandidates(s, seg.chrom) }, [seg.labels.candidates, seg.chrom, seg.sources, meta])
  const nBins = exMeta.dataset.identity.n_bins
  const chromMeta = exMeta.chromosomes.find((c) => c.name === seg.chrom)!
  const cells = meta?.cells || []
  const levels = useMemo(() => (meta ? clusterLevels(cells, 'wnn_leiden_0.3') : []), [meta])
  const clSegs = useMemo(() => clusterSegmentsLocal(meta?.cluster_segments ?? null, exMeta.chromosomes), [meta])
  const metricDef = METRICS.find((m) => m.id === seg.metric)!
  const effectMetric: Metric = metricDef.kind === 'effect' ? seg.metric : seg.effectMetric

  // ---------------- cluster-level karyogram rows (fig 3c)
  const clusterRows: KaryoRow[] = useMemo(() => {
    if (!meta || !meta.cluster_segments) return []
    return meta.sources.map((src) => {
      const a = new Uint8Array(nBins).fill(255)
      for (const g of meta.cluster_segments!) if (g.group === src) a.fill(Math.round(g.cn_call * 2), g.bin_start, g.bin_end + 1)
      const bs = meta.cluster_segments!.find((g) => g.group === src)?.best_s
      return { label: `${src} (n=${meta.source_sizes[src]})`, sub: `best_s ${bs?.toFixed(4)} · ${meta.cluster_segments!.filter((g) => g.group === src).length} segments`, data: a, band: SOURCE_COLORS[src] }
    })
  }, [meta, nBins])
  // ---------------- per-cell karyogram rows (fig 3d)
  const cellOrder = useMemo(() => {
    const idx = cells.map((_, i) => i)
    if (seg.percellOrder === 'coverage') return idx.sort((a, b) => cells[b].total_raw_depth - cells[a].total_raw_depth)
    if (seg.percellOrder === 'chr8') return idx.sort((a, b) => ((cells[b].chr8_score ?? -9) - (cells[a].chr8_score ?? -9)))
    return idx.sort((a, b) => (+clusterOf(cells[a], 'wnn_leiden_0.3') - +clusterOf(cells[b], 'wnn_leiden_0.3')) || a - b)
  }, [cells, seg.percellOrder])
  const cellRows: KaryoRow[] = useMemo(() => {
    if (!ig.karyo) return []
    return cellOrder.map((i) => ({ label: cells[i].cell, sub: `cluster ${clusterOf(cells[i], 'wnn_leiden_0.3')} · best_s ${fmt(cells[i].percell_best_s, 3)}`, data: ig.karyo!.subarray(i * nBins, (i + 1) * nBins), band: clusterColor(clusterOf(cells[i], 'wnn_leiden_0.3'), levels) }))
  }, [ig.karyo, cellOrder, cells, levels, nBins])
  const hlSet = useMemo(() => new Set(ig.highlight), [ig.highlight])
  const selectedCellRow = useMemo(() => (ig.focusCell ? cellOrder.findIndex((i) => cells[i].cell === ig.focusCell) : -1), [ig.focusCell, cellOrder, cells])

  // ---------------- small segments on the chromosome
  const smallHere = useMemo(() => ig.smallSegs.filter((s) => s.chrom === seg.chrom && s.len <= seg.smallMax), [ig.smallSegs, seg.chrom, seg.smallMax])
  const pilot = meta?.pilot_chromosomes?.includes(seg.chrom)

  // ---------------- tracks figure
  const trackFig = useMemo(() => {
    if (!meta || !ig.pb || !meta.breakpoints) return null
    const srcs = seg.sources.filter((s) => meta.sources.includes(s))
    const nS = srcs.length; if (!nS) return null
    const data: any[] = []; const shapes: any[] = []; const annotations: any[] = []; const layout: any = {}
    const mb: number[] = []; for (let i = 0; i < chromMeta.n; i++) mb.push(bins[2 * (chromMeta.offset + i)] / MB)
    const endsMb = Array.from({length:chromMeta.n},(_,i)=>bins[2*(chromMeta.offset+i)+1]/MB)
    const rowH = 1 / nS, gap = 0.06 / nS
    const yr: [number, number] = seg.basis === 'raw' ? [0, 0] : [-4.6, 3.2]
    srcs.forEach((src, k) => {
      const si = meta.sources.indexOf(src)
      const pb = ig.pb!.subarray(si * nBins, (si + 1) * nBins)
      const lt = levelTrack(pb, chromMeta, seg.basis, pb)
      const y = rollingMean(lt.y, seg.smoothing)
      const g = withGaps(mb, y, endsMb)
      const ax = k === 0 ? '' : String(k + 1)
      const yref = `y${ax}`
      const haploAxis = String(nS + k + 1)
      const dom: [number, number] = [1 - (k + 1) * rowH + gap, 1 - k * rowH]
      data.push({ type: 'scatter', mode: 'lines', x: g.x, y: g.y, name: src, line: { color: '#6b7280', width: 0.8 }, connectgaps: false, xaxis: 'x', yaxis: yref, hovertemplate: `${src} · %{x:.2f} Mb · %{y:.2f}<extra></extra>`, showlegend: false })
      if (seg.labels.allelic) {
        const hp = ig.pb!.subarray((6 + si) * nBins, (7 + si) * nBins); const hf = ig.pb!.subarray((12 + si) * nBins, (13 + si) * nBins)
        const hy = Array.from(hp.subarray(chromMeta.offset, chromMeta.offset + chromMeta.n)).map((v) => v)
        const hg = withGaps(mb, rollingMean(hy, Math.max(seg.smoothing, 5)), endsMb)
        const txt = Array.from(hf.subarray(chromMeta.offset, chromMeta.offset + chromMeta.n)).map((f) => `nonzero imbalance in ${(100 * f).toFixed(0)}% of cells`)
        data.push({ type: 'scatter', mode: 'lines', x: hg.x, y: hg.y, text: withGaps(mb, txt as any, endsMb).y, name: `${src} allelic`, line: { color: '#b5651d', width: 1 }, connectgaps: false, xaxis: 'x', yaxis: `y${haploAxis}`, hovertemplate: `${src} mean |dBAF| %{y:.3f} (%{text}; zero = no coverage or balanced)<extra></extra>`, showlegend: false })
        layout[`yaxis${haploAxis}`] = { overlaying: yref, side: 'right', range: [0, 0.5], showgrid: false, tickfont: { size: 8, color: '#b5651d' }, title: { text: '|dBAF|', font: { size: 8, color: '#b5651d' } } }
      }
      layout[`yaxis${ax}`] = { domain: dom, title: { text: `${src}<br><span style="font-size:8px;color:#60717c">${lt.refLabel.replace('chromosome median', 'chr median').replace('genome trimmed mean', 'genome trim. mean').replace('pseudobulk X (mean over member cells)', 'mean X')}</span>`, font: { size: 11, color: SOURCE_COLORS[src] } }, gridcolor: '#f1f5f9', zeroline: false, ...(seg.basis === 'raw' ? { rangemode: 'tozero' } : { range: yr }), tickfont: { size: 9 } }
      if (seg.basis !== 'raw') shapes.push({ type: 'line', xref: 'x', yref, x0: mb[0], x1: mb[mb.length - 1], y0: 0, y1: 0, line: { color: '#000', width: 0.5, dash: 'dot' }, opacity: 0.5 })
      // CN strip from the cluster-level calls
      if (seg.labels.cn) for (const cs of clSegs.filter((c) => c.group === src && c.chromosome === seg.chrom)) shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: mb[cs.s], x1: cs.e < chromMeta.n ? mb[cs.e] : bins[2 * (chromMeta.offset + chromMeta.n - 1) + 1] / MB, y0: dom[1] - 0.045 * rowH, y1: dom[1], fillcolor: (seg.karyoPalette === 'story' ? STORY_CN : CN_COLORS)[cs.cn_call], line: { width: 0 }, layer: 'below' })
      // accepted breakpoints, coloured by recursion depth
      if (seg.labels.breakpoints) for (const bp of meta.breakpoints!.filter((b) => b.source === src && b.chromosome === seg.chrom)) {
        const x = mb[Math.min(bp.absolute_bin, chromMeta.n - 1)]
        shapes.push({ type: 'line', xref: 'x', yref: `${yref} domain`, x0: x, x1: x, y0: 0, y1: 1, line: { color: seg.labels.depthColor ? DEPTH_COLORS[Math.min(bp.recursion_depth, DEPTH_COLORS.length - 1)] : '#440154', width: 1.2 }, opacity: 0.9 })
      }
      // small segments
      if (seg.labels.small) for (const s of smallHere.filter((s) => s.source === src)) {
        const mv = metricValue(s, seg.metric, seg.flank, effectMetric)
        const cc = consistencyClass(s, effectMetric, seg.srdRef)
        // SRD-screen proposal layer: colour by the min-SRD flank's |SRD/√φ| (one of the segment's two flank values) instead of the selected "colour by" metric
        const srdPhiMetric = METRICS.find((m) => m.id === 'srd_phi')!
        const srdPhiMin = metricValue(s, 'srd_phi', 'min_srd', effectMetric)
        const proposalActive = seg.labels.proposal && !!s.proposal
        const color = proposalActive
          ? (seg.scale === 'binned' ? binColor(srdPhiMin.value, srdPhiMetric.bins) : (srdPhiMin.value === null ? '#e5e7eb' : viridisAt(Math.min(1, Math.abs(srdPhiMin.value) / 8))))
          : (seg.scale === 'binned' ? binColor(mv.value, metricDef.bins) : (mv.value === null ? '#e5e7eb' : viridisAt(Math.min(1, Math.abs(mv.value) / (metricDef.kind === 'stat' ? 8 : 2)))))
        const selected = seg.selectedSegment === s.key
        // hatch: SRD-screen proposal type takes priority when its layer is on, else flank-consistency class
        const hatch: 'none' | '/' | 'x' = proposalActive
          ? (s.proposal === 'ambiguous_keep' ? 'x' : s.proposal === 'keep_focal' ? 'none' : '/')
          : (seg.labels.consistency ? (cc.cls === 'weak' ? 'x' : cc.cls === 'inconsistent' ? '/' : 'none') : 'none')
        const x0 = mb[s.s], x1 = s.e < chromMeta.n ? mb[s.e] : bins[2 * (chromMeta.offset + chromMeta.n - 1) + 1] / MB
        shapes.push({ type: 'rect', xref: 'x', yref: `${yref} domain`, x0, x1, y0: 0, y1: 1, fillcolor: color, opacity: selected ? 0.75 : 0.5, line: { width: selected ? 2.5 : 0.6, color: '#111827' }, layer: 'below' })
        shapes.push(...hatchShapes(x0, x1, yref, hatch))
        const labels: string[] = []
        if (seg.labels.proposal && s.proposal) labels.push(PROPOSAL_SHORT[s.proposal] || s.proposal)
        if (seg.labels.consistency) labels.push(cc.cls === 'weak' ? 'x' : cc.cls === 'strong' ? 'o' : cc.cls === 'inconsistent' ? '/' : '')
        if (seg.labels.reason) labels.push(`${REASON_SHORT[s.reasonL] || s.reasonL}|${REASON_SHORT[s.reasonR] || s.reasonR}`)
        if (seg.labels.segLabels) labels.push(`${s.len}b ${mv.value === null ? 'N/A' : mv.value.toFixed(2)}`)
        // vertical, inside the span: adjacent segments then never overlap horizontally
        if (labels.length) annotations.push({ xref: 'x', yref, x: (x0 + x1) / 2, y: seg.basis === 'raw' ? 0 : yr[1] - 0.15, yanchor: 'top', xanchor: 'center', textangle: -90, showarrow: false, font: { size: 9, color: '#111827' }, bgcolor: 'rgba(255,255,255,0.7)', text: labels.filter(Boolean).join(' · ') })
        // hover/click carrier
        data.push({ type: 'scatter', mode: 'markers', x: [(x0 + x1) / 2], y: [seg.basis === 'raw' ? 0 : yr[1] - 0.35], marker: { size: 9, color, line: { color: '#111827', width: 1 }, symbol: 'square' }, customdata: [s.key], xaxis: 'x', yaxis: yref, showlegend: false, name: 'small segment',
          hovertemplate: `<b>${src} ${seg.chrom} segment ${s.segId}</b><br>bins [${s.s}, ${s.e}) · ${s.len} bins · ${s.mbStart}–${s.mbEnd} Mb<br>${metricDef.label}: ${mv.value === null ? 'N/A' : mv.value.toFixed(3)} (${mv.flank} flank, ${seg.flank})<br>Δ median L/R ${fmt(s.left.delta_median, 2)} / ${fmt(s.right.delta_median, 2)} · Δ mean ${fmt(s.left.delta_mean, 2)} / ${fmt(s.right.delta_mean, 2)} · Δ trim ${fmt(s.left.delta_trim, 2)} / ${fmt(s.right.delta_trim, 2)}<br>SRD/√φ L/R ${fmt(s.left.srd_phi, 2)} / ${fmt(s.right.srd_phi, 2)} · SRD ${fmt(s.left.srd, 2)} / ${fmt(s.right.srd, 2)} · log2FC ${fmt(s.left.log2fc, 2)} / ${fmt(s.right.log2fc, 2)}<br>level median ${s.stats.median.toFixed(2)} mean ${s.stats.mean.toFixed(2)} trim ${s.stats.trim.toFixed(2)} · std ${s.stats.std.toFixed(2)} · CV ${s.stats.cv.toFixed(2)}<br>${cc.label}<br>reasons: ${s.reasonL} | ${s.reasonR}${s.proposal ? ` · SRD-screen proposal: ${s.proposal}` : ''}<extra></extra>` })
      }
      // candidate audit
      if (seg.labels.candidates) {
        const rows = ig.candidates[`${src}|${seg.chrom}`] || []
        const rej = rows.filter((r) => !r.accepted && r.selected_for_significance), acc = rows.filter((r) => r.accepted)
        const geo = rows.filter((r) => !r.selected_for_significance && r.candidate_rank <= 1)
        const mk = (rows0: typeof rows, name: string, sym: string, color: string) => { const rs = rows0.filter((r) => r.candidate_absolute_bin !== null && r.candidate_absolute_bin !== undefined); return rs.length && data.push({ type: 'scatter', mode: 'markers', x: rs.map((r) => mb[Math.max(0, Math.min(r.candidate_absolute_bin, chromMeta.n - 1))]), y: rs.map(() => (seg.basis === 'raw' ? 0 : yr[0] + 0.3)), marker: { symbol: sym, size: 7, color }, customdata: rs.map((r) => `${r.rejection_reasons || ''}`), xaxis: 'x', yaxis: yref, showlegend: false, name, hovertemplate: rs.map((r) => `${name} · depth ${r.recursion_depth} · node [${r.node_start_bin},${r.node_end_bin}) · rank ${r.candidate_rank} · AD ${fmt(r.ad_score, 2)} · children (${r.left_child_bins}, ${r.right_child_bins})<br>p_local ${r.p_local ?? 'n/a'} · p_global ${r.p_global ?? 'n/a'}<br>${r.rejection_reasons || 'no rejection reason'}<extra></extra>`) }) }
        mk(acc, 'accepted candidate', 'triangle-up', '#440154'); mk(rej, 'tested, rejected', 'x', '#c1121f'); mk(geo, 'rank-1 candidate blocked by geometry', 'diamond-open', '#f59e0b')
      }
    })
    Object.assign(layout, { height: Math.max(260, 150 * nS + 70), margin: { l: 56, r: seg.labels.allelic ? 46 : 14, t: 30, b: 40 }, paper_bgcolor: 'white', plot_bgcolor: 'white', font: { size: 11 }, hovermode: presentation ? false : 'closest', dragmode: 'zoom', shapes, annotations, showlegend: false, xaxis: { title: { text: `${seg.chrom} position (Mb)` }, range: [mb[0], bins[2 * (chromMeta.offset + chromMeta.n - 1) + 1] / MB], gridcolor: '#f1f5f9', anchor: `y${nS === 1 ? '' : nS}` } })
    layout.title = { text: `${seg.chrom} · pseudobulk depth per WNN cluster${seg.basis === 'chrmedian' ? ' · log2(depth / chromosome median)' : seg.basis === 'genome_trimmed' ? ' · log2(depth / genome trimmed mean)' : ' · mean X'} · small segments (≤ ${seg.smallMax} bins) coloured by ${metricDef.label} (${seg.flank.replace('_', ' ')} flank)`, font: { size: 12 }, x: 0 }
    return { data, layout }
  }, [meta, ig.pb, seg, smallHere, chromMeta, bins, effectMetric, metricDef, clSegs, ig.candidates, presentation, nBins])

  const onTrackClick = (ev: any) => { const p = ev?.points?.[0]; if (p && p.customdata && typeof p.customdata === 'string' && p.customdata.includes('|')) ig.setSeg({ selectedSegment: seg.selectedSegment === p.customdata ? null : p.customdata }) }
  const selectSeg = (s: SmallSeg) => { ig.setSeg({ selectedSegment: s.key, chrom: s.chrom, sources: seg.sources.includes(s.source) ? seg.sources : [...seg.sources, s.source] }); document.getElementById('cluster-tracks')?.scrollIntoView({ behavior: 'smooth', block: 'start' }) }

  // ---------------- diagnostics across all pilot chromosomes
  const allSmall = useMemo(() => ig.smallSegs.filter((s) => s.len <= seg.smallMax), [ig.smallSegs, seg.smallMax])
  const diagCommon = { paper_bgcolor: 'white', plot_bgcolor: 'white', font: { size: 11 }, hovermode: presentation ? false : 'closest', margin: { l: 55, r: 10, t: 36, b: 44 } }
  const volcano = () => {
    const pts = allSmall.map((s) => { const f = pickFlank(s, 'min_delta', effectMetric); if (!f) return null; const d = s[f][effectMetric], sp = s[f].srd_phi; if (d === null || sp === null) return null; return { s, x: Math.abs(d), y: Math.abs(sp), f } }).filter(Boolean) as { s: SmallSeg; x: number; y: number; f: 'left' | 'right' }[]
    const q = { smallStrong: pts.filter((p) => p.x < seg.tau && p.y >= seg.srdRef).length, bigWeak: pts.filter((p) => p.x >= seg.tau && p.y < seg.srdRef).length }
    const data = (['left', 'right'] as const).map((f) => { const sub = pts.filter((p) => p.f === f); return { type: 'scatter', mode: 'markers', name: `weaker flank = ${f}`, x: sub.map((p) => p.x), y: sub.map((p) => p.y), customdata: sub.map((p) => p.s.key), text: sub.map((p) => `${p.s.source} ${p.s.chrom} seg ${p.s.segId} (${p.s.len} bins)`), marker: { size: 7, color: f === 'left' ? '#1f77b4' : '#ff7f0e', opacity: 0.75, line: { width: sub.map((p) => (p.s.key === seg.selectedSegment ? 2 : 0)), color: '#111' } }, hovertemplate: '%{text}<br>|effect| %{x:.3f} · |SRD/√φ| %{y:.2f}<extra></extra>' } })
    return { data, layout: { ...diagCommon, height: 360, title: { text: `weaker-flank |effect| (${METRICS.find((m) => m.id === effectMetric)!.label}) vs |SRD/√φ| of that same flank · ${q.smallStrong} in "small effect, strong SRD", ${q.bigWeak} with effect ≥ τ but SRD < ${seg.srdRef}`, font: { size: 12 }, x: 0 }, xaxis: { title: { text: '|log2 effect| of weaker flank', font: { size: 11 } } }, yaxis: { title: { text: '|SRD/√φ| of that flank', font: { size: 11 } } }, shapes: [{ type: 'line', x0: seg.tau, x1: seg.tau, y0: 0, y1: 1, yref: 'paper', line: { dash: 'dash', color: '#333', width: 1 } }, { type: 'line', x0: 0, x1: 1, xref: 'paper', y0: seg.srdRef, y1: seg.srdRef, line: { dash: 'dash', color: '#333', width: 1 } }], legend: { orientation: 'h', y: -0.2 } } }
  }
  const sizeDep = (stat: 'srd' | 'srd_phi') => {
    const pts = allSmall.map((s) => { const a = s.left[stat], b = s.right[stat]; const v = Math.max(a === null ? -1 : Math.abs(a), b === null ? -1 : Math.abs(b)); return v < 0 ? null : { s, x: s.len, y: v } }).filter(Boolean) as { s: SmallSeg; x: number; y: number }[]
    const rho = pts.length > 3 ? (() => { const rank = (v: number[]) => { const idx = v.map((x, i) => [x, i] as [number, number]).sort((p, q) => p[0] - q[0]); const r = new Array(v.length); let i = 0; while (i < idx.length) { let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++; for (let k = i; k <= j; k++) r[idx[k][1]] = (i + j) / 2; i = j + 1 } return r }; const ra = rank(pts.map((p) => p.x)), rb = rank(pts.map((p) => p.y)); const n = ra.length, ma = ra.reduce((a, b) => a + b, 0) / n, mb2 = rb.reduce((a, b) => a + b, 0) / n; let sab = 0, saa = 0, sbb = 0; for (let i = 0; i < n; i++) { sab += (ra[i] - ma) * (rb[i] - mb2); saa += (ra[i] - ma) ** 2; sbb += (rb[i] - mb2) ** 2 } return sab / Math.sqrt(saa * sbb) })() : NaN
    return { data: [{ type: 'scatter', mode: 'markers', x: pts.map((p) => p.x), y: pts.map((p) => p.y), customdata: pts.map((p) => p.s.key), text: pts.map((p) => `${p.s.source} ${p.s.chrom} seg ${p.s.segId}`), marker: { size: 6, color: '#4c72b0', opacity: 0.7 }, hovertemplate: '%{text}<br>%{x} bins · %{y:.2f}<extra></extra>' }], layout: { ...diagCommon, height: 300, title: { text: `segment length vs max |flank ${stat === 'srd' ? 'SRD' : 'SRD/√φ'}| · Spearman ρ = ${Number.isFinite(rho) ? rho.toFixed(3) : 'N/A'}`, font: { size: 12 }, x: 0 }, xaxis: { title: { text: 'segment length (bins)', font: { size: 11 } } }, yaxis: { title: { text: stat === 'srd' ? '|SRD|' : '|SRD/√φ|', font: { size: 11 } } }, showlegend: false } }
  }
  const noiseFig = (which: 'std' | 'mean' | 'cv') => {
    const rows: { s: SmallSeg; x: number; y: number; c: number }[] = []
    for (const s of allSmall) for (const f of ['left', 'right'] as const) { const d = s[f][effectMetric], sp = s[f].srd_phi; if (d === null || sp === null) continue; const segV = which === 'std' ? s.stats.std : which === 'mean' ? s.stats.mean : s.stats.cv; const flV = which === 'std' ? (f === 'left' ? s.stats.flankStdL : s.stats.flankStdR) : which === 'mean' ? (f === 'left' ? s.stats.flankMeanL : s.stats.flankMeanR) : (f === 'left' ? s.stats.flankCvL : s.stats.flankCvR); const c = Math.max(segV, flV ?? -Infinity); if (!Number.isFinite(c)) continue; rows.push({ s, x: d, y: Math.abs(sp), c }) }
    const color = which === 'cv' ? rows.map((r) => VIRIDIS[Math.min(4, binIndex(r.c, NOISE_CV_BINS))]) : rows.map((r) => r.c)
    return { data: [{ type: 'scatter', mode: 'markers', x: rows.map((r) => r.x), y: rows.map((r) => r.y), customdata: rows.map((r) => r.s.key), text: rows.map((r) => `${r.s.source} ${r.s.chrom} seg ${r.s.segId} · ${which} ${r.c.toFixed(2)}`), marker: { size: 6, color, colorscale: which === 'cv' ? undefined : 'Viridis', showscale: which !== 'cv', colorbar: { thickness: 8, title: { text: `max(seg, flank) ${which}`, side: 'right' } }, opacity: 0.85 }, hovertemplate: '%{text}<br>effect %{x:.2f} · |SRD/√φ| %{y:.2f}<extra></extra>' }], layout: { ...diagCommon, height: 320, title: { text: `colour = max(segment, flank) per-bin ${which}${which === 'cv' ? ' (binned <0.6, 0.6–0.8, 0.8–1.0, 1.0–1.2, ≥1.2)' : ''}`, font: { size: 12 }, x: 0 }, xaxis: { title: { text: `log2 effect (${METRICS.find((m) => m.id === effectMetric)!.label}, signed)`, font: { size: 11 } } }, yaxis: { title: { text: '|SRD/√φ|', font: { size: 11 } } }, shapes: [{ type: 'line', x0: -seg.tau, x1: -seg.tau, y0: 0, y1: 1, yref: 'paper', line: { dash: 'dash', color: '#333', width: 1 } }, { type: 'line', x0: seg.tau, x1: seg.tau, y0: 0, y1: 1, yref: 'paper', line: { dash: 'dash', color: '#333', width: 1 } }, { type: 'line', x0: 0, x1: 1, xref: 'paper', y0: seg.srdRef, y1: seg.srdRef, line: { dash: 'dash', color: '#333', width: 1 } }], showlegend: false } }
  }
  const deltaHist = () => {
    const data: any[] = []; const layout: any = { ...diagCommon, height: 420, grid: { rows: 2, columns: 2, pattern: 'independent' }, title: { text: 'small-segment flank log2 effects: IQR-trimmed-mean vs plain-mean levels, split by the flank with min vs max |SRD/√φ| (dotted = ±τ)', font: { size: 12 }, x: 0 }, showlegend: false }
    const specs: [Metric, string][] = [['delta_trim', 'IQR-trimmed-mean levels'], ['delta_mean', 'mean levels']]
    let k = 1
    for (const [m, ml] of specs) for (const which of ['min', 'max'] as const) {
      const vals = allSmall.map((s) => { const f = pickFlank(s, 'min_srd', effectMetric); if (!f) return null; const side = which === 'min' ? f : (f === 'left' ? 'right' : 'left'); return s[side][m] }).filter((v): v is number => v !== null)
      const mean = vals.reduce((a, b) => a + b, 0) / vals.length; const sorted = [...vals].sort((a, b) => a - b); const med = sorted.length ? (sorted.length%2 ? sorted[sorted.length>>1] : (sorted[(sorted.length>>1)-1]+sorted[sorted.length>>1])/2) : NaN
      const ax = k === 1 ? '' : String(k)
      data.push({ type: 'histogram', x: vals, nbinsx: 40, marker: { color: '#4c72b0' }, opacity: 0.85, xaxis: `x${ax}`, yaxis: `y${ax}`, hovertemplate: '%{x:.2f}: %{y}<extra></extra>' })
      layout[`xaxis${ax}`] = { title: { text: `${ml} — ${which}-SRD flank · mean ${mean.toFixed(2)}, median ${med.toFixed(2)}, n=${vals.length}`, font: { size: 10 } } }
      layout[`yaxis${ax}`] = { title: { text: 'segments', font: { size: 10 } } }
      layout.shapes = [...(layout.shapes || []), { type: 'line', xref: `x${ax}`, yref: 'paper', x0: mean, x1: mean, y0: 0, y1: 1, line: { color: 'red', width: 1.5 } }, { type: 'line', xref: `x${ax}`, yref: 'paper', x0: med, x1: med, y0: 0, y1: 1, line: { color: 'black', width: 1.5, dash: 'dash' } }, { type: 'line', xref: `x${ax}`, yref: 'paper', x0: seg.tau, x1: seg.tau, y0: 0, y1: 1, line: { color: 'grey', width: 1, dash: 'dot' } }, { type: 'line', xref: `x${ax}`, yref: 'paper', x0: -seg.tau, x1: -seg.tau, y0: 0, y1: 1, line: { color: 'grey', width: 1, dash: 'dot' } }]
      k++
    }
    return { data, layout }
  }
  const onDiagClick = (ev: any) => { const p = ev?.points?.[0]; if (p && typeof p.customdata === 'string') { const s = ig.smallSegs.find((x) => x.key === p.customdata); if (s) selectSeg(s) } }

  if (ig.error && !meta) return <div className="integ"><div className="status-box">Integration data unavailable: {ig.error}</div></div>
  if (!meta) return <div className="integ"><div className="status-box">Loading the multimodal run… {ig.status?.index?.step || ''}</div></div>
  const segTable = clSegs.filter((c) => c.chromosome === seg.chrom && seg.sources.includes(c.group))
  const sensitivity = meta.threshold_sensitivity || []
  const sens = (() => { const m = new Map<string, number>(); for (const r of sensitivity) { const k = `${r.tau}|${r.stab}`; m.set(k, (m.get(k) || 0) + (r.n_merged || 0)) } return m })()
  const taus = [...new Set(sensitivity.map((r) => r.tau))].sort(), stabs = [...new Set(sensitivity.map((r) => r.stab))].sort()

  return (
    <div className="integ">
      <div className="section-heading"><div><div className="eyebrow">INTEGRATION · PART 2</div><h2>Pseudobulk segmentation and CN assignment</h2><p>One segmentation per WNN cluster (accepted pseudobulk breakpoints, {meta.breakpoints?.length} rows). Fig 3c calls CN on each cluster's pseudobulk; fig 3d applies the same cluster segmentation to every member cell with per-cell scaling. Tracks below show the pseudobulk depth with the accepted breakpoints and the small internal segments the merge prototype nominated.</p></div>
        <div className="small muted">CN scaling is diploid-anchored and bounded: "2" means high relative to the cluster/cell median, not an absolute copy number.</div></div>

      {!meta.flank_scores?.length && <p className="status-box">Original flank-score and bootstrap exports are unavailable. Empty source diagnostics do not mean zero effects or zero merges. Generate a descriptive alternative below from the available pseudobulk and accepted breakpoints.</p>}
      <AlternativeFlankView />
      <div className="panel">
        <h3>Cluster-level CN (fig 3c) — one row per cluster, plus all cells</h3>
        <p className="lead">Each cluster's pseudobulk X (mean over members) is normalised and called with pyEpi <code>assign_gainloss_new</code> on its own segmentation. Click a row to select that cluster and jump its chromosome to the tracks.</p>
        <div className="controls-row">
          <label>palette<select value={seg.karyoPalette} onChange={(e) => ig.setSeg({ karyoPalette: e.target.value as any })}><option value="story">story figures (purple / green / red)</option><option value="explorer">explorer (blue / white / red)</option></select></label>
          <div className="seg-legend">{[0, 0.5, 1, 1.5, 2].map((v) => <span key={v}><i style={{ background: (seg.karyoPalette === 'story' ? STORY_CN : CN_COLORS)[v] }} />{v} = {STORY_CN_LABELS[v]}</span>)}<span><i style={{ background: '#e6e6e6' }} />not called</span></div>
        </div>
        <ClusterKaryogram rows={clusterRows} palette={seg.karyoPalette} rowHeight={26} labelWidth={150} selectedRow={ig.selectedCluster !== null ? meta.sources.indexOf(`cluster${ig.selectedCluster}`) : null}
          onClick={(ri, _bin, chrom) => { const src = meta.sources[ri]; ig.setSeg({ chrom, sources: src === 'all' ? seg.sources : (seg.sources.includes(src) ? seg.sources : [...seg.sources, src]) }); if (src !== 'all') { const lv = src.replace('cluster', ''); ig.set({ selectedCluster: lv, highlight: cells.filter((c) => clusterOf(c, 'wnn_leiden_0.3') === lv).map((c) => c.cell) }) } document.getElementById('cluster-tracks')?.scrollIntoView({ behavior: 'smooth' }) }}
          hoverText={(ri, bin, chrom) => { const src = meta.sources[ri]; const c = exMeta.chromosomes.find((x) => x.name === chrom)!; const local = bin - c.offset; const g = clSegs.find((x) => x.group === src && x.chromosome === chrom && local >= x.s && local < x.e); return g ? [`segment ${g.segment}: bins [${g.s}, ${g.e}) · ${g.n_bins} bins · continuous ${g.continuous_cn} · integer ${g.integer_cn} · watson ${g.watson}`] : [] }} />
      </div>

      <div className="panel">
        <h3>Per-cell CN under the cluster segmentation (fig 3d) — {cells.length} cells</h3>
        <p className="lead">Same breakpoints as the row above for each cell's cluster, but normalisation, scaling and calling run per cell (<code>assign_gainloss_new(cell_X, cluster_seg_labels)</code>). It shows how consistently member cells support the cluster-level states; a clear pseudobulk event does not imply every member supports it equally.</p>
        <div className="controls-row">
          <label>row order<select value={seg.percellOrder} onChange={(e) => ig.setSeg({ percellOrder: e.target.value as any })}><option value="cluster">WNN cluster (res 0.3)</option><option value="coverage">raw coverage ↓</option><option value="chr8">chr8 score ↓</option></select></label>
          <label>row height<select value={String(seg.percellRowH)} onChange={(e) => ig.setSeg({ percellRowH: e.target.value === 'fit' ? 'fit' : (+e.target.value as 1 | 2) })}><option value="fit">fit (sampled rows)</option><option value="1">1 px (every cell, scroll)</option><option value="2">2 px (every cell, scroll)</option></select></label>
          <div className="seg-legend">{levels.map((lv) => <span key={lv}><i style={{ background: clusterColor(lv, levels), border: 'none' }} />cluster {lv} (n={cells.filter((c) => clusterOf(c, 'wnn_leiden_0.3') === lv).length})</span>)}</div>
          {ig.karyoLoading && <span className="muted small">loading per-cell matrix…</span>}
        </div>
        {cellRows.length ? <ClusterKaryogram rows={cellRows} palette={seg.karyoPalette} rowHeight={seg.percellRowH} maxHeight={760} labelWidth={0} selectedRow={selectedCellRow}
          onClick={(ri) => { const c = cells[cellOrder[ri]]; ig.focus(c.cell) }} hoverText={(ri) => { const c = cells[cellOrder[ri]]; return [`raw ${c.total_raw_depth.toLocaleString()} · chr8 score ${fmt(c.chr8_score, 2)}${hlSet.has(c.cell) ? ' · highlighted' : ''}${c.in_explorer_cohort ? ' · in explorer cohort' : ''}`] }} /> : <div className="muted small">per-cell CN matrix not loaded</div>}
        <CellCard />
      </div>

      <div className="panel" id="cluster-tracks">
        <h3>Chromosome tracks — pseudobulk depth, accepted breakpoints and small segments (figs 9, 10)</h3>
        <p className="lead">Grey line: per-bin pseudobulk depth, broken at genomic gaps. Vertical lines: accepted breakpoints by recursion depth. Shaded spans: internal segments ≤ {seg.smallMax} bins with their flank metrics from the merge prototype. Tick "flank-consistency hatch" to hatch × where the min-delta flank's |SRD/√φ| &lt; {seg.srdRef} (weak) and / where the min-delta flank ≠ the min-SRD flank (inconsistent); plain fill = ≥ {seg.srdRef} and consistent (strong). Tick "SRD-screen proposal" to instead colour and hatch segments by the min-SRD flank's |SRD/√φ| and the screening proposal (/ merge, × ambiguous, plain = keep focal). {pilot ? '' : `Small-segment metrics exist only for the pilot chromosomes (${meta.pilot_chromosomes.join(', ')}); this chromosome shows breakpoints and CN only.`}</p>
        <div className="controls-row">
          <label>chromosome<select value={seg.chrom} onChange={(e) => ig.setSeg({ chrom: e.target.value, selectedSegment: null })}>{exMeta.chromosomes.map((c) => <option key={c.name} value={c.name}>{c.name}{meta.pilot_chromosomes.includes(c.name) ? ' · pilot' : ''}</option>)}</select></label>
          <label>colour by<select value={seg.metric} onChange={(e) => ig.setSeg({ metric: e.target.value as Metric })}>{METRICS.map((m) => <option key={m.id} value={m.id}>{m.label}</option>)}</select></label>
          {metricDef.kind === 'stat' && <label>effect metric for flank choice<select value={seg.effectMetric} onChange={(e) => ig.setSeg({ effectMetric: e.target.value as Metric })}>{METRICS.filter((m) => m.kind === 'effect').map((m) => <option key={m.id} value={m.id}>{m.label}</option>)}</select></label>}
          <label>flank<select value={seg.flank} onChange={(e) => ig.setSeg({ flank: e.target.value as any })}><option value="min_delta">min |effect| flank (figure default)</option><option value="min_srd">min |SRD/√φ| flank</option><option value="max_delta">max |effect| flank</option><option value="left">left flank</option><option value="right">right flank</option></select></label>
          <label>scale<select value={seg.scale} onChange={(e) => ig.setSeg({ scale: e.target.value as any })}><option value="binned">binned bands (as in the figures)</option><option value="continuous">continuous viridis</option></select></label>
          <label>track level<select value={seg.basis} onChange={(e) => ig.setSeg({ basis: e.target.value as any })}><option value="chrmedian">log2(depth / chromosome median)</option><option value="genome_trimmed">log2(depth / genome trimmed mean, pyEpi baseline)</option><option value="raw">raw pseudobulk X</option></select></label>
          <label>smoothing (visual)<select value={seg.smoothing} onChange={(e) => ig.setSeg({ smoothing: +e.target.value })}><option value={1}>none (per bin)</option><option value={5}>5 bins</option><option value={11}>11 bins</option><option value={25}>25 bins</option></select></label>
          <label>small ≤<input type="number" min={5} max={500} value={seg.smallMax} onChange={(e) => ig.setSeg({ smallMax: +e.target.value })} style={{ width: 64 }} /></label>
          <label>SRD ref<input type="number" step={0.1} value={seg.srdRef} onChange={(e) => ig.setSeg({ srdRef: +e.target.value })} style={{ width: 64 }} /></label>
          <label>τ<input type="number" step={0.05} value={seg.tau} onChange={(e) => ig.setSeg({ tau: +e.target.value })} style={{ width: 64 }} /></label>
        </div>
        <div className="controls-row">
          <span className="small muted">rows:</span>{meta.sources.map((s) => <label key={s} className="chk-inline"><input type="checkbox" checked={seg.sources.includes(s)} onChange={(e) => ig.setSeg({ sources: e.target.checked ? meta.sources.filter((x) => x === s || seg.sources.includes(x)) : seg.sources.filter((x) => x !== s) })} /><i className="sw" style={{ background: SOURCE_COLORS[s] }} />{s}</label>)}
          <span className="small muted" style={{ marginLeft: 12 }}>layers:</span>
          {([['breakpoints', 'accepted breakpoints'], ['depthColor', 'colour by recursion depth'], ['small', 'small segments'], ['consistency', 'flank-consistency hatch (×/‒/) + x/o/ label'], ['proposal', 'SRD-screen proposal: hatch + colour by min-SRD flank'], ['reason', 'retained reason (effect ≥ τ / terminal)'], ['segLabels', 'length + value labels'], ['cn', 'cluster CN strip'], ['candidates', 'candidate audit (rejected / blocked)'], ['allelic', 'allelic |dBAF| pseudobulk (right axis)']] as const).map(([k, l]) => <label key={k} className="chk-inline"><input type="checkbox" checked={(seg.labels as any)[k]} onChange={(e) => ig.setLabels({ [k]: e.target.checked } as any)} />{l}</label>)}
        </div>
        <div className="seg-legend">
          <span><b>{metricDef.label}</b> · {metricDef.source}</span>
          {seg.scale === 'binned' ? metricDef.binLabels.map((l, i) => <span key={l}><i style={{ background: metricDef.bins.length === 3 ? [VIRIDIS[0], VIRIDIS[1], VIRIDIS[2], VIRIDIS[4]][i] : VIRIDIS[i] }} />{l}</span>) : <span>continuous 0 → {metricDef.kind === 'stat' ? 8 : 2} (viridis)</span>}
          <span><i style={{ background: '#e5e7eb' }} />N/A</span>
          {seg.labels.breakpoints && seg.labels.depthColor && <span>depth: {DEPTH_COLORS.map((c, i) => <span key={i} className="ln" style={{ borderColor: c }} title={`depth ${i}`} />)}</span>}
          {seg.labels.cn && <span>CN strip: {[0, 0.5, 1, 1.5, 2].map((v) => <i key={v} style={{ background: (seg.karyoPalette === 'story' ? STORY_CN : CN_COLORS)[v] }} />)}</span>}
          {seg.labels.proposal && <span>segments with a proposal: colour = min-SRD flank |SRD/√φ| ({METRICS.find((m) => m.id === 'srd_phi')!.binLabels.join(' / ')}) · hatch / merge, × ambiguous, plain keep-focal</span>}
          {seg.labels.consistency && !seg.labels.proposal && <span>hatch × weak SRD, / inconsistent flank, plain strong/consistent</span>}
        </div>
        {trackFig ? <><Plot data={trackFig.data} layout={trackFig.layout} config={{ displayModeBar: presentation ? false : 'hover' }} onClick={onTrackClick} onReady={setPlotEl} /><button className="btn-xs" onClick={() => plotEl && downloadPlot(plotEl, `tracks_${seg.chrom}`)}>PNG</button></> : <div className="muted small">{ig.pb ? 'select at least one cluster row' : 'loading pseudobulk profiles…'}</div>}
        {seg.selectedSegment && (() => { const s = ig.smallSegs.find((x) => x.key === seg.selectedSegment); if (!s) return null; const cc = consistencyClass(s, effectMetric, seg.srdRef); const bs = (meta.bootstrap || []).find((b: any) => b.source === s.source && b.chromosome === s.chrom && b.bin_start === s.s); return (
          <div className="cell-card"><div className="row wrap"><b>{s.source} {s.chrom} · segment {s.segId}</b> <span className="muted">bins [{s.s}, {s.e}) · {s.len} bins · {s.mbStart}–{s.mbEnd} Mb</span><span className="spacer" /><button className="btn-xs" onClick={() => ig.setSeg({ selectedSegment: null })}>✕</button></div>
            <div className="table-scroll"><table className="tbl"><thead><tr><th></th><th>left flank ({s.flank.n_bins_left_flank ?? '–'} bins)</th><th>right flank ({s.flank.n_bins_right_flank ?? '–'} bins)</th></tr></thead><tbody>
              {METRICS.map((m) => <tr key={m.id} className={m.id === seg.metric ? 'sel' : ''}><td>{m.label}</td><td>{fmt(s.left[m.id])}</td><td>{fmt(s.right[m.id])}</td></tr>)}
              <tr><td>flank level (median / mean / trim)</td><td>{fmt(s.flank.median_left, 2)} / {fmt(s.flank.mean_left_flank_perbin, 2)} / {fmt(s.flank.trimmean_left_flank_perbin, 2)}</td><td>{fmt(s.flank.median_right, 2)} / {fmt(s.flank.mean_right_flank_perbin, 2)} / {fmt(s.flank.trimmean_right_flank_perbin, 2)}</td></tr>
              <tr><td>Original median bootstrap fraction |Δ| &lt; {meta.manifest?.params?.tau_default} ({meta.manifest?.params?.n_bootstrap} cell replicates)</td><td>{fmt(s.flank.frac_under_left, 2)}</td><td>{fmt(s.flank.frac_under_right, 2)}</td></tr>
              <tr><td>eligible / merged / reason</td><td>{String(s.flank.eligible_left)} / {String(s.flank.merged_left)} / {s.reasonL}</td><td>{String(s.flank.eligible_right)} / {String(s.flank.merged_right)} / {s.reasonR}</td></tr>
              {bs && <tr><td>SRD-screen bootstrap SE (log2FC / SRD)</td><td>{fmt(bs.se_log2fc_left, 3)} / {fmt(bs.se_srd_left, 3)} (sign agreement {bs.sign_agree_log2fc_left})</td><td>{fmt(bs.se_log2fc_right, 3)} / {fmt(bs.se_srd_right, 3)} (sign agreement {bs.sign_agree_log2fc_right})</td></tr>}
            </tbody></table></div>
            <div className="small">segment level: median {s.stats.median.toFixed(3)} · mean {s.stats.mean.toFixed(3)} · IQR-trimmed mean {s.stats.trim.toFixed(3)} · per-bin std {s.stats.std.toFixed(3)} · CV {s.stats.cv.toFixed(3)} · {cc.label}{s.proposal ? ` · SRD-screen proposal: ${s.proposal}` : ''}</div>
          </div>) })()}
        <details className="diag" open={chromTable}><summary onClick={(e) => { e.preventDefault(); setChromTable(!chromTable) }}>Segments on {seg.chrom} for the selected rows ({segTable.length} cluster segments · {smallHere.filter((s) => seg.sources.includes(s.source)).length} small segments with metrics)</summary>
          <div className="table-scroll"><table className="tbl"><thead><tr><th>cluster</th><th>segment</th><th>bins [s, e)</th><th>n</th><th>Mb</th><th>CN call</th><th>continuous</th><th>integer</th><th>small-segment metric ({metricDef.label})</th><th>consistency</th><th>proposal</th></tr></thead><tbody>
            {segTable.map((g) => { const sm = smallHere.find((s) => s.source === g.group && s.s === g.s); const mv = sm ? metricValue(sm, seg.metric, seg.flank, effectMetric) : null; return <tr key={`${g.group}-${g.segment}`} className={`${sm ? 'clickable' : ''} ${sm && sm.key === seg.selectedSegment ? 'sel' : ''}`} onClick={() => sm && selectSeg(sm)}><td style={{ color: SOURCE_COLORS[g.group] }}>{g.group}</td><td>{g.segment}</td><td>[{g.s}, {g.e})</td><td>{g.n_bins}</td><td>{(bins[2 * (chromMeta.offset + g.s)] / MB).toFixed(1)}–{(bins[2 * (chromMeta.offset + g.e - 1) + 1] / MB).toFixed(1)}</td><td><i className="sw" style={{ background: (seg.karyoPalette === 'story' ? STORY_CN : CN_COLORS)[g.cn_call] }} />{g.cn_call} {STORY_CN_LABELS[g.cn_call]}</td><td>{g.continuous_cn}</td><td>{g.integer_cn}</td><td>{sm ? `${fmt(mv!.value)} (${mv!.flank})` : ''}</td><td>{sm ? consistencyClass(sm, effectMetric, seg.srdRef).cls : ''}</td><td>{sm?.proposal ? PROPOSAL_LABEL[sm.proposal] || sm.proposal : ''}</td></tr> })}
          </tbody></table></div></details>
      </div>

      {allSmall.length > 0 && <div className="panel">
        <h3>Small-segment diagnostics across the pilot chromosomes ({allSmall.length} segments; figs 4, 6, 7, 8)</h3>
        <p className="lead">Each dot is one small segment (or one segment–flank comparison in the noise panels). Click a dot to open that segment in the tracks. τ and the SRD reference are inspection settings, not calibrated thresholds; SRD/√φ is a dispersion-scaled deviance statistic, not a p-value.</p>
        <div className="two-col">
          {(() => { const f = volcano(); return <Plot data={f.data} layout={f.layout} onClick={onDiagClick} /> })()}
          {(() => { const f = deltaHist(); return <Plot data={f.data} layout={f.layout} /> })()}
          {(() => { const f = sizeDep('srd'); return <Plot data={f.data} layout={f.layout} onClick={onDiagClick} /> })()}
          {(() => { const f = sizeDep('srd_phi'); return <Plot data={f.data} layout={f.layout} onClick={onDiagClick} /> })()}
        </div>
        <div className="emb-grid">
          {(['std', 'mean', 'cv'] as const).map((w) => { const f = noiseFig(w); return <Plot key={w} data={f.data} layout={f.layout} onClick={onDiagClick} /> })}
        </div>
      </div>}

      {!!meta.review_summary?.length && <details className="diag">
        <summary>Merge prototype summary (τ = {meta.manifest?.params?.tau_default}, stability ≥ {meta.manifest?.params?.stab_default})</summary>
        <p className="lead">Original prototype: median segment levels, {meta.manifest?.params?.n_bootstrap} whole-cell bootstrap replicates, seed {meta.manifest?.params?.seed}. A boundary is nominated if it touches an internal segment ≤ {meta.manifest?.params?.max_small_bins} bins. Merging requires observed |log2(m_A/m_B)| &lt; τ, the required bootstrap fraction, all replicates valid, and a chain guard on the original constituent levels. Terminal and gap cases remain protected. Bootstrap fractions are stability scores, not probabilities. The alternative view above computes a separate arithmetic-mean bootstrap.</p>
        <div className="two-col">
          <div className="table-scroll"><table className="tbl"><thead><tr><th>cluster</th><th>chrom</th><th>segments</th><th>boundaries</th><th>nominated</th><th>eligible</th><th>merged</th><th>min |Δ| nominated</th><th>median |Δ|</th><th>max stability</th></tr></thead><tbody>
            {(meta.review_summary || []).map((r: any, i: number) => <tr key={i}><td>{r.source}</td><td>{r.chromosome}</td><td>{r.n_orig_segments}</td><td>{r.n_boundaries}</td><td>{r.n_nominated}</td><td>{r.n_eligible_default}</td><td>{r.n_merged_default}</td><td>{fmt(r.min_abs_delta_nominated, 3)}</td><td>{fmt(r.median_abs_delta_nominated, 3)}</td><td>{fmt(r.max_frac_under_default, 2)}</td></tr>)}
          </tbody></table></div>
          <div><div className="tbl-title">merged boundaries summed over clusters × chromosomes (sensitivity grid)</div><table className="tbl"><thead><tr><th>τ \ stability</th>{stabs.map((s) => <th key={s}>{s}</th>)}</tr></thead><tbody>{taus.map((t) => <tr key={t}><td>{t}</td>{stabs.map((s) => <td key={s}>{sens.get(`${t}|${s}`) ?? '–'}</td>)}</tr>)}</tbody></table>
            <p className="small muted">These are exported results at their original settings; changing display thresholds does not rerun the bootstrap or merge rule. Missing table entries do not imply zero merges.</p></div>
        </div>
      </details>}
    </div>
  )
}
export { setToastUnused as _ }
const setToastUnused = null
