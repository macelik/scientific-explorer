import { useEffect, useMemo, useState } from 'react'
import Plot from './Plot'
import { useIntegration, type ColorBy } from '../integrationStore'
import { useStore } from '../store'
import { median, boxStats, clusterColor, clusterLevels, clusterOf, jitter, pearson, skewness, spearman } from '../integration'
import type { IntegrationCell } from '../types'

interface ColorOpt { id: ColorBy; label: string; get: (c: IntegrationCell) => number | null; continuous: boolean; unit?: string }
const COLOR_OPTIONS: ColorOpt[] = [
  { id: 'cluster', label: 'cluster (selected labelling)', get: () => null, continuous: false },
  { id: 'log10_raw', label: 'log10 raw read depth', get: (c) => Math.log10(c.total_raw_depth), continuous: true, unit: 'log10 raw counts' },
  { id: 'sumX', label: 'ΣX (GC-corrected, normalised)', get: (c) => c.gc_norm_depth, continuous: true, unit: 'ΣX' },
  { id: 'HA_HB', label: 'HA+HB allelic counts', get: (c) => c.snp_coverage, continuous: true, unit: 'HA+HB' },
  { id: 'n_bins_snp', label: 'bins with SNP coverage', get: (c) => c.n_bins_snp_covered, continuous: true, unit: 'bins' },
  { id: 'depth_weight', label: 'WNN depth-modality weight', get: (c) => c.wnn_depth_weight, continuous: true, unit: 'depth weight' },
  { id: 'haplo_weight', label: 'WNN haplotype-modality weight', get: (c) => c.wnn_haplo_weight, continuous: true, unit: 'haplo weight' },
  { id: 'chr8_score', label: 'chr8 log2 regional depth ratio', get: (c) => c.chr8_score ?? null, continuous: true, unit: 'log2(46–120 Mb / 10–35 Mb)' },
  { id: 'zero_frac', label: 'fraction of zero bins', get: (c) => c.zero_frac ?? null, continuous: true, unit: 'zero fraction' },
  { id: 'in_cohort', label: 'in the 300-cell explorer cohort', get: (c) => (c.in_explorer_cohort ? 1 : 0), continuous: false },
  { id: 'selected', label: 'selected / highlighted cells', get: () => null, continuous: false },
]
const EMB = [
  { key: 'depth', title: 'Depth-only UMAP', sub: 'log1p(X) → PCA(50) → cosine kNN(20)', x: 'umap_depth_x', y: 'umap_depth_y', own: 'depth' },
  { key: 'haplo', title: 'Allelic-imbalance UMAP', sub: '|HA−HB| / (2·(AT+2)) → PCA(50) → cosine kNN(20)', x: 'umap_haplo_x', y: 'umap_haplo_y', own: 'haplo' },
  { key: 'wnn', title: 'Integrated WNN UMAP', sub: 'muon WNN of the two graphs', x: 'umap_wnn_x', y: 'umap_wnn_y', own: 'wnn' },
]
const fmt = (v: number | null | undefined, d = 3) => (v === null || v === undefined || !Number.isFinite(v) ? 'N/A' : v.toFixed(d))

export function CellCard() {
  const meta = useIntegration((s) => s.meta)!
  const focusCell = useIntegration((s) => s.focusCell); const detail = useIntegration((s) => s.cellDetail); const focus = useIntegration((s) => s.focus)
  const clusterKey = useIntegration((s) => s.clusterKey)
  const openCell = useStore((s) => s.openCell); const chrom = useStore((s) => s.chrom)
  const explorerCells = useStore((s) => s.meta?.cells)
  if (!focusCell) return <div className="cell-card muted">Click a cell in an embedding or the per-cell karyogram to inspect its memberships, coverage, weights and graph neighbours.</div>
  const c = meta.cells.find((x) => x.cell === focusCell)
  if (!c) return <div className="cell-card">Cell {focusCell} not in the h5mu.</div>
  const inCohort = !!explorerCells?.some((x) => x.cellID === focusCell)
  const keys = [...meta.cluster_keys.wnn, ...meta.cluster_keys.depth, ...meta.cluster_keys.haplo]
  return (
    <div className="cell-card">
      <div className="row wrap"><b>{c.cell}</b> <span className="muted">· {clusterKey} = {clusterOf(c, clusterKey)}</span><span className="spacer" />
        {inCohort ? <button className="btn-sm" onClick={() => openCell(c.cell, chrom || 'chr8', { tab: 'explore' })}>open in cell explorer</button> : <span className="muted small">not in the 300-cell explorer cohort (cell-level tools unavailable)</span>}
        <button className="btn-xs" onClick={() => focus(null)}>✕</button></div>
      <div className="two-col">
        <table className="kv"><tbody>
          <tr><td>memberships</td><td>{keys.map((k) => `${k.replace('_leiden_', ' r')}: ${clusterOf(c, k)}`).join(' · ')}</td></tr>
          <tr><td>coverage</td><td>raw {c.total_raw_depth.toLocaleString()} · ΣX {c.gc_norm_depth.toFixed(0)} · HA+HB {c.snp_coverage.toFixed(0)} · SNP-covered bins {c.n_bins_snp_covered}</td></tr>
          <tr><td>WNN weights</td><td>depth {c.wnn_depth_weight.toFixed(3)} · haplo {c.wnn_haplo_weight.toFixed(3)}</td></tr>
          <tr><td>chr8 score</td><td>{fmt(c.chr8_score)} · zero-bin fraction {fmt(c.zero_frac)} · per-cell best_s (cluster segmentation) {fmt(c.percell_best_s, 4)}</td></tr>
          <tr><td>PC1 / spectral</td><td>depth PC1 {fmt(c.depth_pc1)} · allelic PC1 {fmt(c.haplo_pc1)} · WNN spectral 1 {fmt(c.wnn_spectral1, 4)}</td></tr>
        </tbody></table>
        <div>
          {detail ? (
            <>
              <div className="small"><b>kNN neighbours</b> (connectivities): WNN {detail.overlap.n_wnn}, depth {detail.overlap.n_depth}, allelic {detail.overlap.n_haplo} · shared depth∩allelic {detail.overlap.depth_vs_haplo}, WNN∩depth {detail.overlap.wnn_vs_depth}, WNN∩allelic {detail.overlap.wnn_vs_haplo}</div>
              {(['wnn', 'depth', 'haplo'] as const).map((g) => <div key={g} className="small muted">{g}: {detail.neighbors[g].slice(0, 8).map((n) => `${n.cell.slice(0, 8)}(${n.cluster})`).join(', ')}{detail.neighbors[g].length > 8 ? ', …' : ''}</div>)}
            </>
          ) : <div className="small muted">loading neighbours…</div>}
        </div>
      </div>
    </div>
  )
}

export default function IntegrationClustersView() {
  const ig = useIntegration()
  const presentation = useStore((s) => s.presentation)
  const [drag, setDrag] = useState<'lasso' | 'zoom' | 'pan'>('lasso')
  const [covKey, setCovKey] = useState('wnn_leiden_0.3'); const [covMetric, setCovMetric] = useState<'total_raw_depth' | 'gc_norm_depth' | 'snp_coverage' | 'n_bins_snp_covered'>('gc_norm_depth')
  const [colorScatter, setColorScatter] = useState(true)
  useEffect(() => { ig.load() }, [])
  const meta = ig.meta
  const cells = meta?.cells || []
  const levels = useMemo(() => (meta ? clusterLevels(cells, ig.clusterKey) : []), [meta, ig.clusterKey])
  const hl = useMemo(() => new Set(ig.highlight), [ig.highlight])
  const opt = COLOR_OPTIONS.find((o) => o.id === ig.colorBy)!
  // chr8 diagnostic with resolution comparison
  const chr8 = useMemo(() => {
    const base = ig.selectedCluster ?? (levels.includes('4') ? '4' : levels[levels.length - 1])
    const members = cells.filter((c) => clusterOf(c, ig.clusterKey) === base)
    const cmpLevels = clusterLevels(cells, ig.compareKey)
    const counts = new Map<string, number>(); for (const c of members) counts.set(clusterOf(c, ig.compareKey), (counts.get(clusterOf(c, ig.compareKey)) || 0) + 1)
    const target = [...counts.entries()].sort((a, b) => b[1] - a[1])[0]?.[0] ?? null
    const added = target === null ? [] : cells.filter((c) => clusterOf(c, ig.compareKey) === target && clusterOf(c, ig.clusterKey) !== base)
    const origin = new Map<string, { n: number; hi: number }>(); for (const c of added) { const o = clusterOf(c, ig.clusterKey); const e = origin.get(o) || { n: 0, hi: 0 }; e.n++; if ((c.chr8_score ?? -9) > 0.5) e.hi++; origin.set(o, e) }
    const retained = target === null ? 0 : members.filter((c) => clusterOf(c, ig.compareKey) === target).length
    return { base, members, target, added, origin, retained, cmpLevels }
  }, [cells, ig.clusterKey, ig.compareKey, ig.selectedCluster, levels])

  if (ig.error && !meta) return <div className="integ"><div className="status-box">Integration data unavailable: {ig.error}<br />Expected inputs: <code>integration-story/data/integrated.h5mu</code> (or the path in <code>data/prototype/run_manifest.json</code>), <code>raw_accepted_breakpoints.tsv</code>, <code>data/derived/</code>, <code>data/prototype/</code>. Set <code>PYEPI_SCIENTIFIC_H5MU</code> to override.</div></div>
  if (!meta) return <div className="integ"><div className="status-box">Loading the multimodal run… {ig.status?.index?.status === 'running' ? `building the pseudobulk index (${ig.status.index.step})` : ig.status?.index?.status || ''}</div></div>

  const embFig = (e: typeof EMB[number]) => {
    const data: any[] = []
    const dim = ig.dimOthers && (hl.size > 0)
    if (opt.continuous) {
      const vals = cells.map((c) => opt.get(c))
      data.push({ type: 'scattergl', mode: 'markers', x: cells.map((c) => c[e.x]), y: cells.map((c) => c[e.y]), text: cells.map((c) => c.cell), customdata: cells.map((c, i) => i),
        marker: { size: ig.pointSize, color: vals, colorscale: 'Viridis', showscale: true, colorbar: { title: { text: opt.unit || '', side: 'right' }, thickness: 10, len: 0.8 }, opacity: 0.8 },
        selectedpoints: hl.size ? cells.map((c, i) => (hl.has(c.cell) ? i : -1)).filter((i) => i >= 0) : undefined, unselected: { marker: { opacity: dim ? 0.08 : 0.3 } }, selected: { marker: { opacity: 1 } },
        hovertemplate: '%{text}<br>' + (opt.unit || '') + ' %{marker.color:.3f}<extra></extra>', showlegend: false })
    } else if (ig.colorBy === 'in_cohort' || ig.colorBy === 'selected') {
      const groups = ig.colorBy === 'in_cohort' ? [['in explorer cohort', (c: IntegrationCell) => c.in_explorer_cohort, '#c1121f'], ['other cells', (c: IntegrationCell) => !c.in_explorer_cohort, '#cbd5e1']] as const
        : [['highlighted', (c: IntegrationCell) => hl.has(c.cell), '#c1121f'], ['other cells', (c: IntegrationCell) => !hl.has(c.cell), '#cbd5e1']] as const
      for (const [name, pred, color] of groups) {
        const idx = cells.map((c, i) => (pred(c) ? i : -1)).filter((i) => i >= 0)
        data.push({ type: 'scattergl', mode: 'markers', name: `${name} (n=${idx.length})`, x: idx.map((i) => cells[i][e.x]), y: idx.map((i) => cells[i][e.y]), text: idx.map((i) => cells[i].cell), customdata: idx, marker: { size: ig.pointSize, color, opacity: 0.8 }, hovertemplate: '%{text}<extra>' + name + '</extra>' })
      }
    } else {
      for (const lv of levels) {
        const idx = cells.map((c, i) => (clusterOf(c, ig.clusterKey) === lv ? i : -1)).filter((i) => i >= 0)
        const sel = hl.size ? idx.map((i, k) => (hl.has(cells[i].cell) ? k : -1)).filter((k) => k >= 0) : undefined
        data.push({ type: 'scattergl', mode: 'markers', name: `${lv} (n=${idx.length})`, x: idx.map((i) => cells[i][e.x]), y: idx.map((i) => cells[i][e.y]), text: idx.map((i) => cells[i].cell), customdata: idx,
          marker: { size: ig.pointSize, color: clusterColor(lv, levels), opacity: 0.8 }, selectedpoints: sel, unselected: { marker: { opacity: dim ? 0.08 : 0.3 } }, selected: { marker: { opacity: 1 } },
          hovertemplate: '%{text}<br>' + ig.clusterKey + ' = ' + lv + '<extra></extra>' })
      }
    }
    const layout: any = { height: 400, margin: { l: 10, r: 10, t: 40, b: 10 }, title: { text: `${e.title}<br><span style="font-size:10px;color:#60717c">${e.sub}</span>`, font: { size: 13 }, x: 0.02 }, paper_bgcolor: 'white', plot_bgcolor: 'white', xaxis: { visible: false }, yaxis: { visible: false, scaleanchor: 'x' }, dragmode: drag, hovermode: presentation ? false : 'closest', showlegend: !opt.continuous, legend: { orientation: 'h', y: -0.02, font: { size: 10 } }, clickmode: 'event+select' }
    return { data, layout }
  }
  const onSelected = (ev: any) => { if (!ev || !ev.points) { return } const ids = ev.points.map((p: any) => p.text as string); if (ids.length) ig.setHighlight(ids) }
  const onClick = (ev: any) => { const p = ev?.points?.[0]; if (p && p.text) ig.focus(p.text) }
  const clusterSummary = (lv: string) => { const m = cells.filter((c) => clusterOf(c, ig.clusterKey) === lv); const med = (f: (c: IntegrationCell) => number) => { const v = m.map(f).sort((a, b) => a - b); return median(v) }; return { n: m.length, raw: med((c) => c.total_raw_depth), x: med((c) => c.gc_norm_depth), at: med((c) => c.snp_coverage), dw: med((c) => c.wnn_depth_weight) } }

  // ---- diagnostics
  const covLevels = clusterLevels(cells, covKey)
  const covFig = () => {
    const lab: Record<string, string> = { total_raw_depth: 'raw ATAC counts per cell', gc_norm_depth: 'GC-corrected normalised ΣX per cell', snp_coverage: 'HA+HB allelic counts per cell', n_bins_snp_covered: 'bins with SNP coverage' }
    const data = covLevels.map((lv) => { const v = cells.filter((c) => clusterOf(c, covKey) === lv).map((c) => c[covMetric] as number); return { type: 'box', name: `${lv} (n=${v.length})`, y: v, boxpoints: false, marker: { color: clusterColor(lv, covLevels) }, line: { color: '#334155' }, fillcolor: clusterColor(lv, covLevels) + '66', hoverinfo: 'y+name' } })
    return { data, layout: { height: 320, margin: { l: 60, r: 10, t: 30, b: 60 }, yaxis: { type: 'log', title: { text: lab[covMetric], font: { size: 11 } } }, xaxis: { title: { text: covKey, font: { size: 11 } } }, showlegend: false, paper_bgcolor: 'white', plot_bgcolor: 'white', title: { text: `${lab[covMetric]} by ${covKey} (box = IQR, whiskers = non-outlier range; log y)`, font: { size: 12 }, x: 0 }, hovermode: presentation ? false : 'closest' } }
  }
  const scatterFig = (x: (c: IntegrationCell) => number | null, y: (c: IntegrationCell) => number | null, xl: string, yl: string, logx = false) => {
    const pts = cells.map((c, i) => ({ x: x(c), y: y(c), i })).filter((p) => p.x !== null && p.y !== null && Number.isFinite(p.x!) && Number.isFinite(p.y!)) as { x: number; y: number; i: number }[]
    const r = pts.length > 2 ? pearson(pts.map((p) => p.x), pts.map((p) => p.y)) : NaN
    const data: any[] = [{ type: 'scatter', mode: 'markers', x: pts.map((p) => p.x), y: pts.map((p) => p.y), text: pts.map((p) => cells[p.i].cell), marker: { size: 3, opacity: 0.45, color: colorScatter ? pts.map((p) => clusterColor(clusterOf(cells[p.i], ig.clusterKey), levels)) : '#4c72b0' }, hovertemplate: '%{text}<br>' + xl + ' %{x:.4g}<br>' + yl + ' %{y:.4g}<extra></extra>' }]
    return { data, layout: { height: 300, margin: { l: 60, r: 10, t: 36, b: 44 }, xaxis: { title: { text: xl, font: { size: 11 } }, type: logx ? 'log' : 'linear' }, yaxis: { title: { text: yl, font: { size: 11 } } }, title: { text: `${xl} vs ${yl} · Pearson r = ${Number.isFinite(r) ? (r >= 0 ? '+' : '') + r.toFixed(3) : 'N/A'} (n=${pts.length})`, font: { size: 12 }, x: 0 }, paper_bgcolor: 'white', plot_bgcolor: 'white', showlegend: false, hovermode: presentation ? false : 'closest' } }
  }
  const weightsFig = () => {
    const d = cells.map((c) => c.wnn_depth_weight); const m = d.reduce((a, b) => a + b, 0) / d.length; const s = [...d].sort((a, b) => a - b); const med = median(s)
    return { data: [{ type: 'histogram', x: d, nbinsx: 50, marker: { color: '#1f77b4' }, opacity: 0.85, hovertemplate: 'weight %{x:.3f}: %{y} cells<extra></extra>' }], layout: { height: 300, margin: { l: 50, r: 10, t: 36, b: 44 }, shapes: [{ type: 'line', x0: m, x1: m, y0: 0, y1: 1, yref: 'paper', line: { color: 'red', width: 2 } }, { type: 'line', x0: med, x1: med, y0: 0, y1: 1, yref: 'paper', line: { color: 'black', width: 2, dash: 'dash' } }, { type: 'line', x0: 0.5, x1: 0.5, y0: 0, y1: 1, yref: 'paper', line: { color: 'grey', width: 1, dash: 'dot' } }], title: { text: `Depth-weight distribution · mean ${m.toFixed(3)} (red), median ${med.toFixed(3)} (dashed), skew ${skewness(d) >= 0 ? '+' : ''}${skewness(d).toFixed(2)} · haplo weight = 1 − depth (exact mirror)`, font: { size: 12 }, x: 0 }, xaxis: { title: { text: 'per-cell depth-modality weight', font: { size: 11 } } }, yaxis: { title: { text: 'cells', font: { size: 11 } } }, paper_bgcolor: 'white', plot_bgcolor: 'white', hovermode: presentation ? false : 'closest' } }
  }
  const tercileFig = () => {
    const at = cells.map((c) => c.snp_coverage); const s = [...at].sort((a, b) => a - b); const q1 = s[Math.floor(s.length / 3)], q2 = s[Math.floor((2 * s.length) / 3)]
    const groups = [['low AT', (c: IntegrationCell) => c.snp_coverage <= q1], ['mid AT', (c: IntegrationCell) => c.snp_coverage > q1 && c.snp_coverage <= q2], ['high AT', (c: IntegrationCell) => c.snp_coverage > q2]] as const
    const data = groups.map(([name, pred]) => { const v = cells.filter(pred).map((c) => c.wnn_depth_weight); return { type: 'box', name: `${name} (n=${v.length}, skew ${skewness(v) >= 0 ? '+' : ''}${skewness(v).toFixed(2)})`, y: v, boxpoints: false, fillcolor: '#d9f0d3', line: { color: '#334155' } } })
    return { data, layout: { height: 300, margin: { l: 50, r: 10, t: 36, b: 60 }, title: { text: 'depth weight by HA+HB coverage tercile', font: { size: 12 }, x: 0 }, yaxis: { title: { text: 'depth-modality weight', font: { size: 11 } } }, shapes: [{ type: 'line', x0: 0, x1: 1, xref: 'paper', y0: 0.5, y1: 0.5, line: { color: 'grey', dash: 'dot', width: 1 } }], showlegend: false, paper_bgcolor: 'white', plot_bgcolor: 'white', hovermode: presentation ? false : 'closest' } }
  }
  const chr8Strip = () => {
    const jit = jitter(cells.length, 42, 0.32)
    const data: any[] = levels.map((lv, li) => { const idx = cells.map((c, i) => (clusterOf(c, ig.clusterKey) === lv ? i : -1)).filter((i) => i >= 0); return { type: 'scatter', mode: 'markers', name: `${lv} (n=${idx.length})`, x: idx.map((i) => li + jit[i]), y: idx.map((i) => cells[i].chr8_score), text: idx.map((i) => cells[i].cell), marker: { size: 3.5, color: clusterColor(lv, levels), opacity: 0.5 }, hovertemplate: '%{text}<br>score %{y:.3f}<extra>' + lv + '</extra>' } })
    const meds = levels.map((lv, li) => { const v = cells.filter((c) => clusterOf(c, ig.clusterKey) === lv).map((c) => c.chr8_score).filter((x): x is number => x !== null && x !== undefined && Number.isFinite(x)).sort((a, b) => a - b); return { li, med: v.length ? v[v.length >> 1] : NaN } })
    const shapes = meds.map((m) => ({ type: 'line', x0: m.li - 0.4, x1: m.li + 0.4, y0: m.med, y1: m.med, line: { color: 'black', width: 3 } }))
    shapes.push({ type: 'line', x0: -0.5, x1: levels.length - 0.5, y0: 0.5, y1: 0.5, line: { color: 'grey', width: 1, dash: 'dash' } } as any)
    return { data, layout: { height: 340, margin: { l: 55, r: 10, t: 36, b: 60 }, title: { text: `chr8 log2 regional depth ratio by ${ig.clusterKey} · black = cluster median; dashed 0.5 = inspection reference, not a calling threshold`, font: { size: 12 }, x: 0 }, xaxis: { tickvals: levels.map((_, i) => i), ticktext: levels.map((lv, i) => `${lv}<br>med ${Number.isFinite(meds[i].med) ? meds[i].med.toFixed(2) : 'N/A'}`), title: { text: 'cluster (jitter seed 42, display only)', font: { size: 11 } } }, yaxis: { title: { text: 'log2(mean X 46–120 Mb / mean X 10–35 Mb)', font: { size: 11 } } }, shapes, showlegend: false, paper_bgcolor: 'white', plot_bgcolor: 'white', hovermode: presentation ? false : 'closest' } }
  }
  const chr8Scatter = () => {
    const baseSet = new Set(chr8.members.map((c) => c.cell)), addSet = new Set(chr8.added.map((c) => c.cell))
    const grp = (pred: (c: IntegrationCell) => boolean, name: string, color: string, size: number) => { const m = cells.filter(pred); return { type: 'scatter', mode: 'markers', name: `${name} (n=${m.length})`, x: m.map((c) => c.total_raw_depth / 1000), y: m.map((c) => c.chr8_score), text: m.map((c) => c.cell), marker: { size, color, opacity: name === 'other cells' ? 0.35 : 0.9 }, hovertemplate: '%{text}<br>raw %{x:.1f}k · score %{y:.3f}<extra>' + name + '</extra>' } }
    const data = [grp((c) => !baseSet.has(c.cell) && !addSet.has(c.cell), 'other cells', '#cbd5e1', 3), grp((c) => addSet.has(c.cell), `added in ${ig.compareKey} group ${chr8.target}`, '#e6550d', 5), grp((c) => baseSet.has(c.cell), `baseline cluster ${chr8.base}`, '#9467bd', 5)]
    const rho = spearman(cells.filter((c) => c.chr8_score != null).map((c) => c.total_raw_depth), cells.filter((c) => c.chr8_score != null).map((c) => c.chr8_score as number))
    return { data, layout: { height: 340, margin: { l: 55, r: 10, t: 36, b: 44 }, title: { text: `score vs raw counts · Spearman ρ = ${rho.toFixed(3)} (all cells)`, font: { size: 12 }, x: 0 }, xaxis: { title: { text: 'raw count total per cell (thousands)', font: { size: 11 } } }, yaxis: { title: { text: 'chr8 score', font: { size: 11 } } }, shapes: [{ type: 'line', x0: 0, x1: 1, xref: 'paper', y0: 0.5, y1: 0.5, line: { color: 'grey', dash: 'dash', width: 1 } }], legend: { orientation: 'h', y: -0.2, font: { size: 10 } }, paper_bgcolor: 'white', plot_bgcolor: 'white', hovermode: presentation ? false : 'closest' } }
  }
  const allKeys = [...meta.cluster_keys.wnn, ...meta.cluster_keys.depth, ...meta.cluster_keys.haplo]

  return (
    <div className="integ">
      <div className="section-heading"><div><div className="eyebrow">INTEGRATION · PART 1</div><h2>Cells, modalities and clusters</h2><p>Depth (log1p X) and phased allelic imbalance are clustered separately, fused with a WNN graph, and Leiden-clustered. Selections link the three embeddings; a selected cluster carries over to the segmentation tab. Nothing here recomputes clustering.</p></div>
        <div className="small muted">h5mu: {meta.h5mu}<br />{meta.notes.map((n, i) => <div key={i}>{n}</div>)}{meta.problems.map((n, i) => <div key={i} className="warn">⚠ {n}</div>)}</div></div>

      <div className="panel">
        <div className="controls-row">
          <label>cluster labelling<select value={ig.clusterKey} onChange={(e) => { ig.set({ clusterKey: e.target.value, selectedCluster: null, highlight: [] }) }}>
            {(['wnn', 'depth', 'haplo'] as const).map((m) => <optgroup key={m} label={m === 'wnn' ? 'integrated WNN' : m === 'depth' ? 'depth-only' : 'allelic-only'}>{meta.cluster_keys[m].map((k) => <option key={k} value={k}>{k}</option>)}</optgroup>)}</select></label>
          <label>colour by<select value={ig.colorBy} onChange={(e) => ig.set({ colorBy: e.target.value as ColorBy })}>{COLOR_OPTIONS.map((o) => <option key={o.id} value={o.id}>{o.label}</option>)}</select></label>
          <label>drag<select value={drag} onChange={(e) => setDrag(e.target.value as any)}><option value="lasso">lasso = select cells</option><option value="zoom">zoom</option><option value="pan">pan</option></select></label>
          <label>point size<input type="number" min={1} max={8} value={ig.pointSize} onChange={(e) => ig.set({ pointSize: +e.target.value })} style={{ width: 60 }} /></label>
          <label className="chk-inline"><input type="checkbox" checked={ig.dimOthers} onChange={(e) => ig.set({ dimOthers: e.target.checked })} />dim non-highlighted</label>
          <span className="spacer" />
          <span className="small muted">{hl.size ? `${hl.size} cells highlighted` : 'no highlight'}{ig.selectedCluster ? ` · selected cluster ${ig.selectedCluster}` : ''}</span>
          <button className="btn-sm" onClick={() => ig.set({ highlight: [], selectedCluster: null })} disabled={!hl.size && !ig.selectedCluster}>clear</button>
        </div>
        <div className="legend-chips">
          {levels.map((lv) => { const s = clusterSummary(lv); return <button key={lv} className={`chip ${ig.selectedCluster === lv ? 'on' : ''}`} title={`median raw ${s.raw.toLocaleString()} · median ΣX ${s.x.toFixed(0)} · median HA+HB ${s.at.toFixed(0)} · median depth weight ${s.dw.toFixed(3)}`} onClick={() => { if (ig.selectedCluster === lv) ig.set({ selectedCluster: null, highlight: [] }); else ig.set({ selectedCluster: lv, highlight: cells.filter((c) => clusterOf(c, ig.clusterKey) === lv).map((c) => c.cell) }) }}><i style={{ background: clusterColor(lv, levels) }} />{lv} · n={s.n}</button> })}
          <span className="small muted">click a chip to select that cluster (highlights its members in every embedding and selects it for the segmentation tab)</span>
        </div>
        <div className="emb-grid">{EMB.map((e) => { const f = embFig(e); return <Plot key={e.key} data={f.data} layout={f.layout} config={{ modeBarButtonsToRemove: ['autoScale2d'], displayModeBar: presentation ? false : 'hover' }} onSelected={onSelected} onClick={onClick} onDoubleClick={() => ig.set({ highlight: [] })} /> })}</div>
        <p className="lead">Embeddings from different runs have arbitrary orientation; cells are linked by ID, not by position. UMAP separation is a layout property, not validation of a subclone or CN event.</p>
        <CellCard />
      </div>

      <details className="diag" open>
        <summary>Coverage by cluster (fig 1b)</summary>
        <p className="lead">Coverage differs systematically across clusters; any cluster-specific depth contrast must be read with this in mind. This is a diagnostic, not proof that clusters are defined by coverage.</p>
        <div className="controls-row"><label>grouping<select value={covKey} onChange={(e) => setCovKey(e.target.value)}>{allKeys.map((k) => <option key={k}>{k}</option>)}</select></label>
          <label>metric<select value={covMetric} onChange={(e) => setCovMetric(e.target.value as any)}><option value="total_raw_depth">raw ATAC counts</option><option value="gc_norm_depth">GC-corrected normalised ΣX</option><option value="snp_coverage">HA+HB allelic counts</option><option value="n_bins_snp_covered">bins with SNP coverage</option></select></label></div>
        {(() => { const f = covFig(); return <Plot data={f.data} layout={f.layout} /> })()}
      </details>
      <details className="diag">
        <summary>PC1 and WNN spectral component vs coverage (fig 1a)</summary>
        <p className="lead">PC1 is a per-cell coordinate along the largest-variance axis (sign arbitrary, so |r| matters); the WNN "spectral component 1" is the first non-trivial diffusion component of the joint graph, because WNN produces a graph rather than a joint PCA. {cells.every((c) => c.wnn_spectral1 == null) ? 'Spectral component not available in this run.' : ''}</p>
        <label className="chk-inline"><input type="checkbox" checked={colorScatter} onChange={(e) => setColorScatter(e.target.checked)} />colour points by {ig.clusterKey}</label>
        <div className="two-col">
          {(() => { const f = scatterFig((c) => c.depth_pc1, (c) => c.gc_norm_depth, 'read-depth-modality PC1', 'ΣX (GC-corrected)'); return <Plot data={f.data} layout={f.layout} /> })()}
          {(() => { const f = scatterFig((c) => c.haplo_pc1, (c) => c.snp_coverage, 'allelic-modality PC1', 'HA+HB'); return <Plot data={f.data} layout={f.layout} /> })()}
          {(() => { const f = scatterFig((c) => c.wnn_spectral1 ?? null, (c) => c.gc_norm_depth, 'WNN spectral component 1', 'ΣX (GC-corrected)'); return <Plot data={f.data} layout={f.layout} /> })()}
          {(() => { const f = scatterFig((c) => c.wnn_spectral1 ?? null, (c) => c.snp_coverage, 'WNN spectral component 1', 'HA+HB'); return <Plot data={f.data} layout={f.layout} /> })()}
        </div>
      </details>
      <details className="diag">
        <summary>WNN modality weights (fig 3, 3b)</summary>
        <p className="lead">Per-cell integration weights (depth + haplo = 1): how much of a cell's fused neighbourhood is carried by each graph. Not a coverage measure. Because the two weights sum to one, their histograms are exact mirrors.</p>
        <div className="two-col">
          {(() => { const f = weightsFig(); return <Plot data={f.data} layout={f.layout} /> })()}
          {(() => { const f = tercileFig(); return <Plot data={f.data} layout={f.layout} /> })()}
          {(() => { const f = scatterFig((c) => c.gc_norm_depth, (c) => c.wnn_depth_weight, 'ΣX (GC-corrected)', 'depth-modality weight', true); return <Plot data={f.data} layout={f.layout} /> })()}
          {(() => { const f = scatterFig((c) => c.snp_coverage, (c) => c.wnn_depth_weight, 'HA+HB', 'depth-modality weight', true); return <Plot data={f.data} layout={f.layout} /> })()}
        </div>
      </details>
      <details className="diag">
        <summary>chr8 regional depth diagnostic and resolution comparison (fig 5)</summary>
        <p className="lead">Score = log2(mean X at chr8 {meta.chr8_score.num_mb[0]}–{meta.chr8_score.num_mb[1]} Mb / mean X at {meta.chr8_score.den_mb[0]}–{meta.chr8_score.den_mb[1]} Mb), computed from the explorer's X for every cell. Shared chr8 depth alone is not clonal identity.</p>
        <div className="controls-row"><label>compare against<select value={ig.compareKey} onChange={(e) => ig.set({ compareKey: e.target.value })}>{allKeys.filter((k) => k !== ig.clusterKey).map((k) => <option key={k}>{k}</option>)}</select></label>
          <span className="small">baseline cluster <b>{chr8.base}</b> ({ig.clusterKey}, n={chr8.members.length}) → largest overlapping {ig.compareKey} group <b>{chr8.target ?? 'N/A'}</b> keeps {chr8.retained} of them and adds <b>{chr8.added.length}</b> cells: {[...chr8.origin.entries()].map(([o, v]) => `${v.n} from cluster ${o} (${v.hi} with score > 0.5)`).join(', ') || 'none'}</span></div>
        <div className="two-col">
          {(() => { const f = chr8Strip(); return <Plot data={f.data} layout={f.layout} onClick={onClick} /> })()}
          {(() => { const f = chr8Scatter(); return <Plot data={f.data} layout={f.layout} onClick={onClick} /> })()}
        </div>
      </details>
      <details className="diag">
        <summary>Run provenance</summary>
        <table className="kv"><tbody>
          <tr><td>neighbours</td><td>{JSON.stringify(meta.neighbors_params)}</td></tr>
          <tr><td>cluster sizes ({ig.clusterKey})</td><td>{levels.map((lv) => `${lv}: ${cells.filter((c) => clusterOf(c, ig.clusterKey) === lv).length}`).join(' · ')}</td></tr>
          <tr><td>prototype manifest</td><td>{meta.manifest?.description} · created {meta.manifest?.created} · chromosomes {(meta.pilot_chromosomes || []).join(', ')} · max small segment {meta.manifest?.params?.max_small_bins} bins · τ {meta.manifest?.params?.tau_default} · stability {meta.manifest?.params?.stab_default} · bootstrap {meta.manifest?.params?.n_bootstrap} (seed {meta.manifest?.params?.seed})</td></tr>
          <tr><td>identity</td><td><code>{meta.identity.hash}</code> · {meta.identity.files.map((f: any) => f.path.split('/').slice(-1)[0]).join(', ')}</td></tr>
        </tbody></table>
      </details>
    </div>
  )
}
export { boxStats as _bs }
