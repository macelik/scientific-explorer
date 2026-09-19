import Plot from './Plot'
import { useState } from 'react'
import type { ExperimentResult } from '../experimentStore'
import { useStore } from '../store'

import { experimentOverlays, nodeDepth, ORIGINAL, EDITED } from '../experimentOverlays'
export { ORIGINAL, EDITED } from '../experimentOverlays'

function line(starts: number[], ends: number[], ys: number[]) {
  const x: (number | null)[] = [], y: (number | null)[] = []
  starts.forEach((v, i) => {
    if (i && v > ends[i - 1] + 1) { x.push(null); y.push(null) }
    x.push(v / 1e6); y.push(ys[i])
  })
  return { x, y }
}

export function ExperimentTracks({ result, starts, ends, original, edited, start, end, onRegion }: {
  result?: ExperimentResult | null; starts: number[]; ends: number[]; original: number[]; edited: number[]
  start: number; end: number; onRegion?: (a: number, b: number) => void
}) {
  const presentation = useStore(s => s.presentation)
  const hasAD = !!result
  const data: any[] = []
  const add = (name: string, values: number[], color: string, axis: string, a = starts, b = ends, dash?: string) =>
    data.push({ ...line(a, b, values), name, type: 'scatter', mode: 'lines', yaxis: axis, line: { color, width: 1.5, dash }, hovertemplate: '%{x:.3f} Mb · %{y:.4g}<extra>%{fullData.name}</extra>' })
  if (result) {
    add('Original AD', result.original.nodes[0].ad, ORIGINAL, 'y', starts.slice(1), ends.slice(1))
    add('Edited AD', result.edited.nodes[0].ad, EDITED, 'y', starts.slice(1), ends.slice(1))
  }
  add('Original X', original, ORIGINAL, hasAD ? 'y2' : 'y')
  add('Edited X', edited, EDITED, hasAD ? 'y2' : 'y', starts, ends, 'dot')
  add('Change in X', edited.map((v, i) => v - original[i]), EDITED, hasAD ? 'y3' : 'y2')
  const axis = { gridcolor: '#e8edf0', zerolinecolor: '#bcc7ce', fixedrange: false, automargin: true }
  const shapes: any[] = result ? experimentOverlays(result, 'root') : start >= 0 && end <= starts.length && start < end ? [{ type: 'rect', x0: starts[start] / 1e6, x1: ends[end-1] / 1e6, y0: 0, y1: 1, yref: 'paper', fillcolor: EDITED, opacity: 0.12, line: { width: 0 }, layer:'below' }] : []
  return <Plot data={data} layout={{ height: hasAD ? 510 : 345, margin: { l: 64, r: 22, t: 40, b: 45 }, paper_bgcolor: '#fff', plot_bgcolor: '#fff', font: { family: 'system-ui', size: 12, color: '#36444d' },
    xaxis: { ...axis, title: 'Genomic position (Mb)', anchor: hasAD ? 'y3' : 'y2', showspikes: true, spikemode: 'across', spikesnap: 'cursor' },
    yaxis: { ...axis, title: hasAD ? 'AD statistic' : 'X', domain: hasAD ? [0.68, 1] : [0.43, 1] },
    yaxis2: { ...axis, title: hasAD ? 'X' : 'Δ X', domain: hasAD ? [0.31, 0.6] : [0, 0.29], anchor: 'x' },
    ...(hasAD ? { yaxis3: { ...axis, title: 'Δ X', domain: [0, 0.23], anchor: 'x' } } : {}),
    legend: { orientation: 'h', x: 0, y: 1.12 }, hovermode: presentation ? false : 'x unified', dragmode: onRegion ? 'select' : 'zoom', selectdirection: 'h',
    shapes, uirevision: result ? `${result.request.cell}|${result.request.chrom}` : starts[0],
  }} config={{ displayModeBar: !presentation }} onSelected={e => {
    if (!e?.range?.x || !onRegion) return
    const lo = Math.min(...e.range.x)*1e6, hi = Math.max(...e.range.x)*1e6
    const members = starts.map((_, i) => i).filter(i => ends[i] >= lo && starts[i] <= hi)
    if (members.length) onRegion(members[0], members[members.length-1]+1)
  }} />
}

export function ChildExperimentPlots({ result }: { result: ExperimentResult }) {
  const presentation = useStore(s => s.presentation)
  const [selectedDepth, setSelectedDepth] = useState(1)
  const nodes = [...result.original.nodes, ...result.edited.nodes]
  const depths = [...new Set(nodes.map(nodeDepth).filter(d => d > 0))].sort((a,b) => a-b)
  const depth = depths.includes(selectedDepth) ? selectedDepth : depths[0]
  const sides = [...new Set(nodes.filter(n => nodeDepth(n) === depth).map(n => n.side))].sort()
  if (!depths.length) return <p className="muted">No child nodes were tested at this recursion limit. Increase k and rerun to inspect children.</p>
  return <div className="experiment-children"><label>Inspection depth <select aria-label="Child inspection depth" value={depth} onChange={e=>setSelectedDepth(+e.target.value)}>{depths.map(d=><option key={d} value={d}>{d}{d===1?' · immediate children':''}</option>)}</select></label>
    <p className="muted">Orange shading: modified interval, clipped to the displayed nodes. Solid lines: accepted breakpoints in each node, including accepted descendants. Dashed lines: unaccepted candidates. Blue = original; brown = edited. The shaded interval marks the edit target; total-preserving normalization can also change X elsewhere.</p>
    <div className="two-col">{sides.map(side => {
    const data = [result.original, result.edited].flatMap((tree, i) => {
      const node = tree.nodes.find(v => v.side === side)
      if (!node) return []
      return [{ ...line(result.start_bp.slice(node.start+1, node.end), result.end_bp.slice(node.start+1, node.end), node.ad), type: 'scatter', mode: 'lines',
        name: `${i ? 'Edited' : 'Original'} [${node.start}, ${node.end})${node.hypothetical?' · hypothetical':''}`, line: { color: i ? EDITED : ORIGINAL, width: 1.5 },
        hovertemplate: '%{x:.3f} Mb · AD %{y:.4g}<extra>%{fullData.name}</extra>' }]
    })
    return <section key={side} data-node-path={side}><h3>{side.split('.').map(s=>s[0].toUpperCase()+s.slice(1)).join(' → ')} child AD · depth {depth}</h3><Plot data={data} layout={{ height: 280, margin: { l: 50, r: 15, b: 45, t: 40 },
      xaxis: { title: 'Genomic position (Mb)', gridcolor: '#e8edf0' }, yaxis: { title: 'AD', gridcolor: '#e8edf0', rangemode: 'tozero' },
      font: { family: 'system-ui', size: 11 }, legend: { orientation: 'h', y: 1.18 }, hovermode: presentation ? false : 'x unified',
      shapes:experimentOverlays(result,side),
    }} config={{ displayModeBar: !presentation }} /></section>
  })}</div></div>
}
