import Plot from './Plot'
import type { ExperimentResult } from '../experimentStore'
import { useStore } from '../store'

export const ORIGINAL = '#2a6f97', EDITED = '#b25424'

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
  const shapes: any[] = start >= 0 && end <= starts.length && start < end ? [{ type: 'rect', x0: starts[start] / 1e6, x1: ends[end-1] / 1e6, y0: 0, y1: 1, yref: 'paper', fillcolor: EDITED, opacity: 0.07, line: { width: 0 } }] : []
  if (result) for (const [which, color] of [[result.original, ORIGINAL], [result.edited, EDITED]] as const) {
    const b = which.nodes[0].argmax_b
    if (b !== null) shapes.push({ type: 'line', x0: starts[b]/1e6, x1: starts[b]/1e6, y0: 0, y1: 1, yref: 'paper', line: { color, width: 1.2, dash: 'dash' } })
  }
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
  return <div className="two-col">{['left', 'right'].map(side => {
    const data = [result.original, result.edited].flatMap((tree, i) => {
      const node = tree.nodes.find(v => v.side === side)
      if (!node) return []
      return [{ ...line(result.start_bp.slice(node.start+1, node.end), result.end_bp.slice(node.start+1, node.end), node.ad), type: 'scatter', mode: 'lines',
        name: `${i ? 'Edited' : 'Original'} [${node.start}, ${node.end})`, line: { color: i ? EDITED : ORIGINAL, width: 1.5 },
        hovertemplate: '%{x:.3f} Mb · AD %{y:.4g}<extra>%{fullData.name}</extra>' }]
    })
    return <section key={side}><h3>{side === 'left' ? 'Left' : 'Right'} child AD</h3><Plot data={data} layout={{ height: 260, margin: { l: 50, r: 15, b: 45, t: 40 },
      xaxis: { title: 'Genomic position (Mb)', gridcolor: '#e8edf0' }, yaxis: { title: 'AD', gridcolor: '#e8edf0', rangemode: 'tozero' },
      font: { family: 'system-ui', size: 11 }, legend: { orientation: 'h', y: 1.18 }, hovermode: presentation ? false : 'x unified',
    }} config={{ displayModeBar: !presentation }} /></section>
  })}</div>
}
