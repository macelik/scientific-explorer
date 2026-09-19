import type { ExperimentNode, ExperimentResult } from './experimentStore'

export const ORIGINAL = '#2a6f97', EDITED = '#b25424'
export const nodeDepth = (node: ExperimentNode) => node.depth ?? (node.side === 'root' ? 0 : node.side.split('.').length)

/** All coordinates are chromosome-local, half-open bins until mapped to Mb. */
export function experimentOverlays(result: ExperimentResult, side: string): any[] {
  const shapes: any[] = [], intervals: [number, number][] = []
  const start = result.intervention.target_start ?? result.intervention.event_start
  const end = result.intervention.target_end ?? result.intervention.event_end
  for (const [tree, color, label] of [[result.original, ORIGINAL, 'Original'], [result.edited, EDITED, 'Edited']] as const) {
    const node = tree.nodes.find(n => n.side === side)
    if (!node) continue
    const a = Math.max(start, node.start), b = Math.min(end, node.end)
    if (a < b) intervals.push([a, b])
    const marker = (boundary: number, accepted: boolean) => {
      const x = result.start_bp[boundary] / 1e6
      shapes.push({type:'line', name:`${label} ${accepted ? 'accepted' : 'candidate'} · bin ${boundary}`,
        xref:'x', x0:x, x1:x, yref:'paper', y0:0, y1:1,
        line:{color, width:accepted ? (label === 'Original' ? 2.5 : 1.5) : 1, dash:accepted ? 'solid' : 'dash'}})
    }
    for (const boundary of tree.boundaries) if (node.start < boundary && boundary < node.end) marker(boundary, true)
    if (node.argmax_b !== null && !node.accepted) marker(node.start + node.argmax_b, false)
  }
  // Shade the union once: differing scenario bounds must not double-darken it,
  // extend the plot beyond its nodes, or shade an untouched donor interval.
  const merged: [number, number][] = []
  for (const interval of intervals.sort((a,b) => a[0]-b[0])) {
    const prev = merged[merged.length-1]
    if (prev && interval[0] <= prev[1]) prev[1] = Math.max(prev[1], interval[1])
    else merged.push([...interval])
  }
  for (const [a,b] of merged) shapes.unshift({type:'rect', name:'Modified region', xref:'x',
    x0:result.start_bp[a]/1e6, x1:result.end_bp[b-1]/1e6,
    yref:'paper', y0:0, y1:1, fillcolor:EDITED, opacity:0.12, line:{width:0}, layer:'below'})
  return shapes
}
