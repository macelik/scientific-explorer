import { useEffect, useRef } from 'react'
import Plotly from 'plotly.js-dist-min'

export interface PlotProps {
  data: any[]; layout: any; config?: any; style?: React.CSSProperties; className?: string
  onClick?: (e: any) => void; onSelected?: (e: any) => void; onRelayout?: (e: any) => void; onHover?: (e: any) => void; onDoubleClick?: () => void
  onReady?: (el: any) => void
  /** plain click (no drag) anywhere in the main plot area, in x data coordinates; works in every dragmode */
  onPlainClick?: (x: number, row: number) => void
}

/** Thin Plotly wrapper: Plotly.react on every prop change, event handlers rebound. */
export default function Plot(p: PlotProps) {
  const ref = useRef<HTMLDivElement>(null)
  const handlers = useRef(p)
  handlers.current = p
  useEffect(() => {
    const el = ref.current as any
    if (!el) return
    Plotly.react(el, p.data, p.layout, { responsive: true, displaylogo: false, scrollZoom: false, ...(p.config || {}) })
    // presentation mode: also drop any hover label that was still showing when the mode changed
    if (p.layout && p.layout.hovermode === false) { try { (Plotly as any).Fx.hover(el, []) } catch { /* ignore */ } }
    p.onReady?.(el)
  }, [p.data, p.layout, p.config])
  useEffect(() => {
    const el = ref.current as any
    if (!el) return
    const click = (e: any) => handlers.current.onClick?.(e)
    const sel = (e: any) => handlers.current.onSelected?.(e)
    const rel = (e: any) => handlers.current.onRelayout?.(e)
    const hov = (e: any) => handlers.current.onHover?.(e)
    const dbl = () => handlers.current.onDoubleClick?.()
    el.on('plotly_click', click); el.on('plotly_selected', sel); el.on('plotly_relayout', rel); el.on('plotly_hover', hov); el.on('plotly_doubleclick', dbl)
    return () => { try { el.removeAllListeners?.('plotly_click'); el.removeAllListeners?.('plotly_selected'); el.removeAllListeners?.('plotly_relayout'); el.removeAllListeners?.('plotly_hover'); el.removeAllListeners?.('plotly_doubleclick') } catch { /* ignore */ } }
  }, [])
  useEffect(() => {
    const el = ref.current as HTMLDivElement | null
    if (!el) return
    let down: { x: number; y: number } | null = null
    const onDown = (e: PointerEvent) => {
      const target = e.target as HTMLElement
      // only clicks that start on the main plot drag area (Plotly then covers the window with a dragcover)
      if (e.button === 0 && target.closest('.main-svg') && !target.closest('.modebar') && !target.closest('.legend')) down = { x: e.clientX, y: e.clientY }
    }
    const onUp = (e: PointerEvent) => {
      const d = down; down = null
      if (!d || Math.abs(e.clientX - d.x) > 3 || Math.abs(e.clientY - d.y) > 3) return
      const anyEl = el as any
      const fl = anyEl._fullLayout
      if (!fl || !fl.xaxis) return
      const bb = el.getBoundingClientRect()
      const px = e.clientX - bb.left - fl.margin.l
      const py = e.clientY - bb.top - fl.margin.t
      const x = fl.xaxis.p2d(px)
      // which y row (0 = top) by paper fraction
      const frac = 1 - py / (fl.height - fl.margin.t - fl.margin.b)
      const rows = ['yaxis', 'yaxis2', 'yaxis3'].map((k) => fl[k]).filter(Boolean)
      let row = 0
      rows.forEach((ya: any, i: number) => { if (ya.domain && frac >= ya.domain[0] && frac <= ya.domain[1]) row = i })
      if (Number.isFinite(x)) handlers.current.onPlainClick?.(x, row)
    }
    el.addEventListener('pointerdown', onDown, true)
    window.addEventListener('pointerup', onUp, true)
    return () => { el.removeEventListener('pointerdown', onDown, true); window.removeEventListener('pointerup', onUp, true) }
  }, [])
  useEffect(() => {
    const el = ref.current as any
    if (!el) return
    let last = el.clientWidth
    const ro = new ResizeObserver(() => { if (el.clientWidth > 0 && Math.abs(el.clientWidth - last) > 2) { last = el.clientWidth; try { Plotly.Plots.resize(el) } catch { /* ignore */ } } })
    ro.observe(el)
    return () => ro.disconnect()
  }, [])
  useEffect(() => () => { try { Plotly.purge(ref.current as any) } catch { /* ignore */ } }, [])
  return <div ref={ref} className={p.className} style={{ width: '100%', height: p.layout?.height, ...p.style }} />
}

export function downloadPlot(el: any, filename: string, format: 'png' | 'svg' = 'png') {
  return Plotly.downloadImage(el, { format, filename, width: el.clientWidth * 2, height: el.clientHeight * 2, scale: 1 } as any)
}
