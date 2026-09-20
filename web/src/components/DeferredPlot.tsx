import { useEffect, useRef, useState } from 'react'
import Plot, { type PlotProps } from './Plot'

/** Reserve the chart's space, but initialize Plotly only when it is visible.
 * Once mounted it stays mounted, preserving zoom and selected points on scroll.
 */
export default function DeferredPlot(props: PlotProps) {
  const container = useRef<HTMLDivElement>(null)
  const [visible, setVisible] = useState(false)
  useEffect(() => {
    const el = container.current
    if (!el) return
    if (typeof IntersectionObserver === 'undefined') { setVisible(true); return }
    const observer = new IntersectionObserver(entries => {
      if (entries.some(entry => entry.isIntersecting)) {
        observer.disconnect()
        setVisible(true)
      }
    })
    observer.observe(el)
    return () => observer.disconnect()
  }, [])
  return <div ref={container} className="deferred-plot" style={{height:props.layout?.height, minWidth:0}}>
    {visible ? <Plot {...props} /> : <div className="muted small" style={{padding:12}}>Chart renders when scrolled into view.</div>}
  </div>
}
