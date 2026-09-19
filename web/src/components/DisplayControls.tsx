import { useStore } from '../store'

export default function DisplayControls({ manual }: { manual?: boolean }) {
  const display = useStore((s) => s.display)
  const setDisplay = useStore((s) => s.setDisplay)
  const mode = useStore((s) => s.mode)
  const setMode = useStore((s) => s.setMode)
  const meta = useStore((s) => s.meta)!
  const chrom = useStore((s) => s.chrom)
  const regions = meta.regions.filter((r) => r.chrom === chrom)
  return (
    <div className="card controls">
      <div className="row wrap">
        <span className="seg">
          <button className={mode === 'navigate' ? 'on' : ''} onClick={() => setMode('navigate')} title="drag = zoom, double-click = reset; never adds a window">Navigate</button>
          <button className={mode === 'select' ? 'on' : ''} onClick={() => setMode('select')} title="drag on any track = draw a window; click a point on the AD curve = 40-bin window around that boundary">Select window</button>
          {manual && <button className={mode === 'place' ? 'on' : ''} onClick={() => setMode('place')} title="click a boundary on the active node's AD curve to insert a manual split">Place breakpoint</button>}
        </span>
        <label>baseline <select value={display.baseline} onChange={(e) => setDisplay({ baseline: e.target.value as any })}><option value="trimmed">figure-compatible: trimmed_mean_iqr(lb=False)</option><option value="median">literal median (zeros kept)</option></select></label>
        <label>rolling mean <input type="range" min={5} max={75} step={2} value={display.rolling} onChange={(e) => setDisplay({ rolling: +e.target.value, minPeriods: Math.min(display.minPeriods, +e.target.value) })} /> {display.rolling} bins (min {display.minPeriods}) <span className="muted" title="Smoothing is visual only: it never changes AD, selection statistics, tests or CN assignment.">ⓘ visual only</span></label>
        <label className="chk-inline"><input type="checkbox" checked={display.showGenome} onChange={(e) => setDisplay({ showGenome: e.target.checked })} />vs genome</label>
        <label className="chk-inline"><input type="checkbox" checked={display.showLocal} onChange={(e) => setDisplay({ showLocal: e.target.checked })} />vs node</label>
        <label className="chk-inline"><input type="checkbox" checked={display.showMean} onChange={(e) => setDisplay({ showMean: e.target.checked })} />segment mean</label>
        <label className="chk-inline"><input type="checkbox" checked={display.showMedian} onChange={(e) => setDisplay({ showMedian: e.target.checked })} />segment median</label>
        <label className="chk-inline"><input type="checkbox" checked={display.sharedAdScale} onChange={(e) => setDisplay({ sharedAdScale: e.target.checked })} />shared AD y-scale</label>
        <label className="chk-inline"><input type="checkbox" checked={display.showProduction} onChange={(e) => setDisplay({ showProduction: e.target.checked })} />production boundaries</label>
        <label className="chk-inline"><input type="checkbox" checked={display.showRegions} onChange={(e) => setDisplay({ showRegions: e.target.checked })} />regions</label>
        <label>active event region <select value={display.activeRegionId ?? ''} onChange={(e) => setDisplay({ activeRegionId: e.target.value || null })}><option value="">none</option>{regions.map((r) => <option key={r.id} value={r.id}>{r.label} ({r.start_mb}–{r.end_mb} Mb)</option>)}</select></label>
      </div>
    </div>
  )
}
