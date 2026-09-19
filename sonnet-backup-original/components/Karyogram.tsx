import { useEffect, useMemo, useRef, useState } from 'react'
import { CN_COLORS, CN_LABELS, TERCILE_COLORS } from '../colors'
import { MB } from '../coords'
import { useStore } from '../store'

interface Props { rows: number[]; ruleMatch: Set<number>; onOpen: (cell: string, chrom: string) => void }

const CN_RGB = [0, 0.5, 1, 1.5, 2].map((v) => { const h = CN_COLORS[v].slice(1); return [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16), parseInt(h.slice(4, 6), 16)] })
const LEFT_W = 26, TOP_H = 22

export default function Karyogram({ rows, ruleMatch, onOpen }: Props) {
  const meta = useStore((s) => s.meta)!
  const karyo = useStore((s) => s.karyo)!
  const bins = useStore((s) => s.bins)!
  const view = useStore((s) => s.karyoView)
  const setView = useStore((s) => s.setKaryoView)
  const selectedCell = useStore((s) => s.cell)
  const presentation = useStore((s) => s.presentation)
  const wrap = useRef<HTMLDivElement>(null)
  const canvas = useRef<HTMLCanvasElement>(null)
  const [size, setSize] = useState({ w: 900, h: 480 })
  const [hover, setHover] = useState<{ x: number; y: number; text: string[] } | null>(null)
  const drag = useRef<{ x0: number; binLo: number; binHi: number; moved: boolean } | null>(null)
  const nBins = meta.dataset.identity.n_bins
  const binLo = view?.binLo ?? 0, binHi = view?.binHi ?? nBins
  const rowH = Math.max(2, Math.floor((size.h - TOP_H) / Math.max(1, rows.length)))
  const plotH = rowH * rows.length
  const plotW = size.w - LEFT_W

  useEffect(() => {
    const el = wrap.current; if (!el) return
    const ro = new ResizeObserver(() => { if (el.clientWidth > 100) setSize({ w: el.clientWidth, h: Math.min(720, Math.max(320, Math.round(el.clientWidth * 0.5))) }) })
    ro.observe(el); return () => ro.disconnect()
  }, [])

  const chromSpans = useMemo(() => meta.chromosomes.map((c) => ({ name: c.name, lo: c.offset, hi: c.offset + c.n })), [meta])
  const pxToBin = (px: number) => binLo + Math.floor(((px - LEFT_W) / plotW) * (binHi - binLo))
  const binToPx = (b: number) => LEFT_W + ((b - binLo) / (binHi - binLo)) * plotW

  useEffect(() => {
    const cv = canvas.current; if (!cv || plotW < 10 || rows.length === 0) return
    cv.width = size.w; cv.height = TOP_H + plotH
    const ctx = cv.getContext('2d')!
    ctx.fillStyle = '#ffffff'; ctx.fillRect(0, 0, cv.width, cv.height)
    const img = ctx.createImageData(plotW, plotH)
    const d = img.data
    const span = binHi - binLo
    // nearest-bin per pixel column; categorical values are never averaged
    const colBin = new Int32Array(plotW)
    for (let px = 0; px < plotW; px++) colBin[px] = Math.min(nBins - 1, binLo + Math.floor((px / plotW) * span))
    for (let r = 0; r < rows.length; r++) {
      const base = rows[r] * nBins
      for (let py = 0; py < rowH; py++) {
        const rowOff = ((r * rowH + py) * plotW) * 4
        for (let px = 0; px < plotW; px++) {
          const v = karyo[base + colBin[px]]
          const c = CN_RGB[v]
          const o = rowOff + px * 4
          d[o] = c[0]; d[o + 1] = c[1]; d[o + 2] = c[2]; d[o + 3] = 255
        }
      }
    }
    ctx.putImageData(img, LEFT_W, TOP_H)
    // annotation strips: tercile (left 12px) and rule match (next 12px)
    for (let r = 0; r < rows.length; r++) {
      const c = meta.cells[rows[r]]
      ctx.fillStyle = TERCILE_COLORS[c.tercile]; ctx.fillRect(0, TOP_H + r * rowH, 11, rowH)
      ctx.fillStyle = ruleMatch.has(rows[r]) ? '#f59e0b' : '#f1f5f9'; ctx.fillRect(13, TOP_H + r * rowH, 11, rowH)
    }
    // chromosome separators + labels
    ctx.font = '10px Inter, system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.fillStyle = '#334155'
    for (const c of chromSpans) {
      if (c.hi <= binLo || c.lo >= binHi) continue
      const x0 = Math.max(LEFT_W, binToPx(c.lo)), x1 = Math.min(size.w, binToPx(c.hi))
      ctx.strokeStyle = 'rgba(15,23,42,0.35)'; ctx.lineWidth = 1
      ctx.beginPath(); ctx.moveTo(x0 + 0.5, TOP_H); ctx.lineTo(x0 + 0.5, TOP_H + plotH); ctx.stroke()
      if (x1 - x0 > 18) ctx.fillText(c.name.replace('chr', ''), (x0 + x1) / 2, 15)
    }
    // selected cell outline
    if (selectedCell) {
      const r = rows.indexOf(meta.cells.findIndex((c) => c.cellID === selectedCell))
      if (r >= 0) { ctx.strokeStyle = '#111827'; ctx.lineWidth = 1.5; ctx.strokeRect(LEFT_W + 0.5, TOP_H + r * rowH + 0.5, plotW - 1, Math.max(1, rowH - 1)) }
    }
  }, [rows, karyo, size, binLo, binHi, rowH, plotH, plotW, ruleMatch, selectedCell, meta, chromSpans, nBins])

  const locate = (e: React.MouseEvent) => {
    const rect = canvas.current!.getBoundingClientRect()
    const x = e.clientX - rect.left, y = e.clientY - rect.top
    if (x < LEFT_W || y < TOP_H) return null
    const r = Math.floor((y - TOP_H) / rowH)
    if (r < 0 || r >= rows.length) return null
    const b = Math.min(nBins - 1, Math.max(0, pxToBin(x)))
    const chrom = chromSpans.find((c) => b >= c.lo && b < c.hi)!
    return { x, y, r, b, chrom }
  }
  const onMove = (e: React.MouseEvent) => {
    if (drag.current) {
      const dx = e.clientX - drag.current.x0
      if (Math.abs(dx) > 3) drag.current.moved = true
      const span = drag.current.binHi - drag.current.binLo
      let lo = Math.round(drag.current.binLo - (dx / plotW) * span)
      lo = Math.max(0, Math.min(nBins - span, lo))
      setView({ binLo: lo, binHi: lo + span })
      return
    }
    const l = locate(e)
    if (!l || presentation) { setHover(null); return }
    const cell = meta.cells[rows[l.r]]
    const v = karyo[rows[l.r] * nBins + l.b] / 2
    const s = bins[2 * l.b] / MB, en = bins[2 * l.b + 1] / MB
    setHover({ x: l.x, y: l.y, text: [`${cell.cellID}  (row ${l.r + 1} of ${rows.length})`, `${l.chrom.name}:${s.toFixed(2)}–${en.toFixed(2)} Mb · chrom bin ${l.b - l.chrom.lo} · genome bin ${l.b}`, `CN ${v} = ${CN_LABELS[v]}`, `coverage ${cell.raw_coverage.toLocaleString()} (${cell.tercile}) · ${cell.n_segments} segments · best_s ${cell.best_s.toFixed(3)}`] })
  }
  const onDown = (e: React.MouseEvent) => { drag.current = { x0: e.clientX, binLo, binHi, moved: false } }
  const onUp = (e: React.MouseEvent) => {
    const d = drag.current; drag.current = null
    if (d && d.moved) return
    const l = locate(e)
    if (l) onOpen(meta.cells[rows[l.r]].cellID, l.chrom.name)
  }
  const onWheel = (e: React.WheelEvent) => {
    const rect = canvas.current!.getBoundingClientRect()
    const x = e.clientX - rect.left
    if (x < LEFT_W) return
    const at = pxToBin(x)
    const f = e.deltaY < 0 ? 0.8 : 1.25
    let span = Math.max(50, Math.min(nBins, Math.round((binHi - binLo) * f)))
    let lo = Math.round(at - (at - binLo) * (span / (binHi - binLo)))
    lo = Math.max(0, Math.min(nBins - span, lo))
    setView(span >= nBins ? null : { binLo: lo, binHi: lo + span })
  }
  useEffect(() => { const cv = canvas.current; if (!cv) return; const h = (e: WheelEvent) => e.preventDefault(); cv.addEventListener('wheel', h, { passive: false }); return () => cv.removeEventListener('wheel', h) }, [])

  return (
    <div ref={wrap} className="karyo-wrap">
      <div className="karyo-toolbar">
        <span className="muted">{rows.length} of {meta.cells.length} cells · {binHi - binLo === nBins ? 'whole genome' : `bins ${binLo}–${binHi}`}</span>
        <span className="spacer" />
        <span className="muted">wheel = zoom · drag = pan · click = open cell/chromosome</span>
        {meta.chromosomes.map((c) => <button key={c.name} className="btn-xs" onClick={() => setView({ binLo: c.offset, binHi: c.offset + c.n })}>{c.name.replace('chr', '')}</button>)}
        <button className="btn-xs" onClick={() => setView(null)}>reset</button>
      </div>
      <div style={{ position: 'relative' }}>
        <canvas ref={canvas} style={{ display: 'block', width: size.w, height: TOP_H + plotH, cursor: 'crosshair' }}
          onMouseMove={onMove} onMouseLeave={() => { setHover(null); drag.current = null }} onMouseDown={onDown} onMouseUp={onUp} onWheel={onWheel} onDoubleClick={() => setView(null)} />
        {hover && <div className="tooltip" style={{ left: Math.min(hover.x + 14, size.w - 330), top: hover.y + 12 }}>{hover.text.map((t, i) => <div key={i}>{t}</div>)}</div>}
      </div>
      <div className="karyo-legend">
        {[0, 0.5, 1, 1.5, 2].map((v) => <span key={v}><i style={{ background: CN_COLORS[v] }} /> {v} = {CN_LABELS[v]}</span>)}
        <span className="sep" />
        {(['low', 'mid', 'high'] as const).map((t) => <span key={t}><i style={{ background: TERCILE_COLORS[t] }} /> {t} coverage</span>)}
        <span><i style={{ background: '#f59e0b' }} /> matches enabled event rules</span>
      </div>
    </div>
  )
}
