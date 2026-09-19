import { useEffect, useMemo, useRef, useState } from 'react'
import { CN_COLORS, CN_LABELS } from '../colors'
import { MB } from '../coords'
import { STORY_CN, STORY_CN_LABELS } from '../integration'
import { useStore } from '../store'

export interface KaryoRow { label: string; band?: string; data: Uint8Array; sub?: string }
interface Props {
  rows: KaryoRow[]; palette: 'story' | 'explorer'; rowHeight: 'fit' | number; maxHeight?: number; labelWidth?: number
  onClick?: (rowIndex: number, bin: number, chrom: string) => void
  hoverText?: (rowIndex: number, bin: number, chrom: string, value: number | null) => string[]
  selectedRow?: number | null
}
const TOP_H = 20
const hex = (h: string) => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)]

/** Categorical CN karyogram on a canvas: nearest-bin columns and nearest-row sampling, values are never averaged. */
export default function ClusterKaryogram({ rows, palette, rowHeight, maxHeight = 720, labelWidth = 120, onClick, hoverText, selectedRow }: Props) {
  const meta = useStore((s) => s.meta)!
  const bins = useStore((s) => s.bins)!
  const presentation = useStore((s) => s.presentation)
  const wrap = useRef<HTMLDivElement>(null); const canvas = useRef<HTMLCanvasElement>(null)
  const [w, setW] = useState(900)
  const [hover, setHover] = useState<{ x: number; y: number; text: string[] } | null>(null)
  const nBins = meta.dataset.identity.n_bins
  const bandW = rows.some((r) => r.band) ? 14 : 0
  const left = labelWidth + bandW
  const plotW = Math.max(10, w - left)
  const fit = rowHeight === 'fit'
  const rowH = fit ? Math.max(1, Math.min(24, Math.floor(Math.min(maxHeight, Math.max(160, rows.length * 24)) / Math.max(1, rows.length)))) : rowHeight
  const nPixRows = fit ? Math.min(rows.length, Math.max(1, Math.floor(Math.min(maxHeight, rows.length * rowH) / rowH))) : rows.length
  const plotH = nPixRows * rowH
  const rowOfPix = (py: number) => Math.min(rows.length - 1, Math.floor((py / rowH) * (rows.length / nPixRows)))
  const pal = useMemo(() => (palette === 'story' ? [0, 0.5, 1, 1.5, 2].map((v) => hex(STORY_CN[v])) : [0, 0.5, 1, 1.5, 2].map((v) => hex(CN_COLORS[v]))), [palette])
  useEffect(() => { const el = wrap.current; if (!el) return; const ro = new ResizeObserver(() => { if (el.clientWidth > 100) setW(el.clientWidth) }); ro.observe(el); return () => ro.disconnect() }, [])
  const chromSpans = useMemo(() => meta.chromosomes.map((c) => ({ name: c.name, lo: c.offset, hi: c.offset + c.n })), [meta])

  useEffect(() => {
    const cv = canvas.current; if (!cv || rows.length === 0 || plotW < 10) return
    cv.width = w; cv.height = TOP_H + plotH
    const ctx = cv.getContext('2d')!
    ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, cv.width, cv.height)
    const img = ctx.createImageData(plotW, plotH); const d = img.data
    const colBin = new Int32Array(plotW); for (let px = 0; px < plotW; px++) colBin[px] = Math.min(nBins - 1, Math.floor((px / plotW) * nBins))
    for (let pr = 0; pr < nPixRows; pr++) {
      const r = rows[Math.min(rows.length - 1, Math.floor(pr * rows.length / nPixRows))]
      for (let py = 0; py < rowH; py++) {
        const off = ((pr * rowH + py) * plotW) * 4
        for (let px = 0; px < plotW; px++) {
          const v = r.data[colBin[px]]; const o = off + px * 4
          if (v === 255 || v === undefined) { d[o] = 230; d[o + 1] = 230; d[o + 2] = 230 } else { const c = pal[v]; d[o] = c[0]; d[o + 1] = c[1]; d[o + 2] = c[2] }
          d[o + 3] = 255
        }
      }
    }
    ctx.putImageData(img, left, TOP_H)
    // bands + labels
    ctx.font = '11px Inter, system-ui, sans-serif'; ctx.textBaseline = 'middle'
    for (let pr = 0; pr < nPixRows; pr++) {
      const ri = Math.min(rows.length - 1, Math.floor(pr * rows.length / nPixRows)); const r = rows[ri]
      if (r.band) { ctx.fillStyle = r.band; ctx.fillRect(labelWidth, TOP_H + pr * rowH, bandW - 2, rowH) }
      if (rowH >= 12 && labelWidth > 0) { ctx.fillStyle = '#253640'; ctx.textAlign = 'right'; ctx.fillText(r.label, labelWidth - 6, TOP_H + pr * rowH + rowH / 2) }
    }
    ctx.textAlign = 'center'; ctx.fillStyle = '#334155'; ctx.textBaseline = 'alphabetic'
    for (const c of chromSpans) {
      const x0 = left + (c.lo / nBins) * plotW, x1 = left + (c.hi / nBins) * plotW
      ctx.strokeStyle = 'rgba(15,23,42,0.4)'; ctx.beginPath(); ctx.moveTo(Math.round(x0) + 0.5, TOP_H); ctx.lineTo(Math.round(x0) + 0.5, TOP_H + plotH); ctx.stroke()
      if (x1 - x0 > 16) ctx.fillText(c.name.replace('chr', ''), (x0 + x1) / 2, 14)
    }
    if (selectedRow !== null && selectedRow !== undefined && selectedRow >= 0) { const pr = Math.floor(selectedRow * nPixRows / rows.length); ctx.strokeStyle = '#111827'; ctx.lineWidth = 2; ctx.strokeRect(left + 1, TOP_H + pr * rowH + 1, plotW - 2, Math.max(2, rowH - 2)) }
  }, [rows, w, plotW, plotH, rowH, nPixRows, pal, left, labelWidth, bandW, chromSpans, nBins, selectedRow])

  const locate = (e: React.MouseEvent) => {
    const rect = canvas.current!.getBoundingClientRect(); const x = e.clientX - rect.left, y = e.clientY - rect.top
    if (x < left || y < TOP_H || y >= TOP_H + plotH) return null
    const bin = Math.min(nBins - 1, Math.floor(((x - left) / plotW) * nBins)); const ri = rowOfPix(y - TOP_H)
    const chrom = chromSpans.find((c) => bin >= c.lo && bin < c.hi)!
    return { x, y, bin, ri, chrom }
  }
  const onMove = (e: React.MouseEvent) => {
    const l = locate(e); if (!l || presentation) { setHover(null); return }
    const v = rows[l.ri].data[l.bin]; const val = v === 255 ? null : v / 2
    const base = [`${rows[l.ri].label}${rows[l.ri].sub ? ' · ' + rows[l.ri].sub : ''}`, `${l.chrom.name}:${(bins[2 * l.bin] / MB).toFixed(2)}–${(bins[2 * l.bin + 1] / MB).toFixed(2)} Mb · chrom bin ${l.bin - l.chrom.lo} · genome bin ${l.bin}`, val === null ? 'CN: not called' : `CN ${val} = ${palette === 'story' ? STORY_CN_LABELS[val] : CN_LABELS[val]}`]
    setHover({ x: l.x, y: l.y, text: hoverText ? [...base, ...hoverText(l.ri, l.bin, l.chrom.name, val)] : base })
  }
  return (
    <div ref={wrap} className="ckaryo-wrap">
      <div style={{ position: 'relative', maxHeight: fit ? undefined : maxHeight, overflowY: fit ? undefined : 'auto' }}>
        <canvas ref={canvas} style={{ display: 'block', width: w, height: TOP_H + plotH, cursor: onClick ? 'pointer' : 'crosshair' }} onMouseMove={onMove} onMouseLeave={() => setHover(null)} onClick={(e) => { const l = locate(e); if (l && onClick) onClick(l.ri, l.bin, l.chrom.name) }} />
        {hover && <div className="tooltip" style={{ left: Math.min(hover.x + 14, w - 340), top: hover.y + 12 }}>{hover.text.map((t, i) => <div key={i}>{t}</div>)}</div>}
      </div>
      {fit && nPixRows < rows.length && <div className="muted small">Showing {nPixRows} of {rows.length} rows (nearest-row sampling at {rowH}px per row; no averaging). Switch row height to 1 px to see every cell.</div>}
    </div>
  )
}
