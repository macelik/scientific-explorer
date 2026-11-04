// Integration · allelic segmentation (spec docs/superpowers/specs/2026-10-28-allelic-segmentation-design.md §4–5)
import { useEffect, useMemo, useRef } from 'react'
import Plot from './DeferredPlot'
import { downloadPlot } from './Plot'
import { useAllelic } from '../allelicStore'
import { useIntegration } from '../integrationStore'
import { useStore } from '../store'
import { MB } from '../coords'
import { DEPTH_COLORS, SOURCE_COLORS, STORY_CN, clusterSegmentsLocal, levelTrack, withGaps } from '../integration'
import { wireRounds } from '../mergeRounds'
import { downloadText } from '../cnRecall'
import { COMPARISON_DEFINITIONS, POOL_NOTE, SETTING_FIELDS, SIGNALS, agreementRows, allelicSegmentRows, comparisonTsv, depthSegmentRows, runHeader, sameRun, signalLabel, type AdSettings, type LocalDepthSeg } from '../allelic'

const fmt = (v: number | null | undefined, d = 3) => (v === null || v === undefined || !Number.isFinite(v) ? 'N/A' : v.toFixed(d))

export default function AllelicSegmentationView() {
  const al = useAllelic()
  const ig = useIntegration()
  const exMeta = useStore((s) => s.meta)!
  const bins = useStore((s) => s.bins)!
  const presentation = useStore((s) => s.presentation)
  const plotEl = useRef<any>(null)
  const meta = ig.meta
  const roundsSig = JSON.stringify(wireRounds(ig.rounds))
  const nRounds = ig.rounds.length
  useEffect(() => { ig.load(); al.loadStatus() }, [])
  useEffect(() => { if (meta) ig.ensurePseudobulk() }, [meta])
  useEffect(() => { if (meta && nRounds) ig.refreshBase() }, [meta, roundsSig])
  useEffect(() => { al.ensureTracks() }, [al.status?.state, al.chrom, al.sources.join(','), al.signal, al.settings.min_at])

  const useOriginal = al.depthSource === 'original' || nRounds === 0
  const depthBps = useOriginal ? (meta?.breakpoints ?? null) : (ig.base?.breakpoints ?? null)
  const depthSegsAll = useOriginal ? (meta?.cluster_segments ?? null) : (ig.base?.cluster_segments ?? null)
  const chromMeta = exMeta.chromosomes.find((c) => c.name === al.chrom) ?? exMeta.chromosomes[0]
  const localSegs: LocalDepthSeg[] = useMemo(() => clusterSegmentsLocal(depthSegsAll, exMeta.chromosomes).filter((g) => g.chromosome === chromMeta.name).map((g) => ({ group: g.group, s: g.s, e: g.e, cn_call: g.cn_call })), [depthSegsAll, chromMeta])
  const shown = al.shown
  const stale = !!shown && !sameRun(shown, al.signal, al.settings)
  const signed = SIGNALS.find((s) => s.id === al.signal)!.signed
  const sources = (meta?.sources ?? []).filter((s) => al.sources.includes(s))
  const tracksMatch = !!al.tracks && al.tracks.chrom === chromMeta.name && al.tracks.signal === al.signal && al.tracks.min_at === al.settings.min_at
  const trackValues = useMemo(() => (tracksMatch ? Object.fromEntries(al.tracks!.tracks.map((t) => [t.source, t.values])) : {}), [al.tracks, tracksMatch])

  // ---------------- tracks
  const fig = useMemo(() => {
    if (!meta || !ig.pb || !al.tracks || !tracksMatch || !sources.length) return null
    const nBins = exMeta.dataset.identity.n_bins
    const mb: number[] = []; for (let i = 0; i < chromMeta.n; i++) mb.push(bins[2 * (chromMeta.offset + i)] / MB)
    const endsMb = Array.from({ length: chromMeta.n }, (_, i) => bins[2 * (chromMeta.offset + i) + 1] / MB)
    const xEnd = endsMb[endsMb.length - 1]
    const xAt = (b: number) => (b < chromMeta.n ? mb[b] : xEnd)
    const data: any[] = [], shapes: any[] = [], layout: any = {}
    const nP = sources.length * 2, rowH = 1 / nP, gap = 0.05 / nP
    sources.forEach((src, k) => {
      const si = meta.sources.indexOf(src)
      const dAx = 2 * k + 1, aAx = 2 * k + 2
      const dRef = `y${dAx === 1 ? '' : dAx}`, aRef = `y${aAx}`
      const dDom: [number, number] = [1 - (2 * k + 1) * rowH + gap, 1 - 2 * k * rowH]
      const aDom: [number, number] = [1 - (2 * k + 2) * rowH + gap, 1 - (2 * k + 1) * rowH]
      // depth panel: pseudobulk, CN strip, depth breakpoints
      const lt = levelTrack(ig.pb!.subarray(si * nBins, (si + 1) * nBins), chromMeta, 'chrmedian', ig.pb!.subarray(si * nBins, (si + 1) * nBins))
      const g = withGaps(mb, lt.y, endsMb)
      data.push({ type: 'scatter', mode: 'lines', x: g.x, y: g.y, line: { color: '#6b7280', width: 0.8 }, connectgaps: false, xaxis: 'x', yaxis: dRef, hovertemplate: `${src} depth · %{x:.2f} Mb · %{y:.2f}<extra></extra>`, showlegend: false })
      layout[`yaxis${dAx === 1 ? '' : dAx}`] = { domain: dDom, range: [-4.6, 3.2], title: { text: `${src}<br><span style="font-size:8px">depth log2</span>`, font: { size: 10, color: SOURCE_COLORS[src] } }, gridcolor: '#f1f5f9', zeroline: false, tickfont: { size: 8 } }
      for (const cs of localSegs.filter((c) => c.group === src)) shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: xAt(cs.s), x1: xAt(cs.e), y0: dDom[1] - 0.04 * rowH, y1: dDom[1], fillcolor: STORY_CN[cs.cn_call], line: { width: 0 }, layer: 'below' })
      const dB = (depthBps ?? []).filter((b) => b.source === src && b.chromosome === chromMeta.name)
      for (const b of dB) {
        const x = xAt(b.absolute_bin)
        shapes.push({ type: 'line', xref: 'x', yref: `${dRef} domain`, x0: x, x1: x, y0: 0, y1: 1, line: { color: DEPTH_COLORS[Math.min(b.recursion_depth, DEPTH_COLORS.length - 1)], width: 1.2 }, name: 'depth breakpoint' })
        shapes.push({ type: 'line', xref: 'x', yref: `${aRef} domain`, x0: x, x1: x, y0: 0, y1: 1, line: { color: '#9ca3af', width: 1, dash: 'dash' }, opacity: 0.6, name: 'depth breakpoint (echo)' })
      }
      // allelic panel: dots by AT, skipped ticks, segment means, breakpoints and gap bands
      const t = al.tracks!.tracks.find((x) => x.source === src)
      if (t) {
        const ats = t.at.filter((a, i) => t.values[i] !== null).sort((a, b) => a - b)
        const p90 = ats.length ? ats[Math.floor(0.9 * (ats.length - 1))] : al.settings.min_at
        const opa = (a: number) => Math.max(0.15, Math.min(1, p90 > al.settings.min_at ? 0.15 + 0.85 * (a - al.settings.min_at) / (p90 - al.settings.min_at) : 1))
        const xs: number[] = [], ys: number[] = [], os: number[] = [], txt: string[] = [], skipped: number[] = []
        t.values.forEach((v, i) => { if (v === null) skipped.push(mb[i]); else { xs.push(mb[i]); ys.push(v); os.push(opa(t.at[i])); txt.push(`AT ${t.at[i]}`) } })
        data.push({ type: 'scatter', mode: 'markers', x: xs, y: ys, text: txt, marker: { size: 3, color: SOURCE_COLORS[src], opacity: os }, xaxis: 'x', yaxis: aRef, hovertemplate: `${src} ${signalLabel(al.signal)} %{y:.3f} · %{text}<extra></extra>`, showlegend: false, name: 'allelic bins' })
        const lo = signed ? -0.5 : 0
        if (skipped.length) data.push({ type: 'scatter', mode: 'markers', x: skipped, y: skipped.map(() => lo), marker: { symbol: 'line-ns-open', size: 4, color: '#9ca3af' }, xaxis: 'x', yaxis: aRef, hoverinfo: 'skip', showlegend: false, name: 'skipped bins' })
      }
      layout[`yaxis${aAx}`] = { domain: aDom, range: signed ? [-0.5, 0.5] : [0, 0.5], title: { text: `${signalLabel(al.signal)}`, font: { size: 9, color: '#b5651d' } }, gridcolor: '#f8fafc', zeroline: signed, zerolinecolor: '#9ca3af', tickfont: { size: 8 } }
      if (shown && !stale) {
        for (const s of shown.result!.segments.filter((x) => x.source === src && x.chromosome === chromMeta.name && x.mean !== null))
          shapes.push({ type: 'line', xref: 'x', yref: aRef, x0: xAt(s.bin_start), x1: xAt(s.bin_end + 1), y0: s.mean, y1: s.mean, line: { color: '#111827', width: 2 }, name: 'allelic segment mean' })
        for (const b of shown.result!.breakpoints.filter((x) => x.source === src && x.chromosome === chromMeta.name)) {
          if (b.gap_start !== null && b.gap_end !== null) shapes.push({ type: 'rect', xref: 'x', yref: `${aRef} domain`, x0: xAt(b.gap_start), x1: xAt(b.gap_end), y0: 0, y1: 1, fillcolor: '#b5651d', opacity: 0.18, line: { width: 0 }, name: 'allelic gap breakpoint' })
          else shapes.push({ type: 'line', xref: 'x', yref: `${aRef} domain`, x0: xAt(b.absolute_bin), x1: xAt(b.absolute_bin), y0: 0, y1: 1, line: { color: '#b5651d', width: 1.6 }, name: 'allelic breakpoint' })
        }
      }
    })
    Object.assign(layout, { height: Math.max(300, 120 * nP + 60), margin: { l: 70, r: 14, t: 20, b: 40 }, paper_bgcolor: 'white', plot_bgcolor: 'white', font: { size: 11 }, hovermode: presentation ? false : 'closest', shapes, showlegend: false, uirevision: chromMeta.name, xaxis: { title: { text: `${chromMeta.name} position (Mb)` }, range: [mb[0], xEnd], anchor: `y${nP}`, gridcolor: '#f1f5f9' } })
    return { data, layout }
  }, [meta, ig.pb, al.tracks, tracksMatch, chromMeta, sources.join(','), localSegs, depthBps, shown, stale, signed, al.signal, al.settings.min_at, presentation, bins])

  // ---------------- tables
  const chroms = al.scope === 'genome' ? exMeta.chromosomes.map((c) => c.name) : [chromMeta.name]
  const t1 = useMemo(() => (shown && depthBps ? agreementRows(sources, chroms, depthBps, shown.result!.breakpoints, al.matchN) : []), [shown, depthBps, sources.join(','), chroms.join(','), al.matchN])
  const t2 = useMemo(() => (shown && tracksMatch ? depthSegmentRows(localSegs.filter((g) => sources.includes(g.group)), trackValues, shown.result!.breakpoints, chromMeta.name, al.matchN) : []), [shown, tracksMatch, localSegs, trackValues, sources.join(','), chromMeta, al.matchN])
  const t3 = useMemo(() => (shown ? allelicSegmentRows(shown.result!.segments, localSegs, sources, chromMeta.name) : []), [shown, localSegs, sources.join(','), chromMeta])

  const setSetting = (k: Exclude<keyof AdSettings, 'global_pool'>, v: string) => al.setView({ settings: { ...al.settings, [k]: v === '' ? NaN : +v } })
  const st = al.status
  const disabled = !st || st.state !== 'ready'
  const progress = al.active?.progress
  return (
    <div className="allelic-view">
      <div className="section-heading"><div><div className="eyebrow">INTEGRATION · ALLELIC</div><h2>Allelic segmentation</h2>
        <p>The depth run's constrained AD scan, run on a per-cluster allelic pseudobulk signal. Depth tracks and segments are shown above each allelic track for comparison.</p></div></div>
      {al.notice && <p className="notice small">{al.notice}</p>}
      {st && st.state === 'error' && <div className="notice small" role="alert">Allelic data unavailable: {st.problems.join('; ')}</div>}
      {st && st.state === 'building' && <p className="muted small">Building allelic pseudobulks from {st.allelic_dir}…</p>}
      {al.statusError && <div className="notice small" role="alert">{al.statusError}</div>}

      <div className="panel">
        <div className="controls-row">
          <label>signal<select aria-label="Allelic signal" value={al.signal} onChange={(e) => al.setView({ signal: e.target.value as any })}>{SIGNALS.map((s) => <option key={s.id} value={s.id}>{s.label}</option>)}</select></label>
          <span className="muted small">{SIGNALS.find((s) => s.id === al.signal)!.formula}</span>
          <label>min AT<input aria-label="min AT" type="number" value={Number.isFinite(al.settings.min_at) ? al.settings.min_at : ''} onChange={(e) => setSetting('min_at', e.target.value)} style={{ width: 70 }} /></label>
          <button className="btn-sm" disabled={disabled || !!al.active} onClick={() => al.run()}>Run segmentation</button>
          {al.active && <span className="small" aria-live="polite">{al.active.status} · {progress ? `${progress.source} · ${progress.chromosome} · ${progress.done}/${progress.total}` : 'starting'}</span>}
          {al.runError && <span className="notice small" role="alert">{al.runError}</span>}
        </div>
        <details className="allelic-ad-settings"><summary>AD settings</summary>
          <div className="controls-row">{SETTING_FIELDS.filter((f) => f.key !== 'min_at').map((f) => <label key={f.key}>{f.label}<input aria-label={f.label} type="number" step={f.integer ? 1 : 0.0001} value={Number.isFinite(al.settings[f.key]) ? al.settings[f.key] : ''} onChange={(e) => setSetting(f.key, e.target.value)} style={{ width: 80 }} /></label>)}
            <label>global test pool<select aria-label="Global test pool" value={al.settings.global_pool} onChange={(e) => al.setView({ settings: { ...al.settings, global_pool: e.target.value as any } })}><option value="genome">genome (original pyEpi)</option><option value="chromosome">chromosome (earlier depth run)</option></select></label></div>
          <p className="muted small">{POOL_NOTE}</p>
        </details>
        {stale && <p className="notice small">Segments shown are from an earlier run (signal {signalLabel(shown!.signal)}, settings {JSON.stringify(shown!.settings)}). Run again to update.</p>}
        <div className="controls-row">
          <label>chromosome<select aria-label="Allelic chromosome" value={chromMeta.name} onChange={(e) => al.setView({ chrom: e.target.value })}>{exMeta.chromosomes.map((c) => <option key={c.name}>{c.name}</option>)}</select></label>
          {(meta?.sources ?? []).map((s) => <label key={s} className="chk"><input type="checkbox" checked={al.sources.includes(s)} onChange={(e) => al.setView({ sources: e.target.checked ? [...al.sources, s] : al.sources.filter((x) => x !== s) })} /><span style={{ color: SOURCE_COLORS[s] }}>{s}</span></label>)}
          <div className="segmented" role="group" aria-label="Depth segments">
            <button className={`btn-sm${al.depthSource === 'current' ? ' on' : ''}`} aria-pressed={al.depthSource === 'current'} onClick={() => al.setView({ depthSource: 'current' })}>current round{nRounds ? ` (${nRounds})` : ''}</button>
            <button className={`btn-sm${al.depthSource === 'original' ? ' on' : ''}`} aria-pressed={al.depthSource === 'original'} onClick={() => al.setView({ depthSource: 'original' })}>unmerged</button>
          </div>
        </div>
        {nRounds > 0 && !useOriginal && !ig.base && !ig.baseError && <p className="muted small">loading merged depth segments…</p>}
        {nRounds > 0 && !useOriginal && ig.baseError && <div className="notice small" role="alert">{ig.baseError}</div>}
        {!shown && !disabled && <p className="muted small">Run segmentation to see allelic segments.</p>}
        {al.tracksError && <div className="notice small" role="alert">{al.tracksError}</div>}
        {fig ? <><Plot data={fig.data} layout={fig.layout} config={{ displayModeBar: presentation ? false : 'hover' }} onReady={(el) => { plotEl.current = el }} /><button className="btn-xs" onClick={() => plotEl.current && downloadPlot(plotEl.current, `allelic_${chromMeta.name}`)}>PNG</button></>
          : <div className="muted small">{sources.length ? 'loading tracks…' : 'select at least one cluster row'}</div>}
      </div>

      <div className="panel allelic-compare" role="region" aria-label="Depth vs allelic segments">
        <h3>Depth vs allelic segments</h3>
        {!shown ? <p className="muted small">The comparison needs a finished run.</p> : <>
          <div className="controls-row">
            <div className="segmented" role="group" aria-label="Comparison scope">
              <button className={`btn-sm${al.scope === 'chrom' ? ' on' : ''}`} aria-pressed={al.scope === 'chrom'} onClick={() => al.setView({ scope: 'chrom' })}>this chromosome</button>
              <button className={`btn-sm${al.scope === 'genome' ? ' on' : ''}`} aria-pressed={al.scope === 'genome'} onClick={() => al.setView({ scope: 'genome' })}>genome-wide</button>
            </div>
            <label>match N (bins)<input aria-label="Match N" type="number" min={0} max={50} value={al.matchN} onChange={(e) => al.setView({ matchN: Math.max(0, Math.min(50, Math.round(+e.target.value || 0))) })} style={{ width: 60 }} /></label>
            <button className="btn-sm" disabled={!tracksMatch} onClick={() => downloadText(comparisonTsv(runHeader(shown, al.matchN, al.depthSource, nRounds, al.scope === 'genome' ? 'genome-wide' : chromMeta.name, meta?.identity.hash), t1, t2, t3), 'allelic_vs_depth.tsv')}>allelic_vs_depth.tsv</button>
          </div>
          <h4>Breakpoint agreement · {al.scope === 'genome' ? 'genome-wide' : chromMeta.name}</h4>
          <div className="table-scroll"><table className="tbl small" aria-label="Breakpoint agreement"><thead><tr><th>source</th><th>depth bps</th><th>allelic bps</th><th>matched</th><th>depth only</th><th>allelic only</th><th>median offset</th></tr></thead>
            <tbody>{t1.map((r) => <tr key={r.source}><td style={{ color: SOURCE_COLORS[r.source] }}>{r.source}</td><td>{r.depth}</td><td>{r.allelic}</td><td>{r.matched}</td><td>{r.depthOnly}</td><td>{r.allelicOnly}</td><td>{fmt(r.medianOffset, 1)}</td></tr>)}</tbody></table></div>
          <h4>Depth segments seen by the allelic signal · {chromMeta.name}</h4>
          <div className="table-scroll"><table className="tbl small" aria-label="Depth segments seen by the allelic signal"><thead><tr><th>source</th><th>bins</th><th>CN call</th><th>kept bins</th><th>mean</th><th>sd</th><th>allelic bps inside</th></tr></thead>
            <tbody>{!tracksMatch && <tr><td colSpan={7}>loading tracks…</td></tr>}{t2.map((r) => <tr key={`${r.source}${r.bin_start}`}><td style={{ color: SOURCE_COLORS[r.source] }}>{r.source}</td><td>{r.bin_start}–{r.bin_end}</td><td>{r.cn_call.toFixed(1)}</td><td>{r.kept}</td><td>{fmt(r.mean)}</td><td>{fmt(r.sd)}</td><td>{r.allelicInside || ''}</td></tr>)}</tbody></table></div>
          <h4>Allelic segments seen by depth · {chromMeta.name}</h4>
          <div className="table-scroll"><table className="tbl small" aria-label="Allelic segments seen by depth"><thead><tr><th>source</th><th>bins</th><th>kept bins</th><th>mean</th><th>depth calls spanned</th></tr></thead>
            <tbody>{t3.map((r) => <tr key={`${r.source}${r.bin_start}`}><td style={{ color: SOURCE_COLORS[r.source] }}>{r.source}</td><td>{r.bin_start}–{r.bin_end}</td><td>{r.kept}</td><td>{fmt(r.mean)}</td><td>{r.spanned}</td></tr>)}</tbody></table></div>
          <ul className="small muted" aria-label="Comparison definitions">{COMPARISON_DEFINITIONS.map((d) => <li key={d}>{d}</li>)}</ul>
        </>}
      </div>
    </div>
  )
}
