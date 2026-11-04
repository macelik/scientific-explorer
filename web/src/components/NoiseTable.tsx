import { useState } from 'react'
import { useIntegration } from '../integrationStore'
import { downloadText } from '../cnRecall'
import { NOISE_DEFINITIONS, noiseCellsN, noiseRows, noiseTsv, roundsHeader, sameRounds, type NoiseColumn } from '../mergeRounds'
import { SOURCE_COLORS } from '../integration'

/** Call noise by round, under the chromosome tracks (spec §4). */
export default function NoiseTable() {
  const ig = useIntegration()
  const [scope, setScope] = useState<'chrom' | 'genome'>('genome')
  const key = scope === 'genome' ? 'genome' : ig.seg.chrom
  const st = ig.recall
  const pending = !!st?.diff?.noise && sameRounds(st.settings.rounds ?? [], ig.rounds)
  const cols: NoiseColumn[] = [
    { name: 'Original', noise: ig.base?.noise_original ?? null },
    ...ig.rounds.map((r, i) => ({ name: `Round ${i + 1}`, noise: r.noise })),
    ...(pending ? [{ name: 'After', noise: st!.diff!.noise!, cellsPending: !st!.percell_ready }] : []),
  ]
  const tsvCols = cols.map((c, i) => ({ ...c, name: i === 0 ? 'original' : c.name === 'After' ? 'after' : `round${i}` }))
  const header = [...roundsHeader(ig.rounds), ...(pending ? [`after rules: ${JSON.stringify(st!.settings.rules)} · allow_gap_crossing: ${st!.settings.allow_gap_crossing}`] : [])]
  return <div className="noise-table" role="region" aria-label="Call noise by round">
    <h4>Call noise by round</h4>
    <div className="controls-row">
      <div className="segmented" role="group" aria-label="Noise scope">
        <button className={`btn-sm${scope === 'chrom' ? ' on' : ''}`} aria-pressed={scope === 'chrom'} onClick={() => setScope('chrom')}>this chromosome</button>
        <button className={`btn-sm${scope === 'genome' ? ' on' : ''}`} aria-pressed={scope === 'genome'} onClick={() => setScope('genome')}>genome-wide</button>
      </div>
      <span className="muted small">{scope === 'genome' ? 'all chromosomes' : ig.seg.chrom} · sources shown on the tracks</span>
    </div>
    {!ig.base ? <p className="muted small">loading…</p> : <div className="table-scroll"><table className="tbl small">
      <thead><tr><th>source</th><th>level</th><th>measure</th>{cols.map((c) => <th key={c.name}>{c.name}</th>)}</tr></thead>
      <tbody>{ig.seg.sources.flatMap((src) => noiseRows(cols, src, key).map((r, i) => <tr key={`${src}-${i}`}>
        {i === 0 && <td rowSpan={8} style={{ color: SOURCE_COLORS[src] }}>{src}</td>}
        {i % 4 === 0 && <td rowSpan={4}>{r.level === 'cluster' ? 'cluster row' : `cells (n = ${noiseCellsN(cols, src, key) ?? '…'})`}</td>}
        <td>{r.measure}</td>{r.values.map((v, j) => <td key={j}>{v}</td>)}
      </tr>))}</tbody>
    </table></div>}
    <ul className="small muted" aria-label="Noise definitions">{NOISE_DEFINITIONS.map((d) => <li key={d}>{d}</li>)}</ul>
    <button className="btn-sm" disabled={!ig.base} onClick={() => downloadText(noiseTsv(tsvCols, ig.seg.sources, header), 'noise_by_round.tsv')}>noise_by_round.tsv</button>
  </div>
}
