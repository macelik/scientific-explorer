import { useMemo } from 'react'
import { computeCohort } from '../cohort'
import { useStore } from '../store'
import EventRules from './EventRules'
import Karyogram from './Karyogram'

export default function CohortView() {
  const meta = useStore((s) => s.meta)!
  const bins = useStore((s) => s.bins)!
  const karyo = useStore((s) => s.karyo)!
  const filters = useStore((s) => s.filters)
  const setFilters = useStore((s) => s.setFilters)
  const rules = useStore((s) => s.rules)
  const events = useStore((s) => s.events)
  const openCell = useStore((s) => s.openCell)
  const chrom = useStore((s) => s.chrom) || 'chr9'
  const cohort = useMemo(() => computeCohort(meta, bins, karyo, filters, rules, events.rows), [meta, bins, karyo, filters, rules, events.rows])
  const ruleMatch = useMemo(() => { const s = new Set<number>(); for (const e of cohort.ruleEvals) if (e.status === 'ready') e.matches.forEach((r) => s.add(r)); return s }, [cohort])
  const covs = meta.cells.map((c) => c.raw_coverage)
  const covMin = Math.min(...covs), covMax = Math.max(...covs)
  const open = (cell: string, chr: string) => openCell(cell, chr, { tab: 'explore' })
  return (
    <div className="cohort">
      <div className="card">
        <div className="row wrap">
          <label>Cell ID <input value={filters.search} placeholder="search…" onChange={(e) => setFilters({ search: e.target.value })} /></label>
          <label>Chromosome <select value={chrom} onChange={(e) => useStore.setState({ chrom: e.target.value })}>{meta.chromosomes.map((c) => <option key={c.name}>{c.name}</option>)}</select></label>
          <span className="muted">coverage tercile:</span>
          {(['low', 'mid', 'high'] as const).map((t) => <label key={t} className="chk-inline"><input type="checkbox" checked={filters.terciles.includes(t)} onChange={(e) => setFilters({ terciles: e.target.checked ? [...filters.terciles, t] : filters.terciles.filter((x) => x !== t) })} />{t}</label>)}
          <label>raw coverage <input type="number" placeholder={String(covMin)} value={filters.covMin ?? ''} onChange={(e) => setFilters({ covMin: e.target.value === '' ? null : +e.target.value })} style={{ width: 80 }} /> – <input type="number" placeholder={String(covMax)} value={filters.covMax ?? ''} onChange={(e) => setFilters({ covMax: e.target.value === '' ? null : +e.target.value })} style={{ width: 80 }} /></label>
          <label>order <select value={filters.order} onChange={(e) => setFilters({ order: e.target.value as any })}><option value="clustered">clustered (Ward / Euclidean, computed once)</option><option value="cov_desc">coverage ↓</option><option value="cov_asc">coverage ↑</option><option value="original">original cohort order</option></select></label>
          <button className="btn-xs" onClick={() => setFilters({ search: '', terciles: ['low', 'mid', 'high'], covMin: null, covMax: null })}>reset filters</button>
        </div>
        <div className="row wrap">
          <strong>{cohort.rows.length} matching cells</strong>
          <span className="muted">{cohort.explanation.length ? cohort.explanation.join(' AND ') : 'no filters active'}</span>
          <span className="spacer" />
          <span className="muted">reference cells:</span>
          {meta.reference_cells.map((r) => <button key={r.cell + r.chrom} className="btn-xs" title={r.label} onClick={() => open(r.cell, r.chrom)}>{r.cell} / {r.chrom}</button>)}
        </div>
      </div>
      <div className="card"><EventRules evals={cohort.ruleEvals} /></div>
      <div className="card">
        {cohort.rows.length === 0 ? (
          <div className="empty">
            <p>No cells match the current filters ({cohort.explanation.join(' AND ')}).</p>
            <button className="btn" onClick={() => { setFilters({ search: '', terciles: ['low', 'mid', 'high'], covMin: null, covMax: null }); useStore.setState({ rules: rules.map((r) => ({ ...r, enabled: false })) }) }}>Reset filters and rules</button>
          </div>
        ) : <Karyogram rows={cohort.rows} ruleMatch={ruleMatch} onOpen={open} />}
      </div>
      <div className="card">
        <div className="row"><strong>Matching cells</strong><span className="muted">displayed row → cell ID mapping (click to open on {chrom})</span></div>
        <div className="cell-list">
          {cohort.rows.map((r, i) => { const c = meta.cells[r]; return <button key={c.cellID} className={`cell-chip ${ruleMatch.has(r) ? 'match' : ''}`} onClick={() => open(c.cellID, chrom)}><span className="muted">{i + 1}</span> {c.cellID} <span className="muted">{c.tercile} · {c.raw_coverage.toLocaleString()}</span></button> })}
        </div>
      </div>
    </div>
  )
}
