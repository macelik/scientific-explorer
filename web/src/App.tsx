import React, { useEffect, useMemo, useState } from 'react'
import { computeCohort } from './cohort'
import { useStore } from './store'
import CohortView from './components/CohortView'
import ExplorerView from './components/ExplorerView'
import ManualView from './components/ManualView'
import Inspector from './components/Inspector'
import SessionBar from './components/SessionBar'
import ExperimentView from './components/ExperimentView'
import IntegrationClustersView from './components/IntegrationClustersView'
import IntegrationSegmentationView from './components/IntegrationSegmentationView'
import './components/integration.css'

class ErrorBoundary extends React.Component<{ children: React.ReactNode; label: string }, { error: string | null }> {
  state = { error: null as string | null }
  static getDerivedStateFromError(e: any) { return { error: String(e?.message || e) } }
  render() { return this.state.error ? <div className="notice">{this.props.label} failed to render: {this.state.error} <button className="btn-xs" onClick={() => this.setState({ error: null })}>retry</button></div> : this.props.children }
}

function BootScreen({ error, onRetry }: { error: string | null; onRetry: () => void }) {
  const steps = useStore((s) => s.bootSteps)
  const icon = (st: string) => (st === 'done' ? '✓' : st === 'running' ? '…' : st === 'error' ? '✗' : '·')
  return (
    <div className="boot">
      <h2>pyEpiAneufinder explorer</h2>
      <p className="muted">Starting up. After launching <code>run.py</code> the service needs a few seconds to import and index the data; this page retries automatically.</p>
      <ul className="boot-steps">{steps.map((st, i) => <li key={i} className={st.status}><span className="ico">{icon(st.status)}</span> {st.label}{st.detail ? <span className="muted"> — {st.detail}</span> : null}</li>)}</ul>
      {error && <div className="notice">Failed to load: {error} <button className="btn-sm" onClick={onRetry}>retry</button></div>}
    </div>
  )
}

export default function App() {
  const booted = useStore((s) => s.booted); const bootError = useStore((s) => s.bootError); const boot = useStore((s) => s.boot)
  const meta = useStore((s) => s.meta)
  const tab = useStore((s) => s.tab); const setTab = useStore((s) => s.setTab)
  const cell = useStore((s) => s.cell); const chrom = useStore((s) => s.chrom)
  const openCell = useStore((s) => s.openCell)
  const toast = useStore((s) => s.toast); const setToast = useStore((s) => s.setToast)
  const filters = useStore((s) => s.filters); const rules = useStore((s) => s.rules); const events = useStore((s) => s.events)
  const bins = useStore((s) => s.bins); const karyo = useStore((s) => s.karyo)
  const drawer = useStore((s) => s.cohortDrawerOpen); const setDrawer = useStore((s) => s.setCohortDrawer)
  const [showInspector, setShowInspector] = useState(false)
  const presentation = useStore((s) => s.presentation); const setPresentation = useStore((s) => s.setPresentation)
  useEffect(() => { boot() }, [])
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      const t = e.target as HTMLElement
      if (t && (t.tagName === 'INPUT' || t.tagName === 'SELECT' || t.tagName === 'TEXTAREA' || t.isContentEditable)) return
      if ((e.key === 'p' || e.key === 'P') && !e.ctrlKey && !e.metaKey && !e.altKey) { e.preventDefault(); setPresentation(!useStore.getState().presentation) }
      if (e.key === 'Escape' && useStore.getState().presentation) setPresentation(false)
    }
    window.addEventListener('keydown', onKey)
    return () => window.removeEventListener('keydown', onKey)
  }, [])
  useEffect(() => { if (toast) { const t = setTimeout(() => setToast(null), 6000); return () => clearTimeout(t) } }, [toast])
  const cohort = useMemo(() => (meta && bins && karyo ? computeCohort(meta, bins, karyo, filters, rules, events.rows) : null), [meta, bins, karyo, filters, rules, events.rows])
  if (!booted || !meta) return <BootScreen error={bootError} onRetry={() => boot()} />
  const cm = cell ? meta.cells.find((c) => c.cellID === cell) : null
  const idx = cohort && cm ? cohort.rows.indexOf(cm.row) : -1
  const step = (d: number) => { if (!cohort || !chrom) return; const i = idx < 0 ? 0 : (idx + d + cohort.rows.length) % cohort.rows.length; openCell(meta.cells[cohort.rows[i]].cellID, chrom) }
  return (
    <div className="app">
      <header>
        <div className="brand">pyEpi <span>Scientific explorer</span></div>
        <nav className="tabs">
          <button className={tab === 'cohort' ? 'on' : ''} onClick={() => setTab('cohort')}>Cohort karyogram</button>
          <button className={tab === 'explore' ? 'on' : ''} onClick={() => setTab('explore')}>Explore output</button>
          <button className={tab === 'manual' ? 'on' : ''} onClick={() => setTab('manual')}>Manual segmentation</button>
          <button className={tab === 'experiment' ? 'on' : ''} onClick={() => setTab('experiment')}>X experiments</button>
          <span className="tab-sep" aria-hidden="true" />
          <button className={tab === 'clusters' ? 'on' : ''} onClick={() => setTab('clusters')}>Integration · clusters</button>
          <button className={tab === 'clusterseg' ? 'on' : ''} onClick={() => setTab('clusterseg')}>Integration · segmentation &amp; CN</button>
        </nav>
        <span className="session-bar"><SessionBar /><button className="btn-sm" onClick={() => setShowInspector(!showInspector)}>{showInspector ? 'hide' : 'show'} inspector</button>
          <button className={`btn-sm ${presentation ? 'on' : ''}`} title="Presentation mode: no hover labels, tooltips or toasts so a pointer stays visible (shortcut: P, Esc to leave)" onClick={() => setPresentation(!presentation)}>{presentation ? '● presenting (P)' : 'present (P)'}</button></span>
      </header>
      {cell && chrom && cm && (
        <div className="cellbar">
          <b>{cell}</b> · <select value={chrom} onChange={(e) => openCell(cell, e.target.value, { tab })}>{meta.chromosomes.map((c) => <option key={c.name}>{c.name}</option>)}</select>
          <span className="muted">coverage {cm.raw_coverage.toLocaleString()} ({cm.tercile}) · {cm.n_segments} production segments · best_s {cm.best_s.toFixed(3)}</span>
          <span className="spacer" />
          {cohort && <span className="muted">{idx >= 0 ? `cell ${idx + 1} of ${cohort.rows.length} matching` : `not in the current filtered cohort (${cohort.rows.length} matching)`}</span>}
          <button className="btn-sm" onClick={() => step(-1)} disabled={!cohort?.rows.length}>◀ previous</button>
          <button className="btn-sm" onClick={() => step(1)} disabled={!cohort?.rows.length}>next ▶</button>
          <button className="btn-sm" onClick={() => setDrawer(!drawer)}>{drawer ? 'hide' : 'show'} cohort drawer</button>
          <button className="btn-sm" onClick={() => setTab('cohort')}>back to cohort</button>
        </div>
      )}
      {drawer && cohort && tab !== 'cohort' && (
        <div className="drawer">
          <span className="muted">filters: {cohort.explanation.length ? cohort.explanation.join(' AND ') : 'none'} · {cohort.rows.length} cells</span>
          <div className="cell-list compact">{cohort.rows.map((r, i) => { const c = meta.cells[r]; return <button key={c.cellID} className={`cell-chip ${c.cellID === cell ? 'on' : ''}`} onClick={() => openCell(c.cellID, chrom!)}>{i + 1} {c.cellID}</button> })}</div>
        </div>
      )}
      <div className={`body ${showInspector ? 'with-inspector' : ''}`}>
        <main>
          <div style={{ display: tab === 'cohort' ? 'block' : 'none' }}><ErrorBoundary label="Cohort view"><CohortView /></ErrorBoundary></div>
          <div style={{ display: tab === 'explore' ? 'block' : 'none' }}>{tab === 'explore' && <ErrorBoundary label="Explorer"><ExplorerView /></ErrorBoundary>}</div>
          <div style={{ display: tab === 'manual' ? 'block' : 'none' }}>{tab === 'manual' && <ErrorBoundary label="Manual segmentation"><ManualView /></ErrorBoundary>}</div>
          <div style={{ display: tab === 'experiment' ? 'block' : 'none' }}>{tab === 'experiment' && <ErrorBoundary label="X experiments"><ExperimentView /></ErrorBoundary>}</div>
          <div style={{ display: tab === 'clusters' ? 'block' : 'none' }}>{tab === 'clusters' && <ErrorBoundary label="Integration · clusters"><IntegrationClustersView /></ErrorBoundary>}</div>
          <div style={{ display: tab === 'clusterseg' ? 'block' : 'none' }}>{tab === 'clusterseg' && <ErrorBoundary label="Integration · segmentation"><IntegrationSegmentationView /></ErrorBoundary>}</div>
        </main>
        {showInspector && <ErrorBoundary label="Inspector"><Inspector /></ErrorBoundary>}
      </div>
      {toast && <div className="toast" onClick={() => setToast(null)}>{toast}</div>}
    </div>
  )
}
