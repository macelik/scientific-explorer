import { newId, useStore } from '../store'
import type { EventRule } from '../types'
import type { RuleEval } from '../cohort'

const STATES = [0, 0.5, 1, 1.5, 2]

export default function EventRules({ evals }: { evals: RuleEval[] }) {
  const rules = useStore((s) => s.rules)
  const setRules = useStore((s) => s.setRules)
  const meta = useStore((s) => s.meta)!
  const filters = useStore((s) => s.filters)
  const setFilters = useStore((s) => s.setFilters)
  const events = useStore((s) => s.events)
  const upd = (id: string, patch: Partial<EventRule>) => setRules(rules.map((r) => (r.id === id ? { ...r, ...patch } : r)))
  const add = (family: 'cn' | 'bp') => setRules([...rules, {
    id: newId('rule'), name: family === 'cn' ? 'New regional CN rule' : 'New breakpoint rule', enabled: true, family, chrom: 'chr9', startMb: 19, endMb: 26,
    states: family === 'cn' ? [0, 0.5] : [], minFraction: 0, node: 'root', acceptance: 'candidate',
  }])
  const evalOf = (id: string) => evals.find((e) => e.ruleId === id)
  return (
    <div className="rules">
      <div className="row">
        <strong>Editable event rules</strong>
        <span className="muted">transparent, overlapping exploration tags — not validated biological labels</span>
        <span className="spacer" />
        <label className="muted">combine: <select value={filters.rulesMode} onChange={(e) => setFilters({ rulesMode: e.target.value as any })}><option value="any">match ANY</option><option value="all">match ALL</option></select></label>
        <button className="btn-xs" onClick={() => add('cn')}>+ CN rule</button>
        <button className="btn-xs" onClick={() => add('bp')}>+ breakpoint rule</button>
      </div>
      {events.status !== 'done' && <div className="notice small">Breakpoint-behaviour rules need the cohort candidate index ({events.status}: {events.done}/{events.total} cells). They stay pending until it finishes and are not counted as false.</div>}
      {rules.map((r) => {
        const ev = evalOf(r.id)
        return (
          <div key={r.id} className={`rule ${r.enabled ? 'on' : ''}`}>
            <label className="chk"><input type="checkbox" checked={r.enabled} onChange={(e) => upd(r.id, { enabled: e.target.checked })} /></label>
            <input className="rule-name" value={r.name} onChange={(e) => upd(r.id, { name: e.target.value })} />
            <span className={`pill ${r.family}`}>{r.family === 'cn' ? 'regional CN' : 'breakpoint'}</span>
            <select value={r.chrom} onChange={(e) => upd(r.id, { chrom: e.target.value })}>{meta.chromosomes.map((c) => <option key={c.name}>{c.name}</option>)}</select>
            <input type="number" step="0.1" value={r.startMb} onChange={(e) => upd(r.id, { startMb: +e.target.value })} style={{ width: 62 }} />–
            <input type="number" step="0.1" value={r.endMb} onChange={(e) => upd(r.id, { endMb: +e.target.value })} style={{ width: 62 }} /> Mb
            {r.family === 'cn' ? (
              <>
                <span className="muted">states:</span>
                {STATES.map((s) => <label key={s} className="chk-inline"><input type="checkbox" checked={r.states.includes(s)} onChange={(e) => upd(r.id, { states: e.target.checked ? [...r.states, s] : r.states.filter((x) => x !== s) })} />{s}</label>)}
                <span className="muted">min fraction of retained bins:</span>
                <input type="number" min={0} max={1} step={0.05} value={r.minFraction} onChange={(e) => upd(r.id, { minFraction: +e.target.value })} style={{ width: 56 }} />
              </>
            ) : (
              <>
                <select value={r.node} onChange={(e) => upd(r.id, { node: e.target.value as any })}><option value="root">depth-0 (root) node</option><option value="left">left child</option><option value="right">right child</option><option value="any">any node</option></select>
                <select value={r.acceptance} onChange={(e) => upd(r.id, { acceptance: e.target.value as any })}><option value="candidate">candidate (argmax, tested or not)</option><option value="accepted">accepted production breakpoint</option></select>
              </>
            )}
            <span className="spacer" />
            {ev && ev.status === 'ready' && <span className="count">{ev.matches.size} cells</span>}
            {ev && ev.status === 'pending' && <span className="count pending">pending</span>}
            {ev && ev.status === 'unavailable' && <span className="count warn">unavailable: {ev.note}</span>}
            <button className="btn-xs" title={ev?.note || ''} onClick={() => alert(`Rule: ${r.name}\nfamily: ${r.family}\nregion: ${r.chrom}:${r.startMb}-${r.endMb} Mb (overlap with retained-bin intervals)\n${ev?.note || ''}\nsource: ${r.family === 'cn' ? 'default_300cells/cn_matrix.npy (five-state CN per bin)' : 'node AD argmax computed by the service for root/left/right nodes; accepted = present in production breakpoints.csv with an accepted parent'}`)}>ⓘ</button>
            <button className="btn-xs" onClick={() => setRules(rules.filter((x) => x.id !== r.id))}>✕</button>
          </div>
        )
      })}
    </div>
  )
}
