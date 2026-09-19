/** Cohort filtering / ordering / editable event-rule matching (pure functions). */
import { MB } from './coords'
import type { Filters } from './store'
import type { EventCandidate, EventRule, Meta } from './types'

export interface RuleEval { ruleId: string; status: 'ready' | 'pending' | 'unavailable' | 'disabled'; matches: Set<number>; nBinsInRegion: number; note: string }
export interface CohortResult { rows: number[]; matchAll: number; ruleEvals: RuleEval[]; explanation: string[]; eventsPending: boolean }

export function regionBins(meta: Meta, bins: Int32Array, chrom: string, startMb: number, endMb: number): number[] {
  const c = meta.chromosomes.find((x) => x.name === chrom)
  if (!c) return []
  const out: number[] = []
  const s = startMb * MB, e = endMb * MB
  for (let i = c.offset; i < c.offset + c.n; i++) {
    const bs = bins[2 * i], be = bins[2 * i + 1]
    if (be >= s && bs <= e) out.push(i)
  }
  return out
}

export function evalRule(rule: EventRule, meta: Meta, bins: Int32Array, karyo: Uint8Array, events: EventCandidate[] | null): RuleEval {
  const nBins = meta.dataset.identity.n_bins
  if (!rule.enabled) return { ruleId: rule.id, status: 'disabled', matches: new Set(), nBinsInRegion: 0, note: '' }
  const rb = regionBins(meta, bins, rule.chrom, rule.startMb, rule.endMb)
  if (rule.family === 'cn') {
    if (!rb.length) return { ruleId: rule.id, status: 'unavailable', matches: new Set(), nBinsInRegion: 0, note: 'no retained bins overlap this region' }
    const stateCodes = new Set(rule.states.map((s) => Math.round(s * 2)))
    const m = new Set<number>()
    for (let r = 0; r < meta.cells.length; r++) {
      let cnt = 0
      for (const b of rb) if (stateCodes.has(karyo[r * nBins + b])) cnt++
      const frac = cnt / rb.length
      if (cnt > 0 && frac >= rule.minFraction) m.add(r)
    }
    return { ruleId: rule.id, status: 'ready', matches: m, nBinsInRegion: rb.length, note: `${rb.length} retained bins in region; match = ≥1 bin in states {${rule.states.join(', ')}}${rule.minFraction > 0 ? ` and fraction ≥ ${rule.minFraction}` : ''}` }
  }
  if (!rb.length) return { ruleId: rule.id, status: 'unavailable', matches: new Set(), nBinsInRegion: 0, note: 'no retained bins overlap this region' }
  if (!events) return { ruleId: rule.id, status: 'pending', matches: new Set(), nBinsInRegion: rb.length, note: 'cohort breakpoint index still computing; not counted as false' }
  const c = meta.chromosomes.find((x) => x.name === rule.chrom)!
  const rowOf = new Map(meta.cells.map((x, i) => [x.cellID, i]))
  const m = new Set<number>()
  for (const ev of events) {
    if (ev.chrom !== rule.chrom) continue
    if (rule.node !== 'any' && ev.node !== rule.node) continue
    if (rule.acceptance === 'accepted' && !ev.accepted) continue
    const mb = bins[2 * (c.offset + ev.bp)] / MB  // boundary drawn at the start of its right-hand bin
    if (mb >= rule.startMb && mb <= rule.endMb) { const r = rowOf.get(ev.cell); if (r !== undefined) m.add(r) }
  }
  return { ruleId: rule.id, status: 'ready', matches: m, nBinsInRegion: rb.length, note: `candidate positions from actual node AD argmax (${rule.node} node); ${rule.acceptance === 'accepted' ? 'accepted = production breakpoint with accepted parent' : 'any candidate, accepted or not'}` }
}

export function computeCohort(meta: Meta, bins: Int32Array, karyo: Uint8Array, filters: Filters, rules: EventRule[], events: EventCandidate[] | null): CohortResult {
  const explanation: string[] = []
  let rows: number[]
  switch (filters.order) {
    case 'clustered': rows = [...meta.cluster_order]; break
    case 'cov_asc': rows = meta.cells.map((_, i) => i).sort((a, b) => meta.cells[a].raw_coverage - meta.cells[b].raw_coverage); break
    case 'cov_desc': rows = meta.cells.map((_, i) => i).sort((a, b) => meta.cells[b].raw_coverage - meta.cells[a].raw_coverage); break
    default: rows = meta.cells.map((_, i) => i)
  }
  const search = filters.search.trim().toUpperCase()
  const terc = new Set(filters.terciles)
  rows = rows.filter((r) => {
    const c = meta.cells[r]
    if (search && !c.cellID.toUpperCase().includes(search)) return false
    if (!terc.has(c.tercile)) return false
    if (filters.covMin !== null && c.raw_coverage < filters.covMin) return false
    if (filters.covMax !== null && c.raw_coverage > filters.covMax) return false
    return true
  })
  if (search) explanation.push(`cell ID contains "${filters.search.trim()}"`)
  if (terc.size < 3) explanation.push(`tercile ∈ {${[...terc].join(', ')}}`)
  if (filters.covMin !== null || filters.covMax !== null) explanation.push(`coverage in [${filters.covMin ?? '−∞'}, ${filters.covMax ?? '∞'}]`)
  const ruleEvals = rules.map((r) => evalRule(r, meta, bins, karyo, events))
  const active = ruleEvals.filter((e) => e.status === 'ready')
  const pending = ruleEvals.filter((e) => e.status === 'pending')
  if (active.length) {
    rows = rows.filter((r) => filters.rulesMode === 'any' ? active.some((e) => e.matches.has(r)) : active.every((e) => e.matches.has(r)))
    explanation.push(`${active.length} event rule(s), match ${filters.rulesMode === 'any' ? 'ANY' : 'ALL'}`)
  }
  if (pending.length) explanation.push(`${pending.length} rule(s) pending (not applied yet)`)
  return { rows, matchAll: meta.cells.length, ruleEvals, explanation, eventsPending: pending.length > 0 }
}
