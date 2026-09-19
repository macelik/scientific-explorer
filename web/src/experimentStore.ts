import { create } from 'zustand'
import type { TestResult } from './types'

export interface ExperimentDraft { start: number; end: number; operation: 'multiply' | 'add' | 'set' | 'simulate' | 'extend' | 'duplicate'; value: number; preserve_total: boolean; cn_state?: 'loss' | 'base' | 'gain'; direction?: 'left' | 'right'; length?: number; seed?: number; k?: number }
export interface ExperimentNode {
  start: number; end: number; side: string; depth?: number; n: number; ad: number[]; argmax_b: number | null; argmax_ad: number | null
  accepted: boolean; hypothetical: boolean; passes?: boolean; tests: TestResult | null
}
export interface ExperimentCN {
  status: string; reason?: string; best_s?: number; baseline?: number; cn5: number[] | null; cn_int?: number[]; cn_cont?: number[]
  changed_bins?: number; changed_elsewhere?: { chrom: string; changed_bins: number }[]
  holmes?: number[]; watson?: number[]; event_calls?: Record<string, Record<string, number>>
}
export interface ExperimentResult {
  request: ExperimentDraft & { cell: string; chrom: string; dataset_hash: string; science_hash: string }
  parameters: { k: number; n_permutations: number; seed: number; alpha: number }
  intervention: { clipped_bins: number; normalization_factor: number; original_total: number; edited_total: number; changed_bins_genome: number; event_start: number; event_end: number; target_start: number; target_end: number; actual_extension?: number; requested_extension?: number; donor_count?: number; donor_mean?: number; donor_zero_fraction?: number; sampling?: string }
  original: { nodes: ExperimentNode[]; boundaries: number[] }; edited: { nodes: ExperimentNode[]; boundaries: number[] }
  edited_at_original_boundary: TestResult | null
  start_bp: number[]; end_bp: number[]; x_original: number[]; x_edited: number[]
  cn: { original: ExperimentCN; fixed: ExperimentCN; resegmented: ExperimentCN; event_isolated: ExperimentCN }; compute_ms: number
  event_edges: { edge: string; boundary: number; chromosome_edge: boolean; accepted: boolean; nodes: (TestResult & { side: string; rank: number; is_winner: boolean })[] }[]
  provenance: Record<string, string>
}
export interface ExperimentEntry { draft: ExperimentDraft; result: ExperimentResult | null; status: 'idle' | 'running' | 'done' | 'error'; error?: string; token?: string; history?: ExperimentResult[]; series?: { done: number; total: number; running: boolean; token: string } }

export function sameDraft(a: ExperimentDraft, b: ExperimentDraft) {
  return a.start === b.start && a.end === b.end && a.operation === b.operation && a.value === b.value && a.preserve_total === b.preserve_total &&
    (a.cn_state || 'loss') === (b.cn_state || 'loss') && (a.direction || 'right') === (b.direction || 'right') && (a.length ?? 40) === (b.length ?? 40) && (a.seed ?? 42) === (b.seed ?? 42) && (a.k ?? 2) === (b.k ?? 2)
}
async function json(url: string, body?: unknown): Promise<any> {
  const res = await fetch(url, body === undefined ? undefined : { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) })
  const data = await res.json()
  if (!res.ok) throw new Error(typeof data.detail === 'string' ? data.detail : JSON.stringify(data.detail))
  return data
}

export const experimentKey = (cell: string, chrom: string) => `${cell}|${chrom}`
export const useExperimentStore = create<{
  entries: Record<string, ExperimentEntry>
  initialize: (key: string, draft: ExperimentDraft) => void
  update: (key: string, patch: Partial<ExperimentDraft>) => void
  run: (cell: string, chrom: string) => Promise<void>
  runSeries: (cell: string, chrom: string, n: number) => Promise<void>
  cancelSeries: (key: string) => void
  selectResult: (key: string, result: ExperimentResult) => void
  restore: (entries: Record<string, ExperimentEntry>) => void
}>((set, get) => ({
  entries: {},
  initialize: (key, draft) => { if (!get().entries[key]) set({ entries: { ...get().entries, [key]: { draft, result: null, status: 'idle' } } }) },
  update: (key, patch) => { const e = get().entries[key]; if (e) set({ entries: { ...get().entries, [key]: { ...e, draft: { ...e.draft, ...patch } } } }) },
  run: async (cell, chrom) => {
    const key = experimentKey(cell, chrom), e = get().entries[key]
    if (!e) return
    const token = crypto.randomUUID(), draft = { ...e.draft }
    const change = (patch: Partial<ExperimentEntry>) => {
      const latest = get().entries[key]
      if (latest?.token === token) set({ entries: { ...get().entries, [key]: { ...latest, ...patch } } })
    }
    set({ entries: { ...get().entries, [key]: { ...e, token, status: 'running', error: undefined } } })
    try {
      let job = await json('/api/experiment', { cell, chrom, ...draft })
      while (job.status === 'running') {
        await new Promise(r => setTimeout(r, 250))
        if (get().entries[key]?.token !== token) return
        job = await json(`/api/jobs/${job.id}`)
      }
      if (job.status !== 'done') throw new Error(job.error || 'Experiment failed')
      change({ result: job.result, status: 'done', history: [...(get().entries[key]?.history || []), job.result].slice(-12) })
    } catch (err: any) { change({ status: 'error', error: err.message }) }
  },
  runSeries: async (cell, chrom, n) => {
    const key = experimentKey(cell, chrom), entry = get().entries[key]
    if (!entry || entry.series?.running) return
    const initial = { ...entry.draft }, widths = [10, 25, 50, 100, 150], token = crypto.randomUUID()
    const mark = (done: number, running: boolean) => set({ entries: { ...get().entries, [key]: { ...get().entries[key], series: { done, total: widths.length, running, token } } } })
    mark(0, true)
    for (let i=0;i<widths.length;i++) {
      if (!get().entries[key].series?.running || get().entries[key].series?.token !== token) return
      const extension = ['extend', 'duplicate'].includes(initial.operation), center = Math.floor((initial.start+initial.end)/2)
      const width = widths[i]
      const draft = extension ? { ...initial, length: width } : { ...initial, start: Math.max(0,center-Math.floor(width/2)), end: Math.min(n,center+Math.ceil(width/2)) }
      get().update(key,draft)
      await get().run(cell,chrom)
      if (!get().entries[key].series?.running || get().entries[key].series?.token !== token) return
      if (get().entries[key].status === 'error') { mark(i, false); return }
      mark(i+1, i<widths.length-1)
    }
  },
  cancelSeries: key => { const e=get().entries[key]; if(e?.series) set({entries:{...get().entries,[key]:{...e,series:{...e.series,running:false}}}}) },
  selectResult: (key,result) => { const e=get().entries[key]; if(e) set({entries:{...get().entries,[key]:{...e,draft:{...result.request},result,status:'done'}}}) },
  restore: (entries) => set({ entries: Object.fromEntries(Object.entries(entries || {}).map(([k, v]) => [k, { ...v, status: v.result ? 'done' : 'idle', token: undefined, series: undefined }])) }),
}))
