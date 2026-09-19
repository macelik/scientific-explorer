import type { AssignResult, CandidateRow, CellData, CellDetail, EventCandidate, IntegrationMeta, IntegrationStatus, Meta, NodeData, SegmentEstimator, TestResult } from './types'

async function j<T>(url: string, init?: RequestInit): Promise<T> {
  const r = await fetch(url, init)
  if (!r.ok) {
    let msg = `${r.status} ${r.statusText}`
    try { const b = await r.json(); if (b.detail) msg = typeof b.detail === 'string' ? b.detail : JSON.stringify(b.detail) } catch { /* ignore */ }
    throw new Error(msg)
  }
  return r.json() as Promise<T>
}
const post = (body: unknown): RequestInit => ({ method: 'POST', headers: { 'content-type': 'application/json' }, body: JSON.stringify(body) })

export const api = {
  meta: () => j<Meta>('/api/meta'),
  bins: async () => { const b = await (await fetch('/api/bins.bin')).arrayBuffer(); return new Int32Array(b) },
  karyogram: async () => { const r = await fetch('/api/karyogram.bin'); const b = await r.arrayBuffer(); return new Uint8Array(b) },
  cell: (cell: string) => j<CellData>(`/api/cell/${encodeURIComponent(cell)}`),
  node: (cell: string, chrom: string, start: number, end: number) => j<NodeData>('/api/node', post({ cell, chrom, start, end })),
  test: (cell: string, chrom: string, start: number, end: number, b: number, n_permutations: number, seed: number) =>
    j<{ id: string | null; status: string; result: TestResult | null; error?: string; cached: boolean; elapsed?: number }>('/api/test', post({ cell, chrom, start, end, b, n_permutations, seed })),
  job: (id: string) => j<{ id: string; status: string; result: TestResult | null; error: string | null; elapsed: number }>(`/api/jobs/${id}`),
  assign: (cell: string, chrom: string, boundaries: number[], segment_estimator: SegmentEstimator = 'arithmetic') => j<AssignResult>('/api/assign_cn', post({ cell, chrom, boundaries, segment_estimator })),
  eventsStart: () => j<any>('/api/events/index', { method: 'POST' }),
  eventsStatus: () => j<{ status: string; done: number; total: number; error: string | null }>('/api/events/status'),
  eventsCandidates: () => j<{ status: string; rows: EventCandidate[] | null; done?: number; total?: number }>('/api/events/candidates'),
  sessions: () => j<{ name: string; saved_at?: string; dataset_hash?: string; cell?: string; chrom?: string }[]>('/api/sessions'),
  saveSession: (name: string, payload: unknown) => j<{ ok: boolean }>('/api/sessions', post({ name, payload })),
  loadSession: (name: string) => j<{ format: string; version: number; saved_at: string; dataset: { hash: string }; state: any }>(`/api/sessions/${encodeURIComponent(name)}`),
  integrationStatus: () => j<IntegrationStatus>('/api/integration/status'),
  integrationIndex: () => j<IntegrationStatus>('/api/integration/index', { method: 'POST' }),
  integrationMeta: () => j<IntegrationMeta>('/api/integration/meta'),
  integrationPseudobulk: async () => { const b = await (await fetch('/api/integration/pseudobulk.bin')).arrayBuffer(); return new Float32Array(b) },
  integrationKaryogram: async () => { const r = await fetch('/api/integration/percell_karyogram.bin'); if (!r.ok) throw new Error(`${r.status}`); const b = await r.arrayBuffer(); return new Uint8Array(b) },
  integrationCell: (cell: string) => j<CellDetail>(`/api/integration/cell/${encodeURIComponent(cell)}`),
  integrationCandidates: (source: string, chrom: string) => j<{ rows: CandidateRow[] }>(`/api/integration/candidates?source=${encodeURIComponent(source)}&chrom=${encodeURIComponent(chrom)}`),
}
