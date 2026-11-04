import { create } from 'zustand'
import { api } from './api'
import { buildSmallSegments, type SmallSeg } from './integration'
import { migrateSegSettings, type TrackMetric } from './trackMetrics'
import type { CandidateRow, CellDetail, IntegrationMeta, IntegrationStatus } from './types'
import { fetchHatchScores, type HatchScoresResult } from './hatchMerge'
import { fetchRecall, fetchRecallPercell, startRecallApi, type KaryoMode, type RecallBody, type RecallStatus } from './cnRecall'
import { fetchBase, pinRecallApi, roundLabel, roundsFromSession, sameRounds, wireRounds, type BaseResponse, type Round } from './mergeRounds'
import { migrateHatchLayers, type UserRule, type Rule0Criterion } from './mergeRules'

export type ColorBy = 'cluster' | 'log10_raw' | 'sumX' | 'HA_HB' | 'n_bins_snp' | 'depth_weight' | 'haplo_weight' | 'chr8_score' | 'zero_frac' | 'in_cohort' | 'selected'
export type LevelBasis = 'chrmedian' | 'genome_trimmed' | 'raw'
export interface LabelLayers { small: boolean; signs: boolean; reason: boolean; cn: boolean; candidates: boolean; allelic: boolean; segLabels: boolean }
export interface SegSettings {
  chrom: string; sources: string[]; colorMetric: TrackMetric; criteria: Rule0Criterion[]; combine: 'all' | 'any'; mergeSignMismatch: boolean; mergeSignsOnly: boolean; scale: 'binned' | 'continuous'; basis: LevelBasis; smoothing: number
  labels: LabelLayers; smallMax: number; srdRef: number; tau: number; selectedSegment: string | null
  karyoPalette: 'story' | 'explorer'; percellOrder: 'cluster' | 'coverage' | 'chr8'; percellRowH: 'fit' | 1 | 2
}

export interface IntegrationState {
  status: IntegrationStatus | null; meta: IntegrationMeta | null; loading: boolean; error: string | null
  pb: Float32Array | null; karyo: Uint8Array | null; karyoLoading: boolean
  smallSegs: SmallSeg[]
  clusterKey: string; colorBy: ColorBy; selectedCluster: string | null; highlight: string[]; hoverCell: string | null; focusCell: string | null
  cellDetail: CellDetail | null; dimOthers: boolean; pointSize: number; compareKey: string
  seg: SegSettings
  candidates: Record<string, CandidateRow[]>
  mergeRules: UserRule[]; allowGapCrossing: boolean
  hatchScores: Record<string, HatchScoresResult>
  hatchScoresError: string | null
  setMergeRules: (rules: UserRule[]) => void
  setAllowGapCrossing: (on: boolean) => void
  ensureHatchScores: (source: string, chrom: string) => Promise<void>
  load: () => Promise<void>
  ensurePseudobulk: () => Promise<void>
  ensureKaryo: () => Promise<void>
  set: (p: Partial<IntegrationState>) => void
  setSeg: (p: Partial<SegSettings>) => void
  setLabels: (p: Partial<LabelLayers>) => void
  selectCluster: (c: string | null) => void
  setHighlight: (ids: string[]) => void
  focus: (cell: string | null) => Promise<void>
  ensureCandidates: (source: string, chrom: string) => Promise<void>
  serialize: () => any
  restore: (p: any) => void
  recall: RecallStatus | null; recallKaryo: Uint8Array | null; recallSignature: string | null; recallError: string | null; karyoMode: KaryoMode
  startRecall: (body: RecallBody) => Promise<void>
  rounds: Round[]; base: BaseResponse | null; baseSig: string | null; baseError: string | null; roundsNotice: string | null; savedDepthIdentity: string | null
  baseKaryo: Uint8Array | null; baseCalls: RecallStatus | null
  setRounds: (rounds: Round[]) => void
  refreshBase: () => Promise<void>
  adoptRecall: () => Promise<void>
  undoRound: () => void
  resetRounds: () => void
}

const defaultSeg: SegSettings = {
  chrom: 'chr8', sources: ['cluster0', 'cluster1', 'cluster2', 'cluster3', 'cluster4'], colorMetric: 'delta_mean', criteria: ['delta_mean', 'srd_phi'], combine: 'all', mergeSignMismatch: false, mergeSignsOnly: false, scale: 'binned', basis: 'chrmedian', smoothing: 1,
  labels: { small: true, signs: true, reason: false, cn: false, candidates: false, allelic: false, segLabels: false },
  smallMax: 100, srdRef: 3.3, tau: 0.2, selectedSegment: null, karyoPalette: 'story', percellOrder: 'cluster', percellRowH: 'fit',
}

export const DEPTH_CHANGED_NOTICE = 'These merge rounds were saved against a different depth segmentation; they were re-applied to the current one.'
// compares the identity stored in a saved session with the loaded one; runs once meta is available
export function checkSavedIdentity(get: () => IntegrationState, set: (p: Partial<IntegrationState>) => void) {
  const s = get()
  if (!s.meta || s.savedDepthIdentity === null) return
  if (s.rounds.length && s.savedDepthIdentity !== s.meta.identity.hash) set({ roundsNotice: DEPTH_CHANGED_NOTICE })
  set({ savedDepthIdentity: null })
}

export const useIntegration = create<IntegrationState>((set, get) => ({
  status: null, meta: null, loading: false, error: null, pb: null, karyo: null, karyoLoading: false, smallSegs: [],
  clusterKey: 'wnn_leiden_0.3', colorBy: 'cluster', selectedCluster: null, highlight: [], hoverCell: null, focusCell: null, cellDetail: null, dimOthers: false, pointSize: 3, compareKey: 'wnn_leiden_1.0',
  seg: defaultSeg, candidates: {},
  mergeRules: [], allowGapCrossing: false, hatchScores: {}, hatchScoresError: null,
  recall: null, recallKaryo: null, recallSignature: null, recallError: null, karyoMode: 'before',
  rounds: [], base: null, baseSig: null, baseError: null, roundsNotice: null, savedDepthIdentity: null, baseKaryo: null, baseCalls: null,

  load: async () => {
    if (get().meta || get().loading) return
    set({ loading: true, error: null })
    try {
      let st = await api.integrationStatus()
      set({ status: st })
      if (!st.available) { set({ loading: false, error: st.problems.join('; ') }); return }
      // wait for the background index so meta carries derived per-cell fields
      let tries = 0
      while (st.index.status !== 'done' && st.index.status !== 'error' && tries < 400) { await new Promise((r) => setTimeout(r, 1500)); st = await api.integrationStatus(); set({ status: st }); tries++ }
      const meta = await api.integrationMeta()
      set({ meta, smallSegs: buildSmallSegments(meta.flank_scores, meta.screening), loading: false, clusterKey: meta.default_cluster_key })
      checkSavedIdentity(get, set)
    } catch (e: any) { set({ loading: false, error: String(e.message || e) }) }
  },
  ensurePseudobulk: async () => { if (get().pb) return; try { const pb = await api.integrationPseudobulk(); set({ pb }) } catch (e: any) { set({ error: `pseudobulk: ${e.message || e}` }) } },
  ensureKaryo: async () => { if (get().karyo || get().karyoLoading) return; set({ karyoLoading: true }); try { const k = await api.integrationKaryogram(); set({ karyo: k, karyoLoading: false }) } catch (e: any) { set({ karyoLoading: false, error: `per-cell karyogram: ${e.message || e}` }) } },
  set: (p) => set(p),
  startRecall: async (body) => {
    const signature = JSON.stringify(body)
    set({ recallError: null })
    let id: string
    try { id = (await startRecallApi(body)).id } catch (e: any) { set({ recallError: String(e.message || e) }); return }
    // the rounds changed (undo/adopt/restore) while the POST was in flight: this re-call sits on a stale base
    if (JSON.stringify(wireRounds(get().rounds)) !== JSON.stringify(body.rounds ?? [])) return
    if (get().recall?.id !== id) set({ recall: null, recallKaryo: null, karyoMode: 'before' })
    set({ recallSignature: signature })
    const poll = async () => {
      if (get().recallSignature !== signature) return
      let st: RecallStatus
      try { st = await fetchRecall(id) } catch (e: any) { if (get().recallSignature === signature) { const msg = `re-call ${id}: ${e.message || e}`; const cur = get().recall; set({ recallError: msg, recall: cur ? { ...cur, status: 'error', error: msg } : null }) } return }
      if (get().recallSignature !== signature) return
      set({ recall: st })
      if (st.percell_ready && !get().recallKaryo) {
        try { const k = await fetchRecallPercell(id); if (get().recallSignature === signature) set({ recallKaryo: k }) }
        catch (e: any) { if (get().recallSignature === signature) set({ recallError: `per-cell re-call: ${e.message || e}` }) }
      }
      if (st.status === 'running') setTimeout(poll, 1000)
    }
    await poll()
  },
  setRounds: (rounds) => set({ rounds, base: null, baseSig: null, baseError: null, baseKaryo: null, baseCalls: null, hatchScores: {}, hatchScoresError: null,
    recall: null, recallKaryo: null, recallSignature: null, recallError: null, karyoMode: 'before' }),
  undoRound: () => { set({ roundsNotice: null }); get().setRounds(get().rounds.slice(0, -1)) },
  resetRounds: () => { set({ roundsNotice: null }); get().setRounds([]) },
  refreshBase: async () => {
    const wire = wireRounds(get().rounds), sig = JSON.stringify(wire)
    set({ baseSig: sig, baseError: null })
    let b: BaseResponse
    try { b = await fetchBase(wire) } catch (e: any) {
      if (get().baseSig !== sig) return
      const msg = String(e.message || e)
      if (wire.length && (/^round \d+: /.test(msg) || e?.status === 422)) { get().setRounds([]); set({ roundsNotice: `Saved merge rounds could not be applied: ${msg}` }); return }
      set({ baseError: `segmentation base: ${msg}` }); return
    }
    if (get().baseSig !== sig) return
    set({ base: b })
    const job = b.percell.job
    if (b.percell.status === 'ready' && job && !get().baseKaryo) {
      try { const k = await fetchRecallPercell(job); if (get().baseSig === sig) set({ baseKaryo: k }) }
      catch (e: any) { if (get().baseSig === sig) set({ baseError: `base per-cell: ${e.message || e}` }) }
    } else if (b.percell.status === 'running' && job) {
      const poll = async () => {
        if (get().baseSig !== sig) return
        let st: RecallStatus
        try { st = await fetchRecall(job) } catch (e: any) { if (get().baseSig === sig) set({ baseError: `base per-cell ${job}: ${e.message || e}` }); return }
        if (get().baseSig !== sig) return
        set({ baseCalls: st })
        if (st.status === 'running') { setTimeout(poll, 1000); return }
        if (st.status === 'done') await get().refreshBase()   // now 'ready': loads the matrix and the cells noise
      }
      await poll()
    }
  },
  adoptRecall: async () => {
    const st = get().recall, cur = get().rounds
    if (!st || st.status !== 'done' || !sameRounds(st.settings.rounds ?? [], cur)) return
    const round: Round = { rules: st.settings.rules, allow_gap_crossing: st.settings.allow_gap_crossing, label: roundLabel(st.settings.rules as any[]),
      merges: (st.merge?.pairs ?? []).reduce((a, p) => a + p.merges, 0), noise: st.diff?.noise ?? null }
    const next = [...cur, round]
    const k = get().recall === st ? get().recallKaryo : null
    try { await pinRecallApi(st.id, wireRounds(next)) } catch (e: any) { if (get().recall === st && sameRounds(get().rounds, cur)) set({ recallError: `adopt: ${e.message || e}` }); return }
    if (get().recall !== st || !sameRounds(get().rounds, cur)) return
    set({ roundsNotice: null })
    get().setRounds(next)
    if (k) set({ baseKaryo: k })   // the adopted re-call's per-cell calls are the new base's
  },
  setSeg: (p) => set({ seg: { ...get().seg, ...p } }),
  setLabels: (p) => set({ seg: { ...get().seg, labels: { ...get().seg.labels, ...p } } }),
  selectCluster: (c) => set({ selectedCluster: c }),
  setHighlight: (ids) => set({ highlight: ids }),
  focus: async (cell) => {
    set({ focusCell: cell, cellDetail: null })
    if (!cell) return
    try { const d = await api.integrationCell(cell); if (get().focusCell === cell) set({ cellDetail: d }) } catch { /* ignore */ }
  },
  ensureCandidates: async (source, chrom) => {
    const k = `${source}|${chrom}`; if (get().candidates[k]) return
    try { const r = await api.integrationCandidates(source, chrom); set({ candidates: { ...get().candidates, [k]: r.rows } }) } catch { /* ignore */ }
  },
  setMergeRules: (rules) => set({ mergeRules: rules }),
  setAllowGapCrossing: (on) => set({ allowGapCrossing: on }),
  ensureHatchScores: async (source, chrom) => {
    const k = `${source}|${chrom}`
    if (get().hatchScores[k]) return
    const wire = wireRounds(get().rounds), sig = JSON.stringify(wire)
    try { const r = await fetchHatchScores(chrom, source, wire); if (JSON.stringify(wireRounds(get().rounds)) === sig) set({ hatchScores: { ...get().hatchScores, [k]: r }, hatchScoresError: null }) }
    catch (e: any) { if (JSON.stringify(wireRounds(get().rounds)) === sig) set({ hatchScoresError: `hatch-scores: ${e.message || e}` }) }
  },
  serialize: () => { const s = get(); return { clusterKey: s.clusterKey, colorBy: s.colorBy, selectedCluster: s.selectedCluster, highlight: s.highlight, focusCell: s.focusCell, dimOthers: s.dimOthers, pointSize: s.pointSize, compareKey: s.compareKey, seg: s.seg, mergeRules: s.mergeRules, allowGapCrossing: s.allowGapCrossing, rounds: s.rounds, depthIdentity: s.meta?.identity.hash ?? s.savedDepthIdentity ?? undefined } },
  restore: (p) => { if (!p) return; set({ clusterKey: p.clusterKey || 'wnn_leiden_0.3', colorBy: p.colorBy || 'cluster', selectedCluster: p.selectedCluster ?? null, highlight: p.highlight || [], focusCell: p.focusCell ?? null, dimOthers: !!p.dimOthers, pointSize: p.pointSize || 3, compareKey: p.compareKey || 'wnn_leiden_1.0', seg: (() => { const ms = migrateSegSettings(p.seg); return { ...defaultSeg, ...ms, labels: { ...defaultSeg.labels, ...(ms.labels || {}) } } })(), mergeRules: Array.isArray(p.mergeRules) ? p.mergeRules : migrateHatchLayers(p.hatchLayers), allowGapCrossing: !!p.allowGapCrossing, rounds: roundsFromSession(p), base: null, baseSig: null, baseError: null, baseKaryo: null, baseCalls: null, hatchScores: {}, recall: null, recallKaryo: null, recallSignature: null, recallError: null, roundsNotice: null, savedDepthIdentity: typeof p.depthIdentity === 'string' ? p.depthIdentity : null, karyoMode: 'before' }); checkSavedIdentity(get, set) },
}))
