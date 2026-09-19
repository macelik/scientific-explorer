import { create } from 'zustand'
import { api } from './api'
import { buildSmallSegments, type FlankRule, type Metric, type SmallSeg } from './integration'
import type { CandidateRow, CellDetail, IntegrationMeta, IntegrationStatus } from './types'

export type ColorBy = 'cluster' | 'log10_raw' | 'sumX' | 'HA_HB' | 'n_bins_snp' | 'depth_weight' | 'haplo_weight' | 'chr8_score' | 'zero_frac' | 'in_cohort' | 'selected'
export type LevelBasis = 'chrmedian' | 'genome_trimmed' | 'raw'
export interface LabelLayers { breakpoints: boolean; depthColor: boolean; small: boolean; proposal: boolean; consistency: boolean; reason: boolean; cn: boolean; candidates: boolean; allelic: boolean; segLabels: boolean }

export interface SegSettings {
  chrom: string; sources: string[]; metric: Metric; effectMetric: Metric; flank: FlankRule; scale: 'binned' | 'continuous'; basis: LevelBasis; smoothing: number
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
}

const defaultSeg: SegSettings = {
  chrom: 'chr8', sources: ['cluster0', 'cluster1', 'cluster2', 'cluster3', 'cluster4'], metric: 'delta_mean', effectMetric: 'delta_mean', flank: 'min_delta', scale: 'binned', basis: 'chrmedian', smoothing: 1,
  labels: { breakpoints: true, depthColor: true, small: true, proposal: false, consistency: true, reason: false, cn: false, candidates: false, allelic: false, segLabels: false },
  smallMax: 100, srdRef: 3.3, tau: 0.2, selectedSegment: null, karyoPalette: 'story', percellOrder: 'cluster', percellRowH: 'fit',
}

export const useIntegration = create<IntegrationState>((set, get) => ({
  status: null, meta: null, loading: false, error: null, pb: null, karyo: null, karyoLoading: false, smallSegs: [],
  clusterKey: 'wnn_leiden_0.3', colorBy: 'cluster', selectedCluster: null, highlight: [], hoverCell: null, focusCell: null, cellDetail: null, dimOthers: false, pointSize: 3, compareKey: 'wnn_leiden_1.0',
  seg: defaultSeg, candidates: {},

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
    } catch (e: any) { set({ loading: false, error: String(e.message || e) }) }
  },
  ensurePseudobulk: async () => { if (get().pb) return; try { const pb = await api.integrationPseudobulk(); set({ pb }) } catch (e: any) { set({ error: `pseudobulk: ${e.message || e}` }) } },
  ensureKaryo: async () => { if (get().karyo || get().karyoLoading) return; set({ karyoLoading: true }); try { const k = await api.integrationKaryogram(); set({ karyo: k, karyoLoading: false }) } catch (e: any) { set({ karyoLoading: false, error: `per-cell karyogram: ${e.message || e}` }) } },
  set: (p) => set(p),
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
  serialize: () => { const s = get(); return { clusterKey: s.clusterKey, colorBy: s.colorBy, selectedCluster: s.selectedCluster, highlight: s.highlight, focusCell: s.focusCell, dimOthers: s.dimOthers, pointSize: s.pointSize, compareKey: s.compareKey, seg: s.seg } },
  restore: (p) => { if (!p) return; set({ clusterKey: p.clusterKey || 'wnn_leiden_0.3', colorBy: p.colorBy || 'cluster', selectedCluster: p.selectedCluster ?? null, highlight: p.highlight || [], focusCell: p.focusCell ?? null, dimOthers: !!p.dimOthers, pointSize: p.pointSize || 3, compareKey: p.compareKey || 'wnn_leiden_1.0', seg: { ...defaultSeg, ...(p.seg || {}), labels: { ...defaultSeg.labels, ...((p.seg || {}).labels || {}) } } }) },
}))
