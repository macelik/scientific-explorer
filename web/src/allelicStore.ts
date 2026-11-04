// Allelic segmentation tab state (spec §4). Saved with the session under `allelic`.
import { create } from 'zustand'
import { DEFAULT_VIEW, allelicApi, allelicFromSession, settingsProblems, type AllelicJob, type AllelicStatus, type AllelicTracks, type AllelicView } from './allelic'

export const timing = { pollMs: 1000 }
const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms))
const tracksKeyOf = (v: AllelicView) => `${v.chrom}|${v.sources.join(',')}|${v.signal}|${v.settings.min_at}`
const pickView = (s: AllelicView): AllelicView => ({ signal: s.signal, settings: s.settings, chrom: s.chrom, sources: s.sources, depthSource: s.depthSource, matchN: s.matchN, scope: s.scope })
let runToken = 0

interface AllelicState extends AllelicView {
  status: AllelicStatus | null; statusError: string | null
  tracks: AllelicTracks | null; tracksKey: string | null; tracksError: string | null
  shown: AllelicJob | null; active: AllelicJob | null; runError: string | null; notice: string | null
  loadStatus: () => Promise<void>
  setView: (p: Partial<AllelicView>) => void
  ensureTracks: () => Promise<void>
  run: () => Promise<void>
  serialize: () => AllelicView
  restore: (p: any) => void
}

export const useAllelic = create<AllelicState>((set, get) => ({
  ...DEFAULT_VIEW,
  status: null, statusError: null, tracks: null, tracksKey: null, tracksError: null, shown: null, active: null, runError: null, notice: null,
  loadStatus: async () => {
    try {
      const status = await allelicApi.status()
      set({ status, statusError: null })
      if (status.state === 'idle' || status.state === 'building') setTimeout(() => { get().loadStatus() }, 2000)
    } catch (e: any) { set({ statusError: e.message || String(e) }) }
  },
  setView: (p) => set(p),
  ensureTracks: async () => {
    const s = get()
    if (s.status?.state !== 'ready' || !s.sources.length) return
    const key = tracksKeyOf(s)
    if (key === s.tracksKey) return
    set({ tracksKey: key, tracksError: null })
    try {
      const tracks = await allelicApi.tracks(s.chrom, s.sources, s.signal, s.settings.min_at)
      if (get().tracksKey === key) set({ tracks })
    } catch (e: any) { if (get().tracksKey === key) set({ tracksError: e.message || String(e), tracksKey: null }) }
  },
  run: async () => {
    const s = get()
    const problems = settingsProblems(s.settings)
    if (problems.length) { set({ runError: problems.join('; ') }); return }
    const token = ++runToken
    set({ runError: null })
    try {
      let job = await allelicApi.start(s.signal, s.settings)
      if (token !== runToken) return
      set({ active: job })
      while (job.status === 'queued' || job.status === 'running') {
        await sleep(timing.pollMs)
        if (token !== runToken) return
        job = await allelicApi.job(job.id)
        if (token !== runToken) return
        set({ active: job })
      }
      if (job.status === 'done') set({ shown: job, active: null })
      else set({ active: null, runError: job.error || 'segmentation failed' })
    } catch (e: any) { if (token === runToken) set({ active: null, runError: e.message || String(e) }) }
  },
  serialize: () => pickView(get()),
  restore: (p) => {
    runToken++
    const { view, notice } = allelicFromSession(p)
    set({ ...view, notice, shown: null, active: null, runError: null, tracks: null, tracksKey: null, tracksError: null })
  },
}))
