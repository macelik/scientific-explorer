import { create } from 'zustand'
import { api } from './api'
import { useExperimentStore } from './experimentStore'
import { useIntegration } from './integrationStore'
import { clampWindow } from './coords'
import { windowColor } from './colors'
import type {
  AssignResult, BaselineMethod, CellData, EventCandidate, EventRule, ManualNode, ManualSession, ManualTree, Meta,
  NodeData, Selection, Snapshot, SnapshotWindow, TestState, SegmentEstimator,
} from './types'

export type Tab = 'cohort' | 'explore' | 'manual' | 'experiment' | 'clusters' | 'clusterseg'
export type Mode = 'navigate' | 'select' | 'place'
export type OrderMode = 'clustered' | 'cov_asc' | 'cov_desc' | 'original'

export interface Filters {
  search: string; terciles: string[]; covMin: number | null; covMax: number | null; order: OrderMode; rulesMode: 'any' | 'all'
}
export interface Display {
  baseline: BaselineMethod; rolling: number; minPeriods: number; sharedAdScale: boolean; showMean: boolean; showMedian: boolean
  showLocal: boolean; showGenome: boolean; clip: number; activeRegionId: string | null; showProduction: boolean; showRegions: boolean
}

export const nodeKey = (cell: string, chrom: string, start: number, end: number) => `${cell}|${chrom}|${start}|${end}`
export const testKey = (cell: string, chrom: string, start: number, end: number, b: number) => `${cell}|${chrom}|${start}|${end}|${b}`
export const manualKey = (cell: string, chrom: string) => `${cell}|${chrom}`

let uid = 0
export const newId = (p = 'id') => `${p}${Date.now().toString(36)}${(uid++).toString(36)}`

export interface State {
  meta: Meta | null; bins: Int32Array | null; karyo: Uint8Array | null; bootError: string | null; booted: boolean
  bootSteps: { label: string; status: 'pending' | 'running' | 'done' | 'error'; detail?: string }[]
  presentation: boolean
  tab: Tab
  filters: Filters; rules: EventRule[]; events: { status: string; rows: EventCandidate[] | null; done: number; total: number }
  karyoView: { binLo: number; binHi: number } | null
  cohortDrawerOpen: boolean
  cell: string | null; chrom: string | null
  cells: Record<string, CellData>; nodes: Record<string, NodeData>; nodeLoading: Record<string, boolean>; nodeErrors: Record<string, string>
  focus: 'root' | 'left' | 'right' | null
  display: Display; mode: Mode
  selections: Selection[]; snapshots: Snapshot[]; savedWindows: Selection[]; activeSelectionId: string | null
  compareView: 'overlay' | 'small'; compareCoords: 'offset' | 'absolute'; showLeftRight: boolean
  tests: Record<string, TestState>
  manual: Record<string, ManualSession>
  toast: string | null
  timings: { label: string; ms: number; at: number }[]

  boot: () => Promise<void>
  setPresentation: (on: boolean) => void
  setTab: (t: Tab) => void
  setFilters: (f: Partial<Filters>) => void
  setRules: (r: EventRule[]) => void
  setKaryoView: (v: { binLo: number; binHi: number } | null) => void
  setCohortDrawer: (open: boolean) => void
  openCell: (cell: string, chrom: string, opts?: { tab?: Tab }) => Promise<void>
  ensureCell: (cell: string) => Promise<CellData>
  ensureNode: (cell: string, chrom: string, start: number, end: number) => Promise<NodeData>
  runTest: (cell: string, chrom: string, start: number, end: number, b: number) => Promise<void>
  setFocus: (f: State['focus']) => void
  setDisplay: (d: Partial<Display>) => void
  setMode: (m: Mode) => void
  addSelection: (sel: Omit<Selection, 'id' | 'color' | 'createdAt' | 'visible' | 'name'> & Partial<Pick<Selection, 'name' | 'color'>>) => Selection | null
  updateSelection: (id: string, patch: Partial<Selection>) => void
  removeSelection: (id: string) => void
  clearSelections: (nodeId?: string) => void
  setActiveSelection: (id: string | null) => void
  setCompareView: (v: 'overlay' | 'small') => void
  setCompareCoords: (v: 'offset' | 'absolute') => void
  setShowLeftRight: (v: boolean) => void
  addSnapshot: (s: Snapshot) => void
  removeSnapshot: (id: string) => void
  saveWindow: (sel: Selection) => void
  removeSavedWindow: (id: string) => void
  reapplyWindow: (saved: Selection, node: NodeData, nodeRef: { id: string; version: number; side: string }) => { added: Selection | null; clipped: number; report: string }
  initManual: (cell: string, chrom: string) => Promise<void>
  manualPlace: (cell: string, chrom: string, nodeId: string, b: number) => void
  manualRemove: (cell: string, chrom: string, nodeId: string) => void
  manualUndo: (cell: string, chrom: string) => void
  manualRedo: (cell: string, chrom: string) => void
  manualReset: (cell: string, chrom: string) => void
  manualSetActive: (cell: string, chrom: string, nodeId: string) => void
  manualClearNotice: (cell: string, chrom: string) => void
  manualFinalize: (cell: string, chrom: string) => Promise<void>
  manualSetEstimator: (cell: string, chrom: string, estimator: SegmentEstimator) => void
  manualContinue: (cell: string, chrom: string) => void
  serialize: () => any
  restore: (payload: any, datasetHash?: string) => Promise<string | null>
  setToast: (t: string | null) => void
  pushTiming: (label: string, ms: number) => void
}

const defaultFilters: Filters = { search: '', terciles: ['low', 'mid', 'high'], covMin: null, covMax: null, order: 'clustered', rulesMode: 'any' }
const defaultDisplay: Display = {
  baseline: 'trimmed', rolling: 25, minPeriods: 8, sharedAdScale: false, showMean: true, showMedian: true, showLocal: true,
  showGenome: true, clip: 4, activeRegionId: null, showProduction: true, showRegions: true,
}

function seedRules(meta: Meta): EventRule[] {
  const r = (id: string) => meta.regions.find((x) => x.id === id)!
  const cn = (id: string, name: string, states: number[]): EventRule => {
    const reg = r(id)
    return { id: newId('rule'), name, enabled: false, family: 'cn', chrom: reg.chrom, startMb: reg.start_mb, endMb: reg.end_mb, states, minFraction: 0, node: 'any', acceptance: 'any' }
  }
  const bp = (id: string, name: string, node: EventRule['node'], acc: EventRule['acceptance']): EventRule => {
    const reg = r(id)
    return { id: newId('rule'), name, enabled: false, family: 'bp', chrom: reg.chrom, startMb: reg.start_mb, endMb: reg.end_mb, states: [], minFraction: 0, node, acceptance: acc }
  }
  return [
    cn('9p21', 'Loss calls overlap 9p21', [0, 0.5]),
    cn('chr13-gain', 'Gain calls overlap chr13 94–98 Mb', [1.5, 2]),
    cn('chr1-distal-gain', 'Gain calls overlap chr1 distal 195–249 Mb', [1.5, 2]),
    bp('chr9-transition', 'Depth-0 candidate near chr9 transition (27–33 Mb)', 'root', 'candidate'),
    bp('9p21', 'Any accepted breakpoint inside 9p21', 'any', 'accepted'),
    bp('chr9-shoulder', 'Left-child candidate at chr9 shoulder (6.8–7.4 Mb)', 'left', 'candidate'),
  ]
}

function newTree(cell: string, chrom: string, n: number): ManualTree {
  const root: ManualNode = { id: 'n0', cell, chrom, start: 0, end: n, depth: 0, version: 1, parentId: null, side: 'root', splitB: null, leftId: null, rightId: null, label: 'root' }
  return { rootId: 'n0', nodes: { n0: root }, nextId: 1, version: 1 }
}
function cloneTree(t: ManualTree): ManualTree {
  const nodes: Record<string, ManualNode> = {}
  for (const k in t.nodes) nodes[k] = { ...t.nodes[k] }
  return { rootId: t.rootId, nodes, nextId: t.nextId, version: t.version + 1 }
}
function removeSubtree(t: ManualTree, id: string) {
  const n = t.nodes[id]
  if (!n) return
  if (n.leftId) { removeSubtree(t, n.leftId); delete t.nodes[n.leftId] }
  if (n.rightId) { removeSubtree(t, n.rightId); delete t.nodes[n.rightId] }
  n.leftId = null; n.rightId = null; n.splitB = null
}
export function leafBoundaries(t: ManualTree): number[] {
  const out: number[] = []
  const walk = (id: string) => { const n = t.nodes[id]; if (n.splitB !== null) { out.push(n.start + n.splitB); walk(n.leftId!); walk(n.rightId!) } }
  walk(t.rootId)
  return out.sort((a, b) => a - b)
}
export function leaves(t: ManualTree): ManualNode[] {
  const out: ManualNode[] = []
  const walk = (id: string) => { const n = t.nodes[id]; if (n.splitB === null) out.push(n); else { walk(n.leftId!); walk(n.rightId!) } }
  walk(t.rootId)
  return out
}

export const useStore = create<State>((set, get) => ({
  meta: null, bins: null, karyo: null, bootError: null, booted: false, bootSteps: [], presentation: false,
  tab: 'cohort',
  filters: defaultFilters, rules: [], events: { status: 'idle', rows: null, done: 0, total: 0 },
  karyoView: null, cohortDrawerOpen: false,
  cell: null, chrom: null,
  cells: {}, nodes: {}, nodeLoading: {}, nodeErrors: {},
  focus: null,
  display: defaultDisplay, mode: 'navigate',
  selections: [], snapshots: [], savedWindows: [], activeSelectionId: null,
  compareView: 'overlay', compareCoords: 'offset', showLeftRight: false,
  tests: {},
  manual: {},
  toast: null,
  timings: [],

  setPresentation: (presentation) => { set({ presentation, toast: null }); document.body.classList.toggle('presentation', presentation) },
  boot: async () => {
    const steps: State['bootSteps'] = [
      { label: 'server reachable', status: 'pending' },
      { label: 'cohort metadata (cells, chromosomes, cluster order)', status: 'pending' },
      { label: 'bin coordinates (25k bins)', status: 'pending' },
      { label: 'categorical CN matrix (7.6 MB karyogram)', status: 'pending' },
    ]
    const setStep = (i: number, status: 'pending' | 'running' | 'done' | 'error', detail?: string) => { steps[i] = { ...steps[i], status, detail }; set({ bootSteps: [...steps] }) }
    set({ bootSteps: [...steps], bootError: null })
    // wait for the service (it needs a few seconds after launch to import and load the data); retry until ready
    setStep(0, 'running')
    let attempt = 0
    for (;;) {
      try {
        const r = await fetch('/api/health', { cache: 'no-store' })
        if (r.ok) { const h = await r.json(); if (h.loaded) break }
        throw new Error(`server answered ${r.status}`)
      } catch (e: any) {
        attempt++
        setStep(0, 'running', `waiting for the server on this port (attempt ${attempt}: ${e.message || e})`)
        await new Promise((res) => setTimeout(res, 1500))
      }
    }
    setStep(0, 'done')
    try {
      const t0 = performance.now()
      const timed = async <T,>(i: number, fn: () => Promise<T>): Promise<T> => { setStep(i, 'running'); const t = performance.now(); try { const v = await fn(); setStep(i, 'done', `${Math.round(performance.now() - t)} ms`); return v } catch (e: any) { setStep(i, 'error', String(e.message || e)); throw e } }
      const [meta, bins, karyo] = await Promise.all([timed(1, api.meta), timed(2, api.bins), timed(3, api.karyogram)])
      set({ meta, bins, karyo, rules: seedRules(meta), booted: true })
      get().pushTiming('boot: meta + bins + karyogram', performance.now() - t0)
      // background event indexing (pending classifications stay pending)
      const poll = async () => {
        try {
          const st = await api.eventsStatus()
          if (st.status === 'done') {
            const c = await api.eventsCandidates()
            set({ events: { status: 'done', rows: c.rows, done: st.done, total: st.total } })
            return
          }
          set({ events: { status: st.status, rows: null, done: st.done, total: st.total } })
          if (st.status === 'idle') await api.eventsStart()
          if (st.status !== 'error') setTimeout(poll, 2000)
        } catch (e) { set({ events: { status: 'error', rows: null, done: 0, total: 0 } }) }
      }
      poll()
    } catch (e: any) { set({ bootError: String(e.message || e) }) }
  },
  setTab: (tab) => set({ tab }),
  setFilters: (f) => set({ filters: { ...get().filters, ...f } }),
  setRules: (rules) => set({ rules }),
  setKaryoView: (karyoView) => set({ karyoView }),
  setCohortDrawer: (cohortDrawerOpen) => set({ cohortDrawerOpen }),

  ensureCell: async (cell) => {
    const have = get().cells[cell]
    if (have) return have
    const t0 = performance.now()
    const d = await api.cell(cell)
    get().pushTiming(`cell ${cell}`, performance.now() - t0)
    set({ cells: { ...get().cells, [cell]: d } })
    return d
  },
  ensureNode: async (cell, chrom, start, end) => {
    const k = nodeKey(cell, chrom, start, end)
    const have = get().nodes[k]
    if (have) return have
    if (get().nodeLoading[k]) {
      // wait for in-flight
      return new Promise((resolve, reject) => {
        const iv = setInterval(() => {
          const n = get().nodes[k]; const err = get().nodeErrors[k]
          if (n) { clearInterval(iv); resolve(n) } else if (err) { clearInterval(iv); reject(new Error(err)) }
        }, 30)
      })
    }
    set({ nodeLoading: { ...get().nodeLoading, [k]: true } })
    try {
      const t0 = performance.now()
      const d = await api.node(cell, chrom, start, end)
      get().pushTiming(`node ${chrom}[${start},${end})${d.cached ? ' (server cache)' : ''}`, performance.now() - t0)
      set({ nodes: { ...get().nodes, [k]: d }, nodeLoading: { ...get().nodeLoading, [k]: false } })
      return d
    } catch (e: any) {
      set({ nodeLoading: { ...get().nodeLoading, [k]: false }, nodeErrors: { ...get().nodeErrors, [k]: String(e.message || e) } })
      throw e
    }
  },
  openCell: async (cell, chrom, opts) => {
    const t0 = performance.now()
    set({ cell, chrom, focus: null, tab: opts?.tab ?? (get().tab === 'cohort' ? 'explore' : get().tab), activeSelectionId: null })
    const meta = get().meta!
    const c = meta.chromosomes.find((x) => x.name === chrom)!
    await get().ensureCell(cell)
    const root = await get().ensureNode(cell, chrom, 0, c.n)
    if (root.argmax_b !== null) {
      await Promise.all([get().ensureNode(cell, chrom, 0, root.argmax_b), get().ensureNode(cell, chrom, root.argmax_b, c.n)])
    }
    get().pushTiming(`open ${cell}/${chrom} (cell + root + both children)`, performance.now() - t0)
  },
  runTest: async (cell, chrom, start, end, b) => {
    const k = testKey(cell, chrom, start, end, b)
    if (get().tests[k] && get().tests[k].status !== 'error') return
    const meta = get().meta!
    set({ tests: { ...get().tests, [k]: { status: 'pending' } } })
    try {
      const t0 = performance.now()
      let r = await api.test(cell, chrom, start, end, b, meta.defaults.n_permutations, meta.defaults.seed)
      while (r.status === 'running' && r.id) {
        await new Promise((res) => setTimeout(res, 250))
        const j = await api.job(r.id)
        r = { ...r, status: j.status, result: j.result, error: j.error ?? undefined, elapsed: j.elapsed }
      }
      get().pushTiming(`test ${chrom}[${start},${end}) b=${b}${r.cached ? ' (cached)' : ''}`, performance.now() - t0)
      if (r.status === 'done' && r.result) set({ tests: { ...get().tests, [k]: { status: 'done', result: r.result, elapsed: r.elapsed } } })
      else set({ tests: { ...get().tests, [k]: { status: 'error', error: r.error || 'unknown error' } } })
    } catch (e: any) {
      set({ tests: { ...get().tests, [k]: { status: 'error', error: String(e.message || e) } } })
    }
  },
  setFocus: (focus) => set({ focus }),
  setDisplay: (d) => set({ display: { ...get().display, ...d } }),
  setMode: (mode) => set({ mode }),

  addSelection: (sel) => {
    const n = sel.nodeEnd - sel.nodeStart
    const cw = clampWindow(sel.s, sel.e, n)
    if (cw.e <= cw.s) return null
    const dup = get().selections.find((x) => x.nodeId === sel.nodeId && x.nodeVersion === sel.nodeVersion && x.s === cw.s && x.e === cw.e && x.anchor === sel.anchor && x.cell === sel.cell && x.chrom === sel.chrom)
    if (dup) return dup
    const idx = get().selections.length
    const s: Selection = {
      ...sel,
      id: newId('w'), name: sel.name ?? `W${idx + 1}`, color: sel.color ?? windowColor(idx), visible: true, createdAt: Date.now(),
      s: cw.s, e: cw.e, clipped: cw.clipped.left || cw.clipped.right ? cw.clipped : undefined,
    }
    set({ selections: [...get().selections, s], activeSelectionId: s.id })
    return s
  },
  updateSelection: (id, patch) => set({ selections: get().selections.map((s) => (s.id === id ? { ...s, ...patch } : s)) }),
  removeSelection: (id) => set({ selections: get().selections.filter((s) => s.id !== id), activeSelectionId: get().activeSelectionId === id ? null : get().activeSelectionId }),
  clearSelections: (nodeId) => set({ selections: nodeId ? get().selections.filter((s) => s.nodeId !== nodeId) : [] }),
  setActiveSelection: (activeSelectionId) => set({ activeSelectionId }),
  setCompareView: (compareView) => set({ compareView }),
  setCompareCoords: (compareCoords) => set({ compareCoords }),
  setShowLeftRight: (showLeftRight) => set({ showLeftRight }),
  addSnapshot: (s) => set({ snapshots: [...get().snapshots, s] }),
  removeSnapshot: (id) => set({ snapshots: get().snapshots.filter((s) => s.id !== id) }),
  saveWindow: (sel) => set({ savedWindows: [...get().savedWindows, { ...sel, id: newId('saved'), source: 'saved' }] }),
  removeSavedWindow: (id) => set({ savedWindows: get().savedWindows.filter((s) => s.id !== id) }),
  reapplyWindow: (saved, node, ref) => {
    // genomic window (Mb from the saved selection's chromosome-local bounds) applied to the new node by bin overlap
    const savedChromStart = saved.nodeStart + saved.s, savedChromEnd = saved.nodeStart + saved.e
    if (saved.chrom !== node.chrom) return { added: null, clipped: 0, report: `saved window is on ${saved.chrom}, target node is on ${node.chrom}` }
    const s = savedChromStart - node.start, e = savedChromEnd - node.start
    const cw = clampWindow(s, e, node.n)
    const clipped = cw.clipped.left + cw.clipped.right
    if (cw.e <= cw.s) return { added: null, clipped, report: 'saved window lies entirely outside the target node' }
    const anchor = saved.anchor !== null ? saved.anchor + saved.nodeStart - node.start : null
    const added = get().addSelection({
      cell: node.cell, chrom: node.chrom, nodeId: ref.id, nodeVersion: ref.version, nodeStart: node.start, nodeEnd: node.end, nodeSide: ref.side,
      s: cw.s, e: cw.e, anchor: anchor !== null && anchor > 0 && anchor < node.n ? anchor : null, source: 'reapplied',
      name: `${saved.name} (reapplied)`, color: saved.color,
      note: `reapplied from ${saved.cell}/${saved.chrom} ${saved.nodeSide} node; ${clipped ? `${clipped} bins clipped` : 'no clipping'}`,
    })
    return { added, clipped, report: clipped ? `${clipped} bin(s) clipped at the node boundary` : 'no clipping' }
  },

  initManual: async (cell, chrom) => {
    const k = manualKey(cell, chrom)
    if (get().manual[k]) return
    const meta = get().meta!
    const c = meta.chromosomes.find((x) => x.name === chrom)!
    await get().ensureCell(cell)
    await get().ensureNode(cell, chrom, 0, c.n)
    const tree = newTree(cell, chrom, c.n)
    set({ manual: { ...get().manual, [k]: { cell, chrom, tree, history: [], future: [], activeNodeId: 'n0', notice: null, cnResult: null, cnOutdated: false, cnBoundaries: null } } })
  },
  manualPlace: (cell, chrom, nodeId, b) => {
    const k = manualKey(cell, chrom)
    const ms = get().manual[k]; if (!ms) return
    const t = cloneTree(ms.tree)
    const n = t.nodes[nodeId]
    if (!n || b <= 0 || b >= n.end - n.start) return
    let notice: string | null = null
    if (n.splitB !== null) {
      const removed = Object.keys(t.nodes).length
      removeSubtree(t, nodeId)
      notice = `Replaced boundary of ${n.label}: ${removed - Object.keys(t.nodes).length} descendant node(s) discarded and rebuilt. Undo restores the previous tree.`
    }
    n.splitB = b; n.version += 1
    const li = `n${t.nextId++}`, ri = `n${t.nextId++}`
    t.nodes[li] = { id: li, cell, chrom, start: n.start, end: n.start + b, depth: n.depth + 1, version: 1, parentId: n.id, side: 'left', splitB: null, leftId: null, rightId: null, label: `${n.label === 'root' ? '' : n.label + '·'}L` }
    t.nodes[ri] = { id: ri, cell, chrom, start: n.start + b, end: n.end, depth: n.depth + 1, version: 1, parentId: n.id, side: 'right', splitB: null, leftId: null, rightId: null, label: `${n.label === 'root' ? '' : n.label + '·'}R` }
    n.leftId = li; n.rightId = ri
    set({ manual: { ...get().manual, [k]: { ...ms, tree: t, history: [...ms.history, ms.tree], future: [], notice, cnOutdated: ms.cnResult !== null, activeNodeId: nodeId } } })
    get().ensureNode(cell, chrom, n.start, n.start + b); get().ensureNode(cell, chrom, n.start + b, n.end)
    get().runTest(cell, chrom, n.start, n.end, b)
  },
  manualRemove: (cell, chrom, nodeId) => {
    const k = manualKey(cell, chrom)
    const ms = get().manual[k]; if (!ms) return
    const t = cloneTree(ms.tree)
    const n = t.nodes[nodeId]; if (!n || n.splitB === null) return
    removeSubtree(t, nodeId); n.version += 1
    set({ manual: { ...get().manual, [k]: { ...ms, tree: t, history: [...ms.history, ms.tree], future: [], notice: `Removed split and subtree of ${n.label}.`, cnOutdated: ms.cnResult !== null, activeNodeId: nodeId } } })
  },
  manualUndo: (cell, chrom) => {
    const k = manualKey(cell, chrom); const ms = get().manual[k]; if (!ms || !ms.history.length) return
    const prev = ms.history[ms.history.length - 1]
    const active = prev.nodes[ms.activeNodeId] ? ms.activeNodeId : prev.rootId
    set({ manual: { ...get().manual, [k]: { ...ms, tree: prev, history: ms.history.slice(0, -1), future: [ms.tree, ...ms.future], notice: null, cnOutdated: ms.cnResult !== null, activeNodeId: active } } })
  },
  manualRedo: (cell, chrom) => {
    const k = manualKey(cell, chrom); const ms = get().manual[k]; if (!ms || !ms.future.length) return
    const next = ms.future[0]
    const active = next.nodes[ms.activeNodeId] ? ms.activeNodeId : next.rootId
    set({ manual: { ...get().manual, [k]: { ...ms, tree: next, history: [...ms.history, ms.tree], future: ms.future.slice(1), notice: null, cnOutdated: ms.cnResult !== null, activeNodeId: active } } })
  },
  manualReset: (cell, chrom) => {
    const k = manualKey(cell, chrom); const ms = get().manual[k]; if (!ms) return
    const meta = get().meta!; const c = meta.chromosomes.find((x) => x.name === chrom)!
    set({ manual: { ...get().manual, [k]: { ...ms, tree: newTree(cell, chrom, c.n), history: [...ms.history, ms.tree], future: [], notice: 'Tree reset to one unsplit chromosome.', cnOutdated: ms.cnResult !== null, activeNodeId: 'n0' } } })
  },
  manualSetActive: (cell, chrom, nodeId) => {
    const k = manualKey(cell, chrom); const ms = get().manual[k]; if (!ms || !ms.tree.nodes[nodeId]) return
    set({ manual: { ...get().manual, [k]: { ...ms, activeNodeId: nodeId } } })
  },
  manualClearNotice: (cell, chrom) => {
    const k = manualKey(cell, chrom); const ms = get().manual[k]; if (!ms) return
    set({ manual: { ...get().manual, [k]: { ...ms, notice: null } } })
  },
  manualSetEstimator: (cell, chrom, estimator) => {
    const k=manualKey(cell,chrom), ms=get().manual[k]; if(!ms) return
    set({manual:{...get().manual,[k]:{...ms,cnEstimator:estimator,cnOutdated:ms.cnResult!==null,notice:null}}})
  },
  manualFinalize: async (cell, chrom) => {
    const k = manualKey(cell, chrom); const ms = get().manual[k]; if (!ms) return
    const bounds = leafBoundaries(ms.tree)
    const estimator=ms.cnEstimator??'arithmetic', token=newId('cn')
    set({manual:{...get().manual,[k]:{...ms,cnBusy:true,cnToken:token,notice:null}}})
    const t0 = performance.now()
    try {
      const res: AssignResult = await api.assign(cell, chrom, bounds, estimator)
      get().pushTiming(`assign CN ${chrom} (${estimator}, ${bounds.length} boundaries)`, performance.now() - t0)
      const cur = get().manual[k]
      if(cur?.cnToken!==token) return
      set({ manual: { ...get().manual, [k]: { ...cur, cnResult: res, cnBusy:false, cnOutdated: (cur.cnEstimator??'arithmetic')!==estimator || JSON.stringify(leafBoundaries(cur.tree)) !== JSON.stringify(bounds), cnBoundaries: bounds } } })
    } catch(e:any) {
      const cur=get().manual[k]; if(cur?.cnToken!==token)return
      set({manual:{...get().manual,[k]:{...cur,cnBusy:false,notice:`CN assignment could not be computed: ${e.message}`}}})
    }
  },
  manualContinue: (cell, chrom) => { const k = manualKey(cell, chrom); const ms = get().manual[k]; if (!ms) return; set({ manual: { ...get().manual, [k]: { ...ms } } }) },

  serialize: () => {
    const s = get()
    const manual: Record<string, any> = {}
    for (const k in s.manual) { const m = s.manual[k]; manual[k] = { cell: m.cell, chrom: m.chrom, tree: m.tree, history: m.history, future: m.future, activeNodeId: m.activeNodeId, cnBoundaries: m.cnBoundaries, cnOutdated: m.cnOutdated, cnResult: m.cnResult, cnEstimator: m.cnEstimator ?? 'arithmetic' } }
    return {
      experiments: useExperimentStore.getState().entries,
      integration: useIntegration.getState().serialize(),
      tab: s.tab, cell: s.cell, chrom: s.chrom, filters: s.filters, rules: s.rules, display: s.display, karyoView: s.karyoView,
      selections: s.selections, snapshots: s.snapshots, savedWindows: s.savedWindows, focus: s.focus, compareView: s.compareView, compareCoords: s.compareCoords,
      manual, scientific: s.meta ? { n_permutations: s.meta.defaults.n_permutations, seed: s.meta.defaults.seed, alpha: s.meta.defaults.alpha, k_production: s.meta.defaults.production_k } : null,
      definitions: {
        selection: 'node-local half-open [s,e) of retained bins; AD boundary j included when max(s,1)<=j<e, value d[j-1]',
        stats: 'mean, median, population std (ddof=0), CV=std/mean (N/A if mean 0), lag-k = Pearson(v[:-k], v[k:]), descriptive only',
        cn: 'five-state 0=loss 0.5=putative loss 1=base/diploid 1.5=putative gain 2=gain',
      },
    }
  },
  restore: async (p, datasetHash) => {
    const s = get()
    if (datasetHash && s.meta && datasetHash !== s.meta.dataset.identity.hash) {
      return `Session was saved against dataset ${datasetHash}, but the loaded dataset is ${s.meta.dataset.identity.hash}. Not restoring (indices would not be comparable).`
    }
    useExperimentStore.getState().restore(p.experiments || {})
    useIntegration.getState().restore(p.integration)
    set({
      filters: { ...defaultFilters, ...(p.filters || {}) }, rules: p.rules || s.rules, display: { ...defaultDisplay, ...(p.display || {}) }, karyoView: p.karyoView ?? null,
      selections: p.selections || [], snapshots: p.snapshots || [], savedWindows: p.savedWindows || [], focus: p.focus ?? null,
      compareView: p.compareView || 'overlay', compareCoords: p.compareCoords || 'offset',
      manual: p.manual ? Object.fromEntries(Object.entries(p.manual).map(([k, m]: [string, any]) => [k, { ...m, notice: null, cnBusy: false, cnToken: undefined, cnEstimator: m.cnEstimator ?? 'arithmetic', cnResult: m.cnResult ?? null, cnOutdated: m.cnOutdated ?? false, cnBoundaries: m.cnBoundaries ?? null }])) : {},
    })
    if (p.cell && p.chrom) {
      await get().openCell(p.cell, p.chrom, { tab: p.tab || 'explore' })
      // preload manual nodes
      for (const k in get().manual) { const m = get().manual[k]; for (const id in m.tree.nodes) { const n = m.tree.nodes[id]; get().ensureNode(n.cell, n.chrom, n.start, n.end) } }
    } else set({ tab: p.tab || 'cohort' })
    return null
  },
  setToast: (toast) => set({ toast: get().presentation ? null : toast }),
  pushTiming: (label, ms) => set({ timings: [{ label, ms: Math.round(ms * 10) / 10, at: Date.now() }, ...get().timings].slice(0, 60) }),
}))

/** Explorer node references derived from the loaded root node. */
export function explorerNodes(state: State): { root: { ref: NodeRefLike; data: NodeData } | null; left: { ref: NodeRefLike; data: NodeData } | null; right: { ref: NodeRefLike; data: NodeData } | null } {
  const { cell, chrom, meta, nodes } = state
  if (!cell || !chrom || !meta) return { root: null, left: null, right: null }
  const c = meta.chromosomes.find((x) => x.name === chrom)!
  const root = nodes[nodeKey(cell, chrom, 0, c.n)]
  if (!root) return { root: null, left: null, right: null }
  const rootRef: NodeRefLike = { id: nodeKey(cell, chrom, 0, c.n), version: 0, side: 'root', depth: 0 }
  const out: any = { root: { ref: rootRef, data: root }, left: null, right: null }
  if (root.argmax_b !== null) {
    const l = nodes[nodeKey(cell, chrom, 0, root.argmax_b)], r = nodes[nodeKey(cell, chrom, root.argmax_b, c.n)]
    if (l) out.left = { ref: { id: nodeKey(cell, chrom, 0, root.argmax_b), version: 0, side: 'left', depth: 1 }, data: l }
    if (r) out.right = { ref: { id: nodeKey(cell, chrom, root.argmax_b, c.n), version: 0, side: 'right', depth: 1 }, data: r }
  }
  return out
}
export interface NodeRefLike { id: string; version: number; side: 'root' | 'left' | 'right'; depth: number }

export function makeSnapshotWindow(w: SnapshotWindow): SnapshotWindow { return w }
