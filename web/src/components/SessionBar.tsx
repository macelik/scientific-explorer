import { useEffect, useState } from 'react'
import { api } from '../api'
import { useStore } from '../store'

export default function SessionBar() {
  const serialize = useStore((s) => s.serialize)
  const restore = useStore((s) => s.restore)
  const meta = useStore((s) => s.meta)!
  const setToast = useStore((s) => s.setToast)
  const [list, setList] = useState<{ name: string; saved_at?: string; dataset_hash?: string; cell?: string; chrom?: string }[]>([])
  const [pick, setPick] = useState('')
  const refresh = () => api.sessions().then(setList).catch(() => undefined)
  useEffect(() => { refresh() }, [])
  const save = async () => {
    const name = prompt('Session name (letters, digits, . _ -):', `session-${new Date().toISOString().slice(0, 16).replace(/[:T]/g, '-')}`)
    if (!name) return
    try { await api.saveSession(name, serialize()); setToast(`Saved session "${name}"`); refresh() } catch (e: any) { setToast(`Save failed: ${e.message}`) }
  }
  const load = async () => {
    if (!pick) return
    try { const d = await api.loadSession(pick); const err = await restore(d.state, d.dataset?.hash); setToast(err ? err : `Loaded session "${pick}" (saved ${d.saved_at})`) } catch (e: any) { setToast(`Load failed: ${e.message}`) }
  }
  const download = () => {
    const doc = { format: 'pyepi-explorer-session', version: 1, saved_at: new Date().toISOString(), dataset: meta.dataset.identity, state: serialize() }
    const blob = new Blob([JSON.stringify(doc)], { type: 'application/json' }); const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = 'pyepi-explorer-session.json'; a.click()
  }
  const upload = (f: File | null) => {
    if (!f) return
    f.text().then(async (t) => { const d = JSON.parse(t); const err = await restore(d.state ?? d, d.dataset?.hash); setToast(err ? err : `Loaded session from ${f.name}`) }).catch((e) => setToast(`Upload failed: ${e.message}`))
  }
  return (
    <span className="session-inner">
      <button className="btn-sm" onClick={save}>save session</button>
      <select value={pick} onChange={(e) => setPick(e.target.value)}><option value="">saved sessions…</option>{list.map((s) => <option key={s.name} value={s.name}>{s.name}{s.cell ? ` · ${s.cell}/${s.chrom}` : ''}{s.dataset_hash && s.dataset_hash !== meta.dataset.identity.hash ? ' (other dataset)' : ''}</option>)}</select>
      <button className="btn-sm" onClick={load} disabled={!pick}>load</button>
      <button className="btn-sm" onClick={download}>download JSON</button>
      <label className="btn-sm file">upload JSON<input type="file" accept="application/json" onChange={(e) => upload(e.target.files?.[0] ?? null)} /></label>
    </span>
  )
}
