// tests/allelic_store.ts
// Run with: web/node_modules/.bin/esbuild tests/allelic_store.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/allelic_store.mjs
//           node tests/_bundle/allelic_store.mjs
import { useAllelic, timing } from '../web/src/allelicStore'
import { DEFAULT_SETTINGS, DEFAULT_VIEW } from '../web/src/allelic'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}
timing.pollMs = 1
const calls: string[] = []
const jobs: Record<string, any[]> = {}
const reply = (body: unknown, status = 200) => ({ ok: status < 400, status, json: async () => body }) as any
;(globalThis as any).fetch = async (url: string, init?: any) => {
  calls.push(`${init?.method ?? 'GET'} ${url}`)
  if (url.startsWith('/api/allelic/status')) return reply({ available: true, state: 'ready', problems: [], notes: [], allelic_dir: '/x', sources: ['cluster0', 'all'] })
  if (url.startsWith('/api/allelic/tracks')) return reply({ chrom: 'chr8', signal: 'abs_dbaf', min_at: 10, n: 2, tracks: [] })
  if (url === '/api/allelic/segment') { const sig = JSON.parse(init.body).signal; return reply(jobs[sig].shift()) }
  const id = url.split('/').pop()!
  const sig = Object.keys(jobs).find((k) => k.startsWith(id.slice(0, 1)))!
  return reply(jobs[sig].shift())
}
const job = (id: string, signal: string, status: string) => ({ id, status, signal, settings: DEFAULT_SETTINGS, progress: null, error: null, result: status === 'done' ? { breakpoints: [], audit: [], segments: [] } : null })

;(async () => {
  const s = useAllelic.getState()
  await s.loadStatus()
  eq(useAllelic.getState().status?.state, 'ready', 'status loaded')
  await useAllelic.getState().ensureTracks()
  await useAllelic.getState().ensureTracks()
  eq(calls.filter((c) => c.includes('/tracks')).length, 1, 'tracks fetched once per key')

  jobs.abs_dbaf = [job('a1', 'abs_dbaf', 'queued'), job('a1', 'abs_dbaf', 'running'), job('a1', 'abs_dbaf', 'done')]
  await useAllelic.getState().run()
  eq([useAllelic.getState().shown?.id, useAllelic.getState().active], ['a1', null], 'run polls until done')

  // run supersedes: a second run started while the first polls; the first must not overwrite shown
  jobs.dbaf = [job('d1', 'dbaf', 'running'), job('d1', 'dbaf', 'running'), job('d1', 'dbaf', 'done')]
  jobs.shrunk_abs_dbaf = [job('s1', 'shrunk_abs_dbaf', 'done')]
  useAllelic.getState().setView({ signal: 'dbaf' })
  const first = useAllelic.getState().run()
  useAllelic.getState().setView({ signal: 'shrunk_abs_dbaf' })
  await useAllelic.getState().run()
  await first
  eq(useAllelic.getState().shown?.id, 's1', 'run supersedes')

  jobs.abs_dbaf = [{ ...job('a2', 'abs_dbaf', 'error'), error: 'kaput' }]
  useAllelic.getState().setView({ signal: 'abs_dbaf', settings: { ...DEFAULT_SETTINGS, alpha: 0.01 } })
  await useAllelic.getState().run()
  eq([useAllelic.getState().runError, useAllelic.getState().shown?.id], ['kaput', 's1'], 'error keeps the previous result')

  useAllelic.getState().setView({ settings: { ...DEFAULT_SETTINGS, alpha: 0 } })
  await useAllelic.getState().run()
  eq(useAllelic.getState().runError, 'alpha: must be > 0 and ≤ 0.1', 'invalid settings never reach the server')

  const saved = useAllelic.getState().serialize()
  eq(Object.keys(saved).sort(), Object.keys(DEFAULT_VIEW).sort(), 'serialize writes the view only')
  useAllelic.getState().restore({ signal: 'nope' })
  eq([useAllelic.getState().signal, useAllelic.getState().shown, useAllelic.getState().notice?.startsWith('Saved allelic settings could not be applied')], ['abs_dbaf', null, true], 'restore: defaults, result cleared, notice')

  if (failures) { console.error(`${failures} failure(s)`); process.exit(1) } else console.log('all allelic store checks passed')
})()
