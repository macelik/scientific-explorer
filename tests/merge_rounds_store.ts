// Run with: web/node_modules/.bin/esbuild tests/merge_rounds_store.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rounds_store.mjs
//           node tests/_bundle/merge_rounds_store.mjs
import { useIntegration } from '../web/src/integrationStore'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}
const r0 = { id: 'rule0', type: 'flank_threshold', tau: 0.5, srdRef: 8, combine: 'any' }
const round = { rules: [r0], allow_gap_crossing: false, label: 'x', merges: 1, noise: null }
let release: () => void = () => {}
;(globalThis as any).fetch = () => new Promise((res) => { release = () => res({ ok: true, status: 200, json: async () => ({ pinned: 'p' }) }) })

const st: any = { id: 'job', status: 'done', settings: { rules: [r0], allow_gap_crossing: false, rounds: [round] }, merge: { pairs: [{ merges: 2 }] }, diff: null }
useIntegration.setState({ rounds: [round], recall: st, recallKaryo: new Uint8Array(1) })
const p = useIntegration.getState().adoptRecall()
useIntegration.getState().undoRound()
release()
await p
eq(useIntegration.getState().rounds, [], 'undo during pin is not resurrected')
eq(useIntegration.getState().baseKaryo, null, 'no stale karyo adopted')
if (failures) { console.error(`${failures} failure(s)`); process.exit(1) } else console.log('all merge_rounds_store checks passed')
