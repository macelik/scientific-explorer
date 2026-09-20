// tests/hatch_merge_helpers.ts
// Run with: web/node_modules/.bin/esbuild tests/hatch_merge_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/hatch_merge_helpers.mjs
//           node tests/_bundle/hatch_merge_helpers.mjs
import { classifyTransition, classifyMergeProposal, classifyAmbiguous, weakerSide } from '../web/src/hatchMerge'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}

eq(classifyTransition(1, -1), true, 'transition opposite sign')
eq(classifyTransition(1, 1), false, 'transition same sign')
eq(classifyTransition(0, -1), false, 'transition zero is not opposite')
eq(classifyTransition(null, -1), false, 'transition missing left')
eq(classifyTransition(NaN, 1), false, 'transition NaN')

eq(weakerSide(-1, -3), 'left', 'weaker side both negative picks smaller magnitude')
eq(weakerSide(-2, 2), 'tie', 'weaker side exact tie')
eq(weakerSide(null, 1), null, 'weaker side missing')

eq(classifyMergeProposal(0.5, 5, 1), { eligibleLeft: true, eligibleRight: false, weakerSide: 'left' }, 'merge proposal one side')
eq(classifyMergeProposal(1, 5, 1), { eligibleLeft: false, eligibleRight: false, weakerSide: null }, 'merge proposal strict threshold')
eq(classifyMergeProposal(0.5, 0.5, 1), { eligibleLeft: true, eligibleRight: true, weakerSide: 'both' }, 'merge proposal tie -> both')
eq(classifyMergeProposal(0.2, null, 1), { eligibleLeft: true, eligibleRight: false, weakerSide: 'left' }, 'merge proposal terminal segment')

eq(classifyAmbiguous(1, 3, 3, 1), 'ambiguous', 'ambiguous differing weaker flanks')
eq(classifyAmbiguous(1, 3, 1, 3), 'consistent', 'ambiguous consistent weaker flanks')
eq(classifyAmbiguous(2, 2, 1, 3), 'no_unique_preference', 'ambiguous exact SRD tie')
eq(classifyAmbiguous(1, 3, null, 3), 'undefined', 'ambiguous missing input')

if (failures) { console.error(`${failures} failure(s)`); process.exit(1) }
console.log('all hatch_merge_helpers checks passed')
