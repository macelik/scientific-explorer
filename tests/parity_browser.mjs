// Runs the browser statistics module (bundled by esbuild) against Python reference values.
import fs from 'node:fs'
import path from 'node:path'
import { fileURLToPath } from 'node:url'
const here = path.dirname(fileURLToPath(import.meta.url))
const mod = await import(path.join(here, '_bundle', 'stats_bundle.mjs'))
const { summarize, acf, adSamples, xSamples, rollingMean, log2Track, findCompetingPeaks, selectionCoords, boundaryMb } = mod
const ref = JSON.parse(fs.readFileSync(path.join(here, 'parity_reference.json'), 'utf8'))
const base = process.env.API || 'http://127.0.0.1:8765'
let fails = 0, checks = 0
const close = (a, b, tol = 1e-9, label = '') => {
  checks++
  if (a === null && b === null) return
  if (a === null || b === null || Math.abs(a - b) > tol * Math.max(1, Math.abs(b))) { fails++; console.log('FAIL', label, a, b) }
}
for (const c of ref) {
  const r = await fetch(`${base}/api/node`, { method: 'POST', headers: { 'content-type': 'application/json' }, body: JSON.stringify({ cell: c.cell, chrom: c.chrom, start: c.start, end: c.end }) })
  const node = await r.json()
  close(node.argmax_b, c.argmax_b, 0, 'argmax_b')
  for (const w of c.windows) {
    const ad = adSamples(node, w.s, w.e), xs = xSamples(node, w.s, w.e)
    if (JSON.stringify(ad.j) !== JSON.stringify(w.ad_j)) { fails++; console.log('FAIL ad membership', c.cell, w.s, w.e, ad.j.slice(0, 3), w.ad_j.slice(0, 3)) }
    checks++
    const sa = summarize(ad.values), sx = summarize(xs.values)
    for (const k of ['n', 'mean', 'median', 'std', 'cv', 'lag1']) { close(sa[k], w.ad_stats[k], 1e-6, `ad.${k} [${w.s},${w.e})`); close(sx[k], w.x_stats[k], 1e-6, `x.${k} [${w.s},${w.e})`) }
    const a = acf(ad.values, 10)
    a.values.forEach((v, i) => close(v, w.ad_acf[i], 1e-6, `acf lag${i + 1} [${w.s},${w.e})`))
    const co = selectionCoords(node, w.s, w.e)
    close(co.spanBp, w.span_bp, 0, 'span'); close(co.retainedBp, w.retained_bp, 0, 'retained')
  }
  const roll = rollingMean(node.x, 25, 8)
  roll.forEach((v, i) => close(v, c.rolling[i], 1e-5, `rolling[${i}]`))
  const lfc = log2Track(roll, c.baseline_trimmed_genome, 4)
  lfc.values.forEach((v, i) => close(v, c.log2fc[i], 1e-5, `log2fc[${i}]`))
  const peaks = findCompetingPeaks(node, c.argmax_b, 5, 6.0, null, (b) => boundaryMb(node, b))
  const pj = peaks.map((p) => p.b), rj = c.peaks.map((p) => p.b)
  checks++
  if (JSON.stringify(pj) !== JSON.stringify(rj)) { fails++; console.log('FAIL peaks', c.cell, c.start, pj, rj) }
  console.log(`${c.cell} ${c.chrom} [${c.start},${c.end}) checked`)
}
console.log(`parity: ${checks - fails}/${checks} checks passed`)
process.exit(fails ? 1 : 0)
