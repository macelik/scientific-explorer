// web/src/coords.ts
var MB = 1e6;
function clampWindow(s, e, n) {
  const cs = Math.max(0, Math.min(n, s));
  const ce = Math.max(0, Math.min(n, e));
  return { s: cs, e: Math.max(ce, cs), clipped: { left: cs - s, right: e - ce } };
}
function peakWindow(b, n, half = 20) {
  return clampWindow(b - half, b + half, n);
}
function adIndices(s, e, n) {
  const out = [];
  for (let j = Math.max(s, 1); j < Math.min(e, n); j++) out.push(j);
  return out;
}
function selectionCoords(node, s, e) {
  const n = node.n;
  const nBins = Math.max(0, e - s);
  const startBp = nBins ? node.start_bp[s] : NaN;
  const endBp = nBins ? node.end_bp[e - 1] : NaN;
  let retained = 0;
  let gaps = 0;
  for (let i = s; i < e; i++) {
    retained += node.end_bp[i] - node.start_bp[i] + 1;
    if (i > s && node.start_bp[i] !== node.end_bp[i - 1] + 1) gaps++;
  }
  const ad = adIndices(s, e, n);
  return {
    s,
    e,
    nBins,
    chromStart: node.start + s,
    chromEnd: node.start + e,
    genomeStart: node.offset_genome + node.start + s,
    genomeEnd: node.offset_genome + node.start + e,
    startBp,
    endBp,
    spanBp: nBins ? endBp - startBp + 1 : 0,
    retainedBp: retained,
    startMb: startBp / MB,
    endMb: endBp / MB,
    spanMb: nBins ? (endBp - startBp + 1) / MB : 0,
    retainedMb: retained / MB,
    gaps,
    adCount: ad.length,
    adFirst: ad.length ? ad[0] : null,
    adLast: ad.length ? ad[ad.length - 1] : null
  };
}
function boundaryMb(node, b) {
  if (b >= node.n) return node.end_bp[node.n - 1] / MB;
  return node.start_bp[b] / MB;
}
function mbToBin(node, mb) {
  const bp = mb * MB;
  let lo = 0, hi = node.n - 1;
  if (bp <= node.start_bp[0]) return 0;
  if (bp >= node.end_bp[hi]) return hi;
  while (lo < hi) {
    const mid = lo + hi >> 1;
    if (node.end_bp[mid] < bp) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}
function mbToBoundary(node, mb) {
  const bp = mb * MB;
  const i = mbToBin(node, mb);
  const mid = (node.start_bp[i] + node.end_bp[i]) / 2;
  return bp < mid ? i : i + 1;
}
function gapMask(startBp) {
  const n = startBp.length;
  const diffs = [];
  for (let i = 1; i < n; i++) diffs.push(startBp[i] - startBp[i - 1]);
  const sorted = [...diffs].sort((a, b) => a - b);
  const med = sorted.length ? sorted[Math.floor(sorted.length / 2)] : 0;
  const out = new Array(n).fill(false);
  for (let i = 1; i < n; i++) if (diffs[i - 1] > 3 * med) out[i] = true;
  return out;
}
function breakGaps(xs, ys, gapAt) {
  const x = [];
  const y = [];
  for (let i = 0; i < xs.length; i++) {
    if (gapAt[i] && i > 0) {
      x.push(xs[i] - 1e-9);
      y.push(null);
    }
    x.push(xs[i]);
    y.push(ys[i]);
  }
  return { x, y };
}

// web/src/stats.ts
function mean(v) {
  return v.length ? v.reduce((a, b) => a + b, 0) / v.length : null;
}
function median(v) {
  if (!v.length) return null;
  const s = [...v].sort((a, b) => a - b);
  const m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}
function stdPop(v) {
  if (!v.length) return null;
  const mu = mean(v);
  return Math.sqrt(v.reduce((a, b) => a + (b - mu) * (b - mu), 0) / v.length);
}
function pearson(a, b) {
  const n = Math.min(a.length, b.length);
  if (n < 3) return { r: null, reason: `needs >= 3 paired observations (have ${n})` };
  const ma = mean(a.slice(0, n)), mb = mean(b.slice(0, n));
  let sab = 0, saa = 0, sbb = 0;
  for (let i = 0; i < n; i++) {
    const da = a[i] - ma, db = b[i] - mb;
    sab += da * db;
    saa += da * da;
    sbb += db * db;
  }
  if (saa === 0 || sbb === 0) return { r: null, reason: "constant paired array" };
  return { r: sab / Math.sqrt(saa * sbb) };
}
function lagCorr(v, k) {
  if (v.length - k < 3) return { r: null, reason: `needs >= 3 pairs at lag ${k} (have ${Math.max(0, v.length - k)})` };
  return pearson(v.slice(0, v.length - k), v.slice(k));
}
function summarize(values) {
  const finite = values.filter((x) => Number.isFinite(x));
  const nonfinite = values.length - finite.length;
  const mu = mean(finite);
  const sd = stdPop(finite);
  let cv = null, cvReason;
  if (mu === null) cvReason = "empty";
  else if (mu === 0) cvReason = "mean is zero";
  else cv = sd / mu;
  const l1 = lagCorr(finite, 1);
  return { n: finite.length, mean: mu, median: median(finite), std: sd, cv, cvReason, lag1: l1.r, lag1Reason: l1.reason, nonfinite };
}
function acf(values, maxLag = 10) {
  const finite = values.filter((x) => Number.isFinite(x));
  const lags = [], vals = [], reasons = [];
  for (let k = 1; k <= maxLag; k++) {
    const r = lagCorr(finite, k);
    lags.push(k);
    vals.push(r.r);
    reasons.push(r.reason ?? null);
  }
  return { lags, values: vals, reasons };
}
function ecdf(values) {
  const s = [...values].filter(Number.isFinite).sort((a, b) => a - b);
  return { x: s, y: s.map((_, i) => (i + 1) / s.length) };
}
function adSamples(node, s, e) {
  const j = adIndices(s, e, node.n);
  return { j, values: j.map((jj) => node.ad[jj - 1]) };
}
function xSamples(node, s, e) {
  const i = [];
  for (let k = Math.max(0, s); k < Math.min(e, node.n); k++) i.push(k);
  return { i, values: i.map((k) => node.x[k]) };
}
function rollingMean(x, window, minPeriods) {
  const n = x.length;
  const half = Math.floor((window - 1) / 2);
  const out = new Array(n);
  const ps = new Float64Array(n + 1), pc = new Int32Array(n + 1);
  for (let i = 0; i < n; i++) {
    const f = Number.isFinite(x[i]);
    ps[i + 1] = ps[i] + (f ? x[i] : 0);
    pc[i + 1] = pc[i] + (f ? 1 : 0);
  }
  for (let i = 0; i < n; i++) {
    const lo = Math.max(0, i - half), hi = Math.min(n - 1, i + (window - 1 - half));
    const c = pc[hi + 1] - pc[lo];
    out[i] = c >= minPeriods ? (ps[hi + 1] - ps[lo]) / c : null;
  }
  return out;
}
function log2Track(roll, baseline, clip) {
  const n = roll.length;
  const values = new Array(n), status = new Array(n);
  if (baseline === null || !(baseline > 0)) {
    for (let i = 0; i < n; i++) {
      values[i] = null;
      status[i] = "undefined";
    }
    return { values, status, clip, undefinedReason: baseline === 0 ? "reference baseline is exactly zero: log2 ratio undefined (no substitution applied)" : "reference baseline unavailable" };
  }
  for (let i = 0; i < n; i++) {
    const r = roll[i];
    if (r === null || !Number.isFinite(r)) {
      values[i] = null;
      status[i] = "nan";
      continue;
    }
    if (r <= 0) {
      values[i] = -clip;
      status[i] = "zero";
      continue;
    }
    const v = Math.log2(r / baseline);
    values[i] = Math.max(-clip, Math.min(clip, v));
    status[i] = "ok";
  }
  return { values, status, clip, undefinedReason: null };
}
function findCompetingPeaks(node, winnerB, nPeaks, minSepMb, eventRange, mbOfBoundary) {
  const d = node.ad, n = d.length;
  if (n === 0) return [];
  const isPeak = new Array(n).fill(false);
  for (let k = 1; k < n - 1; k++) isPeak[k] = d[k] > d[k - 1] && d[k] > d[k + 1];
  isPeak[0] = n > 1 ? d[0] > d[1] : true;
  isPeak[n - 1] = n > 1 ? d[n - 1] > d[n - 2] : true;
  const cand = /* @__PURE__ */ new Set();
  for (let k = 0; k < n; k++) if (isPeak[k]) cand.add(k);
  const w0 = winnerB - 1;
  cand.add(w0);
  const sorted = [...cand].sort((a, b) => d[b] - d[a]);
  const chosen = [{ b: winnerB, ad: d[w0], isWinner: true, isEvent: false, relaxed: false, mb: mbOfBoundary(winnerB) }];
  const chosenMb = [chosen[0].mb];
  if (eventRange) {
    let best = -1;
    for (let k = 0; k < n; k++) {
      const mb = mbOfBoundary(k + 1);
      if (mb >= eventRange.startMb && mb <= eventRange.endMb && (best < 0 || d[k] > d[best])) best = k;
    }
    if (best >= 0 && best !== w0) {
      const mb = mbOfBoundary(best + 1);
      const relaxed = chosenMb.some((m) => Math.abs(mb - m) < minSepMb);
      chosen.push({ b: best + 1, ad: d[best], isWinner: false, isEvent: true, relaxed, mb });
      chosenMb.push(mb);
    } else if (best === w0) {
      chosen[0].isEvent = true;
    }
  }
  for (const k of sorted) {
    if (chosen.length >= nPeaks) break;
    if (chosen.some((p) => p.b === k + 1)) continue;
    const mb = mbOfBoundary(k + 1);
    if (chosenMb.every((m) => Math.abs(mb - m) >= minSepMb)) {
      chosen.push({ b: k + 1, ad: d[k], isWinner: false, isEvent: false, relaxed: false, mb });
      chosenMb.push(mb);
    }
  }
  return chosen;
}
function fmt(v, digits = 3) {
  if (v === null || v === void 0 || Number.isNaN(v)) return "N/A";
  if (!Number.isFinite(v)) return v > 0 ? "+\u221E" : "\u2212\u221E";
  return v.toFixed(digits);
}
function fmtP(p, exact) {
  return `${p.toFixed(6)}${exact ? ` (${exact})` : ""}`;
}
export {
  MB,
  acf,
  adIndices,
  adSamples,
  boundaryMb,
  breakGaps,
  clampWindow,
  ecdf,
  findCompetingPeaks,
  fmt,
  fmtP,
  gapMask,
  lagCorr,
  log2Track,
  mbToBin,
  mbToBoundary,
  mean,
  median,
  peakWindow,
  pearson,
  rollingMean,
  selectionCoords,
  stdPop,
  summarize,
  xSamples
};
