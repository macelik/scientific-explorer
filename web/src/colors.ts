export const AD0 = '#2a6f97'
export const BP0 = '#440154'
export const LEFT = '#2a9d8f'
export const RIGHT = '#e76f51'
export const SEGSTAT = '#c1121f'
export const OVERLAP = '#7209b7'
export const LOCALBASE = '#495057'
export const EVENT = '#c1121f'
export const PEAKS = ['#264653', '#e9c46a', '#8ab17d', '#6a4c93', '#f4a261', '#219ebc', '#bc4749', '#8338ec']
export const CN_COLORS: Record<number, string> = { 0: '#1d4e89', 0.5: '#7fa7d6', 1: '#f3f4f6', 1.5: '#f2a8a0', 2: '#c1121f' }
export const CN_LABELS: Record<number, string> = { 0: 'loss', 0.5: 'putative loss', 1: 'base / diploid', 1.5: 'putative gain', 2: 'gain' }
export const TERCILE_COLORS: Record<string, string> = { low: '#cbd5e1', mid: '#64748b', high: '#0f172a' }
export const NODE_COLOR: Record<string, string> = { root: AD0, left: LEFT, right: RIGHT }
export function windowColor(i: number): string { return PEAKS[i % PEAKS.length] }
