export interface ChromMeta { name: string; offset: number; n: number; start_bp: number; end_bp: number }
export interface CellMeta {
  cellID: string; row: number; raw_coverage: number; tercile: 'low' | 'mid' | 'high';
  n_segments: number; best_s: number; frac_loss: number; frac_gain: number; frac_base: number
}
export interface RegionSeed { id: string; label: string; chrom: string; start_mb: number; end_mb: number; note: string }
export interface Meta {
  dataset: { identity: { hash: string; files: any[]; n_cells: number; n_bins: number; n_cells_total_h5ad: number };
    h5ad: string; cohort_dir: string; package_dir: string; bin_width_bp: number | null; signal: string }
  chromosomes: ChromMeta[]
  cells: CellMeta[]
  cluster_order: number[]
  cn_levels: { value: number; label: string }[]
  regions: RegionSeed[]
  reference_cells: { cell: string; chrom: string; label: string }[]
  defaults: {
    n_permutations: number; seed: number; alpha: number; rolling_window: number; rolling_min_periods: number;
    half_window_bins: number; n_competing_peaks: number; min_peak_separation_mb: number; log2_clip: number; production_k: number
  }
}
export interface ProdBreakpoint { bp: number; ad: number; p_local: number; p_global: number }
export interface CellData {
  cellID: string; row: number; raw_coverage: number; tercile: string; best_s_export: number; n_segments: number
  genome: { baseline_trimmed: number; baseline_median: number; zero_fraction: number; x_nonfinite: number }
  x: number[]
  production_breakpoints: Record<string, ProdBreakpoint[]>
  production_clusters: number[]
  production_cn5: number[]
  load_ms: number
}
export interface NodeData {
  cell: string; chrom: string; start: number; end: number; n: number; offset_genome: number
  ad: number[]; argmax_b: number | null; argmax_ad: number | null
  baseline_trimmed: number | null; baseline_median: number | null; x_nonfinite: number
  x: number[]; start_bp: number[]; end_bp: number[]
  argmax_chrom: number | null; argmax_in_production: boolean
  production_breakpoints_in_node: ProdBreakpoint[]
  compute_ms: number; cached: boolean
}
export interface TestResult {
  observed_ad: number; p_local: number; p_global: number; p_local_exact: string; p_global_exact: string
  n_permutations: number; seed: number; n_left: number; n_right: number
}
export interface TestState { status: 'pending' | 'done' | 'error'; jobId?: string; result?: TestResult; error?: string; elapsed?: number }
export interface SegmentRow {
  segment_id: number; start_bin: number; end_bin: number; n_bins: number; start_bp: number; end_bp: number
  mean_x: number; median_x: number; mean_norm: number | null; cn5: number; cn_int: number; cn_cont: number
  summary_x?: number; summary_norm?: number | null; retained_bins?: number
}
export type SegmentEstimator = 'arithmetic' | 'iqr_upper' | 'iqr_two_sided'
export interface AssignResult {
  cell: string; chrom: string; boundaries: number[]
  production: { best_s: number; best_s_export: number; matches_export: boolean; cn5: number[]; cn_int: number[]; cn_cont: number[]; segments: SegmentRow[] }
  manual: { best_s: number; cn5: number[]; cn_int: number[]; cn_cont: number[]; segments: SegmentRow[]; n_segments_genome: number }
  changed_bins_this_chrom: number; changed_elsewhere: { chrom: string; changed_bins: number; n: number }[]
  genome_baseline_trimmed: number; compute_ms: number
  segment_estimator?: SegmentEstimator; estimator_scope?: string
  arithmetic_manual?: { best_s: number; cn5: number[]; cn_int: number[]; cn_cont: number[]; segments: SegmentRow[] }
  estimator_changed_bins?: number; estimator_changed_elsewhere?: { chrom: string; changed_bins: number }[]
}
export interface EventCandidate {
  cell: string; chrom: string; node: 'root' | 'left' | 'right'; bp: number; ad: number; accepted: boolean; parent_accepted: boolean
  n_left: number; n_right: number
}

/** A node = one contiguous chromosome-local half-open interval [start, end) of one cell. */
export interface NodeRef {
  id: string; cell: string; chrom: string; start: number; end: number; depth: number; version: number
  parentId: string | null; side: 'root' | 'left' | 'right'
}

export type BaselineMethod = 'trimmed' | 'median'

/** Selection window: node-local half-open [s, e), optional anchor boundary b (node-local). */
export interface Selection {
  id: string; name: string; color: string; visible: boolean
  cell: string; chrom: string; nodeId: string; nodeVersion: number; nodeStart: number; nodeEnd: number; nodeSide: string
  s: number; e: number; anchor: number | null
  source: 'drawn' | 'peak' | 'winner' | 'event' | 'reapplied' | 'saved'
  note?: string
  clipped?: { left: number; right: number }
  createdAt: number
}

export interface EventRule {
  id: string; name: string; enabled: boolean
  family: 'cn' | 'bp'
  chrom: string; startMb: number; endMb: number
  states: number[]; minFraction: number
  node: 'root' | 'left' | 'right' | 'any'; acceptance: 'candidate' | 'accepted' | 'any'
}

export interface ManualNode extends NodeRef { splitB: number | null; leftId: string | null; rightId: string | null; label: string }
export interface ManualTree { rootId: string; nodes: Record<string, ManualNode>; nextId: number; version: number }
export interface ManualSession {
  cell: string; chrom: string; tree: ManualTree; history: ManualTree[]; future: ManualTree[]
  activeNodeId: string; notice: string | null
  cnResult: AssignResult | null; cnOutdated: boolean; cnBoundaries: number[] | null
  cnEstimator?: SegmentEstimator; cnBusy?: boolean; cnToken?: string
}

export interface SnapshotWindow {
  selection: Selection; adOffsets: number[]; adValues: number[]; xOffsets: number[]; xValues: number[]
  adStats: SummaryStats; xStats: SummaryStats; coords: SelectionCoords; adAcf: AcfResult; xAcf: AcfResult
}
export interface Snapshot { id: string; name: string; createdAt: number; cell: string; chrom: string; nodeLabel: string; windows: SnapshotWindow[]; provenance: string }

export interface SummaryStats {
  n: number; mean: number | null; median: number | null; std: number | null; cv: number | null; cvReason?: string
  lag1: number | null; lag1Reason?: string; nonfinite: number
}
export interface AcfResult { lags: number[]; values: (number | null)[]; reasons: (string | null)[] }
export interface SelectionCoords {
  s: number; e: number; nBins: number; chromStart: number; chromEnd: number; genomeStart: number; genomeEnd: number
  startBp: number; endBp: number; spanBp: number; retainedBp: number; startMb: number; endMb: number; spanMb: number; retainedMb: number
  gaps: number; adCount: number; adFirst: number | null; adLast: number | null
}

// ---------------------------------------------------------------- integration story
export interface IntegrationCell {
  cell: string; gc_norm_depth: number; snp_coverage: number; total_raw_depth: number; n_bins_snp_covered: number
  wnn_depth_weight: number; wnn_haplo_weight: number; in_explorer_cohort: boolean; percell_best_s: number | null
  umap_wnn_x: number; umap_wnn_y: number; umap_depth_x: number; umap_depth_y: number; umap_haplo_x: number; umap_haplo_y: number
  depth_pc1: number; haplo_pc1: number; depth_pc2: number; haplo_pc2: number
  chr8_score?: number | null; zero_frac?: number | null; wnn_spectral1?: number | null
  [clusterKey: string]: any
}
export interface ClusterBreakpoint { source: string; chromosome: string; absolute_bin: number; bp_genomic_start: number; bp_genomic_end: number; recursion_depth: number; ad_score: number; p_local: number; p_global: number; n_cells: number }
export interface ClusterSegment { group: string; n_cells: number; best_s: number; segment: number; chromosome: string; bin_start: number; bin_end: number; n_bins: number; continuous_cn: number; integer_cn: number; watson: number; cn_call: number }
export interface FlankRow {
  source: string; chromosome: string; orig_seg_id: number; bin_start: number; bin_end: number; seg_length: number; Mb_start: number; Mb_end: number
  median_seg: number; mean_seg_perbin: number; trimmean_seg_perbin: number; var_seg_perbin: number; std_seg_perbin: number; cv_seg_perbin: number; n_bins_seg: number
  median_left: number | null; delta_left: number | null; srd_phi_left: number | null; frac_under_left: number | null; eligible_left: boolean; merged_left: boolean; reason_left: string
  mean_left_flank_perbin: number | null; trimmean_left_flank_perbin: number | null; std_left_flank_perbin: number | null; cv_left_flank_perbin: number | null; n_bins_left_flank: number | null; delta_mean_left: number | null; delta_trim_left: number | null
  median_right: number | null; delta_right: number | null; srd_phi_right: number | null; frac_under_right: number | null; eligible_right: boolean; merged_right: boolean; reason_right: string
  mean_right_flank_perbin: number | null; trimmean_right_flank_perbin: number | null; std_right_flank_perbin: number | null; cv_right_flank_perbin: number | null; n_bins_right_flank: number | null; delta_mean_right: number | null; delta_trim_right: number | null
}
export interface ScreeningRow { source: string; chromosome: string; n_cells: number; seg_id: number; bin_start: number; bin_end: number; seg_length: number; mean_rate: number; median_level: number; log2fc_left: number | null; log2fc_right: number | null; srd_left: number | null; srd_right: number | null; srd_phi_left: number | null; srd_phi_right: number | null; pooled_phi: number; proposal: string; proposed_boundary_bin: number | null }
export interface IntegrationMeta {
  identity: { hash: string; files: any[] }; h5mu: string; notes: string[]; problems: string[]
  sources: string[]; source_sizes: Record<string, number>
  cluster_keys: Record<string, string[]>; default_cluster_key: string
  chromosomes: string[]; pilot_chromosomes: string[]; manifest: any; neighbors_params: any
  chr8_score: { chrom: string; num_mb: [number, number]; den_mb: [number, number] }
  cells: IntegrationCell[]
  breakpoints: ClusterBreakpoint[] | null; cluster_segments: ClusterSegment[] | null
  flank_scores: FlankRow[] | null; screening: ScreeningRow[] | null; bootstrap: any[] | null; srd_descriptive: any[] | null
  retained_boundaries: any[] | null; simplified_segments: any[] | null; initial_boundary_scores: any[] | null; unresolved_segments: any[] | null
  threshold_sensitivity: any[] | null; review_summary: any[] | null
  index_done: boolean
}
export interface IntegrationStatus { available: boolean; h5mu?: string; problems: string[]; notes?: string[]; index: { status: string; step?: string; done?: number; total?: number; error?: string | null } }
export interface CellDetail { cell: string; index: number; neighbors: Record<string, { cell: string; index: number; w: number; cluster: string }[]>; overlap: Record<string, number> }
export interface CandidateRow { source: string; chromosome: string; recursion_depth: number; node_start_bin: number; node_end_bin: number; candidate_rank: number; candidate_local_cut: number; candidate_absolute_bin: number; ad_score: number; left_child_bins: number; right_child_bins: number; geometry_valid: boolean; rejection_reasons: string | null; selected_for_significance: boolean; p_local: number | null; p_global: number | null; accepted: boolean }
