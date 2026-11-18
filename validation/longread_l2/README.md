# Long-read L2 gate

Checks that `relocaTE3 lr-align` -> `lr-signatures` -> `lr-cluster` finds the
true insertion sites of the small simulated Chr3 long-read panel.

**Gate:** at 20x, on both `ont-hq` and `hifi`, at least 95% of the 200 truth
sites have a cluster whose start is within 50 bp. 5x recall and false clusters
are reported, not gated.

## Run

From the repo root, on a checkout with the pixi `default` env installed:

```bash
# Link the panel once (built in simulate_data / make_simulation_new; gitignored)
ln -s /rhome/nmath020/wessler_bigdata/rice/make_simulated_genome/make_simulation_new/results/chr3_lr_small \
  tests/data/longread_small

mkdir -p logs && sbatch validation/longread_l2/run_gate.slurm
```

Overridable environment variables: `PANEL` (default
`tests/data/longread_small`), `GENOME` (`tests/data/sim_genome/MSU7.Chr3_2M.fa`),
`OUT` (`validation/longread_l2/output`), `MAX_DIST` (50), `MIN_RECALL` (0.95).

The verdict is the `GATE: PASS` / `GATE: FAIL` line in
`logs/longread_l2_gate.<jobid>.log`; the job exits 0 either way.

Rerun-safe, per stage: `lr-align` is skipped when `<set>.lr.bam` and its
`.bai` exist, `lr-signatures` when `<set>.signatures.tsv.gz` exists,
`lr-cluster` when `<set>.clusters.tsv.gz` exists (each stage writes its output
atomically). A BAM without its `.bai` stops the job: remove that set directory
and resubmit.

Expected run time: about 3 minutes on 8 CPUs for all eight read sets (job 32185663, 2026-11-18); peak memory about 2 GB (minimap2 on ONT-HQ 20x).

## Outputs

`$OUT/<platform>_<sample>/` and `$OUT/<platform>_<sample>_control/`, one per
read set: `<set>.lr.bam`, `<set>.signatures.tsv.gz`, `<set>.clusters.tsv.gz`,
`timing.tsv`. The control reads come from the reference with no inserts.

The log's timing summary shows the last `status == ok` row per stage of each
set's `timing.tsv` (reruns append `FileExistsError` rows, which are dropped).
`records_in` is 0 for `lr-align` with FASTQ input (only `records_out`, the
mapped records, is counted).

`$OUT/gate.tsv`, one row per sample (`score_sites.py`):

| Column | Meaning |
|---|---|
| `label` | `<platform>_<sample>` |
| `n_truth` | truth sites in `truth_events.tsv` |
| `found` | truth sites with a cluster on the same chrom whose start is within `MAX_DIST` bp |
| `recall` | `found / n_truth` (gated at 20x) |
| `found_min2`, `recall_min2` | the same, using only clusters with `n_reads >= 2` |
| `n_clusters` | clusters in the sample |
| `unmatched` | sample clusters not within `MAX_DIST` bp of any truth site |
| `unmatched_min2` | the same, among clusters with `n_reads >= 2` |
| `control_clusters` | clusters in the control: all false positives |
| `control_clusters_min2` | control clusters with `n_reads >= 2` |
| `chance_recall` | recall reached by the control clusters alone: how often noise "finds" a site |
| `chance_recall_min2` | the same, control clusters with `n_reads >= 2` |

Score one sample by hand:

```bash
pixi run --frozen python validation/longread_l2/score_sites.py \
  --truth tests/data/longread_small/truth_events.tsv \
  --clusters validation/longread_l2/output/hifi_cov20x_rep1/hifi_cov20x_rep1.clusters.tsv.gz \
  --control validation/longread_l2/output/hifi_cov20x_rep1_control/hifi_cov20x_rep1_control.clusters.tsv.gz \
  --label hifi_cov20x_rep1 --header
```
