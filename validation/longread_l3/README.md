# Long-read L3 gate

Checks that `relocaTE3 lr-align` -> `lr-signatures` -> `lr-cluster` ->
`lr-annotate` -> `lr-call` calls the true insertions of the small simulated
Chr3 long-read panel, with the right TE family.

**Gate:** at 20x, on both `ont-hq` and `hifi`, recall >= 0.95 and precision
>= 0.95 over the 200 truth sites. A match is a PASS spanning call within 50 bp
of a truth site that names the truth family; matching is greedy one-to-one
(`score_calls.py`). Strand accuracy, the 5x rows, clip-only calls and the
control calls are reported, not gated.

## Run

From the repo root, on a checkout with the pixi `default` env installed:

```bash
# Link the panel once at tests/data/longread_small (gitignored); the ln -s
# command is in docs/long-read-test-panel.md ("Use from RelocaTE3").
mkdir -p logs && sbatch validation/longread_l3/run_gate.slurm
```

Overridable environment variables:

| Variable | Default |
|---|---|
| `PANEL` | `tests/data/longread_small` |
| `GENOME` | `tests/data/sim_genome/MSU7.Chr3_2M.fa` |
| `TE_LIBRARY` | `$PANEL/config/te_library.fa` |
| `OUT` | `validation/longread_l3/output` |
| `MAX_DIST` | 50 (bp) |
| `MIN_RECALL` | 0.95 |
| `MIN_PRECISION` | 0.95 |

## Read the log

`logs/longread_l3_gate.<jobid>.log`, in order:

1. Run parameters: commit (and dirty file count), panel, genome, TE library,
   thresholds.
2. One line per stage per read set: run, or skipped because its output exists.
3. `gate.tsv`, then each sample's per-family table.
4. The timing summary: the last `status == ok` row per stage of each set's
   `timing.tsv` (reruns append `FileExistsError` rows, which are dropped).
   `records_in` is 0 for `lr-align` with FASTQ input.
5. One check line per gated row and metric, e.g.
   `hifi_cov20x_rep1: precision 0.9850 >= 0.95`, then the verdict line
   `GATE: PASS` or `GATE: FAIL`. The job exits 0 either way.

## Rerun behaviour

Each stage is skipped when its output exists (each stage writes its output
atomically): `<set>.lr.bam` + `.bai` (`lr-align`), `<set>.signatures.tsv.gz`
(`lr-signatures`), `<set>.clusters.tsv.gz` (`lr-cluster`),
`<set>.te_hits.tsv.gz` (`lr-annotate`), `<set>.calls.tsv.gz` (`lr-call`).

The job stops with an error, asking you to remove that set directory and
resubmit, when:

- a BAM exists without its `.bai` (an interrupted `lr-align`);
- a stage's existing output is older than an input it was made from, i.e. an
  earlier stage was rebuilt after it. Signatures are checked against the BAM,
  clusters and TE hits against the signatures, calls against the clusters,
  signatures and TE hits.

Scoring always reruns; `gate.tsv` and the per-family tables are written to a
`.tmp` file and moved into place.

## Outputs

`$OUT/<platform>_<sample>/` and `$OUT/<platform>_<sample>_control/`, one per
read set: `<set>.lr.bam` (+ `.bai`), `<set>.lr.meta.tsv`,
`<set>.signatures.tsv.gz`, `<set>.clusters.tsv.gz`, `<set>.te_hits.tsv.gz`,
`<set>.calls.tsv.gz`, `timing.tsv`. The control reads come from the reference
with no inserts.

`$OUT/gate.tsv`, one row per sample (`score_calls.py`); scored calls are the
`spanning` calls that PASS:

| Column | Meaning |
|---|---|
| `label` | `<platform>_<sample>` |
| `n_truth` | truth sites in `truth_events.tsv` |
| `n_pass` | PASS spanning calls |
| `tp` | matched (truth, call) pairs: same chrom, within `MAX_DIST` bp, same family |
| `recall` | `tp / n_truth` (gated at 20x) |
| `precision` | `tp / n_pass` (gated at 20x) |
| `f1` | `2 * tp / (n_truth + n_pass)` |
| `tp_position_only` | the same matching ignoring family |
| `strand_correct`, `strand_accuracy` | matched pairs with the truth strand, and that over `tp` |
| `n_filtered` | spanning calls that do not PASS |
| `clip_only_pass` | PASS clip-only calls |
| `clip_only_hits` | truth sites left unmatched that a PASS clip-only call matches |
| `control_pass` | PASS spanning calls on the control: all false positives |
| `control_clip_only_pass` | PASS clip-only calls on the control |

`$OUT/<platform>_<sample>.per_family.tsv`: `family n_truth tp recall n_pass
precision`, one row per family in the truth or in the scored calls.

Score one sample by hand:

```bash
pixi run --frozen python validation/longread_l3/score_calls.py \
  --truth tests/data/longread_small/truth_events.tsv \
  --calls validation/longread_l3/output/hifi_cov20x_rep1/hifi_cov20x_rep1.calls.tsv.gz \
  --control validation/longread_l3/output/hifi_cov20x_rep1_control/hifi_cov20x_rep1_control.calls.tsv.gz \
  --label hifi_cov20x_rep1 --header
```
