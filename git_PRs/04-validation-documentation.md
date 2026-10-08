# Add reproducible parity audits and checkpoint validation

Add tools to audit truth coordinates, paired detections, residual false
positives, read evidence, and frozen downstream replays. Record the finalized
checkpoint's accuracy/resource evidence and remaining release requirements.

The local scoring copy makes the TSD-interval policy explicit: compare insertion
positions against the truth TSD interval plus the stated tolerance, keeping
genotype and TSD accuracy separate. Benchmark integration is supplied as a
patch; it has not been applied to the sibling benchmark repository or dashboard.

SLURM runners use submission-directory paths, explicit resources, frozen
sources/input metadata, protected outputs, and completion markers. The full
checkpoint runner records executable paths and source/test-control hashes and
fails if those inputs change while the suite runs. Large generated benchmark
outputs and frozen run directories are excluded from this PR.

Validation: **331 tests passed, zero failures/errors/skips**, job 29419250.
The final 72-sample local replay removed 28 FPs without sacrificing a baseline
detection. Historical rejected experiments are retained as diagnostic context,
not instructions to enable their behavior.

Document release-first priorities and future Rust, Nextflow, long-read sequence
reconstruction, and multi-reference/pangenome evidence requirements. These
future capabilities are planning notes, not implemented features.

Depends on the local-filter PR. This completes a development checkpoint.
A fresh end-to-end benchmark of the merged revision, full resource measurements,
remaining precision/detection differences, and broader R2 utility coverage
remain prerequisites for a release promise.
