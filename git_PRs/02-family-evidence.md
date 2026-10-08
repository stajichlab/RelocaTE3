# Preserve tied TE-family evidence and resolve supported assignments

Alphabetical alignment tie-breaking could turn ambiguous reads into confident
votes for the wrong TE family. Preserve equally best family alternatives while
keeping the selected alignment geometry, and resolve the insertion's primary
family using distinguishing junction reads when support is sufficient.

Resolution requires at least two distinct distinguishing reads, a strict
majority among distinguishing reads, membership in a strict majority of all
informative best-family sets, and compatible trim geometry. Unresolved ties
retain their previous selected-vote primary, including the pre-deduplication
fallback; repeated read names cannot satisfy the distinguishing-read quorum.

Outputs retain one primary family and append ambiguity, candidate-support, and
resolution metadata. Three-column legacy read mappings remain supported;
malformed new evidence fails explicitly. Metadata survives insertion and
characterization TXT/GFF round trips. Confidence is an evidence fraction, not
a calibrated probability.

In the 72-sample stored-alignment replay, this recovered **three truth
detections and removed three false positives without losing a baseline
detection**. The gains are three benchmark conditions of the same Os3328 LINE
locus, not three independent biological loci. Original evidence reproduced
historical calls in all 72 samples; geometry and genotype were unchanged.

Portable real-read fixtures and tests cover tie ordering, compatibility,
deduplication, sparse/conflicting evidence, and output contracts. The assembled
checkpoint passed **331 tests with zero failures or skips** (job 29419250).

Depends on the memory-streaming PR. Full trimming resource costs of retaining
ties still require measurement on the final committed revision.
