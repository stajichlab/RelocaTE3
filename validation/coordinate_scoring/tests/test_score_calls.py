import csv, sys, tempfile, unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scoring.score_calls import score, main


def _write(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        w.writeheader(); w.writerows(rows)


class TestScore(unittest.TestCase):
    def setUp(self):
        self.d = Path(tempfile.mkdtemp())
        self.truth = self.d / "truth.tsv"
        self.calls = self.d / "calls.tsv"
        _write(self.truth,
               ["event_id","chrom","position","te_family","biological_class","tsd"],
               [{"event_id":"E1","chrom":"Chr1","position":1000,"te_family":"mPing","biological_class":"homozygous","tsd":"TTA"},
                {"event_id":"E2","chrom":"Chr1","position":5000,"te_family":"mPing","biological_class":"somatic_insertion","tsd":"TAA"}])
        _write(self.calls,
               ["chrom","position","te_family","tsd","strand","status","caller","sample"],
               [{"chrom":"Chr1","position":1003,"te_family":"mPing","tsd":"TTA","strand":"+","status":"homozygous","caller":"relocate3","sample":"S1"},
                {"chrom":"Chr1","position":9999,"te_family":"mPing","tsd":"XX","strand":"+","status":"homozygous","caller":"relocate3","sample":"S1"}])

    def test_score(self):
        summary, matches, fps, precision_row = score(self.truth, self.calls, sample="S1", caller="relocate3", window=10)
        by_class = {r["biological_class"]: r for r in summary}
        self.assertEqual(by_class["homozygous"]["detected_events"], 1)
        self.assertEqual(by_class["homozygous"]["status_correct_events"], 1)
        self.assertEqual(by_class["homozygous"]["tsd_exact_events"], 1)
        self.assertEqual(by_class["somatic_insertion"]["detected_events"], 0)
        self.assertEqual(len(fps), 1)
        # class summary column renamed precision -> class_call_share (same formula)
        self.assertIn("class_call_share", by_class["homozygous"])
        self.assertNotIn("precision", by_class["homozygous"])
        self.assertEqual(by_class["homozygous"]["class_call_share"], 1 / 2)
        # per-sample precision row: 1 matched of 2 calls -> 0.5
        self.assertEqual(precision_row["caller"], "relocate3")
        self.assertEqual(precision_row["sample"], "S1")
        self.assertEqual(precision_row["total_calls"], 2)
        self.assertEqual(precision_row["matched_calls"], 1)
        self.assertEqual(precision_row["false_positive_calls"], 1)
        self.assertEqual(precision_row["overall_precision"], 0.5)
        self.assertEqual(precision_row["false_discovery_rate"], 0.5)

    def test_somatic_cellular_fraction_split(self):
        # Two somatic truth events at cellular_fraction 0.1 and 0.4 must produce
        # SEPARATE summary rows keyed by cellular_fraction with independent recall.
        # Only the 0.4 event is detected here.
        truth = self.d / "truth_cf.tsv"
        calls = self.d / "calls_cf.tsv"
        _write(truth,
               ["event_id","chrom","position","te_family","biological_class","tsd","cellular_fraction","expected_vaf"],
               [{"event_id":"C1","chrom":"Chr1","position":2000,"te_family":"mPing","biological_class":"somatic_insertion","tsd":"TTA","cellular_fraction":"0.1","expected_vaf":"0.05"},
                {"event_id":"C2","chrom":"Chr1","position":8000,"te_family":"mPing","biological_class":"somatic_insertion","tsd":"TAA","cellular_fraction":"0.4","expected_vaf":"0.2"}])
        _write(calls,
               ["chrom","position","te_family","tsd","strand","status","caller","sample"],
               [{"chrom":"Chr1","position":8002,"te_family":"mPing","tsd":"TAA","strand":"+","status":"somatic","caller":"relocate3","sample":"S1"}])
        summary, matches, fps, precision_row = score(truth, calls, sample="S1", caller="relocate3", window=10)
        by_cf = {r["cellular_fraction"]: r for r in summary}
        self.assertEqual(set(by_cf), {"0.1", "0.4"})
        self.assertEqual(by_cf["0.1"]["biological_class"], "somatic_insertion")
        self.assertEqual(by_cf["0.1"]["truth_events"], 1)
        self.assertEqual(by_cf["0.1"]["detected_events"], 0)
        self.assertEqual(by_cf["0.1"]["detection_recall"], 0.0)
        self.assertEqual(by_cf["0.1"]["expected_vaf"], "0.05")
        self.assertEqual(by_cf["0.4"]["detected_events"], 1)
        self.assertEqual(by_cf["0.4"]["detection_recall"], 1.0)
        self.assertEqual(by_cf["0.4"]["expected_vaf"], "0.2")

    def test_main_writes_precision(self):
        # 2 of 3 calls match a truth event -> overall_precision 0.666..., FDR 0.333...
        truth = self.d / "truth_p.tsv"
        calls = self.d / "calls_p.tsv"
        outdir = self.d / "report_p"
        _write(truth,
               ["event_id","chrom","position","te_family","biological_class","tsd"],
               [{"event_id":"P1","chrom":"Chr1","position":1000,"te_family":"mPing","biological_class":"homozygous","tsd":"TTA"},
                {"event_id":"P2","chrom":"Chr1","position":2000,"te_family":"mPing","biological_class":"homozygous","tsd":"TAA"}])
        _write(calls,
               ["chrom","position","te_family","tsd","strand","status","caller","sample"],
               [{"chrom":"Chr1","position":1001,"te_family":"mPing","tsd":"TTA","strand":"+","status":"homozygous","caller":"relocate3","sample":"S1"},
                {"chrom":"Chr1","position":2001,"te_family":"mPing","tsd":"TAA","strand":"+","status":"homozygous","caller":"relocate3","sample":"S1"},
                {"chrom":"Chr1","position":9999,"te_family":"mPing","tsd":"XX","strand":"+","status":"homozygous","caller":"relocate3","sample":"S1"}])
        argv = ["--truth", str(truth), "--calls", str(calls),
                "--sample", "S1", "--caller", "relocate3",
                "--window", "10", "--outdir", str(outdir)]
        self.assertEqual(main(argv), 0)
        self.assertTrue((outdir / "precision.tsv").exists())
        with open(outdir / "precision.tsv") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["total_calls"], "3")
        self.assertEqual(rows[0]["matched_calls"], "2")
        self.assertAlmostEqual(float(rows[0]["overall_precision"]), 2 / 3, places=6)
        self.assertAlmostEqual(float(rows[0]["false_discovery_rate"]), 1 / 3, places=6)

    def test_somatic_status_mapping(self):
        # A matched somatic event: truth biological_class="somatic_insertion",
        # call status="somatic". _norm maps somatic_insertion -> somatic, so this
        # must count as detected AND status-correct. Fails if that mapping is removed.
        truth = self.d / "truth_som.tsv"
        calls = self.d / "calls_som.tsv"
        _write(truth,
               ["event_id","chrom","position","te_family","biological_class","tsd"],
               [{"event_id":"S1","chrom":"Chr1","position":3000,"te_family":"mPing","biological_class":"somatic_insertion","tsd":"TTA"}])
        _write(calls,
               ["chrom","position","te_family","tsd","strand","status","caller","sample"],
               [{"chrom":"Chr1","position":3004,"te_family":"mPing","tsd":"TTA","strand":"+","status":"somatic","caller":"relocate3","sample":"S1"}])
        summary, matches, fps, precision_row = score(truth, calls, sample="S1", caller="relocate3", window=10)
        by_class = {r["biological_class"]: r for r in summary}
        self.assertEqual(by_class["somatic_insertion"]["detected_events"], 1)
        self.assertEqual(by_class["somatic_insertion"]["status_correct_events"], 1)
        self.assertEqual(len(fps), 0)

    def test_te_family_repeatmasker_suffix_matches(self):
        # Real callers emit RepeatMasker-style TE names (e.g. RelocaTE3 emits
        # "mPing#DNA/Harbinger") while truth is the bare family "mPing". _norm
        # strips the "#class/family" suffix so these still match. Fails (0 recall)
        # if the suffix stripping is removed.
        truth = self.d / "truth_rm.tsv"
        calls = self.d / "calls_rm.tsv"
        _write(truth,
               ["event_id","chrom","position","te_family","biological_class","tsd"],
               [{"event_id":"R1","chrom":"Chr1","position":247385,"te_family":"mPing","biological_class":"homozygous","tsd":"AAG"}])
        _write(calls,
               ["chrom","position","te_family","tsd","strand","status","caller","sample"],
               [{"chrom":"Chr1","position":247381,"te_family":"mPing#DNA/Harbinger","tsd":"AAG","strand":"+","status":"homozygous","caller":"relocate3","sample":"S1"}])
        summary, matches, fps, precision_row = score(truth, calls, sample="S1", caller="relocate3", window=10)
        by_class = {r["biological_class"]: r for r in summary}
        self.assertEqual(by_class["homozygous"]["detected_events"], 1)
        self.assertEqual(len(fps), 0)

    def test_multi_te_truth_is_stratified_by_curated_group(self):
        truth = self.d / "truth_multi.tsv"
        calls = self.d / "calls_multi.tsv"
        _write(
            truth,
            [
                "event_id", "chrom", "position", "te_family", "te_group",
                "te_class", "te_order", "te_superfamily",
                "biological_class", "tsd",
            ],
            [
                {
                    "event_id": "T1", "chrom": "Chr1", "position": 1000,
                    "te_family": "TE_A", "te_group": "LINE",
                    "te_class": "Class_I", "te_order": "LINE",
                    "te_superfamily": "LINE", "biological_class": "homozygous",
                    "tsd": "AAAAAAA",
                },
                {
                    "event_id": "T2", "chrom": "Chr1", "position": 2000,
                    "te_family": "TE_B", "te_group": "SINE",
                    "te_class": "Class_I", "te_order": "SINE",
                    "te_superfamily": "SINE", "biological_class": "homozygous",
                    "tsd": "TTTTTTT",
                },
            ],
        )
        _write(
            calls,
            ["chrom", "position", "te_family", "tsd", "strand", "status", "caller", "sample"],
            [{
                "chrom": "Chr1", "position": 1001, "te_family": "TE_A#LINE/L1",
                "tsd": "AAAAAAA", "strand": "+", "status": "homozygous",
                "caller": "relocate3", "sample": "S1",
            }],
        )
        summary, _, _, _ = score(
            truth, calls, sample="S1", caller="relocate3", window=10
        )
        by_group = {row["te_group"]: row for row in summary}
        self.assertEqual(set(by_group), {"LINE", "SINE"})
        self.assertEqual(by_group["LINE"]["detected_events"], 1)
        self.assertEqual(by_group["SINE"]["detected_events"], 0)

    def test_main_writes_reports(self):
        # Truth/calls where every call matches a truth event -> zero false positives,
        # exercising the empty-file branch of _write for false_positive_calls.tsv.
        truth = self.d / "truth_main.tsv"
        calls = self.d / "calls_main.tsv"
        outdir = self.d / "report"
        _write(truth,
               ["event_id","chrom","position","te_family","biological_class","tsd"],
               [{"event_id":"M1","chrom":"Chr1","position":1000,"te_family":"mPing","biological_class":"homozygous","tsd":"TTA"}])
        _write(calls,
               ["chrom","position","te_family","tsd","strand","status","caller","sample"],
               [{"chrom":"Chr1","position":1002,"te_family":"mPing","tsd":"TTA","strand":"+","status":"homozygous","caller":"relocate3","sample":"S1"}])
        argv = ["--truth", str(truth), "--calls", str(calls),
                "--sample", "S1", "--caller", "relocate3",
                "--window", "10", "--outdir", str(outdir)]
        rc = main(argv)
        self.assertEqual(rc, 0)
        for name in ("matches.tsv", "false_positive_calls.tsv", "correctness.tsv", ".complete"):
            self.assertTrue((outdir / name).exists(), f"missing {name}")
        # zero false positives -> empty file
        self.assertEqual((outdir / "false_positive_calls.tsv").read_text(), "")
        # A completed report dir (.complete present) is idempotently re-scorable:
        # a second run overwrites cleanly and returns 0 instead of refusing.
        self.assertEqual(main(argv), 0)
        self.assertTrue((outdir / ".complete").exists())
        # But a non-empty dir WITHOUT the .complete marker is foreign/partial and
        # must still be refused, so we never clobber unrelated data.
        (outdir / ".complete").unlink()
        with self.assertRaises(FileExistsError):
            main(argv)


if __name__ == "__main__":
    unittest.main()
