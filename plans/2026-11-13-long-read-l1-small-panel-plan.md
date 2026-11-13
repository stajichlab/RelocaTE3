# Long-read L1a: Small Chr3 Long-Read Panel — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a small, reproducible ONT-HQ + HiFi read panel over the
repository's Chr3 2 Mb fixture (200 known insertions), so long-read code can
be tested in minutes.

**Architecture:** Three small scripts under `validation/longread_small/`.
(1) Turn the fixture's truth GFF into a config for the `simulate_data`
toolkit's `te-panel` command, which builds the sample genome. (2) Simulate
reads from that genome with the toolkit's long-read functions, the same ones
that built the full `riceTElib_lr_benchmark` panel. (3) Check the finished
panel. The panel is written outside git and has the same manifest layout as
the full panel.

**Tech stack:** Python ≥ 3.10 standard library; `simulate_data` toolkit
(PBSIM3, `ccs`) through its own pixi environment; minimap2 + samtools from
RelocaTE3's pixi environment; SLURM.

Design: `plans/2026-11-09-long-read-algorithm-proposal.md`, "Test tiers",
tier 1. This plan covers milestone L1 part (a) only. Parts (b) benchmark
long-read manifests and (c) HiFi 30× repair live in other repositories and
get their own plans; see "Not in this plan".

---

## Facts checked on 2026-11-13 (do not re-derive)

- **Truth GFF convention.** `tests/data/sim_genome/MSU7.Chr3_2M.ALL.gff` has
  200 rows; column 3 is the TE family, column 4 the position `s`, column 7
  the strand, attribute `TSD=` the target-site duplication. For all 200 rows
  the TSD equals the reference bases `seq[s-L:s]` (0-based slice, `L` = TSD
  length). In `te-panel` terms the insertion position is `P = s`.
- **TSD length is constant per family:** mPing, mGing, Gaijin, spmlike 3;
  Bajie, Copia2_LTR, Dasheng, Retro1, RIRE2, RIRE3, Truncator 5; nDaiz,
  nDart 8; karma 12.
- **All 14 families are in `tests/data/TE_lib/RiceTE.fa`.** Lengths run from
  145 bp (mGing) to 10,989 bp (spmlike).
- **`RiceTE.fa` has three ambiguity codes** (one each of `S`, `M`, `R`).
  `te-panel` accepts only A/C/G/T, so the simulation uses a normalized copy.
  Callers still get the original `RiceTE.fa`.
- **Closest two events are 26 bp apart,** which `te-panel` accepts
  (footprint is `[P-max(L,5), P+5)`).
- **The mutated genome of the original fixture does not exist** (only its
  reads and truth, from RelocaTE2's `test_data`). This plan rebuilds it from
  the truth. Task 5 checks the rebuild against the original short reads.
- **`te-panel`** (`simulate-data te-panel --config X.json --output DIR`):
  JSON paths resolve relative to the JSON file; it refuses an existing output
  directory; it writes `genomes/<genome_id>.fa`, `genomes.tsv`
  (`genome_id sample haplotype role fasta`), `events.tsv`,
  `genome_truth.tsv`
  (`event_id genome_id path state start end duplicated_tsd_start duplicated_tsd_end`)
  and `.complete`. It needs at least one `sample` and one `reference` genome.
  FASTA headers are `sample#haplotype#chromosome`.
- **Toolkit state.** `/rhome/nmath020/bigdata/github/github_tools/data_sim/simulate_data`
  is on branch `feat/refactor-replan` with uncommitted work, and
  `te_panel.py` / `panel.py` are untracked there. **Run the toolkit; never
  edit, commit, stash or switch branches in it.** Use
  `pixi run --as-is --manifest-path <toolkit>/pyproject.toml ...`; a plain
  `pixi run` once tried to refresh the install and removed the launcher.
- **Long-read functions** in `simulate_data.longread`: `PLATFORMS["ont-hq"]`,
  `PLATFORMS["hifi"]`, `build_pbsim_command(platform, genome, depth, prefix, seed)`,
  `run_ccs(bam, output_fastq_gz, threads, min_yield)`,
  `pool_fastq_gz(source, destination, prefix)`, `GZIP_LEVEL`.
  `simulate_data.modules.te_benchmark_panel._pbsim_contig_outputs(prefix)`
  returns `(ref, maf, fq.gz-or-bam)` per contig. It is private; the rice
  panel driver imports from that module the same way.
- **Never run computation on the login node.** Everything past unit tests
  goes through `sbatch`.

## Panel layout (what this plan produces)

Default root: `validation/longread_small/panel_v1/` (gitignored). The
environment variable `RELOCATE3_LR_SMALL_PANEL` overrides it.

```text
panel_v1/
├── config/                  backbone.fa, te_library.fa, panel.json,
│                            te_library_changes.tsv
├── truth/                   the te-panel bundle (genomes/, events.tsv, ...)
├── truth_events.tsv         one row per insertion, reference coordinates
├── panel_manifest.tsv       sample platform coverage replicate reads control_reads
├── reads/<platform>/<sample>/reads.fastq.gz + .complete
├── reads/<platform>/<sample>_reference_control/reads.fastq.gz + .complete
├── qc/panel_check.tsv
└── .complete
```

- Platforms: `ont-hq`, `hifi`. Coverage: 5× and 20×. One replicate. Seed 916
  (the full panel's seed).
- All insertions are homozygous (one sample haplotype). Heterozygous and
  low-frequency sites are added as `panel_v2` when genotyping (L5) needs
  them.
- Expected size: about 130 MB.

---

### Task 0: Worktree

**Step 1: Create the worktree from the design branch**

Run from `RelocaTE3/` (the main checkout):

```bash
git worktree add ../RelocaTE3-lr-l1-small-panel -b feat/long-read-l1-small-panel design/long-read-algorithm
cd ../RelocaTE3-lr-l1-small-panel
pixi install
```

Expected: a new directory on branch `feat/long-read-l1-small-panel`
containing `plans/2026-11-09-long-read-algorithm-proposal.md`. If
`design/long-read-algorithm` has been merged, branch from `origin/main`
instead.

**Step 2: Baseline**

Run: `pixi run test`
Expected: all pass. Record the count in the Task 7 note.

---

### Task 1: Read the truth GFF

**Files:**
- Create: `validation/longread_small/__init__.py` (empty)
- Create: `validation/longread_small/make_panel_config.py`
- Test: `tests/longread_small_config_test.py`

**Step 1: Write the failing test**

```python
"""Tests for the small long-read panel config builder."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation" / "longread_small"))

import make_panel_config as mpc  # noqa: E402

DATA = ROOT / "tests" / "data"
GFF = DATA / "sim_genome" / "MSU7.Chr3_2M.ALL.gff"
GENOME = DATA / "sim_genome" / "MSU7.Chr3_2M.fa"
TE_LIB = DATA / "TE_lib" / "RiceTE.fa"


def test_read_truth_gff_parses_one_row(tmp_path):
    gff = tmp_path / "t.gff"
    gff.write_text("Chr3\tMSU7\tmPing\t100\t101\t.\t-\t.\tID=mPing.Chr3.100;TSD=TTA;\n")
    events = mpc.read_truth_gff(gff)
    assert events == [
        mpc.TruthEvent(
            event_id="mPing.Chr3.100",
            chrom="Chr3",
            position=100,
            family="mPing",
            strand="-",
            tsd="TTA",
        )
    ]


def test_read_truth_gff_rejects_missing_tsd(tmp_path):
    gff = tmp_path / "t.gff"
    gff.write_text("Chr3\tMSU7\tmPing\t100\t101\t.\t+\t.\tID=mPing.Chr3.100;\n")
    with pytest.raises(ValueError, match="TSD"):
        mpc.read_truth_gff(gff)


def test_fixture_has_200_events():
    assert len(mpc.read_truth_gff(GFF)) == 200
```

**Step 2: Run it to see it fail**

Run: `pixi run pytest tests/longread_small_config_test.py -q`
Expected: FAIL, `ModuleNotFoundError: No module named 'make_panel_config'`.

**Step 3: Minimal implementation**

`validation/longread_small/make_panel_config.py`:

```python
#!/usr/bin/env python3
"""Build the te-panel config for the small Chr3 long-read panel.

Input: the Chr3 2 Mb fixture (reference FASTA, truth GFF, TE library).
Output: a directory holding everything `simulate-data te-panel` needs.
Standard library only, so it runs in any Python >= 3.10.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class TruthEvent:
    """One simulated insertion, in reference coordinates."""

    event_id: str
    chrom: str
    position: int  # GFF column 4; the TSD is reference[position-L:position]
    family: str
    strand: str
    tsd: str


def read_truth_gff(path: Path) -> list[TruthEvent]:
    """Read the fixture's truth GFF."""
    events = []
    with open(path) as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(f"{path}:{number}: expected 9 columns")
            ident = re.search(r"ID=([^;]+)", fields[8])
            tsd = re.search(r"TSD=([ACGT]+)", fields[8])
            if ident is None or tsd is None:
                raise ValueError(f"{path}:{number}: missing ID or TSD attribute")
            if fields[6] not in "+-":
                raise ValueError(f"{path}:{number}: bad strand {fields[6]!r}")
            events.append(
                TruthEvent(
                    event_id=ident.group(1),
                    chrom=fields[0],
                    position=int(fields[3]),
                    family=fields[2],
                    strand=fields[6],
                    tsd=tsd.group(1),
                )
            )
    return events
```

**Step 4: Run tests**

Run: `pixi run pytest tests/longread_small_config_test.py -q`
Expected: 3 passed.

**Step 5: Commit**

```bash
git add validation/longread_small tests/longread_small_config_test.py
git commit -m "longread_small: read the Chr3 truth GFF"
```

---

### Task 2: Normalize the TE library and build the config

**Files:**
- Modify: `validation/longread_small/make_panel_config.py`
- Test: `tests/longread_small_config_test.py`

**Step 1: Write the failing tests** (append)

```python
def test_normalize_te_library_replaces_ambiguity_codes(tmp_path):
    src = tmp_path / "te.fa"
    src.write_text(">A desc\nacgS\nMR\n>B\nACGT\n")
    dst = tmp_path / "out.fa"
    changes = mpc.normalize_te_library(src, dst)
    assert dst.read_text() == ">A\nACGCAA\n>B\nACGT\n"
    assert changes == [("A", 4, "S", "C"), ("A", 5, "M", "A"), ("A", 6, "R", "A")]


def test_normalize_te_library_rejects_n(tmp_path):
    src = tmp_path / "te.fa"
    src.write_text(">A\nACNT\n")
    with pytest.raises(ValueError, match="N"):
        mpc.normalize_te_library(src, tmp_path / "out.fa")


def test_build_config_shape():
    events = [
        mpc.TruthEvent("mPing.Chr3.100", "Chr3", 100, "mPing", "-", "TTA"),
        mpc.TruthEvent("karma.Chr3.900", "Chr3", 900, "karma", "+", "A" * 12),
    ]
    config = mpc.build_config(events)
    assert config["schema_version"] == 1
    assert config["backbone"] == "backbone.fa"
    assert config["te_library"] == "te_library.fa"
    assert config["chromosomes"] == ["Chr3"]
    assert config["genomes"] == [
        {"id": "reference", "sample": "Reference", "haplotype": 1, "role": "reference"},
        {"id": "sample_h1", "sample": "Sample", "haplotype": 1, "role": "sample"},
    ]
    assert config["te_rules"] == {
        "karma": {"family": "karma", "tsd_length": 12},
        "mPing": {"family": "mPing", "tsd_length": 3},
    }
    assert config["events"][0] == {
        "id": "mPing.Chr3.100",
        "chromosome": "Chr3",
        "position": 100,
        "te_id": "mPing",
        "strand": "-",
        "carriers": ["sample_h1"],
    }


def test_build_config_rejects_mixed_tsd_length():
    events = [
        mpc.TruthEvent("a", "Chr3", 100, "mPing", "+", "TTA"),
        mpc.TruthEvent("b", "Chr3", 900, "mPing", "+", "TTAA"),
    ]
    with pytest.raises(ValueError, match="mPing"):
        mpc.build_config(events)


def test_fixture_tsd_matches_reference():
    genome = mpc.read_fasta(GENOME)
    for event in mpc.read_truth_gff(GFF):
        start = event.position - len(event.tsd)
        assert genome[event.chrom][start : event.position] == event.tsd, event
```

**Step 2: Run to see them fail**

Run: `pixi run pytest tests/longread_small_config_test.py -q`
Expected: FAIL, `AttributeError: module 'make_panel_config' has no attribute 'normalize_te_library'`.

**Step 3: Implementation** (add to `make_panel_config.py`)

```python
# First base of each IUPAC code, so the choice is deterministic and recorded.
AMBIGUITY = {
    "R": "A", "Y": "C", "S": "C", "W": "A", "K": "G", "M": "A",
    "B": "C", "D": "A", "H": "A", "V": "A",
}


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA into {first header word: uppercase sequence}."""
    records: dict[str, list[str]] = {}
    name = None
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError(f"{path}: duplicate record {name!r}")
                records[name] = []
            elif name is None:
                raise ValueError(f"{path}: sequence before first header")
            else:
                records[name].append(line.strip().upper())
    return {key: "".join(parts) for key, parts in records.items()}


def normalize_te_library(src: Path, dst: Path) -> list[tuple[str, int, str, str]]:
    """Write an A/C/G/T-only copy of the TE library.

    Returns the changes as (record, 1-based position, original, replacement).
    `N` is refused: there is no defensible single-base replacement.
    """
    changes = []
    lines = []
    for name, sequence in read_fasta(src).items():
        bases = list(sequence)
        for index, base in enumerate(bases):
            if base in "ACGT":
                continue
            if base not in AMBIGUITY:
                raise ValueError(f"{src}: {name} has unsupported base {base!r} (N?)")
            changes.append((name, index + 1, base, AMBIGUITY[base]))
            bases[index] = AMBIGUITY[base]
        lines.append(f">{name}\n{''.join(bases)}\n")
    dst.write_text("".join(lines))
    return changes


def build_config(events: list[TruthEvent]) -> dict:
    """Build the te-panel JSON: one reference genome, one sample haplotype."""
    tsd_length: dict[str, int] = {}
    for event in events:
        length = len(event.tsd)
        if tsd_length.setdefault(event.family, length) != length:
            raise ValueError(f"family {event.family} has more than one TSD length")
    chromosomes = sorted({event.chrom for event in events})
    return {
        "schema_version": 1,
        "backbone": "backbone.fa",
        "te_library": "te_library.fa",
        "chromosomes": chromosomes,
        "flank_size": 5,
        "genomes": [
            {"id": "reference", "sample": "Reference", "haplotype": 1, "role": "reference"},
            {"id": "sample_h1", "sample": "Sample", "haplotype": 1, "role": "sample"},
        ],
        "te_rules": {
            family: {"family": family, "tsd_length": tsd_length[family]}
            for family in sorted(tsd_length)
        },
        "events": [
            {
                "id": event.event_id,
                "chromosome": event.chrom,
                "position": event.position,
                "te_id": event.family,
                "strand": event.strand,
                "carriers": ["sample_h1"],
            }
            for event in events
        ],
    }
```

**Step 4: Run tests**

Run: `pixi run pytest tests/longread_small_config_test.py -q`
Expected: 8 passed.

**Step 5: Commit**

```bash
git add validation/longread_small/make_panel_config.py tests/longread_small_config_test.py
git commit -m "longread_small: TE library normalization and te-panel config"
```

---

### Task 3: Write the config directory (command line)

**Files:**
- Modify: `validation/longread_small/make_panel_config.py`
- Test: `tests/longread_small_config_test.py`

**Step 1: Write the failing test** (append)

```python
def test_write_config_dir_on_fixture(tmp_path):
    out = tmp_path / "config"
    mpc.write_config_dir(GFF, GENOME, TE_LIB, out)
    assert sorted(p.name for p in out.iterdir()) == [
        "backbone.fa",
        "panel.json",
        "te_library.fa",
        "te_library_changes.tsv",
        "truth_events.tsv",
    ]
    import json

    config = json.loads((out / "panel.json").read_text())
    assert len(config["events"]) == 200
    assert len(config["te_rules"]) == 14
    truth = (out / "truth_events.tsv").read_text().splitlines()
    assert truth[0].split("\t") == [
        "event_id", "chrom", "position", "te_family", "strand", "tsd",
        "tsd_length", "biological_class",
    ]
    assert len(truth) == 201
    changes = (out / "te_library_changes.tsv").read_text().splitlines()
    assert len(changes) == 4  # header + the three ambiguity codes
    assert set("".join(mpc.read_fasta(out / "te_library.fa").values())) <= set("ACGT")


def test_write_config_dir_refuses_existing(tmp_path):
    out = tmp_path / "config"
    out.mkdir()
    with pytest.raises(FileExistsError):
        mpc.write_config_dir(GFF, GENOME, TE_LIB, out)
```

**Step 2: Run to see it fail**

Run: `pixi run pytest tests/longread_small_config_test.py -q`
Expected: FAIL, no attribute `write_config_dir`.

**Step 3: Implementation** (add; put `import argparse, json, shutil, sys` at the top)

```python
TRUTH_HEADER = (
    "event_id", "chrom", "position", "te_family", "strand", "tsd",
    "tsd_length", "biological_class",
)


def write_config_dir(gff: Path, genome: Path, te_library: Path, outdir: Path) -> None:
    """Write backbone, normalized TE library, panel.json and truth table."""
    outdir = Path(outdir)
    if outdir.exists():
        raise FileExistsError(f"refusing to overwrite {outdir}")
    events = read_truth_gff(gff)
    reference = read_fasta(genome)
    for event in events:
        start = event.position - len(event.tsd)
        if reference[event.chrom][start : event.position] != event.tsd:
            raise ValueError(f"{event.event_id}: TSD does not match the reference")
    outdir.mkdir(parents=True)
    shutil.copyfile(genome, outdir / "backbone.fa")
    changes = normalize_te_library(te_library, outdir / "te_library.fa")
    with open(outdir / "te_library_changes.tsv", "w") as handle:
        handle.write("te_id\tposition\toriginal\treplacement\n")
        for row in changes:
            handle.write("\t".join(str(value) for value in row) + "\n")
    (outdir / "panel.json").write_text(
        json.dumps(build_config(events), indent=2, sort_keys=True) + "\n"
    )
    with open(outdir / "truth_events.tsv", "w") as handle:
        handle.write("\t".join(TRUTH_HEADER) + "\n")
        for event in events:
            handle.write(
                "\t".join(
                    (
                        event.event_id, event.chrom, str(event.position),
                        event.family, event.strand, event.tsd,
                        str(len(event.tsd)), "homozygous",
                    )
                )
                + "\n"
            )


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--te-library", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    write_config_dir(args.gff, args.genome, args.te_library, args.outdir)
    print(f"config written: {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**Step 4: Run tests**

Run: `pixi run pytest tests/longread_small_config_test.py -q`
Expected: 10 passed.

**Step 5: Commit**

```bash
git add validation/longread_small/make_panel_config.py tests/longread_small_config_test.py
git commit -m "longread_small: write the te-panel config directory"
```

---

### Task 4: Read simulation driver

This script runs inside the **toolkit's** environment, so RelocaTE3's test
suite cannot import it. Its pure parts (job list, manifest, skip rules) are
tested without the toolkit; the PBSIM3 call is exercised by the real build
in Task 6.

**Files:**
- Create: `validation/longread_small/simulate_reads.py`
- Test: `tests/longread_small_reads_test.py`

**Step 1: Write the failing test**

```python
"""Tests for the pure parts of the small-panel read driver."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation" / "longread_small"))

import simulate_reads as sr  # noqa: E402


def test_jobs_are_platform_by_coverage():
    jobs = sr.jobs(["ont-hq", "hifi"], [5, 20], replicates=1, seed=916)
    assert [(j.platform, j.sample) for j in jobs] == [
        ("ont-hq", "cov5x_rep1"),
        ("ont-hq", "cov20x_rep1"),
        ("hifi", "cov5x_rep1"),
        ("hifi", "cov20x_rep1"),
    ]
    assert len({j.seed for j in jobs} | {j.control_seed for j in jobs}) == 8


def test_manifest_matches_full_panel_columns(tmp_path):
    jobs = sr.jobs(["ont-hq"], [5], replicates=1, seed=916)
    sr.write_manifest(jobs, tmp_path / "panel_manifest.tsv")
    lines = (tmp_path / "panel_manifest.tsv").read_text().splitlines()
    assert lines[0].split("\t") == [
        "sample", "platform", "coverage", "replicate", "reads", "control_reads",
    ]
    assert lines[1].split("\t") == [
        "cov5x_rep1", "ont-hq", "5", "1",
        "reads/ont-hq/cov5x_rep1/reads.fastq.gz",
        "reads/ont-hq/cov5x_rep1_reference_control/reads.fastq.gz",
    ]


def test_needs_run_states(tmp_path):
    target = tmp_path / "reads" / "ont-hq" / "cov5x_rep1"
    assert sr.needs_run(target) is True  # absent
    target.mkdir(parents=True)
    assert sr.needs_run(target) is True  # empty directory
    (target / "reads.fastq.gz").write_bytes(b"x")
    with pytest.raises(RuntimeError, match="incomplete"):
        sr.needs_run(target)  # output without completion marker
    (target / ".complete").touch()
    assert sr.needs_run(target) is False
```

**Step 2: Run to see it fail**

Run: `pixi run pytest tests/longread_small_reads_test.py -q`
Expected: FAIL, `ModuleNotFoundError: No module named 'simulate_reads'`.

**Step 3: Implementation**

`validation/longread_small/simulate_reads.py`:

```python
#!/usr/bin/env python3
"""Simulate ONT-HQ and HiFi reads for the small Chr3 long-read panel.

Run inside the simulate_data toolkit's pixi environment (PBSIM3 + ccs). The
toolkit is imported only inside `simulate`, so the job list and manifest
code can be tested without it.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

MANIFEST_HEADER = ("sample", "platform", "coverage", "replicate", "reads", "control_reads")


@dataclass(frozen=True)
class Job:
    """One platform x coverage x replicate sample, plus its control."""

    platform: str
    coverage: int
    replicate: int
    seed: int
    control_seed: int

    @property
    def sample(self) -> str:
        return f"cov{self.coverage}x_rep{self.replicate}"

    @property
    def reads(self) -> str:
        return f"reads/{self.platform}/{self.sample}/reads.fastq.gz"

    @property
    def control_reads(self) -> str:
        return f"reads/{self.platform}/{self.sample}_reference_control/reads.fastq.gz"


def jobs(platforms: list[str], coverages: list[int], replicates: int, seed: int) -> list[Job]:
    """List jobs in platform, coverage, replicate order with distinct seeds."""
    out = []
    for platform in platforms:
        for coverage in coverages:
            for replicate in range(1, replicates + 1):
                index = len(out)
                out.append(
                    Job(platform, coverage, replicate, seed + 2 * index, seed + 2 * index + 1)
                )
    return out


def write_manifest(job_list: list[Job], path: Path) -> None:
    """Write the manifest with the full long-read panel's columns."""
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(MANIFEST_HEADER)
        for job in job_list:
            writer.writerow(
                (job.sample, job.platform, job.coverage, job.replicate, job.reads, job.control_reads)
            )


def needs_run(directory: Path) -> bool:
    """True if the directory is absent or empty; False if complete.

    Output without a completion marker is never overwritten.
    """
    if not directory.exists() or not any(directory.iterdir()):
        return True
    if (directory / ".complete").is_file():
        return False
    raise RuntimeError(
        f"incomplete output in {directory}: move it aside, do not overwrite it"
    )


def simulate(
    platform_key: str,
    genome: Path,
    coverage: float,
    seed: int,
    directory: Path,
    name_prefix: str,
    threads: int,
    min_ccs_yield: float,
) -> int:
    """Simulate one FASTQ into `directory`; return the read count."""
    from simulate_data.longread import (
        GZIP_LEVEL,
        PLATFORMS,
        build_pbsim_command,
        pool_fastq_gz,
        run_ccs,
    )
    from simulate_data.modules.te_benchmark_panel import _pbsim_contig_outputs
    from simulate_data.utils import run_command

    platform = PLATFORMS[platform_key]
    directory.mkdir(parents=True, exist_ok=True)
    partial = directory / "reads.fastq.gz.partial"
    count = 0
    yields = []
    with tempfile.TemporaryDirectory(dir=directory) as tmp:
        prefix = Path(tmp) / "sim"
        run_command(build_pbsim_command(platform, genome, coverage, prefix, seed))
        with gzip.open(partial, "wt", compresslevel=GZIP_LEVEL) as pooled:
            for _ref, _maf, payload in _pbsim_contig_outputs(prefix):
                if payload.suffix == ".bam":
                    fastq = Path(str(payload)[: -len(".bam")] + ".fastq.gz")
                    stats = run_ccs(payload, fastq, threads, min_ccs_yield)
                    yields.append(stats["yield_fraction"])
                else:
                    fastq = payload
                count += pool_fastq_gz(fastq, pooled, name_prefix)
    if count == 0:
        raise RuntimeError(f"no reads simulated for {directory}")
    partial.rename(directory / "reads.fastq.gz")
    with open(directory / "simulation.tsv", "w") as handle:
        handle.write("platform\tcoverage\tseed\treads\tccs_yield\n")
        mean_yield = f"{sum(yields) / len(yields):.4f}" if yields else "NA"
        handle.write(f"{platform_key}\t{coverage:g}\t{seed}\t{count}\t{mean_yield}\n")
    (directory / ".complete").touch()
    return count


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True, help="panel root")
    parser.add_argument("--sample-genome", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--platform", action="append", default=None)
    parser.add_argument("--coverage", action="append", type=int, default=None)
    parser.add_argument("--replicates", type=int, default=1)
    parser.add_argument("--seed", type=int, default=916)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--min-ccs-yield", type=float, default=0.80)
    args = parser.parse_args(argv)

    job_list = jobs(
        args.platform or ["ont-hq", "hifi"],
        args.coverage or [5, 20],
        args.replicates,
        args.seed,
    )
    for job in job_list:
        for genome, relative, seed, label in (
            (args.sample_genome, job.reads, job.seed, "sample"),
            (args.reference, job.control_reads, job.control_seed, "control"),
        ):
            directory = (args.panel / relative).parent
            if not needs_run(directory):
                print(f"skip (complete): {directory}")
                continue
            count = simulate(
                job.platform, genome, job.coverage, seed, directory,
                f"{job.sample}:{label}:s{seed}", args.threads, args.min_ccs_yield,
            )
            print(f"wrote {count} reads: {directory}")
    write_manifest(job_list, args.panel / "panel_manifest.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**Step 4: Run tests**

Run: `pixi run pytest tests/longread_small_reads_test.py -q`
Expected: 3 passed.

**Step 5: Commit**

```bash
git add validation/longread_small/simulate_reads.py tests/longread_small_reads_test.py
git commit -m "longread_small: read simulation driver (job list, manifest, skip rules)"
```

---

### Task 5: Panel check

The check runs in RelocaTE3's environment. Besides file integrity, it tests
the rebuild of the mutated genome: the fixture's original short reads must
fit the rebuilt sample genome much better than they fit the reference.

**Files:**
- Create: `validation/longread_small/check_panel.py`
- Test: `tests/longread_small_check_test.py`

**Step 1: Write the failing test**

```python
"""Tests for the small-panel checker's pure functions."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation" / "longread_small"))

import check_panel as cp  # noqa: E402


def _fastq(path: Path, lengths: list[int]) -> None:
    with gzip.open(path, "wt") as handle:
        for index, length in enumerate(lengths):
            handle.write(f"@r{index}\n{'A' * length}\n+\n{'I' * length}\n")


def test_fastq_stats(tmp_path):
    path = tmp_path / "r.fastq.gz"
    _fastq(path, [100, 300, 600])
    stats = cp.fastq_stats(path)
    assert stats == {"reads": 3, "bases": 1000, "n50": 600, "max_length": 600}


def test_fastq_stats_rejects_truncated(tmp_path):
    path = tmp_path / "r.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@r0\nACGT\n+\n")
    with pytest.raises(ValueError, match="truncated"):
        cp.fastq_stats(path)


def test_count_clipped_reads():
    sam = [
        "@SQ\tSN:Chr3\tLN:10",
        "r1\t0\tChr3\t1\t60\t100M\t*\t0\t0\t*\t*",
        "r2\t0\tChr3\t1\t60\t30S70M\t*\t0\t0\t*\t*",
        "r3\t0\tChr3\t1\t60\t80M20S\t*\t0\t0\t*\t*",
        "r4\t0\tChr3\t1\t60\t90M10S\t*\t0\t0\t*\t*",
        "r5\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*",
    ]
    assert cp.count_clipped_reads(sam, min_clip=20) == 2


def test_coverage_within_tolerance():
    assert cp.coverage_ok(bases=38_000_000, genome_length=2_000_000, requested=20) is True
    assert cp.coverage_ok(bases=20_000_000, genome_length=2_000_000, requested=20) is False
```

**Step 2: Run to see it fail**

Run: `pixi run pytest tests/longread_small_check_test.py -q`
Expected: FAIL, `ModuleNotFoundError: No module named 'check_panel'`.

**Step 3: Implementation**

`validation/longread_small/check_panel.py`:

```python
#!/usr/bin/env python3
"""Check a finished small long-read panel and write qc/panel_check.tsv.

Checks: manifest paths exist with completion markers; every FASTQ is a
complete gzip stream; delivered coverage is within tolerance; the truth
table has 200 rows; and the fixture's original short reads fit the rebuilt
sample genome (few clipped reads) far better than the reference.
Needs minimap2 on PATH for the last check.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
import sys
from pathlib import Path

COVERAGE_TOLERANCE = 0.25  # HiFi loses about 8% of molecules in ccs
MAX_CLIP_RATIO = 0.05


def fastq_stats(path: Path) -> dict[str, int]:
    """Read a whole gzipped FASTQ; return reads, bases, N50, longest read."""
    lengths = []
    with gzip.open(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not quality or not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"truncated or malformed FASTQ: {path}")
            lengths.append(len(sequence.strip()))
    lengths.sort(reverse=True)
    total = sum(lengths)
    running = 0
    n50 = 0
    for length in lengths:
        running += length
        if running * 2 >= total:
            n50 = length
            break
    return {
        "reads": len(lengths),
        "bases": total,
        "n50": n50,
        "max_length": lengths[0] if lengths else 0,
    }


def count_clipped_reads(sam_lines, min_clip: int) -> int:
    """Count mapped primary records with a soft clip of at least `min_clip`."""
    count = 0
    for line in sam_lines:
        if line.startswith("@"):
            continue
        fields = line.split("\t")
        flag = int(fields[1])
        if flag & 0x904:  # unmapped, secondary or supplementary
            continue
        clips = [int(n) for n in re.findall(r"(\d+)S", fields[5])]
        if clips and max(clips) >= min_clip:
            count += 1
    return count


def coverage_ok(bases: int, genome_length: int, requested: float) -> bool:
    """True if delivered coverage is within tolerance of the request."""
    delivered = bases / genome_length
    return abs(delivered - requested) <= COVERAGE_TOLERANCE * requested


def clipped_short_reads(genome: Path, r1: Path, r2: Path, threads: int) -> int:
    """Map the fixture's short reads to `genome`; count clipped reads."""
    process = subprocess.Popen(
        ["minimap2", "-ax", "sr", "-t", str(threads), str(genome), str(r1), str(r2)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    assert process.stdout is not None
    count = count_clipped_reads((line.rstrip("\n") for line in process.stdout), 20)
    if process.wait() != 0:
        raise RuntimeError(f"minimap2 failed on {genome}")
    return count


def fasta_length(path: Path) -> int:
    """Total sequence length of a FASTA."""
    with open(path) as handle:
        return sum(len(line.strip()) for line in handle if not line.startswith(">"))


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; exit 1 if any check fails."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--short-r1", type=Path, required=True)
    parser.add_argument("--short-r2", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args(argv)

    failures = []
    rows = []
    genome_length = fasta_length(args.reference)

    truth = (args.panel / "truth_events.tsv").read_text().splitlines()
    if len(truth) != 201:
        failures.append(f"truth_events.tsv has {len(truth) - 1} events, expected 200")

    with open(args.panel / "panel_manifest.tsv", newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    for entry in manifest:
        for column in ("reads", "control_reads"):
            fastq = args.panel / entry[column]
            if not (fastq.parent / ".complete").is_file():
                failures.append(f"no completion marker: {fastq.parent}")
                continue
            stats = fastq_stats(fastq)
            ok = coverage_ok(stats["bases"], genome_length, float(entry["coverage"]))
            if not ok:
                failures.append(f"coverage off target: {fastq}")
            rows.append(
                (
                    entry["platform"], entry["sample"], column, stats["reads"],
                    stats["bases"], f"{stats['bases'] / genome_length:.2f}",
                    stats["n50"], stats["max_length"], "ok" if ok else "FAIL",
                )
            )

    sample_genome = args.panel / "truth" / "genomes" / "sample_h1.fa"
    clipped_reference = clipped_short_reads(
        args.reference, args.short_r1, args.short_r2, args.threads
    )
    clipped_sample = clipped_short_reads(
        sample_genome, args.short_r1, args.short_r2, args.threads
    )
    if clipped_sample > MAX_CLIP_RATIO * clipped_reference:
        failures.append(
            f"rebuilt genome does not explain the original short reads: "
            f"{clipped_sample} clipped reads vs {clipped_reference} on the reference"
        )

    qc = args.panel / "qc"
    qc.mkdir(exist_ok=True)
    with open(qc / "panel_check.tsv", "w") as handle:
        handle.write(
            "platform\tsample\tkind\treads\tbases\tcoverage\tn50\tmax_length\tstatus\n"
        )
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")
        handle.write(f"# clipped_short_reads_reference\t{clipped_reference}\n")
        handle.write(f"# clipped_short_reads_sample\t{clipped_sample}\n")

    for failure in failures:
        print(f"FAIL: {failure}", file=sys.stderr)
    if failures:
        return 1
    (args.panel / ".complete").touch()
    print(f"RESULT: PASS ({len(rows)} FASTQs checked)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**Step 4: Run tests**

Run: `pixi run pytest tests/longread_small_check_test.py -q`
Expected: 4 passed.

**Step 5: Commit**

```bash
git add validation/longread_small/check_panel.py tests/longread_small_check_test.py
git commit -m "longread_small: panel checker"
```

---

### Task 6: SLURM build script, then build the panel

**Files:**
- Create: `validation/longread_small/build_panel.slurm`
- Modify: `.gitignore`

**Step 1: Ignore the panel output**

Append to `.gitignore`:

```text
/validation/longread_small/panel_v*/
```

**Step 2: Write the script**

`validation/longread_small/build_panel.slurm`:

```bash
#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --job-name=lr_small_panel
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH -o logs/lr_small_panel.%j.log

# Build the small Chr3 long-read panel. Submit from the repository root:
#   mkdir -p logs && sbatch validation/longread_small/build_panel.slurm
# Rerun-safe: finished steps are skipped; incomplete output is refused.

set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${BASE_DIR}"
date

THREADS="${SLURM_CPUS_PER_TASK:-1}"
SIM_PROJECT_DIR="${SIM_PROJECT_DIR:-/rhome/nmath020/bigdata/github/github_tools/data_sim/simulate_data}"
PANEL="${RELOCATE3_LR_SMALL_PANEL:-validation/longread_small/panel_v1}"
SCRIPTS="validation/longread_small"
GFF="tests/data/sim_genome/MSU7.Chr3_2M.ALL.gff"
GENOME="tests/data/sim_genome/MSU7.Chr3_2M.fa"
TE_LIB="tests/data/TE_lib/RiceTE.fa"
SHORT_R1="tests/data/sim_reads/MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
SHORT_R2="tests/data/sim_reads/MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"

echo "PANEL=${PANEL} THREADS=${THREADS} SIM_PROJECT_DIR=${SIM_PROJECT_DIR}"
echo "toolkit commit: $(git -C "${SIM_PROJECT_DIR}" rev-parse HEAD) (working tree may be dirty)"

for input in "${GFF}" "${GENOME}" "${TE_LIB}" "${SHORT_R1}" "${SHORT_R2}" \
    "${SIM_PROJECT_DIR}/pyproject.toml"; do
  if [[ ! -f "${input}" ]]; then
    echo "ERROR: Missing input: ${input}" >&2
    exit 1
  fi
done
command -v pixi >/dev/null || { echo "ERROR: pixi not on PATH" >&2; exit 1; }

sim() {
  pixi run --as-is --manifest-path "${SIM_PROJECT_DIR}/pyproject.toml" "$@"
}

mkdir -p "${PANEL}"

if [[ ! -d "${PANEL}/config" ]]; then
  python3 "${SCRIPTS}/make_panel_config.py" --gff "${GFF}" --genome "${GENOME}" \
    --te-library "${TE_LIB}" --outdir "${PANEL}/config"
fi
cp -n "${PANEL}/config/truth_events.tsv" "${PANEL}/truth_events.tsv"

if [[ ! -f "${PANEL}/truth/.complete" ]]; then
  if [[ -e "${PANEL}/truth" ]]; then
    echo "ERROR: incomplete ${PANEL}/truth; move it aside and resubmit" >&2
    exit 1
  fi
  sim simulate-data te-panel --config "${PANEL}/config/panel.json" \
    --output "${PANEL}/truth"
fi

sim python "${SCRIPTS}/simulate_reads.py" --panel "${PANEL}" \
  --sample-genome "${PANEL}/truth/genomes/sample_h1.fa" \
  --reference "${GENOME}" --threads "${THREADS}"

pixi run python "${SCRIPTS}/check_panel.py" --panel "${PANEL}" \
  --reference "${GENOME}" --short-r1 "${SHORT_R1}" --short-r2 "${SHORT_R2}" \
  --threads "${THREADS}"

date
```

**Step 3: Check the script parses**

Run: `bash -n validation/longread_small/build_panel.slurm`
Expected: no output.

**Step 4: Submit**

```bash
mkdir -p logs
sbatch validation/longread_small/build_panel.slurm
```

Expected: `Submitted batch job <id>`. Wait for it (`squeue -j <id>`); do not
poll in a tight loop.

**Step 5: Read the result**

Run: `tail -5 logs/lr_small_panel.<id>.log; cat validation/longread_small/panel_v1/qc/panel_check.tsv`
Expected: `RESULT: PASS (8 FASTQs checked)`; eight rows with coverage near 5
or 20; `clipped_short_reads_sample` at most 5% of
`clipped_short_reads_reference`.

**If the job fails, stop and report; do not adjust thresholds to pass.**
Likely causes and what they mean:
- `te-panel` rejects the config (overlapping footprints, unknown field):
  the config builder is wrong for this toolkit version. Fix the builder with
  a failing test first.
- `python3` on the node is older than 3.10: change that one call to
  `pixi run python`.
- The clipped-read check fails: the rebuilt genome is not the genome the
  original fixture was simulated from (for example, the original inserted
  truncated elements). Report the two counts to the maintainer. The fallback
  is to simulate matching short reads from the rebuilt genome so both read
  types share one truth; that is a design change, not a fix to make alone.
- `ccs` yield below 0.80: report the yield; do not lower the floor.

**Step 6: Commit**

```bash
git add .gitignore validation/longread_small/build_panel.slurm
git commit -m "longread_small: SLURM build script; ignore panel output"
```

---

### Task 7: Panel-present test, README, progress note

**Files:**
- Create: `tests/longread_small_panel_test.py`
- Create: `validation/longread_small/README.md`
- Create: `docs/2026-11-13-long-read-l1-small-panel.md`

**Step 1: Write the test** (skips when the panel has not been built)

```python
"""Shape of the built small long-read panel (skipped when it is absent)."""

from __future__ import annotations

import csv
import os
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PANEL = Path(
    os.environ.get(
        "RELOCATE3_LR_SMALL_PANEL", ROOT / "validation" / "longread_small" / "panel_v1"
    )
)

pytestmark = pytest.mark.skipif(
    not (PANEL / ".complete").is_file(),
    reason="small long-read panel not built (sbatch validation/longread_small/build_panel.slurm)",
)


def test_manifest_lists_four_samples_with_reads():
    with open(PANEL / "panel_manifest.tsv", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert sorted((r["platform"], r["coverage"]) for r in rows) == [
        ("hifi", "20"), ("hifi", "5"), ("ont-hq", "20"), ("ont-hq", "5"),
    ]
    for row in rows:
        for column in ("reads", "control_reads"):
            path = PANEL / row[column]
            assert path.is_file() and path.stat().st_size > 0, path


def test_truth_has_200_homozygous_events():
    with open(PANEL / "truth_events.tsv", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 200
    assert {r["biological_class"] for r in rows} == {"homozygous"}
```

**Step 2: Run the whole suite**

Run: `pixi run test`
Expected: the Task 0 baseline count plus 19 new tests (10 + 3 + 4 + 2), all
passing, with the panel built.

**Step 3: Write `validation/longread_small/README.md`**

Sections: Purpose (two sentences); Build (`sbatch` command, expected run
time taken from the job log); Layout (the tree above); Truth conventions
(position = GFF column 4, TSD = `reference[position-L:position]`, TE library
changes file); Use from tests (`RELOCATE3_LR_SMALL_PANEL`, the skip rule);
Limits (homozygous only; simulated error profiles are not yet compared with
the lab's real R10.4.1 reads; toolkit was run from an uncommitted working
tree, commit hash in the job log).

**Step 4: Write `docs/2026-11-13-long-read-l1-small-panel.md`**

Use the repository's note format: date and time, purpose, current status,
commands, the numbers from `qc/panel_check.tsv`, failures met, decisions,
next steps (L2 plan; benchmark manifest plan; HiFi 30× repair plan).

**Step 5: Lint and commit**

```bash
pre-commit run --files validation/longread_small/*.py tests/longread_small_*_test.py
git add tests/longread_small_panel_test.py validation/longread_small/README.md docs/2026-11-13-long-read-l1-small-panel.md
git commit -m "longread_small: panel test, README and progress note"
```

Expected: pre-commit passes (fix what it reports, then re-run).

**Step 6: Handoff**

End with a `HANDOFF: INTENT` paragraph for the reviewer: what was added, the
measured panel numbers, and the fact that no existing RelocaTE3 source file
was changed.

---

## Not in this plan

- **Benchmark long-read manifests** (`relocate-benchmark`: `lib/panel.py`
  requires `r1`/`r2`; sample keys must include the platform). The benchmark
  cannot score anything until RelocaTE3 produces long-read calls, so this is
  planned next to milestone L3, in that repository.
- **HiFi 30× repair** (`make_simulation_new`): archive the three partial
  tasks, split sample and control into separate jobs, rerun. Planned in that
  workspace; it blocks only the 30× HiFi benchmark rows.
- **Heterozygous and low-frequency sites** in the small panel (`panel_v2`,
  with milestone L5).
- **Comparing simulated read profiles with the real R10.4.1 data.**
