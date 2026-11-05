"""GitHub Actions workflows: the image build/publish and the test suites.

``.github/workflows/container.yml`` builds ``containers/Dockerfile``, runs
``check_image.sh`` in the built image, pushes it to ghcr (not from forks), then
runs the pipeline in that image. ``.github/workflows/ci.yml`` runs pytest,
nf-test and ``nextflow lint`` in the locked pixi environments.

PyYAML is not in the pixi environments, so these are text checks on the
workflow files (with a small block splitter for jobs and steps). Syntax and
expressions are checked by actionlint (``pixi exec -s actionlint actionlint``),
not here.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = ROOT / ".github" / "workflows"
CONTAINER = WORKFLOWS / "container.yml"
CI = WORKFLOWS / "ci.yml"
DOCKERFILE = ROOT / "containers" / "Dockerfile"

PIXI_VERSION = "0.70.2"
IMAGE = "ghcr.io/stajichlab/relocate3"
# Push and login only for this repository's own events, never fork PRs.
NOT_FORK = (
    "github.event_name != 'pull_request' || "
    "github.event.pull_request.head.repo.full_name == github.repository"
)


def _read(path: Path) -> str:
    assert path.is_file(), f"missing {path.relative_to(ROOT)}"
    return path.read_text()


def _strip_comments(text: str) -> str:
    return "\n".join(
        line for line in text.splitlines() if not line.lstrip().startswith("#")
    )


def _job(text: str, name: str) -> str:
    """The block of job ``name`` (lines under ``  name:`` up to the next job)."""
    lines = _strip_comments(text).splitlines()
    try:
        start = lines.index(f"  {name}:")
    except ValueError:
        pytest.fail(f"no job {name!r}")
    block = []
    for line in lines[start + 1 :]:
        if re.match(r"^ {0,2}\S", line):
            break
        block.append(line)
    return "\n".join(block)


def _steps(job: str) -> list[str]:
    """Each step of a job block (split on the ``- `` that opens a step)."""
    parts = re.split(r"\n(?=      - )", job)
    return [p for p in parts if p.lstrip().startswith("- ")]


def _step_index(steps: list[str], needle: str) -> int:
    hits = [i for i, s in enumerate(steps) if needle in s]
    assert hits, f"no step containing {needle!r}"
    return hits[0]


# --------------------------------------------------------------------- common


@pytest.mark.parametrize("path", [CONTAINER, CI], ids=["container", "ci"])
def test_actions_pinned_to_released_versions(path):
    uses = re.findall(r"uses:\s*(\S+)", _strip_comments(_read(path)))
    assert uses
    for ref in uses:
        assert re.fullmatch(r"[\w.-]+/[\w.-]+@v\d+(\.\d+){0,2}", ref), ref


@pytest.mark.parametrize("path", [CONTAINER, CI], ids=["container", "ci"])
def test_triggers_pull_request_and_main(path):
    text = _strip_comments(_read(path))
    assert re.search(r"^  pull_request:", text, re.M)
    assert re.search(r"^  push:\n    branches:\s*\[\s*main\s*\]", text, re.M)


@pytest.mark.parametrize("path", [CONTAINER, CI], ids=["container", "ci"])
def test_checkout_fetches_full_history(path):
    # versioningit derives the version from git tags.
    text = _strip_comments(_read(path))
    n_checkout = len(re.findall(r"uses:\s*actions/checkout@", text))
    assert n_checkout >= 1
    assert text.count("fetch-depth: 0") == n_checkout


# ------------------------------------------------------------- container.yml


def test_container_triggers_version_tags():
    text = _strip_comments(_read(CONTAINER))
    assert re.search(r"^    tags:\s*\[\s*'v\*'\s*\]", text, re.M)


def test_container_permissions():
    text = _strip_comments(_read(CONTAINER))
    assert re.search(r"^permissions:\n  contents: read\n  packages: write\n", text, re.M)


def test_container_builds_pr_head_commit():
    # The image version (versioningit, from the copied .git HEAD) and the
    # sha-<short> tag must name the branch commit, not GitHub's merge commit.
    job = _job(_read(CONTAINER), "image")
    assert "ref: ${{ github.event.pull_request.head.sha || github.sha }}" in job
    assert "DOCKER_METADATA_PR_HEAD_SHA: true" in job


def test_image_checkout_does_not_persist_the_token():
    # The Dockerfile copies .git into the build stage; a persisted token
    # (http.extraheader in .git/config) would go with it.
    steps = _steps(_job(_read(CONTAINER), "image"))
    checkout = steps[_step_index(steps, "actions/checkout@")]
    assert "persist-credentials: false" in checkout


def test_container_build_inputs():
    steps = _steps(_job(_read(CONTAINER), "image"))
    builds = [s for s in steps if "docker/build-push-action@" in s]
    assert len(builds) == 2, "one build to load and test, one to push"
    for b in builds:
        assert "file: containers/Dockerfile" in b
        assert "context: ." in b
        assert "platforms: linux/amd64" in b


def test_container_metadata_tags():
    text = _strip_comments(_read(CONTAINER))
    assert re.search(rf"^env:\n  IMAGE: {re.escape(IMAGE)}\n", text, re.M)
    job = _job(text, "image")
    meta = _steps(job)[_step_index(_steps(job), "docker/metadata-action@")]
    assert "images: ${{ env.IMAGE }}" in meta
    assert "type=ref,event=pr" in meta  # pr-<n>
    assert "type=raw,value=dev,enable={{is_default_branch}}" in meta
    # Tags are PEP 440 (v0.2.0rc2): semver would not parse them.
    assert "type=pep440,pattern={{version}}" in meta
    assert "type=sha" in meta
    assert "latest=false" in meta


def test_container_checks_image_before_push():
    steps = _steps(_job(_read(CONTAINER), "image"))
    load = _step_index(steps, "load: true")
    check = _step_index(steps, "bash /opt/relocate3/check_image.sh")
    push = _step_index(steps, "push: true")
    assert load < check < push
    assert "push: true" not in steps[load]
    assert "docker run --rm" in steps[check]


RELEASE_REF = "startsWith(github.ref, 'refs/tags/v')"


def test_container_builds_pass_the_release_tag():
    # build_env.sh builds exactly the tag being released (RELEASE_TAG), even
    # with an rc and a final tag on one commit. Both builds get the same
    # build-arg, so the push build reuses the tested build's cache.
    steps = _steps(_job(_read(CONTAINER), "image"))
    builds = [s for s in steps if "docker/build-push-action@" in s]
    arg = f"RELEASE_TAG=${{{{ {RELEASE_REF} && github.ref_name || '' }}}}"
    for b in builds:
        assert re.search(r"build-args:\s*\|?\s*\n?\s*" + re.escape(arg), b), b


def test_container_release_check_runs_on_tags_before_push():
    steps = _steps(_job(_read(CONTAINER), "image"))
    check = _step_index(steps, "bash /opt/relocate3/check_image.sh")
    release = _step_index(steps, "containers/check_release.sh")
    login = _step_index(steps, "docker/login-action@")
    push = _step_index(steps, "push: true")
    assert check < release < login < push
    step = steps[release]
    assert f"if: {RELEASE_REF}" in step
    assert 'docker run --rm "${TEST_TAG}" relocaTE3 --version' in step
    assert '"${GITHUB_REF_NAME#v}"' in step
    assert "nextflow.config" in step


def test_container_push_and_login_skip_forks():
    steps = _steps(_job(_read(CONTAINER), "image"))
    for needle in ("docker/login-action@", "push: true"):
        step = steps[_step_index(steps, needle)]
        assert f"if: {NOT_FORK}" in step, needle
    login = steps[_step_index(steps, "docker/login-action@")]
    assert "registry: ghcr.io" in login
    assert "password: ${{ secrets.GITHUB_TOKEN }}" in login


def test_pipeline_container_job_runs_after_image():
    job = _job(_read(CONTAINER), "pipeline-container")
    assert re.search(r"^    needs:\s*image\s*$", job, re.M)
    assert f"if: {NOT_FORK}" in job
    # Task 3 adds the docker profile, container_tag and the parity mode.
    assert re.search(r"-profile test,docker\b.*--container_tag", job)
    assert "tests/nextflow/parity_chr3.sh docker" in job


def test_pipeline_container_job_checks_the_recorded_image():
    # pipeline_info/software_environment.txt names the image the run used.
    job = _job(_read(CONTAINER), "pipeline-container")
    assert "pipeline_info/software_environment.txt" in job
    assert 'container: ${IMAGE}:${TAG}' in job


# ------------------------------------------------------------------- ci.yml


def test_ci_setup_pixi_locked_and_pinned():
    job = _job(_read(CI), "tests")
    step = _steps(job)[_step_index(_steps(job), "prefix-dev/setup-pixi@")]
    assert f"pixi-version: v{PIXI_VERSION}" in step
    assert "locked: true" in step
    assert re.search(r"environments:\s*pipeline nf\b", step)


def test_ci_pixi_matches_container_build():
    # The CI pixi and the image's pixi are the one that wrote pixi.lock.
    assert f"ghcr.io/prefix-dev/pixi:{PIXI_VERSION}-" in _read(DOCKERFILE)


def test_ci_runs_all_suites():
    job = _job(_read(CI), "tests")
    assert re.search(r"-m pytest tests\b", job)
    assert ".pixi/envs/pipeline/bin" in job  # the real tools (and multiqc) on PATH
    assert "pixi run -e nf nf-test test" in job
    assert "pixi run -e nf nextflow lint ." in job
