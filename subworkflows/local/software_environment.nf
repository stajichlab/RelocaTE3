/*
 * The software environment of a run, written by main.nf to
 * pipeline_info/software_environment.txt. Tool versions are
 * in software_versions.yml. Nextflow resolves a container tag when it pulls
 * and does not expose the image digest, so this records the image reference
 * and says where the digest is. tests/nextflow_profiles_test.py runs it under
 * each profile.
 */
def softwareEnvironment() {
    def lines = [
        "date: ${workflow.start}",
        "run_name: ${workflow.runName}",
        "session_id: ${workflow.sessionId}",
        "pipeline: ${workflow.manifest.name} ${workflow.manifest.version}",
        "revision: ${workflow.revision ?: 'none'}",
        "commit: ${workflow.commitId ?: 'none'}",
        "profile: ${workflow.profile}",
        "container_engine: ${workflow.containerEngine ?: 'none'}",
        "container: ${workflow.container ?: 'none'}",
    ]
    if (workflow.containerEngine) {
        lines << "digest: not recorded by Nextflow; apptainer/singularity: `apptainer inspect` of the image in NXF_APPTAINER_CACHEDIR (org.opencontainers.image.revision); docker: `docker image inspect --format '{{index .RepoDigests 0}}' <image>`"
    }
    if (workflow.profile.tokenize(',').contains('pixi')) {
        lines << "pixi_env: ${params.pixi_env}"
    }
    return lines.join('\n') + '\n'
}
