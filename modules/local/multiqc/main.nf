process MULTIQC {
    label 'process_single'

    conda "bioconda::multiqc=1.35"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data' :
        'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc' }"

    input:
    path files, stageAs: 'input/?/*'
    path config

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "versions.yml"        , emit: versions

    script:
    def args  = task.ext.args ?: ''
    def title = params.multiqc_title ? "--title \"${params.multiqc_title}\"" : ''
    """
    multiqc --force ${args} ${title} --config ${config} input/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_report.html versions.yml
    """
}
