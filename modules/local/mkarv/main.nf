process MKARV {
    label 'process_single'

    conda "bioconda::ataqv=1.3.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/ataqv:1.3.1--py310ha155cf9_1' :
        'biocontainers/ataqv:1.3.1--py310ha155cf9_1' }"

    input:
    path jsons, stageAs: 'json/*'

    output:
    path "ataqv_report" , emit: html
    path "versions.yml" , emit: versions

    script:
    """
    mkarv --concurrency ${task.cpus} --force ataqv_report json/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ataqv: \$(ataqv --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir ataqv_report
    touch ataqv_report/index.html versions.yml
    """
}
