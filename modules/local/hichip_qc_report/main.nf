/*
 * PEDDEP graded QC report for one HiChIP library (target or INPUT).
 */
process HICHIP_QC_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::deeptools=3.5.6 bioconda::samtools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' :
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta),
          path(raw_fastqc, stageAs: 'fastqc_raw/*'),
          path(dedup_stats),
          path(frip, stageAs: 'frip/*'),
          path(loops, stageAs: 'loops/*')
    path thresholds

    output:
    tuple val(meta), path("*.hichip_qc_report.txt")  , emit: txt
    tuple val(meta), path("*.hichip_qc_report.html") , emit: html
    tuple val(meta), path("*.hichip_qc_metrics.json"), emit: json
    path "versions.yml"                              , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def opts = [
        meta.is_control ? '--control' : '',
        frip ? "--frip ${frip}" : '',
        loops ? "--loops ${loops}" : ''
    ].findAll { o -> o }.join(' ')
    """
    hichip_qc_report.py \\
        --sample ${meta.id} \\
        --prefix ${prefix} \\
        --thresholds ${thresholds} \\
        --fastqc fastqc_raw/* \\
        --dedup-stats ${dedup_stats} \\
        ${opts} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hichip_qc_report.txt ${prefix}.hichip_qc_report.html versions.yml
    echo '{"sample": "${meta.id}", "overall": {"label": "GOOD", "score": 1, "mean": 1.0}, "metrics": [{"name": "Cis_Pct", "value": 80.0, "display": "80.00", "score": 2, "graded": true, "description": "", "source": "", "rule": ""}]}' > ${prefix}.hichip_qc_metrics.json
    """
}
