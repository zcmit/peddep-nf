/*
 * PEDDEP graded QC report for one sample: every metric scored -2..2, the mean
 * score gives POOR / BELOW_AVERAGE / AVERAGE / GOOD / EXCELLENT.
 */
process ATAC_QC_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::deeptools=3.5.6 bioconda::samtools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' :
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta),
          path(raw_fastqc, stageAs: 'fastqc_raw/*'),
          path(trim_logs, stageAs: 'trim/*'),
          path(bowtie2_logs, stageAs: 'bowtie2/*'),
          path(idxstats),
          path(markdup_metrics),
          path(complexity),
          path(filtered_flagstat),
          path(frip),
          path(ataqv_json),
          path(fragment_qc)
    path thresholds
    val  mito_name

    output:
    tuple val(meta), path("*.atac_qc_report.txt") , emit: txt
    tuple val(meta), path("*.atac_qc_report.html"), emit: html
    tuple val(meta), path("*.atac_qc_metrics.json"), emit: json
    path "versions.yml"                           , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    atac_qc_report.py \\
        --sample ${meta.id} \\
        --prefix ${prefix} \\
        --thresholds ${thresholds} \\
        --mito-name '${mito_name}' \\
        --fastqc fastqc_raw/* \\
        --trim-logs trim/* \\
        --bowtie2-logs bowtie2/* \\
        --idxstats ${idxstats} \\
        --markdup-metrics ${markdup_metrics} \\
        --complexity ${complexity} \\
        --filtered-flagstat ${filtered_flagstat} \\
        --frip ${frip} \\
        --ataqv ${ataqv_json} \\
        --fragment-qc ${fragment_qc}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.atac_qc_report.txt ${prefix}.atac_qc_report.html versions.yml
    echo '{"sample": "${meta.id}", "overall": {"label": "GOOD", "score": 1}, "metrics": []}' > ${prefix}.atac_qc_metrics.json
    """
}
