/*
 * PEDDEP graded QC report for one ChIP-seq / CUT&RUN library (case or control).
 */
process CHIP_QC_REPORT {
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
          path(markdup_metrics),
          path(complexity),
          path(filtered_flagstat),
          path(spp),
          path(frip, stageAs: 'frip/*'),
          path(rose, stageAs: 'rose/*')
    path thresholds

    output:
    tuple val(meta), path("*.chip_qc_report.txt") , emit: txt
    tuple val(meta), path("*.chip_qc_report.html"), emit: html
    tuple val(meta), path("*.chip_qc_metrics.json"), emit: json
    path "versions.yml"                           , emit: versions

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def rose_files = [rose].flatten()
    def stitched   = rose_files.find { f -> f.name.endsWith('AllStitched.table.txt') }
    def superse    = rose_files.find { f -> f.name.endsWith('SuperStitched.table.txt') }
    def opts = [
        meta.single_end ? '--single-end' : '',
        meta.is_control ? '--control' : '',
        "--peak-type ${meta.peak_type}",
        frip ? "--frip ${frip}" : '',
        stitched ? "--rose-stitched ${stitched}" : '',
        superse ? "--rose-super ${superse}" : ''
    ].findAll { o -> o }.join(' ')
    """
    chip_qc_report.py \\
        --sample ${meta.id} \\
        --prefix ${prefix} \\
        --thresholds ${thresholds} \\
        --fastqc fastqc_raw/* \\
        --trim-logs trim/* \\
        --bowtie2-logs bowtie2/* \\
        --markdup-metrics ${markdup_metrics} \\
        --complexity ${complexity} \\
        --filtered-flagstat ${filtered_flagstat} \\
        --spp ${spp} \\
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
    touch ${prefix}.chip_qc_report.txt ${prefix}.chip_qc_report.html versions.yml
    echo '{"sample": "${meta.id}", "overall": {"label": "GOOD", "score": 1, "mean": 1.0}, "metrics": [{"name": "NSC", "value": 1.1, "display": "1.1000", "score": 2, "graded": true, "description": "", "source": "", "rule": ""}]}' > ${prefix}.chip_qc_metrics.json
    """
}
