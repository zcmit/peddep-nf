/*
 * One table for the whole run (TSV + coloured HTML) and the MultiQC section.
 */
process ATAC_QC_SUMMARY {
    label 'process_single'

    conda "bioconda::deeptools=3.5.6 bioconda::samtools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' :
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    path metrics, stageAs: 'metrics/*'

    output:
    path "atac_qc_summary.tsv"        , emit: tsv
    path "atac_qc_summary.html"       , emit: html
    path "atac_qc_grades_mqc.html"    , emit: multiqc
    path "versions.yml"               , emit: versions

    script:
    """
    atac_qc_summary.py --out-prefix atac_qc_summary --multiqc atac_qc_grades_mqc.html metrics/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    """
    touch atac_qc_summary.tsv atac_qc_summary.html atac_qc_grades_mqc.html versions.yml
    """
}
