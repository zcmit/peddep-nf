/*
 * ataqv: TSS enrichment, mitochondrial fraction, fragment-length distribution
 * and peak metrics. Needs the duplicate-marked BAM that still contains mito reads.
 */
process ATAQV {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::ataqv=1.3.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/ataqv:1.3.1--py310ha155cf9_1' :
        'biocontainers/ataqv:1.3.1--py310ha155cf9_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(peaks)
    path tss
    path excluded_regions
    path autosomes
    val  mito_name

    output:
    tuple val(meta), path("*.ataqv.json"), emit: json
    tuple val(meta), path("*.ataqv.txt") , emit: txt
    path "versions.yml"                  , emit: versions

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def excluded = excluded_regions ? "--excluded-region-file ${excluded_regions}" : ''
    """
    ataqv \\
        ${args} \\
        --threads ${task.cpus} \\
        --name ${meta.id} \\
        --ignore-read-groups \\
        --peak-file ${peaks} \\
        --tss-file ${tss} \\
        ${excluded} \\
        --autosomal-reference-file ${autosomes} \\
        --mitochondrial-reference-name ${mito_name} \\
        --metrics-file ${prefix}.ataqv.json \\
        NA \\
        ${bam} \\
        > ${prefix}.ataqv.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ataqv: \$(ataqv --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '[{"metrics": {"name": "${meta.id}", "tss_enrichment": 8.5}}]' > ${prefix}.ataqv.json
    touch ${prefix}.ataqv.txt versions.yml
    """
}
