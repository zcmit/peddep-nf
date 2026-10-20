/*
 * Final analysis BAM: proper pairs, MAPQ filter, duplicates removed, reads
 * outside the include regions (mito, blacklist, non-canonical contigs) dropped,
 * then (paired-end) any mate left without its partner removed.
 */
process BAM_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::deeptools=3.5.6 bioconda::samtools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' :
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path include_bed

    output:
    tuple val(meta), path("*.filtered.bam"), path("*.filtered.bam.bai"), emit: bam
    path "versions.yml"                                                , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pairs  = meta.single_end ?
        "samtools view -@ ${task.cpus} -b ${args} -L ${include_bed} -o ${prefix}.filtered.bam ${bam}" :
        [
            "samtools view -@ ${task.cpus} -b ${args} -L ${include_bed} -o ${prefix}.view.bam ${bam}",
            "samtools sort -n -@ ${task.cpus} -T ${prefix}.nsort -o ${prefix}.nsort.bam ${prefix}.view.bam",
            "rm_orphan_pairs.py ${prefix}.nsort.bam ${prefix}.paired.bam",
            "samtools sort -@ ${task.cpus} -T ${prefix}.csort -o ${prefix}.filtered.bam ${prefix}.paired.bam",
            "rm ${prefix}.view.bam ${prefix}.nsort.bam ${prefix}.paired.bam"
        ].join('\n    ')
    """
    ${pairs}
    samtools index ${prefix}.filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(deeptools --version | sed "s/^deeptools //")
        samtools: \$(samtools --version | head -1 | sed "s/^samtools //")
        python: \$(python3 --version | sed "s/^Python //")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.bam ${prefix}.filtered.bam.bai versions.yml
    """
}
