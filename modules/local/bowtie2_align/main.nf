process BOWTIE2_ALIGN {
    tag "${meta.id}_L${meta.lane}"
    label 'process_high'

    conda "bioconda::bowtie2=2.5.4 bioconda::samtools=1.21 bioconda::htslib=1.21 conda-forge::pigz=2.8"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b41b403e81883126c3227fc45840015538e8e2212f13abc9ae84e4b98891d51c/data' :
        'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6' }"

    input:
    tuple val(meta), path(reads)
    path index, stageAs: 'bowtie2_index/*'

    output:
    tuple val(meta), path("*.sorted.bam")  , emit: bam
    tuple val(meta), path("*.bowtie2.log") , emit: log
    path "versions.yml"                    , emit: versions

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def sort_cpus  = Math.max(1, task.cpus.intdiv(4))
    def fastqs     = [reads].flatten()
    def read_args  = meta.single_end ? "-U ${fastqs[0]}" : "-1 ${fastqs[0]} -2 ${fastqs[1]}"
    """
    INDEX=\$(find -L bowtie2_index/ -name "*.rev.1.bt2" -o -name "*.rev.1.bt2l" | head -1 | sed -E 's/\\.rev\\.1\\.bt2l?\$//')
    if [ -z "\$INDEX" ]; then
        echo "ERROR: no bowtie2 index (*.rev.1.bt2) found in the staged --bowtie2_index" >&2
        exit 1
    fi

    bowtie2 \\
        -p ${task.cpus} \\
        ${args} \\
        --rg-id ${prefix} \\
        --rg SM:${meta.id} \\
        --rg PL:ILLUMINA \\
        -x \$INDEX \\
        ${read_args} \\
        2> ${prefix}.bowtie2.log \\
        | samtools sort -@ ${sort_cpus} -T ${prefix}.sorttmp -o ${prefix}.sorted.bam -

    cat ${prefix}.bowtie2.log >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | head -1 | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(samtools --version | head -1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.bam ${prefix}.bowtie2.log versions.yml
    """
}
