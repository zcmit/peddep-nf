/*
 * HiChIP alignment: each read of a pair aligned on its own (-SP) with chimeric
 * ligation reads split into 5'-primary pieces (-5), every alignment reported (-T0)
 * so pairtools can choose. Output is the unsorted BAM pairtools parse reads.
 */
process BWAMEM2_MEM {
    tag "${meta.id}_L${meta.lane}_${meta.chunk}"
    label 'process_high'

    conda "bioconda::bwa-mem2=2.3 bioconda::htslib=1.22.1 bioconda::samtools=1.22.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05ce34b46ad42810eb29f74e4e304c0cb592b2ca15572929ed8bbaee58faf01/data' :
        'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113' }"

    input:
    tuple val(meta), path(reads, stageAs: 'input*/*')
    path index, stageAs: 'bwamem2_index/*'

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_L${meta.lane}_${meta.chunk}"
    def fastqs = [reads].flatten()
    """
    INDEX=\$(find -L bwamem2_index/ -name "*.bwt.2bit.64" | head -1 | sed 's/\\.bwt\\.2bit\\.64\$//')
    if [ -z "\$INDEX" ]; then
        echo "ERROR: no bwa-mem2 index (*.bwt.2bit.64) found in the staged --bwamem2_index" >&2
        exit 1
    fi

    bwa-mem2 mem \\
        -t ${task.cpus} \\
        ${args} \\
        \$INDEX \\
        ${fastqs[0]} ${fastqs[1]} \\
        | samtools view -b -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>/dev/null | tail -1)
        samtools: \$(samtools --version | head -1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_L${meta.lane}_${meta.chunk}"
    """
    touch ${prefix}.bam versions.yml
    """
}
