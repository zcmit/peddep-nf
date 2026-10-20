/*
 * Splits one lane's read pairs into chunks of --split_fastq pairs so the chunks
 * align side by side. Chunks are written with fast gzip (level 1): they live only
 * in the work directory.
 */
process SEQKIT_SPLIT2 {
    tag "${meta.id}_L${meta.lane}"
    label 'process_low'

    conda "bioconda::seqkit=2.13.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/29/297e985dd395d87f9089b314dbe685ddce59fa0770c4d1f5e1009bc7b6ec7a8d/data' :
        'community.wave.seqera.io/library/seqkit:2.13.0--3bacfb68370dd835' }"

    input:
    tuple val(meta), path(reads, stageAs: 'input*/*')
    val chunk_pairs

    output:
    tuple val(meta), path("chunks/*_R{1,2}.part_*.fastq.gz"), emit: reads
    path "versions.yml"                                     , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_L${meta.lane}"
    def fastqs = [reads].flatten()
    """
    ln -s ${fastqs[0]} ${prefix}_R1.fastq.gz
    ln -s ${fastqs[1]} ${prefix}_R2.fastq.gz
    seqkit split2 \\
        --threads ${task.cpus} \\
        --compress-level 1 \\
        --by-size ${chunk_pairs} \\
        --out-dir chunks \\
        -1 ${prefix}_R1.fastq.gz \\
        -2 ${prefix}_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/^seqkit v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_L${meta.lane}"
    """
    mkdir chunks
    for i in 001 002; do
        echo "" | gzip > chunks/${prefix}_R1.part_\$i.fastq.gz
        echo "" | gzip > chunks/${prefix}_R2.part_\$i.fastq.gz
    done
    touch versions.yml
    """
}
