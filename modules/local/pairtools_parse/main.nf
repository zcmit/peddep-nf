/*
 * Read pairs -> sorted contact pairs for one alignment chunk. Ligation junctions
 * inside a read are resolved with --walks-policy 5unique; the SAM fields and read
 * IDs are dropped, which keeps the pairs files about 20x smaller.
 */
process PAIRTOOLS_PARSE {
    tag "${meta.id}_L${meta.lane}_${meta.chunk}"
    label 'process_medium'

    conda "bioconda::pairtools=1.1.3"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cc/cc26df7ecf3557f443f89f4dc83380c68b413809cb6ca2213f14ea405c726080/data' :
        'community.wave.seqera.io/library/pairtools:1.1.3--f7f48429aa47aeaf' }"

    input:
    tuple val(meta), path(bam)
    path chrom_sizes

    output:
    tuple val(meta), path("*.pairs.gz"), emit: pairs
    path "versions.yml"                , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}_L${meta.lane}_${meta.chunk}"
    def sort_mem  = "${Math.max(1, (task.memory.toGiga() * 0.6).intValue())}G"
    """
    mkdir -p sort_tmp
    pairtools parse \\
        --drop-sam \\
        --drop-readid \\
        --chroms-path ${chrom_sizes} \\
        --nproc-in 2 \\
        --nproc-out 2 \\
        ${args} \\
        ${bam} \\
        | pairtools sort \\
            --nproc ${task.cpus} \\
            --memory ${sort_mem} \\
            --tmpdir sort_tmp \\
            -o ${prefix}.pairs.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version 2>/dev/null | sed 's/^pairtools, version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_L${meta.lane}_${meta.chunk}"
    """
    echo "" | gzip > ${prefix}.pairs.gz
    touch versions.yml
    """
}
