/*
 * One library: merge the chunks, remove PCR/optical duplicates, keep contacts
 * between the kept contigs (genome.sizes), and write
 *   <id>.pairs.gz        deduplicated contacts (4DN pairs format, bgzipped)
 *   <id>.dedup.stats     pairtools stats: pair types, duplicates, cis/trans, distances
 *   <id>.ends.bam        each read end of every kept contact as a 1-bp alignment, outside
 *                        the blacklist and mitochondria -- input for peaks, coverage and FRiP
 */
process PAIRTOOLS_DEDUP {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::pairtools=1.1.3"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cc/cc26df7ecf3557f443f89f4dc83380c68b413809cb6ca2213f14ea405c726080/data' :
        'community.wave.seqera.io/library/pairtools:1.1.3--f7f48429aa47aeaf' }"

    input:
    tuple val(meta), path(chunks, stageAs: 'chunks/*')
    path chrom_sizes
    path include_bed

    output:
    tuple val(meta), path("${meta.id}.pairs.gz")               , emit: pairs
    tuple val(meta), path("${meta.id}.dedup.stats")            , emit: stats
    tuple val(meta), path("${meta.id}.ends.bam"), path("${meta.id}.ends.bam.bai"), emit: bam
    path "versions.yml"                                        , emit: versions

    script:
    def args     = task.ext.args ?: ''
    def prefix   = meta.id
    def sort_mem = "${Math.max(1, (task.memory.toGiga() * 0.5 / task.cpus).intValue())}G"
    """
    pairtools merge --nproc ${task.cpus} chunks/* \\
        | pairtools dedup ${args} --output-stats ${prefix}.dedup.stats \\
        | pairtools select --chrom-subset ${chrom_sizes} 'True' --nproc-out 4 -o ${prefix}.pairs.gz

    # read ends -> 1-bp alignments (flag 16 = minus strand); MAPQ was checked by pairtools parse
    (
        awk -v OFS='\\t' '{ print "@SQ", "SN:" \$1, "LN:" \$2 }' ${chrom_sizes}
        bgzip -dc ${prefix}.pairs.gz | awk -v OFS='\\t' '
            !/^#/ {
                n++
                print "c" n "/1", (\$6 == "-") ? 16 : 0, \$2, \$3, 60, "1M", "*", 0, 0, "*", "*"
                print "c" n "/2", (\$7 == "-") ? 16 : 0, \$4, \$5, 60, "1M", "*", 0, 0, "*", "*"
            }'
    ) | samtools view -u -L ${include_bed} - \\
      | samtools sort -@ ${task.cpus} -m ${sort_mem} -T ends_sort -o ${prefix}.ends.bam -
    samtools index ${prefix}.ends.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version 2>/dev/null | sed 's/^pairtools, version //')
        samtools: \$(samtools --version | head -1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    echo "" | gzip > ${prefix}.pairs.gz
    printf 'total\\t1000\\ntotal_unmapped\\t10\\ntotal_single_sided_mapped\\t40\\ntotal_mapped\\t950\\ntotal_dups\\t200\\ntotal_nodups\\t750\\ncis\\t600\\ntrans\\t150\\ncis_1kb+\\t400\\ncis_20kb+\\t250\\n' > ${prefix}.dedup.stats
    touch ${prefix}.ends.bam ${prefix}.ends.bam.bai versions.yml
    """
}
