/*
 * Genome-browser peak track: sorted, de-duplicated 3-column BED, bgzipped and
 * tabix-indexed (the format ProteinPaint "bedj" and IGV tracks load directly).
 */
process BROWSER_PEAKS {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::samtools=1.24 bioconda::htslib=1.24"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd' }"

    input:
    tuple val(meta), path(peaks)

    output:
    tuple val(meta), path("*.bed.gz"), path("*.bed.gz.tbi"), emit: bed
    path "versions.yml"                                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_peaks"
    """
    cut -f1-3 ${peaks} | LC_ALL=C sort -k1,1 -k2,2n -u | bgzip -c > ${prefix}.bed.gz
    tabix -p bed ${prefix}.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htslib: \$(tabix --version | head -1 | sed 's/^tabix (htslib) //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_peaks"
    """
    echo "" | gzip > ${prefix}.bed.gz
    touch ${prefix}.bed.gz.tbi versions.yml
    """
}
