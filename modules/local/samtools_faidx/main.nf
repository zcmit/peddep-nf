process SAMTOOLS_FAIDX {
    tag "${fasta}"
    label 'process_single'

    conda "bioconda::samtools=1.24 bioconda::htslib=1.24"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd' }"

    input:
    path fasta

    output:
    path "${fasta}.fai"  , emit: fai
    path "${fasta}.sizes", emit: sizes
    path "versions.yml"  , emit: versions

    script:
    """
    samtools faidx ${fasta}
    cut -f1,2 ${fasta}.fai > ${fasta}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    printf 'chr1\\t1000000\\t6\\t60\\t61\\nchrM\\t16569\\t1016955\\t60\\t61\\n' > ${fasta}.fai
    cut -f1,2 ${fasta}.fai > ${fasta}.sizes
    touch versions.yml
    """
}
