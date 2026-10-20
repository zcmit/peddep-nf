/*
 * Sizes of the contigs HiChIP contacts are kept on: the canonical contigs
 * (--canonical_chroms) minus the mitochondrial one unless --keep_mito, in FASTA
 * order (contact maps list chromosomes in this order).
 */
process CHROM_SIZES {
    tag "${fai}"
    label 'process_single'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    path fai
    val  mito_name
    val  keep_mito
    val  canonical_chroms

    output:
    path "genome.sizes", emit: sizes
    path "versions.yml", emit: versions

    script:
    def regex     = canonical_chroms ?: '.*'
    def drop_mito = keep_mito ? 'false' : 'true'
    """
    awk -v OFS='\\t' -v re='${regex}' -v m='${mito_name}' -v drop='${drop_mito}' \\
        '\$1 ~ re && !(drop == "true" && \$1 == m) { print \$1, \$2 }' ${fai} > genome.sizes
    if [ ! -s genome.sizes ]; then
        echo "ERROR: no contigs in ${fai} match --canonical_chroms '${regex}'" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/^bedtools v//')
    END_VERSIONS
    """

    stub:
    """
    printf 'chr1\\t1000000\\n' > genome.sizes
    touch versions.yml
    """
}
