/*
 * Builds the regions reads are allowed to fall in (canonical contigs minus the
 * mitochondrial contig minus the blacklist) and the autosome list used by ataqv.
 */
process GENOME_FILTER {
    tag "${fai}"
    label 'process_single'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    path fai
    path blacklist
    val  mito_name
    val  keep_mito
    val  canonical_chroms
    val  sex_chroms

    output:
    path "include_regions.bed", emit: include_bed
    path "autosomes.txt"      , emit: autosomes
    path "blacklist.bed"      , emit: blacklist
    path "versions.yml"       , emit: versions

    script:
    def regex     = canonical_chroms ?: '.*'
    def drop_mito = keep_mito ? 'false' : 'true'
    def blacklist_cmd = blacklist ?
        "cut -f1-3 ${blacklist} | sort -k1,1 -k2,2n > blacklist.bed" :
        "touch blacklist.bed"
    """
    awk -v OFS='\\t' -v re='${regex}' '\$1 ~ re { print \$1, 0, \$2 }' ${fai} | sort -k1,1 -k2,2n > canonical.bed
    if [ ! -s canonical.bed ]; then
        echo "ERROR: no contigs in ${fai} match --canonical_chroms '${regex}'" >&2
        exit 1
    fi

    if [ "${drop_mito}" = "true" ]; then
        awk -v OFS='\\t' -v m='${mito_name}' '\$1 != m' canonical.bed > keep.bed
    else
        cp canonical.bed keep.bed
    fi

    ${blacklist_cmd}
    if [ -s blacklist.bed ]; then
        bedtools subtract -a keep.bed -b blacklist.bed > include_regions.bed
    else
        cp keep.bed include_regions.bed
    fi

    echo '${mito_name}' | tr ',' '\\n' > not_autosomes.txt
    echo '${sex_chroms}' | tr ',' '\\n' >> not_autosomes.txt
    cut -f1 canonical.bed | grep -v -x -F -f not_autosomes.txt > autosomes.txt || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/^bedtools v//')
    END_VERSIONS
    """

    stub:
    """
    printf 'chr1\\t0\\t1000000\\n' > include_regions.bed
    echo chr1 > autosomes.txt
    touch blacklist.bed versions.yml
    """
}
