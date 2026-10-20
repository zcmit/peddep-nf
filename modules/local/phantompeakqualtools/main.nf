/*
 * Strand cross-correlation (NSC, RSC, fragment length) with run_spp.R.
 * Paired-end: read 1 only, with NSC/RSC taken at the median insert size (-speak); SPP's own
 * estimate can lock onto shift 0 for CUT&RUN. Subsampled to at most --spp_max_reads reads
 * (ENCODE uses 15M).
 */
process PHANTOMPEAKQUALTOOLS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::phantompeakqualtools=1.2.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4a/4a1cddfad5b503ee347cc5de17d172e1876c547fca00aa844559c9e764fb400f/data' :
        'community.wave.seqera.io/library/phantompeakqualtools:1.2.2--f8026fe2526a5e18' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.spp.out"), emit: spp
    tuple val(meta), path("*.spp.pdf"), emit: pdf
    path "versions.yml"               , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def max_reads = params.spp_max_reads ?: 15000000
    def select    = meta.single_end ? '' : '-f 64'
    """
    samtools view -@ ${task.cpus} -b ${select} -o ${prefix}.all.bam ${bam}
    N=\$(samtools view -c -@ ${task.cpus} ${prefix}.all.bam)
    if [ "\$N" -gt ${max_reads} ]; then
        FRAC=\$(awk -v n="\$N" 'BEGIN { printf "%.6f", ${max_reads} / n }')
        samtools view -@ ${task.cpus} -b -s "42\${FRAC#0}" -o ${prefix}.bam ${prefix}.all.bam
        rm ${prefix}.all.bam
    else
        mv ${prefix}.all.bam ${prefix}.bam
    fi

    # paired-end: median insert size (|TLEN| of read 1) as the cross-correlation peak
    SPEAK=''
    if [ "${meta.single_end}" != "true" ]; then
        SPEAK=\$(samtools view ${prefix}.bam \\
            | awk -v max=${params.max_fragment_length} '{ t = \$9 < 0 ? -\$9 : \$9 } t > 0 && t <= max { c[t]++; n++ }
                END { if (n) { s = 0; for (i = 1; i <= max; i++) { s += c[i]; if (s >= n / 2) { print i; exit } } } }')
    fi

    RUN_SPP=\$(command -v run_spp.R)
    Rscript -e "library(caTools); source(\\"\$RUN_SPP\\")" \\
        -c=${prefix}.bam \\
        -p=${task.cpus} \\
        -savp=${prefix}.spp.pdf \\
        -out=${prefix}.spp.out \\
        -rf \\
        -tmpdir=. \\
        \${SPEAK:+-speak=\$SPEAK} \\
        ${args}
    rm ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phantompeakqualtools: 1.2.2
        samtools: \$(samtools --version | head -1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '${meta.id}\\t1000\\t180,200,220\\t0.3\\t50\\t0.25\\t1500\\t0.2\\t1.12\\t1.05\\t1\\n' > ${prefix}.spp.out
    touch ${prefix}.spp.pdf versions.yml
    """
}
