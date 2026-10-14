/*
 * Fraction of reads in peaks and the peak count, from the filtered BAM and the
 * narrowPeak file (summit-split duplicate rows and overlaps are merged first).
 */
process FRIP {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::deeptools=3.5.6 bioconda::samtools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' :
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(peaks)

    output:
    tuple val(meta), path("*.frip.txt"), emit: txt
    path "versions.yml"                , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cut -f1-3 ${peaks} | sort -k1,1 -k2,2n -u > peaks.unique.bed
    awk -v OFS='\\t' '
        NR == 1 { c = \$1; s = \$2; e = \$3; next }
        \$1 == c && \$2 <= e { if (\$3 > e) e = \$3; next }
        { print c, s, e; c = \$1; s = \$2; e = \$3 }
        END { if (NR > 0) print c, s, e }
    ' peaks.unique.bed > peaks.merged.bed

    TOTAL=\$(samtools view -c -@ ${task.cpus} ${bam})
    if [ -s peaks.merged.bed ]; then
        IN_PEAKS=\$(samtools view -c -@ ${task.cpus} -L peaks.merged.bed ${bam})
    else
        IN_PEAKS=0
    fi
    N_PEAKS=\$(wc -l < peaks.unique.bed | tr -d ' ')
    FRIP=\$(awk -v a="\$IN_PEAKS" -v b="\$TOTAL" 'BEGIN { if (b > 0) printf "%.4f", a / b; else printf "0" }')

    printf 'Sample\\tReads_In_Peaks\\tTotal_Reads\\tFRiP\\tPeaks\\n' > ${prefix}.frip.txt
    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "${meta.id}" "\$IN_PEAKS" "\$TOTAL" "\$FRIP" "\$N_PEAKS" >> ${prefix}.frip.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(deeptools --version | sed "s/^deeptools //")
        samtools: \$(samtools --version | head -1 | sed "s/^samtools //")
        python: \$(python3 --version | sed "s/^Python //")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'Sample\tReads_In_Peaks\tTotal_Reads\tFRiP\tPeaks\n${meta.id}\t30\t100\t0.3\t10\n' > ${prefix}.frip.txt
    touch versions.yml
    """
}
