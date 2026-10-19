/*
 * ROSE stitched regions / super-enhancer-like regions (ROSE 1.3.2, the version the
 * PEDDEP/SEAseq pipeline uses). Narrow peaks with p <= --rose_peak_pvalue (ROSE's
 * usual MACS p 1e-9 input) are merged into unionpeaks.gff, so output names
 * (unionpeaks_AllStitched.table.txt, ...) match the old STITCHED_peaks folder.
 * ROSE's read mapping (ROSE_bamToGFF.py) is swapped for bin/rose_bamtogff.py, which
 * gives the same table in minutes instead of hours.
 * Needs a container engine: ROSE is not packaged for conda.
 */
process ROSE {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/stjude/abralab/rose:v1.3.2'

    input:
    tuple val(meta), path(bam), path(bai), path(control), path(control_bai), path(peaks)
    path refseq

    output:
    tuple val(meta), path("STITCHED_peaks")                                 , emit: dir
    tuple val(meta), path("STITCHED_peaks/unionpeaks_{All,Super}Stitched.table.txt"), emit: tables, optional: true
    path "versions.yml"                                                     , emit: versions

    script:
    def stitch   = params.rose_stitch ?: 12500
    def tss      = params.rose_tss ?: 2000
    def min_logp = params.rose_peak_pvalue ? -Math.log10(params.rose_peak_pvalue as double) : 0
    def ctrl     = control ? "-c ${control}" : ''
    """
    # narrowPeak column 8 is -log10(p)
    awk -F'\\t' -v t=${min_logp} '\$8 >= t' ${peaks} | cut -f1-3 | sort -k1,1 -k2,2n | mergeBed -i - \\
        | awk -F'\\t' '{ print \$1 "\\t" NR "\\t\\t" \$2 "\\t" \$3 "\\t\\t.\\t\\t" NR }' > unionpeaks.gff

    mkdir -p rose_fast
    ln -s "\$(command -v rose_bamtogff.py)" rose_fast/ROSE_bamToGFF.py
    export PATH="\$PWD/rose_fast:\$PATH"

    if [ -s unionpeaks.gff ]; then
        ROSE_main.py -s ${stitch} -t ${tss} --custom ${refseq} -i unionpeaks.gff -r ${bam} ${ctrl} -o STITCHED_peaks
    else
        mkdir -p STITCHED_peaks
        echo "No peaks with p <= ${params.rose_peak_pvalue} for ${meta.id}: ROSE not run" > STITCHED_peaks/NO_PEAKS.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rose: 1.3.2
    END_VERSIONS
    """

    stub:
    """
    mkdir -p STITCHED_peaks
    touch STITCHED_peaks/unionpeaks_AllStitched.table.txt STITCHED_peaks/unionpeaks_SuperStitched.table.txt versions.yml
    """
}
