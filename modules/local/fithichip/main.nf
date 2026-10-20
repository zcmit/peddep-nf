/*
 * FitHiChIP loops with the HiLOW settings: anchors = the library's 1D peaks (or the
 * samplesheet's `peaks` file), peak-to-all interactions, 5-kb bins, 5 kb - 4 Mb,
 * coverage bias correction, nearby loops merged. FitHiChIP 10 from the HiLOW image,
 * the version behind the PEDDEP loops so far. Needs a container engine.
 *
 * As in HiLOW, contacts with an end in the blacklist and anchors overlapping it are
 * left out. Besides FitHiChIP's own files, loops at each --fithichip_qvalues cutoff
 * are written as BEDPE and as a both-direction BED (bgzip + tabix), the layout of
 * the HiLOW ProteinPaint arc tracks.
 */
process FITHICHIP {
    tag "${meta.id}"
    label 'process_medium'

    container 'ghcr.io/stjude/abralab/hilow:v1.0'

    input:
    tuple val(meta), path(pairs), path(peaks)
    path chrom_sizes
    path blacklist

    output:
    tuple val(meta), path("*.loops.q*.bedpe")                              , emit: bedpe
    tuple val(meta), path("*.loops.q*.bed.gz"), path("*.loops.q*.bed.gz.tbi"), emit: browser
    tuple val(meta), path("*.FitHiChIP.*")                                 , emit: fithichip
    tuple val(meta), path("*.anchors.bed")                                 , emit: anchors
    tuple val(meta), path("*.loops.txt")                                   , emit: counts
    path "versions.yml"                                                    , emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def bin      = params.fithichip_bin_size
    def low      = params.fithichip_low_dist
    def upp      = params.fithichip_upp_dist
    def int_type = params.fithichip_int_type
    def qvalues  = params.fithichip_qvalues.toString().tokenize(',').collect { q -> q.trim() }.join(' ')
    def type_dir = int_type.toString() == '1' ? "FitHiChIP_Peak2Peak_b${bin}_L${low}_U${upp}/Coverage_Bias" :
                                                "FitHiChIP_Peak2ALL_b${bin}_L${low}_U${upp}/P2PBckgr_0/Coverage_Bias"
    """
    # Blacklist, merged so a binary search can look up read ends
    if [ -s "${blacklist}" ]; then
        cut -f1-3 ${blacklist} | sort -k1,1 -k2,2n | bedtools merge -i - > blacklist.merged.bed
    else
        : > blacklist.merged.bed
    fi

    # Anchors: peak intervals outside the blacklist
    awk -v OFS='\\t' '!/^(#|track|browser)/ && NF >= 3 { print \$1, \$2, \$3 }' ${peaks} | sort -k1,1 -k2,2n -u > anchors.all.bed
    if [ -s blacklist.merged.bed ]; then
        bedtools intersect -v -a anchors.all.bed -b blacklist.merged.bed > ${prefix}.anchors.bed
    else
        cp anchors.all.bed ${prefix}.anchors.bed
    fi

    LOOPDIR=fithichip_out/${type_dir}/FitHiC_BiasCorr
    ALL=\$LOOPDIR/FitHiChIP.interactions_FitHiC.bed
    SIG=\$LOOPDIR/FitHiChIP.interactions_FitHiC_Q0.01.bed
    MERGED=\$LOOPDIR/Merge_Nearby_Interactions/FitHiChIP.interactions_FitHiC_Q0.01_MergeNearContacts.bed

    if [ -s ${prefix}.anchors.bed ]; then
        # Contacts in HiC-Pro validPairs layout (readID chr1 pos1 strand1 chr2 pos2 strand2),
        # dropping those with an end in the blacklist
        bgzip -dc ${pairs} | awk -v OFS='\\t' -v bl=blacklist.merged.bed '
            function hit(chr, pos,   lo, hi, mid) {
                lo = 1; hi = n[chr] + 0
                while (lo <= hi) {
                    mid = int((lo + hi) / 2)
                    if (pos <= s[chr, mid]) hi = mid - 1
                    else if (pos > e[chr, mid]) lo = mid + 1
                    else return 1
                }
                return 0
            }
            BEGIN {
                while ((getline line < bl) > 0) {
                    split(line, f, "\\t"); k = ++n[f[1]]; s[f[1], k] = f[2]; e[f[1], k] = f[3]
                }
            }
            /^#/ { next }
            hit(\$2, \$3) || hit(\$4, \$5) { dropped++; next }
            { print "c" NR, \$2, \$3, \$6, \$4, \$5, \$7 }
            END { print dropped + 0 > "blacklist_dropped.txt" }
        ' | gzip -1 > validpairs.txt.gz

        printf '%s\\n' \\
            "ValidPairs=\$PWD/validpairs.txt.gz" \\
            "Interval=" \\
            "Matrix=" \\
            "Bed=" \\
            "CircularGenome=0" \\
            "PeakFile=\$PWD/${prefix}.anchors.bed" \\
            "OutDir=\$PWD/fithichip_out" \\
            "IntType=${int_type}" \\
            "BINSIZE=${bin}" \\
            "LowDistThr=${low}" \\
            "UppDistThr=${upp}" \\
            "UseP2PBackgrnd=0" \\
            "BiasType=1" \\
            "MergeInt=1" \\
            "QVALUE=0.01" \\
            "ChrSizeFile=\$PWD/${chrom_sizes}" \\
            "PREFIX=FitHiChIP" \\
            "OverWrite=1" \\
            > fithichip.config

        bash /FitHiChIP/FitHiChIP_HiCPro.sh -C fithichip.config
        rm -f validpairs.txt.gz
        if [ ! -s "\$ALL" ]; then
            echo "ERROR: FitHiChIP wrote no interaction table (\$ALL); see the log above" >&2
            exit 1
        fi
    else
        echo "WARNING: ${meta.id} has no anchors outside the blacklist, so no loops are called" >&2
        echo 0 > blacklist_dropped.txt
        mkdir -p \$LOOPDIR/Merge_Nearby_Interactions
        printf 'chr1\\ts1\\te1\\tchr2\\ts2\\te2\\tcc\\tp\\tQ\\n' > \$ALL
    fi
    gzip -c \$ALL > ${prefix}.FitHiChIP.interactions.bed.gz
    if [ -f "\$SIG" ]; then cp \$SIG ${prefix}.FitHiChIP.Q0.01.bed; else head -1 \$ALL > ${prefix}.FitHiChIP.Q0.01.bed; fi
    if [ -f "\$MERGED" ]; then cp \$MERGED ${prefix}.FitHiChIP.Q0.01.merged.bed; else head -1 \$ALL > ${prefix}.FitHiChIP.Q0.01.merged.bed; fi

    # Loops per q-value cutoff (last column of the interaction table), HiLOW layout
    for Q in ${qvalues}; do
        awk -F '\\t' -v OFS='\\t' -v q=\$Q 'NR > 1 && \$NF < q { print \$1, \$2, \$3, \$4, \$5, \$6, \$7 }' \$ALL \\
            | sort -k1,1 -k2,2n > ${prefix}.loops.q\$Q.bedpe
        awk -F '\\t' -v OFS='\\t' '{ print \$1, \$2, \$3, \$4, \$5, \$6, \$7; print \$4, \$5, \$6, \$1, \$2, \$3, \$7 }' ${prefix}.loops.q\$Q.bedpe \\
            | sort -k1,1 -k2,2n | bgzip -c > ${prefix}.loops.q\$Q.bed.gz
        tabix -f -p bed ${prefix}.loops.q\$Q.bed.gz
    done

    count() { if [ -s "\$1" ]; then echo \$(( \$(wc -l < "\$1") - \$2 )); else echo 0; fi; }
    printf 'Sample\\tAnchors\\tContacts_Blacklisted\\tInteractions_Tested\\tLoops_Q0.01\\tMerged_Loops_Q0.01\\n' > ${prefix}.loops.txt
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "${meta.id}" "\$(wc -l < ${prefix}.anchors.bed | tr -d ' ')" "\$(cat blacklist_dropped.txt)" \\
        "\$(count \$ALL 1)" "\$(count ${prefix}.FitHiChIP.Q0.01.bed 1)" "\$(count ${prefix}.FitHiChIP.Q0.01.merged.bed 1)" >> ${prefix}.loops.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fithichip: \$(grep -m1 -o 'Version [0-9.]*' /FitHiChIP/Release_Notes.txt | sed 's/Version //')
        bedtools: \$(bedtools --version | sed 's/^bedtools v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for Q in 0.01 0.05 0.1; do
        touch ${prefix}.loops.q\$Q.bedpe ${prefix}.loops.q\$Q.bed.gz ${prefix}.loops.q\$Q.bed.gz.tbi
    done
    touch ${prefix}.FitHiChIP.interactions.bed.gz ${prefix}.FitHiChIP.Q0.01.bed ${prefix}.anchors.bed versions.yml
    printf 'Sample\\tAnchors\\tContacts_Blacklisted\\tInteractions_Tested\\tLoops_Q0.01\\tMerged_Loops_Q0.01\\n${meta.id}\\t10\\t0\\t100\\t5\\t4\\n' > ${prefix}.loops.txt
    """
}
