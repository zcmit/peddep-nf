/*
 * Signal around transcription start sites: matrix, average profile and heatmap.
 */
process DEEPTOOLS_TSS_PLOTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::deeptools=3.5.6 bioconda::samtools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' :
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta), path(bigwig)
    path tss

    output:
    tuple val(meta), path("*.computeMatrix.mat.gz"), emit: matrix
    tuple val(meta), path("*.plotProfile.tab")     , emit: profile_data
    tuple val(meta), path("*.tss_profile.png")     , emit: profile_png
    tuple val(meta), path("*.tss_heatmap.png")     , emit: heatmap_png
    path "versions.yml"                            , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib

    computeMatrix reference-point \\
        --referencePoint TSS \\
        --regionsFileName ${tss} \\
        --scoreFileName ${bigwig} \\
        --samplesLabel ${meta.id} \\
        --missingDataAsZero \\
        --skipZeros \\
        --numberOfProcessors ${task.cpus} \\
        ${args} \\
        --outFileName ${prefix}.computeMatrix.mat.gz

    plotProfile \\
        --matrixFile ${prefix}.computeMatrix.mat.gz \\
        --outFileName ${prefix}.tss_profile.png \\
        --outFileNameData ${prefix}.plotProfile.tab \\
        --plotTitle "${meta.id}: signal around TSS"

    plotHeatmap \\
        --matrixFile ${prefix}.computeMatrix.mat.gz \\
        --outFileName ${prefix}.tss_heatmap.png \\
        --plotTitle "${meta.id}" \\
        --colorMap Blues

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
    echo "" | gzip > ${prefix}.computeMatrix.mat.gz
    touch ${prefix}.plotProfile.tab ${prefix}.tss_profile.png ${prefix}.tss_heatmap.png versions.yml
    """
}
