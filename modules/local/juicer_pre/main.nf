/*
 * Juicebox / ProteinPaint .hic contact map from the deduplicated pairs. Uses the
 * Juicer Tools 1.22.01 build from the HiLOW image, so the .hic version is the one
 * the PEDDEP ProteinPaint tracks already read. Needs a container engine.
 */
process JUICER_PRE {
    tag "${meta.id}"
    label 'process_medium'

    container 'ghcr.io/stjude/abralab/hilow:v1.0'

    input:
    tuple val(meta), path(pairs)
    path chrom_sizes

    output:
    tuple val(meta), path("*.hic"), emit: hic
    path "versions.yml"           , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def heap   = "${Math.max(2, (task.memory.toGiga() * 0.85).intValue())}g"
    """
    mkdir -p juicer_tmp
    java -Djava.awt.headless=true -Xmx${heap} -jar /usr/local/bin/juicer_tools_1.22.01.jar pre \\
        --threads ${task.cpus} \\
        -t juicer_tmp \\
        ${args} \\
        ${pairs} \\
        ${prefix}.hic \\
        ${chrom_sizes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicer_tools: 1.22.01
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hic versions.yml
    """
}
