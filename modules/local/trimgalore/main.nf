process TRIMGALORE {
    tag "${meta.id}_L${meta.lane}"
    label 'process_medium'

    conda "bioconda::trim-galore=0.6.10"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/trim-galore:0.6.10--hdfd78af_1' :
        'biocontainers/trim-galore:0.6.10--hdfd78af_1' }"

    input:
    tuple val(meta), path(reads, stageAs: 'input*/*')

    output:
    tuple val(meta), path("*{_val_1,_val_2,_trimmed}.fq.gz"), emit: reads
    tuple val(meta), path("*_trimming_report.txt")    , emit: log
    tuple val(meta), path("*_fastqc.zip")             , emit: zip , optional: true
    tuple val(meta), path("*_fastqc.html")            , emit: html, optional: true
    path "versions.yml"                               , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastqs = [reads].flatten()
    def ext    = fastqs[0].name.endsWith('.gz') ? 'fastq.gz' : 'fastq'
    def names  = meta.single_end ? ["${prefix}.${ext}"] : ["${prefix}_R1.${ext}", "${prefix}_R2.${ext}"]
    def links  = [fastqs, names].transpose().collect { pair -> "ln -s ${pair[0]} ${pair[1]}" }.join('\n    ')
    def paired = meta.single_end ? '' : '--paired'
    // Trim Galore uses ~3 threads per core plus overhead (see its CHANGELOG for v0.6.0)
    def cores  = Math.min(8, Math.max(1, task.cpus - 4))
    """
    ${links}

    trim_galore --cores ${cores} ${paired} --gzip ${args} ${names.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$(trim_galore --version | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        echo "" | gzip > ${prefix}_trimmed.fq.gz
        touch ${prefix}.fastq.gz_trimming_report.txt ${prefix}_trimmed_fastqc.zip versions.yml
        """
    } else {
        """
        echo "" | gzip > ${prefix}_R1_val_1.fq.gz
        echo "" | gzip > ${prefix}_R2_val_2.fq.gz
        touch ${prefix}_R1.fastq.gz_trimming_report.txt ${prefix}_R2.fastq.gz_trimming_report.txt
        touch ${prefix}_R1_val_1_fastqc.zip ${prefix}_R2_val_2_fastqc.zip versions.yml
        """
    }
}
