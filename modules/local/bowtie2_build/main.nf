process BOWTIE2_BUILD {
    tag "${fasta}"
    label 'process_high'

    conda "bioconda::bowtie2=2.5.4 bioconda::samtools=1.21 bioconda::htslib=1.21 conda-forge::pigz=2.8"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b41b403e81883126c3227fc45840015538e8e2212f13abc9ae84e4b98891d51c/data' :
        'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6' }"

    input:
    path fasta

    output:
    path "bowtie2"      , emit: index
    path "versions.yml" , emit: versions

    script:
    """
    mkdir bowtie2
    bowtie2-build --threads ${task.cpus} ${fasta} bowtie2/${fasta.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | head -1 | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir bowtie2
    for ext in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do touch bowtie2/${fasta.baseName}.\$ext; done
    touch versions.yml
    """
}
