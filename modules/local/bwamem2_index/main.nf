/*
 * bwa-mem2 index of the reference FASTA (HiChIP). Built once per run unless
 * --bwamem2_index is given; keep it with --save_reference. A human genome needs
 * about 64 GB of memory to index.
 */
process BWAMEM2_INDEX {
    tag "${fasta}"
    label 'process_single'
    label 'process_long'

    conda "bioconda::bwa-mem2=2.3 bioconda::htslib=1.22.1 bioconda::samtools=1.22.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05ce34b46ad42810eb29f74e4e304c0cb592b2ca15572929ed8bbaee58faf01/data' :
        'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113' }"

    input:
    path fasta

    output:
    path "bwamem2"      , emit: index
    path "versions.yml" , emit: versions

    script:
    """
    mkdir bwamem2
    bwa-mem2 index -p bwamem2/${fasta.baseName} ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>/dev/null | tail -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir bwamem2
    for ext in 0123 amb ann bwt.2bit.64 pac; do touch bwamem2/${fasta.baseName}.\$ext; done
    touch versions.yml
    """
}
