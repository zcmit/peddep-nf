process GUNZIP {
    tag "${archive}"
    label 'process_single'

    conda "conda-forge::coreutils=9.5 conda-forge::gzip=1.13"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data' :
        'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8' }"

    input:
    path archive

    output:
    path "${archive.name - ~/\.gz$/}", emit: file
    path "versions.yml"              , emit: versions

    script:
    def unzipped = archive.name - ~/\.gz$/
    """
    gzip -cd ${archive} > ${unzipped}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(gzip --version 2>&1 | head -1 | sed 's/^gzip //')
    END_VERSIONS
    """

    stub:
    def unzipped = archive.name - ~/\.gz$/
    """
    touch ${unzipped} versions.yml
    """
}
