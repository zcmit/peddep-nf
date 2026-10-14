process GTF2TSS {
    tag "${gtf}"
    label 'process_single'

    conda "bioconda::deeptools=3.5.6 bioconda::samtools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' :
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    path gtf

    output:
    path "tss.bed"      , emit: tss
    path "versions.yml" , emit: versions

    script:
    """
    gtf2tss.py --gtf ${gtf} --out tss.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'chr1\\t999\\t1000\\tGENE1\\t0\\t+\\n' > tss.bed
    touch versions.yml
    """
}
