nextflow.enable.dsl = 2

process RSEM_PREPARE_REFERENCE {
    publishDir "${params.outdir}/rsem_reference", mode: 'copy'
    conda 'bioconda::rsem=1.3.3'
    cpus 16
    memory '64.GB'

    input:
    path(genome_fasta)
    path(gtf)

    output:
    path "*"

    script:
    """
    rsem-prepare-reference --gtf ${gtf} \\
        --num-threads ${task.cpus} \\
        ${genome_fasta} \\
        rsem_ref 
    """
}

process RSEM_QUANTIFICATION {
    tag "$sample_id"
    publishDir "${params.outdir}/rsem/${sample_id}", mode: 'copy'
    conda 'bioconda::rsem=1.3.3'
    cpus 32
    memory '64.GB'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample_id), path(transcriptome_bam)
    path(rsem_ref_dir)

    output:
    tuple val(sample_id), path("${sample_id}.genes.results"), path("${sample_id}.isoforms.results"), emit: results
    path "${sample_id}.stat", emit: stats
    
    script:
    """
    rsem-calculate-expression --paired-end \\
        --bam \\
        -p ${task.cpus} \\
        ${transcriptome_bam} \\
        rsem_ref \\
        ${sample_id}
    """
} 