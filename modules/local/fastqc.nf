nextflow.enable.dsl = 2

process FASTQC {
    tag "$sample_id - $stage"
    publishDir "${params.outdir}/fastqc_${stage}", mode: 'copy'
    conda 'bioconda::fastqc=0.12.1'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    val stage
    
    output:
    path("*.html", emit: report)
    path("*.zip", emit: zip)
    
    script:
    """
    fastqc -t ${task.cpus} ${read1} ${read2}
    """
} 