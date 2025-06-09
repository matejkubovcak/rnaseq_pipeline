#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW } from './modules/local/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/local/fastqc'
include { FASTQC as FASTQC_FILTERED } from './modules/local/fastqc'
include { RSEM_PREPARE_REFERENCE; RSEM_QUANTIFICATION } from './modules/local/rsem'

// Pipeline parameters
params.input_dir = null
params.outdir = "results"
params.tRNA_sortmerna_db= "/home/jovyan/age/analysis/genome/rna_sequences/human_tRNA.fa"
params.rRNA_sortmerna_db= "/home/jovyan/age/analysis/genome/rna_sequences/human_rRNA.fa"
params.genome_fasta = "/home/jovyan/age/analysis/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"    // Reference genome FASTA file
params.gtf = "/home/jovyan/age/analysis/genome/Homo_sapiens.GRCh38.114.gtf"             // GTF annotation file
params.star_index = "/home/jovyan/age/analysis/genome/star_index"    // Pre-built STAR index directory (optional)
params.rsem_index = null    // Pre-built RSEM index directory (optional)
params.sjdbOverhang = 149 // 150 - 1 for 150bp reads(illumina novaseq 6000)
params.help = false


// Function to create FASTQ pairs channel
def createFastqPairsChannel(input_dir) {
    return Channel
        .fromPath("${input_dir}/*/*.fastq.gz", checkIfExists: true)
        .filter { it.name.contains('_R1_') }
        .map { r1_file ->
            def r2_file = file(r1_file.toString().replaceAll('_R1_', '_R2_'))
            if (r2_file.exists()) {
                def sample_name = r1_file.parent.name
                return tuple(sample_name, r1_file, r2_file)
            } else {
                return null
            }
        }
        .filter { it != null }
}

// Process: Generate STAR genome index
process STAR_INDEX {
    publishDir "${params.outdir}/star_index", mode: 'copy'
    conda 'bioconda::star=2.7.11a'

    cpus 32
    memory '32.GB'

    input:
    path(genome_fasta)
    path(gtf)
    
    output:
    path "star_index"
    
    script:
    """
    mkdir -p star_index
    
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang ${params.sjdbOverhang} \\
        --runThreadN ${task.cpus - 1}
    """
}

// Process: Convert GTF to BED12
process GTF_TO_BED12 {
    publishDir "${params.outdir}/annotations", mode: 'copy'
    conda 'bioconda::ucsc-gtftogenepred=447 bioconda::ucsc-genepredtobed=447'
    
    input:
    path(gtf)
    
    output:
    path "genes.bed12"
    
    script:
    """
    # Convert GTF to genePred format
    gtfToGenePred ${gtf} genes.genePred
    
    # Convert genePred to BED12 format
    genePredToBed genes.genePred genes.bed12
    """
}

// Process: Run fastp for trimming and filtering
process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp", mode: 'copy'
    conda 'bioconda::fastp=0.23.4'

    cpus 8
    memory '16.GB'

    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz")
    path "${sample_id}_fastp.html"
    path "${sample_id}_fastp.json"
    
    script:
    """
    fastp \\
        -i ${read1} \\
        -I ${read2} \\
        -o ${sample_id}_R1_trimmed.fastq.gz \\
        -O ${sample_id}_R2_trimmed.fastq.gz \\
        -h ${sample_id}_fastp.html \\
        -j ${sample_id}_fastp.json \\
        --detect_adapter_for_pe \\
        --thread ${task.cpus}
    """
}

// Process: Run SortMeRNA to remove rRNA
process SORTMERNA {
    tag "$sample_id"
    publishDir "${params.outdir}/sortmerna", mode: 'copy'
    conda 'bioconda::sortmerna=4.3.6'


    cpus 16
    memory '32.GB'
    
    input:
    tuple val(sample_id), path(read1_trimmed), path(read2_trimmed)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1_non_rRNA.fastq.gz"), path("${sample_id}_R2_non_rRNA.fastq.gz")
    path "${sample_id}_sortmerna.log"
    

    
    script:
    """
    # Create working directory
    echo "Creating working directory"
    mkdir -p sortmerna_work
    
    echo "Running SortMeRNA"
    # Merge paired reads for SortMeRNA
    sortmerna \\
        --ref ${params.rRNA_sortmerna_db} \\
        --ref ${params.tRNA_sortmerna_db} \\
        --reads ${read1_trimmed} \\
        --reads ${read2_trimmed} \\
        --fastx \\
        --paired_out \\
        --out2 \\
        --other \\
        --workdir sortmerna_work \\
        -a ${task.cpus -1}

    echo "Renaming files"
    # Rename files to standard naming
    mv sortmerna_work/out/other_fwd.fq.gz ${sample_id}_R1_non_rRNA.fastq.gz
    mv sortmerna_work/out/other_rev.fq.gz ${sample_id}_R2_non_rRNA.fastq.gz
    mv sortmerna_work/out/aligned.log ${sample_id}_sortmerna.log
    """
}

// Process: Align reads with STAR
process STAR_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/star", mode: 'copy'
    conda 'bioconda::star=2.7.11a bioconda::samtools=1.19.2'

    cpus 32
    memory '64.GB'
    
    input:
    tuple val(sample_id), path(read1_filtered), path(read2_filtered)
    path(star_index)
    path(gtf)
    
    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), path("${sample_id}.Aligned.sortedByCoord.out.bam.bai"), path("${sample_id}.Aligned.toTranscriptome.out.bam")
    path "${sample_id}.Log.final.out"
    path "${sample_id}.SJ.out.tab"
    
    script:
    """
    STAR \\
        --genomeDir ${star_index} \\
        --readFilesIn ${read1_filtered} ${read2_filtered} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --quantMode TranscriptomeSAM \\
        --outSAMunmapped Within \\
        --outSAMattributes Standard \\
        --sjdbGTFfile ${gtf} \\
        --runThreadN ${task.cpus} \\
        --outBAMsortingThreadN ${task.cpus} \\
        --limitBAMsortRAM 8000000000
    
    # Index the BAM file
    samtools index ${sample_id}.Aligned.sortedByCoord.out.bam
    """
}

// Process: RSeQC bam_stat
process RSEQC_BAMSTAT {
    tag "$sample_id"
    publishDir "${params.outdir}/rseqc", mode: 'copy'
    conda 'bioconda::rseqc=5.0.3'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    path "${sample_id}.bam_stat.txt"
    
    script:
    """
    bam_stat.py -i ${bam} > ${sample_id}.bam_stat.txt
    """
}

// Process: RSeQC read_distribution
process RSEQC_READ_DISTRIBUTION {
    tag "$sample_id"
    publishDir "${params.outdir}/rseqc", mode: 'copy'
    conda 'bioconda::rseqc=5.0.3'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(bed12)
    
    output:
    path "${sample_id}.read_distribution.txt"
    
    script:
    """
    read_distribution.py -i ${bam} -r ${bed12} > ${sample_id}.read_distribution.txt
    """
}

// Process: RSeQC infer_experiment
process RSEQC_INFER_EXPERIMENT {
    tag "$sample_id"
    publishDir "${params.outdir}/rseqc", mode: 'copy'
    conda 'bioconda::rseqc=5.0.3'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(bed12)
    
    output:
    path "${sample_id}.infer_experiment.txt"
    
    script:
    """
    infer_experiment.py -i ${bam} -r ${bed12} > ${sample_id}.infer_experiment.txt
    """
}

// Process: RSeQC inner_distance
process RSEQC_INNER_DISTANCE {
    tag "$sample_id"
    publishDir "${params.outdir}/rseqc", mode: 'copy'
    conda 'bioconda::rseqc=5.0.3'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(bed12)
    
    output:
    path "${sample_id}.inner_distance.txt"
    path "${sample_id}.inner_distance_freq.txt"
    path "${sample_id}.inner_distance_plot.r"
    path "${sample_id}.inner_distance_plot.pdf"
    
    script:
    """
    inner_distance.py -i ${bam} -r ${bed12} -o ${sample_id}
    """
}

// Process: Preseq library complexity
process PRESEQ {
    tag "$sample_id"
    publishDir "${params.outdir}/preseq", mode: 'copy'
    conda 'bioconda::preseq=3.2.0'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    path "${sample_id}_preseq_lc_extrap.txt"
    
    script:
    """
    preseq lc_extrap -B -o ${sample_id}_preseq_lc_extrap.txt ${bam}
    """
}

// Process: Run MultiQC to aggregate all results
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    conda 'bioconda::multiqc=1.19'
    
    input:
    path('*')
    
    output:
    path "multiqc_report.html"
    path "multiqc_report_plots"
    path "multiqc_report_data"
    
    script:
    """
    multiqc . \\
        --title "RNAseq quantification pipeline QC report" \\
        --filename multiqc_report.html \\
        --dirs \\
        --dirs-depth 3 \\
        --export
    """
}

// Main workflow
workflow {

    // Help message
    if (params.help) {
        println """
        ===========================================
        FASTQ QC Pipeline with RSeQC and Preseq
        ===========================================
        Usage:
        nextflow run main.nf --input_dir <path_to_directory_with_subfolders> [options]
        
        Parameters:
        --input_dir       Directory containing subfolders with paired FASTQ files
        --outdir          Output directory (default: results)
        --rRNA_sortmerna_db    SortMeRNA database file (default: /home/jovyan/age/analysis/genome/rna_sequences/human_rRNA.fa)
        --tRNA_sortmerna_db    SortMeRNA database file (default: /home/jovyan/age/analysis/genome/rna_sequences/human_tRNA.fa)
        --genome_fasta    Path to reference genome FASTA file (required if --star_index not provided)
        --gtf             Path to GTF annotation file
        --star_index      Path to pre-built STAR index directory (optional, will generate if not provided)
        --rsem_index      Path to pre-built RSEM index directory (optional, will generate if not provided)
        --sjdbOverhang    Value for STAR --sjdbOverhang parameter (default: 99)
        --help            Show this help message
        
        Each subfolder should contain paired FASTQ files with naming pattern:
        *_R1_*.fastq.gz and *_R2_*.fastq.gz
        
        Note: This pipeline uses conda environments that will be created automatically.
        If --star_index is not provided, the pipeline will generate STAR index and BED12 files from the provided references.
        If --star_index is provided, you still need to provide --gtf for other analyses.
        """
        exit 0
    }

    if ((!params.star_index || !params.rsem_index) && !params.genome_fasta) {
        error "Please provide --genome_fasta when --star_index or --rsem_index are not provided."
    }

    // Create FASTQ pairs channel directly
    fastq_pairs_ch = createFastqPairsChannel(params.input_dir)

    
    // Create reference channels
    gtf_ch = Channel.value(file(params.gtf, checkIfExists: true))
    genome_fasta_ch = Channel.value(file(params.genome_fasta, checkIfExists: true))
    
    // Handle STAR index - use provided or generate new one
    if (params.star_index) {
        // Use provided STAR index
        star_index_ch = Channel.value(file(params.star_index, checkIfExists: true))
        println "Using provided STAR index: ${params.star_index}"
    } else {
        // Generate STAR index from genome fasta
        star_index_ch = STAR_INDEX(genome_fasta_ch, gtf_ch)
        println "Generating STAR index from genome FASTA"
    }
    
    // Handle RSEM index - use provided or generate new one
    if (params.rsem_index) {
        // Use provided RSEM index
        rsem_index_ch = Channel.value(file(params.rsem_index, checkIfExists: true))
        println "Using provided RSEM index: ${params.rsem_index}"
    } else {
        // Generate RSEM index from genome fasta
        rsem_index_ch = RSEM_PREPARE_REFERENCE(genome_fasta_ch, gtf_ch)
        println "Generating RSEM index from genome FASTA"
    }
    
    // Generate BED12 file
    bed12_file = GTF_TO_BED12(gtf_ch)
    
    // Run FastQC on raw reads
    fastqc_raw = FASTQC_RAW(fastq_pairs_ch, 'raw')
    
    // Run fastp trimming
    fastp_results = FASTP(fastq_pairs_ch)
    
    // Run FastQC on fastp trimmed reads
    fastqc_trimmed = FASTQC_TRIMMED(fastp_results[0], 'trimmed')
    
    // Run SortMeRNA to remove rRNA
    sortmerna_results = SORTMERNA(fastp_results[0])
    
    // Run FastQC on rRNA-filtered reads
    fastqc_filtered = FASTQC_FILTERED(sortmerna_results[0], 'filtered')
    
    // Run STAR alignment
    star_results = STAR_ALIGN(sortmerna_results[0], star_index_ch, gtf_ch)
    
    star_for_downstream = star_results[0]
        .map { sample_id, sorted_bam, sorted_bai, transcriptome_bam -> tuple(sample_id, sorted_bam, sorted_bai) }

    rsem_input = star_results[0]
        .map { sample_id, sorted_bam, sorted_bai, transcriptome_bam -> tuple(sample_id, transcriptome_bam) }

    // Run RSeQC analyses
    bamstat_results = RSEQC_BAMSTAT(star_for_downstream)
    read_dist_results = RSEQC_READ_DISTRIBUTION(star_for_downstream, bed12_file)
    infer_exp_results = RSEQC_INFER_EXPERIMENT(star_for_downstream, bed12_file)
    inner_dist_results = RSEQC_INNER_DISTANCE(star_for_downstream, bed12_file)
    
    // Run Preseq
    preseq_results = PRESEQ(star_for_downstream)

    // Run RSEM
    rsem_results = RSEM_QUANTIFICATION(rsem_input, rsem_index_ch)
    
    // Collect all results for MultiQC
    multiqc_input = fastqc_raw.report
        .mix(fastqc_raw.zip)
        .mix(fastp_results[1])
        .mix(fastp_results[2])
        .mix(fastqc_trimmed.report)
        .mix(fastqc_trimmed.zip)
        .mix(fastqc_filtered.report)
        .mix(fastqc_filtered.zip)
        .mix(star_results[1])
        .mix(star_results[2])
        .mix(bamstat_results)
        .mix(read_dist_results)
        .mix(infer_exp_results)
        .mix(inner_dist_results[0])
        .mix(inner_dist_results[1])
        .mix(preseq_results[0])
        .collect()
    
    // Run MultiQC only after all analyses are complete
    MULTIQC(multiqc_input)

}
