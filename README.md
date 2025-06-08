# RNA-seq Quantification Pipeline

This Nextflow pipeline performs a comprehensive RNA-sequencing analysis, from raw FASTQ files to gene and transcript quantification. It includes steps for quality control, adapter trimming, rRNA removal, alignment, post-alignment QC, and expression quantification.

## Features

- **Automatic package installation**: Uses conda environments to install all required tools.
- **Quality Control**: FastQC on raw, trimmed, and rRNA-filtered reads.
- **Read Trimming**: `fastp` for adapter trimming and quality filtering.
- **rRNA Removal**: `SortMeRNA` to filter out ribosomal RNA.
- **Alignment**: `STAR` for splice-aware alignment to a reference genome.
- **Quantification**: `RSEM` for gene and isoform level quantification.
- **Post-Alignment QC**: A suite of `RSeQC` tools and `Preseq` for library complexity analysis.
- **Aggregate Reporting**: `MultiQC` to generate a comprehensive report for all pipeline steps.
- **Reproducibility**: Conda environments and Nextflow provide a fully reproducible analysis environment.
- **Reference Management**: Can generate STAR and RSEM indices on the fly, or use pre-built ones.

## Requirements

- Nextflow (>=21.04.0)
- Conda or Miniconda/Mambaforge
- Internet connection (for initial package downloads)

## Quick Start

To run the pipeline, use the following command:

```bash
nextflow run main.nf --input_dir /path/to/your/fastq/directory --gtf /path/to/genome.gtf
```

**Important:** You must provide an input directory with FASTQ files and a genome annotation file (`--gtf`). If you do not provide pre-built genome indices (`--star_index`, `--rsem_index`), you must also provide a genome FASTA file (`--genome_fasta`).

## Input Data Structure

Your input directory should contain subfolders, one for each sample, with paired-end FASTQ files inside.

```
input_directory/
├── sample1/
│   ├── sample1_R1_001.fastq.gz
│   └── sample1_R2_001.fastq.gz
├── sample2/
│   ├── sample2_R1_001.fastq.gz
│   └── sample2_R2_001.fastq.gz
└── sample3/
    ├── sample3_R1_001.fastq.gz
    └── sample3_R2_001.fastq.gz
```

**Supported naming patterns:** `*_R1_*.fastq.gz` and `*_R2_*.fastq.gz`.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_dir` | *required* | Directory containing subfolders with paired FASTQ files. |
| `--outdir` | `results` | Output directory for all results. |
| `--gtf` | *required* | Path to the genome annotation file in GTF format. |
| `--genome_fasta`| `null` | Path to the reference genome FASTA file. Required if `--star_index` or `--rsem_index` are not provided. |
| `--star_index` | `null` | Path to a pre-built STAR index directory. If not provided, it will be generated. |
| `--rsem_index` | `null` | Path to a pre-built RSEM index directory. If not provided, it will be generated. |
| `--sortmerna_db`| `/path/in/main.nf` | Path to the SortMeRNA database. |
| `--sjdbOverhang`| `149` | Value for STAR's `--sjdbOverhang` parameter. Should be `read_length - 1`. |
| `--help` | `false` | Show the pipeline's help message. |

## Output Structure

The `outdir` will have the following structure:

```
results/
├── star_index/          # Generated STAR index (if not provided)
├── annotations/         # BED12 annotation file generated from GTF
├── rsem_index/          # Generated RSEM index (if not provided)
├── fastqc_raw/          # FastQC reports for raw reads
├── fastp/               # fastp trimming results and reports
├── fastqc_trimmed/      # FastQC reports for trimmed reads
├── sortmerna/           # SortMeRNA filtered reads and logs
├── fastqc_filtered/     # FastQC reports for rRNA-filtered reads
├── star/                # STAR alignment BAM files and logs
├── rseqc/               # RSeQC quality control reports
├── preseq/              # Preseq library complexity reports
├── rsem/                # RSEM quantification results
├── multiqc/             # MultiQC aggregated report
├── timeline.html        # Pipeline execution timeline
├── report.html          # Nextflow execution report
├── trace.txt            # Process execution trace
└── dag.html             # Pipeline DAG visualization
```

## Pipeline Workflow

1.  **Reference Indexing (Optional)**: If not provided, `STAR_INDEX` and `RSEM_PREPARE_REFERENCE` processes will build indices from the provided `--genome_fasta` and `--gtf`. `GTF_TO_BED12` converts the GTF for RSeQC.
2.  **Raw Read QC**: `FASTQC_RAW` runs on the input FASTQ files.
3.  **Trimming**: `FASTP` performs adapter and quality trimming.
4.  **Trimmed QC**: `FASTQC_TRIMMED` runs on the reads post-`fastp`.
5.  **rRNA Removal**: `SORTMERNA` filters out ribosomal RNA.
6.  **Filtered QC**: `FASTQC_FILTERED` runs on the reads post-`SortMeRNA`.
7.  **Alignment**: `STAR_ALIGN` aligns the clean reads to the genome.
8.  **Post-Alignment QC**: A series of `RSEQC` processes (`BAM_STAT`, `READ_DISTRIBUTION`, `INFER_EXPERIMENT`, `INNER_DISTANCE`) and `PRESEQ` are run on the aligned BAM files.
9.  **Quantification**: `RSEM_QUANTIFICATION` calculates gene and transcript expression levels.
10. **Reporting**: `MULTIQC` aggregates all logs and reports into a single HTML file.

## Execution Profiles

The `nextflow.config` file defines several execution profiles.

### Local (Standard)
This is the default profile.
```bash
nextflow run main.nf --input_dir /path/to/data --gtf /path/to/anno.gtf
```

### Docker / Singularity
The pipeline has `docker` and `singularity` profiles, but they are currently incomplete and only specify containers for a few processes. They can be used as a template for a fully containerized execution.
```bash
nextflow run main.nf --input_dir /path/to/data --gtf /path/to/anno.gtf -profile docker
nextflow run main.nf --input_dir /path/to/data --gtf /path/to/anno.gtf -profile singularity
```

### Cluster (SLURM)
A basic `cluster` profile for SLURM is provided. You will need to customize `clusterOptions` in `nextflow.config` for your specific cluster environment.
```bash
nextflow run main.nf --input_dir /path/to/data --gtf /path/to/anno.gtf -profile cluster
```

## Resource Configuration

You can adjust resource allocation in `nextflow.config`:

```groovy
process {
    withName: 'FASTP' {
        cpus = 8
        memory = '16.GB'
        time = '4.h'
    }
}
```

## Example Run

For the provided SRA data:

```bash
# Run the pipeline on the age RNA-seq data
nextflow run main.nf --input_dir sra/age_rnaseq --outdir age_qc_results

# With Docker for reproducibility
nextflow run main.nf \
    --input_dir sra/age_rnaseq \
    --outdir age_qc_results \
    -profile docker
```

## Troubleshooting

### Common Issues

1. **No FASTQ pairs found**
   - Check that your files follow the naming convention `*_R1_*.fastq.gz` and `*_R2_*.fastq.gz`
   - Ensure files are in subfolders, not the root directory

2. **Out of memory errors**
   - Reduce the number of parallel processes in `nextflow.config`
   - Increase memory allocation for specific processes

3. **Permission errors with Docker**
   - The pipeline includes user mapping (`-u $(id -u):$(id -g)`) to handle permissions

### Getting Help

```bash
# Show help message
nextflow run main.nf --help

# Check pipeline version
nextflow run main.nf --version
```

## Citation

If you use this pipeline in your research, please cite the tools it uses:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
- **FastQC**: Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.
- **fastp**: Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890.
- **SortMeRNA**: Kopylova, E., Noé, L., & Touzet, H. (2012). SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics, 28(24), 3211-3217.
- **STAR**: Dobin, A., et al. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15-21.
- **RSEM**: Li, B., & Dewey, C. N. (2011). RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC bioinformatics, 12(1), 1-16.
- **RSeQC**: Wang, L., et al. (2012). RSeQC: quality control of RNA-seq experiments. Bioinformatics, 28(16), 2184-2185.
- **Preseq**: Daley, T., & Smith, A. D. (2013). Predicting the molecular complexity of sequencing libraries. Nature methods, 10(4), 325-327.
- **MultiQC**: Ewels, P., et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048.

## License

This pipeline is released under the MIT License. 