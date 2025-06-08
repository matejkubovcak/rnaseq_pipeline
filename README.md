# FASTQ QC Pipeline with Conda

This Nextflow pipeline performs quality control analysis on paired FASTQ files using conda environments for automatic package management.

## Features

- **Automatic package installation**: Uses conda environments to install required tools
- **FastQC analysis**: Quality control on raw and trimmed reads
- **Read trimming**: Uses fastp for adapter trimming and quality filtering
- **MultiQC reporting**: Aggregated quality control report
- **Conda integration**: Each process uses its own conda environment

## Requirements

- Nextflow (>=21.04.0)
- Conda or Miniconda/Mambaforge
- Internet connection (for initial package downloads)

## Quick Start

1. **Run the pipeline**:
```bash
nextflow run main.nf --input_dir /path/to/your/fastq/directory
```

2. **With custom output directory**:
```bash
nextflow run main.nf --input_dir /path/to/your/fastq/directory --outdir my_results
```

## Input Data Structure

Your input directory should contain subfolders, each with paired FASTQ files:

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

**Supported naming patterns:**
- `*_R1_*.fastq.gz` and `*_R2_*.fastq.gz`
- `*_1.fastq.gz` and `*_2.fastq.gz`
- Any pattern where R1/R2 or 1/2 distinguish the pairs

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_dir` | *required* | Directory containing subfolders with paired FASTQ files |
| `--outdir` | `results` | Output directory for results |
| `--help` | `false` | Show help message |

## Output Structure

```
results/
├── fastqc/                    # FastQC reports for raw reads
│   ├── sample1_R1_001_fastqc.html
│   ├── sample1_R1_001_fastqc.zip
│   └── ...
├── fastp/                     # fastp results and reports
│   ├── sample1_fastp.html
│   ├── sample1_fastp.json
│   └── ...
├── fastqc_trimmed/           # FastQC reports for trimmed reads
│   ├── sample1_R1_trimmed_fastqc.html
│   └── ...
├── multiqc_report.html       # Aggregated QC report
├── multiqc_data/            # MultiQC data directory
├── timeline.html            # Pipeline execution timeline
├── report.html              # Nextflow execution report
├── trace.txt               # Process execution trace
└── dag.html                # Pipeline DAG visualization
```

## Execution Profiles

### Standard (Local)
```bash
nextflow run main.nf --input_dir /path/to/data
```

### Docker
```bash
nextflow run main.nf --input_dir /path/to/data -profile docker
```

### Singularity
```bash
nextflow run main.nf --input_dir /path/to/data -profile singularity
```

### Cluster (SLURM)
```bash
nextflow run main.nf --input_dir /path/to/data -profile cluster
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

## Pipeline Workflow

1. **Discovery**: Scan input directory for paired FASTQ files
2. **Raw QC**: Run FastQC on original reads
3. **Trimming**: Use fastp for adapter trimming and quality filtering
4. **Trimmed QC**: Run FastQC on processed reads
5. **Aggregation**: Generate comprehensive MultiQC report

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

If you use this pipeline in your research, please cite:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
- **FastQC**: Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.
- **fastp**: Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890.
- **MultiQC**: Ewels, P., et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048.

## License

This pipeline is released under the MIT License. 