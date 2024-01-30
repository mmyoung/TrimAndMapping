# TrimAndMapping
A nextflow workflow for trimming and mapping of pair-end sequencing data

usage: 
```
nextflow run ./ -params-file params.config
```
Get help:

```
nextflow run ./ --help
```

Test:
```
nextflow run ./ -c ./test/params.config
```

Parameters:
```
--sample_sheet sample_sheet    A tab-delimited file storing the samples information, with three columns: sample_id, fq1, fq2. See an example in ./test/sample_sheet.csv.
--prime5_trim_len    Trim the first # bases in 5-end of each read.
--prime3_trim_len    Trim the first # bases in 3-end of each read.
--output_dir OUTDIR   Name for directory for saving the results. Default: results/
--bowtie_idx bowtie_idx    bowtie index prefix for analysis.
--data_dir raw_data    The folder where the raw .fq files are.
--remove_rRNA False  Whether to remove rRNAs from the library. If Ture, provide rRNA sequences. Default: false
--rRNA_fa rRNA sequences.
```

## Workflow Summary
1. Trimming adapter and low-quality bases (trim_galore)
2. Removal of ribosomal RNA (SortMeRNA)
3. Quality control for clean reads (fastQC)
4. Mapping to genome (bowtie2)
5. Aligment filter (samtools, bedtools)
    reads mapping to blacklisted regions (..to be added)
    unmapped & multiply mapped reads (samtools Q>20)
6. Make .bw file from .bam
7. Coverage evaluation

## Results Sturcture
```
├── alignment
│   ├── sample1_sorted.bam
│   ├── sample1_sorted.bam.bai
│   ├── sample2_sorted.bam
│   └── sample2_sorted.bam.bai
├── bw_output
│   ├── sample1_sorted_bam.bw
│   └── sample2_sorted_bam.bw
├── clean_reads_fastQC
│   ├── sample1
│   │   ├── sample1_val_1_fastqc.html
│   │   ├── sample1_val_1_fastqc.zip
│   │   ├── sample1_val_2_fastqc.html
│   │   └── sample1_val_2_fastqc.zip
│   └── sample2
│       ├── sample2_val_1_fastqc.html
│       ├── sample2_val_1_fastqc.zip
│       ├── sample2_val_2_fastqc.html
│       └── sample2_val_2_fastqc.zip
├── coverage_out
│   ├── sample1_base.depth
│   ├── sample1_depth.pdf
│   ├── sample2_base.depth
│   └── sample2_depth.pdf
└── trimm
    ├── sample1_val_1.fq.gz
    ├── sample1_val_2.fq.gz
    ├── sample2_val_1.fq.gz
    └── sample2_val_2.fq.gz
```