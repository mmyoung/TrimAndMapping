nextflow.enable.dsl=2

params.help = false
params.threads = 1
params.output_dir = './results'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to preprocess the genomic DNA library and QC'
    log.info '--------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow Preprocess.nf --sample_sheet sample_sheet --bowtie_idx bowtie_idx  --output_dir results --data_dir raw_data'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--sample_sheet sample_sheet    A tab-delimited file storing the samples information, with three columns: sample_id, fq1, fq2.'
    log.info '	--prime5_trim_len    Trim the first ? number of bases in 5-end of each read.'
    log.info '	--prime3_trim_len    Trim the first ? number of bases in 3-end of each read.'
    log.info '  --output_dir OUTDIR   Name for directory for saving the results. Default: results/'
    log.info '  --bowtie_idx bowtie_idx    bowtie index prefix for analysis.'
    log.info '  --data_dir raw_data    The folder where the raw .fq files are.'
    log.info '  --remove_rRNA False  Whether to remove rRNAs from the library. If Ture, provide rRNA sequences.'
    log.info '  --rRNA_fa rRNA sequences.' 
    exit 1
}

process fastq_trim {

tag "Adapter and low-quality based trimming"

publishDir "${params.output_dir}/trimm", mode: 'copy'

input:
    tuple val(sample_id), path(fq1,stageAs:"*_1.fastq"), path(fq2,stageAs:"*_2.fastq")

output:
    tuple val(sample_id), path("${sample_id}_val_1.fq.gz"), path("${sample_id}_val_2.fq.gz")

script:

"""
    echo "Trimming ${sample_id}"
    /project/zhuzhuzhang/lyang/software/TrimGalore-0.6.10/trim_galore --gzip -o ./ \
    --path_to_cutadapt /project/zhuzhuzhang/lyang/software/miniconda3/envs/cutadaptenv/bin/cutadapt \
    --basename ${sample_id} \
    --fastqc \
    --clip_R1 ${params.prime5_trim_len} \
    --three_prime_clip_R1 ${params.prime3_trim_len} \
    -paired \
    --clip_R2 ${params.prime5_trim_len} \
    --three_prime_clip_R2 ${params.prime3_trim_len} \
    ${fq1} \
    ${fq2}

"""
}

process fastq_QC {
tag "QC for clean .fastq ..."

publishDir "${params.output_dir}/clean_reads_fastQC", mode: 'copy'

input:
    tuple val(sample_id), path(clean_fq1), path(clean_fq2)

output:
    path("${sample_id}/*")

script:

"""
    mkdir ${sample_id}
    /project/zhuzhuzhang/lyang/software/FastQC/fastqc -o ${sample_id} -t 2 ${clean_fq1} ${clean_fq2}

"""
}

process fastq_map {
    
tag "Mapping clean reads to bowtie index"

publishDir "${params.output_dir}/alignment", mode: 'copy'

input:
    tuple val(sample_id), path(clean_fq1), path(clean_fq2)

output:
    tuple val(sample_id), path("${sample_id}_map.bam")

script:

"""
    module load samtools/1.13
    
    /project/zhuzhuzhang/lyang/software/bowtie2-2.4.2-sra-linux-x86_64/bowtie2 \
    -x ${params.bowtie_idx} \
    -1 ${clean_fq1} \
    -2 ${clean_fq2} \
    | samtools view -bS -o ${sample_id}_map.bam

"""
}


process bam_2_bw {

tag "Making .bw files from .bam files"

publishDir "${params.output_dir}/bw_output", mode: 'copy'

input:
    tuple val(sample_id), path(bam_file), path(bai)

output:
    path("${sample_id}_sorted_bam.bw")

script: 
"""

    bamCoverage -b ${bam_file} -o ${sample_id}_sorted_bam.bw

"""
}

process coverage_cal {

tag "Calculate the reads coverage in each base ... "

publishDir "${params.output_dir}/coverage_out", mode: 'copy'

input:
    tuple val(sample_id), path(bam_file), path(bai)

output:
    tuple path("${sample_id}_base.depth"), path("${sample_id}_depth.pdf")

script: 

"""

    module load samtools/1.13
    samtools depth -a ${bam_file} |cut -f 3|sort|uniq -c > ${sample_id}_base.depth
    Rscript ${workflow.projectDir}/script/line_plot.r ${sample_id}_base.depth ${sample_id} ${sample_id}_depth.pdf

"""
}

sample_sheet = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                      .ifEmpty { exit 1, "sample sheet not found" }
                      .splitCsv(header:true, sep: ',')

include {SORTMERNA} from './modules/sortmerna'
include {MARK_DUPLICATES} from './module/markduplicates'

workflow {

    aln_in = sample_sheet.map { row -> row.fq1 = params.data_dir + "/" + row.fq1; row }
                .map { row -> row.fq2 = params.data_dir + "/" + row.fq2; row }
                .map { row -> [row.sample, file(row.fq1), file(row.fq2)] }
    trim_fq = aln_in | fastq_trim
    if(params.remove_rRNA) {
        trim_fq = SORTMERNA(trim_fq.concat(params.rRNA_fa))
    }
    trim_fq |
    fastq_map |
    MARK_DUPLICATES |
    (bam_2_bw & coverage_cal)
    fastq_QC(fastq_trim.out)
}