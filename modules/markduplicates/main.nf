process MARK_DUPLICATES {
    tag "$sample_id"
    label 'Mark_duplication'
    publishDir "${params.output_dir}/deduplicateion_out", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}_dedup_sorted.bam"), path("${sample_id}_dedup_sorted.bam.bai"), emit: txt

    script:

    """
    java -jar /project/gzy8899/softwares/picard.jar \\
        MarkDuplicates \\
        --INPUT $bam_file \\
        --OUTPUT ${sample_id}_dedup.bam \\
        --METRICS_FILE ${sample_id}.MarkDuplicates.metrics.txt \\
        --REMOVE_DUPLICATES
        
    samtools view -bq 20  ${sample_id}_dedup.bam | samtools sort - > ${sample_id}_sorted.bam

    samtools index ${sample_id}_dedup_sorted.bam
    
    """
     

}