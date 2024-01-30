process SORTMERNA {
    tag "$sample_id"
    label "Removing rRNA reads"
    conda "/project/zhuzhuzhang/lyang/software/miniconda3/envs/sortmerna_env"

    input:
    tuple val(sample_id), path(clean_fq1), path(clean_fq2)
    path fastas

    output:
    tuple val(sample_id), path("*non_rNRA.fastq.gz"), emit: reads
    tuple val(sample_id), path("*.log"), emit: log

    script:
    """
    sortmerna \\
            ${'--ref '+fastas.join(' --ref ')} \\
            --reads ${reads[0]} \\
            --reads ${reads[1]} \\
            --threads $task.cpus \\
            --workdir . \\
            --aligned rRNA_reads \\
            --fastx \\
            --other non_rRNA_reads \\
            --paired_in \\
            --out2

    mv non_rRNA_reads_fwd.f*q.gz ${sample_id}_1.non_rRNA.fastq.gz
    mv non_rRNA_reads_rev.f*q.gz ${sample_id}_2.non_rRNA.fastq.gz
    mv rRNA_reads.log ${sample_id}.sortmerna.log
    """

}