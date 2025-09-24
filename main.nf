nextflow.enable.dsl=2
params.input_dir = "/home/khanle/Documents/catfish/rna_seq/PRJNA902092/fastq"

workflow {
    Channel
        .fromFilePairs("${params.input_dir}/*_{1,2}.fastq.gz", flat: true)
        .map { sample_id, r1, r2 -> tuple(sample_id.replaceFirst(/_1$/, ''), [r1, r2]) }
        .view()
        .set { raw_reads }

    trimmed_reads = raw_reads | trimmomatic_trim
    mapped_bams = trimmed_reads | hisat2_map
    mapped_bams | bcftools_call
    mapped_bams | featurecounts
}

process trimmomatic_trim {
    tag "$sample_id"
    publishDir "${params.trim_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1_paired.fastq.gz"), path("${sample_id}_2_paired.fastq.gz")

    script:
    def (r1, r2) = reads
    """
    mkdir -p ${params.trim_dir}
    trimmomatic PE -threads ${task.cpus} \
      $r1 $r2 \
      ${sample_id}_1_paired.fastq.gz ${sample_id}_1_unpaired.fastq.gz \
      ${sample_id}_2_paired.fastq.gz ${sample_id}_2_unpaired.fastq.gz \
      ILLUMINACLIP:${params.trimmomatic_adapters}:2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:70
    """
}

process hisat2_map {
    tag "$sample_id"
    publishDir "${params.map_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    hisat2 -p ${task.cpus} -x ${params.ref_index} \
      -1 $r1 -2 $r2 | \
      samtools view -bS - | \
      samtools sort -@ ${task.cpus} -o ${sample_id}.bam
    """
}

process bcftools_call {
    tag "$sample_id"
    publishDir "${params.vcf_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.csi")

    script:
    """
    # Index BAM nếu chưa có .bai
    if [ ! -f ${bam}.bai ]; then
        samtools index $bam
    fi

    # Gọi biến thể
    bcftools mpileup -Ou -f ${params.ref_genome} --threads ${task.cpus} $bam | \
      bcftools call -mv -Oz -o ${sample_id}.vcf.gz

    # Index VCF
    bcftools index ${sample_id}.vcf.gz
    """
}


process featurecounts {
    tag "$sample_id"
    publishDir "${params.project_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_counts.txt")

    script:
    """
    featureCounts -p -T ${task.cpus} -t exon -g gene_id \
      -a ${params.gtf_file} -o ${sample_id}_counts.txt $bam
    """
}
