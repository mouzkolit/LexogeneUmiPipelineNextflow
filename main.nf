params.project_dir = "/home/data-science/Desktop/Lexogene_Analysis/"
params.reads = "$params.project_dir/data/*"
params.star = "$params.project_dir/Star"
params.genome_alignment = "$params.project_dir/reference/transcriptome.fa"
params.multiqc = "$params.project_dir/multiqc"
params.lexogen = false
params.adapter_seq = "$projectDir/reference/polyA.fa.gz"
params.illumina_seq = "$projectDir/reference/truseq.fa.gz"
params.genome_annotation = "$params.project_dir/reference/annotation.gtf"
params.publish_dir = "$params.project_dir/quality/"
params.aligning_path = "$params.project_dir/alignments/"
params.outdir = "$params.project_dir/differential_expression/"
params.diff_table = "$params.project_dir/Test/diff_table.csv"

// Define your input channel
Channel
    .fromPath(params.reads)
    .set { read_ch }


process fastqc {
    
    container {params.containers.fastqc}
    publishDir params.publish_dir, mode: 'copy'

    input:
    path(reads)

    output:
    file("${reads}_fastqc/*")

    script:
    """
    mkdir ${reads}_fastqc
    fastqc -o ${reads}_fastqc $reads
    """

}


process multiqc {

    publishDir params.publish_dir, mode: 'copy'

    input:
        file('*') 

    output:
        file('multiqc_report.html')

    """
    multiqc $projectDir
    """

}

process INDEX {

    container { params.containers.star }
    cpus 8
    publishDir params.star, mode: 'copy'

    input:
    path transcriptome
    path reference
    
    output:
    path reference

    script:
    """
    mkdir -p $params.star
    STAR --runThreadN $task.cpus --runMode genomeGenerate --genomeDir $reference --genomeFastaFiles $transcriptome 
    """
}

process deduplicate {

    container {params.containers.umi_tools}
    cpus 8
    input:
    path(reads)

    output:
    path "uma_extract_${reads}"


    script:
    """
    umi_tools extract --extract-method=regex --bc-pattern "(?P<umi_1>.{6})(?P<discard_1>.{4}).*"  -I  $reads --stdout  uma_extract_${reads}
    """
}

process bbduk {

    conda 'bbmap'
    input:
    path(reads)

    output:
    path "bbduk_${reads}.fastq"
    

    script:
    """
    bbduk.sh in=$reads out=bbduk_${reads}.fastq ref=$params.adapter_seq,$params.illumina_seq k=13 ktrim=r useshortkmers=t mink=5 qtrim=t trimq=10 minlength=20 
    """
}

process ALIGN {

    container { params.containers.star }
    cpus 8
    publishDir params.aligning_path, mode: 'copy'

    input:
    path transcriptome
    path reads

    output:
    path '*.sortedByCoord.out.bam', emit: samples_bam
    path '*final*', emit: quality

    maxForks 1 
    script:
    """
    STAR --genomeDir $transcriptome\
	     --limitBAMsortRAM 10000000000 \
	     --readFilesIn $reads \
	     --runThreadN $task.cpus \
	     --outFilterMismatchNoverLmax 0.6 \
	     --outFilterMultimapNmax 200 \
	     --alignIntronMax 1000000 \
	     --outFilterType BySJout \
	     --limitOutSJcollapsed 5000000 \
	     --alignMatesGapMax 1000000 \
         --alignSJDBoverhangMin 1 \
         --alignIntronMin 20 --outSAMtype BAM SortedByCoordinate \
         --chimOutType WithinBAM \
         --outFileNamePrefix aligned_${reads} 
    """

}

process indexing {

    conda "samtools"

    input: 
    path read_bam
    
    output:
    path "${read_bam}"
    path "${read_bam}.bai"

    script:
    """
    samtools index $read_bam
    """

}

process dedup {

    container { params.containers.umi_tools }
    input: 
    path read_bam
    path index    

    output:
    path "dedup_${read_bam}"

    script:
    """
    umi_tools dedup -I $read_bam -S dedup_${read_bam}
    """
}

process printBaseDir {
    echo true

    script:
    """
    echo 'Base directory: ${params.project_dir}'
    """
}


process feature_count_files {

    conda "bioconductor-DESeq2 bioconductor-Rsubread"
    publishDir params.outdir, mode: 'copy'

    input: 
    path reads 
    path annotation

    output:
    path "*.rds"

    script:
    """
    Rscript $baseDir/DESeq2.R $reads $annotation
    """
}


workflow {
    printBaseDir()
    index_ch = INDEX(params.genome_alignment, params.star)
    trial = fastqc(read_ch)
    path = deduplicate(read_ch)
    reads = bbduk(path)
    bam = ALIGN(index_ch, reads)
    read = indexing(bam.samples_bam)
    deduped = dedup(read)
    multiqc(deduped.collect())
    feature_count_files(deduped.collect(), params.genome_annotation)
}