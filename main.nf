params.reads = "$projectDir/data/*"
params.genome_alignment = "$projectDir/reference/transcriptome.fa"
params.star_index = "$projectDir/reference/STAR/"
params.multiqc = "$projectDir/multiqc"
params.lexogen = false
params.adapter_seq = "$projectDir/reference/polyA.fa.gz"
params.illumina_seq = "$projectDir/reference/truseq.fa.gz"
params.genome_annotation = "$projectDir/reference/annotation.gtf"
params.publish_dir = "$projectDir/quality/"
params.aligning_path = "$projectDir/alignments/"
params.outdir = "$projectDir/differential_expression/"
params.diff_table = "$projectDir/Test/diff_table.csv"

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

    input:
    path transcriptome
    path reference

    output:
    path reference

    script:
    """
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


process printBaseDir {
    echo true

    script:
    """
    echo 'Base directory: ${baseDir}'
    """
}


workflow {
    baseDir '.'
    printBaseDir()
    index_ch = INDEX(params.genome_alignment, params.star_index)
    trial = fastqc(read_ch)
    path = deduplicate(read_ch)
    reads = bbduk(path)
    bam = ALIGN(params.star_index, reads)
    read = indexing(bam.samples_bam)
    deduped = dedup(read)
    multiqc(deduped.collect())
    feature_count_files(deduped.collect(), params.genome_annotation)
}