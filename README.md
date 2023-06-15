# LexogeneUmiPipelineNextflow

This pipelines is specifically developed for Quantseq lexogene data with UMI deduplications
To use this pipeline:

```
nextflow run ../nextflow_lexogen/main.nf --project_dir ---add_your_project_dir---
```


## Prerequisites
Important within your ProjectDir you need a folder /data containing the reads in fastq.gz format. In addition you need a reference folder holding the transcriptome.fa as well
as the annotation.gtf from you genome that you would like to align against.


## What does this pipeline do?

In general it goes through the basic steps of an alignment process, by first running FastQC as well as the index creation of the Star genome. Followed by deduplication of the UMI counts
to prevent PCR amplification issues. Finally the data will be aligned against the used reference genome with Star and after this quality control will be performed by means of MultiQC.
Finally counting procedure will be performed using FeatureCounts using the R Language which will create a Folder called differential_expression holding and RDS file with the counts that you 
can use than for further downstream analysis such as Differential Gene Expression Analysis.

