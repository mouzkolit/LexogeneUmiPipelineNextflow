library(Rsubread)
library("DESeq2")

# test if there is at least one argument: if not, return an error


file_path <- commandArgs(trailingOnly = TRUE)
files = file_path[-length(file_path)]
annotation = tail(file_path, 1)

reads<-featureCounts(files = files,
                     annot.ext = annotation,
                     isGTFAnnotationFile = T,
                     countMultiMappingReads = T,
                     isPairedEnd = F,
                     nthreads = 8,
                     GTF.attrType = "gene_name",
                     GTF.featureType = "gene",
                     strandSpecific = 1,
                     fraction = T
)



saveRDS(reads, file = "counted_data.rds")