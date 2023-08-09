# /Volumes/texas/spinazzi_bgi/parl_bgi

library (ggplot2)
library (ggrepel)
library (DESeq2)

counts <- read.delim ("parl_testis_subread.counts.txt")

annot <- counts[ ,c("Geneid", "gene_name", "gene_type")]

# remove ribosomal RNAs and miRNAs
counts <- counts[grep ("miRNA|rRNA", counts$gene_type, invert=TRUE), ]
counts.s <- counts[ ,grep ("Enrico", colnames (counts))]
row.names (counts.s) <- counts$Geneid

colnames (counts.s) <- gsub ("Aligned.out.bam", "", colnames (counts.s))
colnames (counts.s) <- gsub (".*([W|K])", "\\1", colnames (counts.s))
counts <- counts.s

coldata <- data.frame (matrix (nrow=dim (counts)[2], ncol=0))
coldata$sample <- row.names (coldata) <- colnames (counts)
coldata$condition <- factor (gsub ("[0-9].*", "", coldata$sample))
coldata

stopifnot (colnames (counts) == row.names (coldata))


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# keep <- rowSums(counts(dds)) >= 30*(dim(counts)[2]/2)
keep <- apply (counts (dds), 1, function (x) sum (x >30) > dim(counts)[2]/2 )
dds <- dds[keep,]
dds


dds$condition <- relevel(dds$condition, ref = "WT")

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, name= "condition_KO_vs_WT")

res <- merge (data.frame (res), round (counts(dds, norm=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid") 

res <- res[order (res$padj), ]
colnames (res)[1] <- "gene_id"

table (res$padj < 0.05)

write.table (res, "parl_star_rnaseq_testis_v3.txt", sep="\t", quote=F, row.names= F)
