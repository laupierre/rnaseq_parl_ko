library (ggplot2)
library (ggrepel)
library (DESeq2)


#### brain

counts <- read.delim ("parl_star_raw_data_annot_v3.txt", row.names=1)
# remove ribosomal RNAs and miRNAs
counts <- counts[grep ("miRNA|rRNA", counts$gene_type, invert=TRUE), ]

annot <- counts[ ,grep ("Parl", colnames (counts), invert=TRUE)]
counts <- counts[ ,grep ("Parl", colnames (counts))]



## brain tissue
counts <- counts[ ,grep ("_C", colnames (counts))]

coldata <- data.frame (matrix (nrow=dim (counts)[2], ncol=0))
coldata$sample <- row.names (coldata) <- colnames (counts)
coldata$sample <- gsub ("Parl_", "", coldata$sample)
coldata$condition <- factor (gsub ("_.*", "", coldata$sample))
coldata

stopifnot (colnames (counts) == row.names (coldata))

counts1 <- counts
coldata1 <- coldata

rm (counts)
rm (coldata)



counts <- read.delim ("ndufs4_star_raw_data_annot_v3.txt", row.names=1)
# remove ribosomal RNAs and miRNAs
counts <- counts[grep ("miRNA|rRNA", counts$gene_type, invert=TRUE), ]

annot <- counts[ ,grep ("NDU", colnames (counts), invert=TRUE)]
counts <- counts[ ,grep ("NDU", colnames (counts))]


## brain tissue
counts <- counts[ ,grep ("\\.C", colnames (counts))]

coldata <- data.frame (matrix (nrow=dim (counts)[2], ncol=0))
coldata$sample <- row.names (coldata) <- colnames (counts)
coldata$sample <- gsub ("NDU.*_", "", coldata$sample)
coldata$condition <- factor (gsub ("\\..*", "", coldata$sample))
coldata

stopifnot (colnames (counts) == row.names (coldata))


coldata <- rbind (coldata1, coldata)

coldata$batch <- factor (c(rep (1, 12), rep (2,12)))

levels (coldata$condition) <- c(levels (coldata$condition), "Parl_KO", "Ndufs4_KO", "Parl_WT", "Ndufs4_WT")

coldata$condition [grep ("Parl.*KO", row.names (coldata))] <- "Parl_KO"
coldata$condition [grep ("Parl.*WT", row.names (coldata))] <- "Parl_WT"
coldata$condition [grep ("NDU.*KO", row.names (coldata))] <- "Ndufs4_KO"
coldata$condition [grep ("NDU.*WT", row.names (coldata))] <- "Ndufs4_WT"
coldata

counts <- cbind (counts1, counts)

stopifnot (colnames (counts) == row.names (coldata))


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# keep <- rowSums(counts(dds)) >= 30*(dim(counts)[2]/2)
keep <- apply (counts (dds), 1, function (x) sum (x >30) > dim(counts)[2]/2 )
dds <- dds[keep,]
dds

dds$condition <- relevel(dds$condition, ref = "Parl_WT")

dds <- DESeq(dds)
resultsNames(dds)

res1 <- data.frame (results(dds, name= "condition_Parl_KO_vs_Parl_WT"))



dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# keep <- rowSums(counts(dds)) >= 30*(dim(counts)[2]/2)
keep <- apply (counts (dds), 1, function (x) sum (x >30) > dim(counts)[2]/2 )
dds <- dds[keep,]
dds

dds$condition <- relevel(dds$condition, ref = "Ndufs4_WT")

dds <- DESeq(dds)
resultsNames(dds)

res2 <- data.frame (results(dds, name= "condition_Ndufs4_KO_vs_Ndufs4_WT"))

res <- merge (res1, res2,by= "row.names")
colnames (res) <- gsub ("\\.x", ".Parl", colnames (res))
colnames (res) <- gsub ("\\.y", ".Ndufs4", colnames (res))
head (res)

res <- merge (res, annot, by.x="Row.names", by.y="row.names")
colnames (res)[1] <- "Ensembl_ID"
head (res)

res <- res[order (res$padj.Parl), ]
head (res)

write.table (res, "brain_deg_of_the_two_genotypes_for_heatmap.txt", sep="\t", quote=F, row.names=F)




###########
#### muscle

library (ggplot2)
library (ggrepel)
library (DESeq2)

counts <- read.delim ("parl_star_raw_data_annot_v3.txt", row.names=1)
# remove ribosomal RNAs and miRNAs
counts <- counts[grep ("miRNA|rRNA", counts$gene_type, invert=TRUE), ]

annot <- counts[ ,grep ("Parl", colnames (counts), invert=TRUE)]
counts <- counts[ ,grep ("Parl", colnames (counts))]



## muscle tissue
counts <- counts[ ,grep ("_M", colnames (counts))]

coldata <- data.frame (matrix (nrow=dim (counts)[2], ncol=0))
coldata$sample <- row.names (coldata) <- colnames (counts)
coldata$sample <- gsub ("Parl_", "", coldata$sample)
coldata$condition <- factor (gsub ("_.*", "", coldata$sample))
coldata

stopifnot (colnames (counts) == row.names (coldata))

counts1 <- counts
coldata1 <- coldata

rm (counts)
rm (coldata)



counts <- read.delim ("ndufs4_star_raw_data_annot_v3.txt", row.names=1)
# remove ribosomal RNAs and miRNAs
counts <- counts[grep ("miRNA|rRNA", counts$gene_type, invert=TRUE), ]

annot <- counts[ ,grep ("NDU", colnames (counts), invert=TRUE)]
counts <- counts[ ,grep ("NDU", colnames (counts))]


## muscle tissue
counts <- counts[ ,grep ("\\.M", colnames (counts))]

coldata <- data.frame (matrix (nrow=dim (counts)[2], ncol=0))
coldata$sample <- row.names (coldata) <- colnames (counts)
coldata$sample <- gsub ("NDU.*_", "", coldata$sample)
coldata$condition <- factor (gsub ("\\..*", "", coldata$sample))
coldata

stopifnot (colnames (counts) == row.names (coldata))


coldata <- rbind (coldata1, coldata)

coldata$batch <- factor (c(rep (1, 12), rep (2,12)))

levels (coldata$condition) <- c(levels (coldata$condition), "Parl_KO", "Ndufs4_KO", "Parl_WT", "Ndufs4_WT")

coldata$condition [grep ("Parl.*KO", row.names (coldata))] <- "Parl_KO"
coldata$condition [grep ("Parl.*WT", row.names (coldata))] <- "Parl_WT"
coldata$condition [grep ("NDU.*KO", row.names (coldata))] <- "Ndufs4_KO"
coldata$condition [grep ("NDU.*WT", row.names (coldata))] <- "Ndufs4_WT"
coldata

counts <- cbind (counts1, counts)

stopifnot (colnames (counts) == row.names (coldata))


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# keep <- rowSums(counts(dds)) >= 30*(dim(counts)[2]/2)
keep <- apply (counts (dds), 1, function (x) sum (x >30) > dim(counts)[2]/2 )
dds <- dds[keep,]
dds

dds$condition <- relevel(dds$condition, ref = "Parl_WT")

dds <- DESeq(dds)
resultsNames(dds)

res1 <- data.frame (results(dds, name= "condition_Parl_KO_vs_Parl_WT"))



dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# keep <- rowSums(counts(dds)) >= 30*(dim(counts)[2]/2)
keep <- apply (counts (dds), 1, function (x) sum (x >30) > dim(counts)[2]/2 )
dds <- dds[keep,]
dds

dds$condition <- relevel(dds$condition, ref = "Ndufs4_WT")

dds <- DESeq(dds)
resultsNames(dds)

res2 <- data.frame (results(dds, name= "condition_Ndufs4_KO_vs_Ndufs4_WT"))

res <- merge (res1, res2,by= "row.names")
colnames (res) <- gsub ("\\.x", ".Parl", colnames (res))
colnames (res) <- gsub ("\\.y", ".Ndufs4", colnames (res))
head (res)

res <- merge (res, annot, by.x="Row.names", by.y="row.names")
colnames (res)[1] <- "Ensembl_ID"
head (res)

res <- res[order (res$padj.Parl), ]
head (res)

write.table (res, "muscle_deg_of_the_two_genotypes_for_heatmap.txt", sep="\t", quote=F, row.names=F)


#### heatmap

library (pheatmap)

a <- read.delim ("muscle_deg_of_the_two_genotypes_for_heatmap.txt")
head (a)
               













