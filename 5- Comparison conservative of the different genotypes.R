library (ggplot2)
library (ggrepel)
library (DESeq2)


#### brain

counts <- read.delim ("../parl_star_results/parl_star_raw_data_annot_v3.txt", row.names=1)
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

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "WT")

dds <- DESeq(dds)
resultsNames(dds)

res1 <- data.frame (results(dds, name= "condition_KO_vs_WT"))



counts <- read.delim ("../ndufs4_star_results/ndufs4_star_raw_data_annot_v3.txt", row.names=1)
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


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "WT")

dds <- DESeq(dds)
resultsNames(dds)

res2 <- data.frame (results(dds, name= "condition_KO_vs_WT"))


res <- merge (res1, res2, by="row.names")
 
colnames (res) <- gsub ("\\.x", ".Parl", colnames (res))
colnames (res) <- gsub ("\\.y", ".Ndufs4", colnames (res))

res <- merge (res, annot, by.x="Row.names", by.y="row.names")
colnames (res)[1] <- "Ensembl_ID"

res <- res[order (res$padj.Parl), ]
head (res)

write.table (res, "brain_deg_of_the_two_genotypes_for_heatmap.txt", sep="\t", quote=F, row.names=F)



#### heatmap
# see: https://www.reneshbedre.com/blog/heatmap-with-pheatmap-package-r.html?utm_content=cmp-true

library (pheatmap)

a <- read.delim ("brain_deg_of_the_two_genotypes_for_heatmap.txt", row.names=1)


a <- a[a$padj.Parl < 0.05 | a$padj.Ndufs4 < 0.05, ] 
a <- a[grep ("ENS", row.names (a)), ]
row.names (a) <- a$external_gene_name
a <- unique (a[ ,grep ("log2", colnames (a))])
colnames (a) <- gsub ("log2FoldChange.", "", colnames (a))
a[is.na (a)] <- 0
head (a)

test <- a
paletteLength <- 100

library (RColorBrewer)
               
myColor <- colorRampPalette(rev(brewer.pal(n = 3, name = "RdYlBu"))) (paletteLength)
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
          
pdf ("heatmap brain.pdf")
pheatmap(a, angle_col = 0, fontsize_row=7, color=myColor, breaks=myBreaks)
dev.off ()





###########
###########
#### muscle

library (ggplot2)
library (ggrepel)
library (DESeq2)


counts <- read.delim ("../parl_star_results/parl_star_raw_data_annot_v3.txt", row.names=1)
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

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "WT")

dds <- DESeq(dds)
resultsNames(dds)

res1 <- data.frame (results(dds, name= "condition_KO_vs_WT"))



counts <- read.delim ("../ndufs4_star_results/ndufs4_star_raw_data_annot_v3.txt", row.names=1)
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


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "WT")

dds <- DESeq(dds)
resultsNames(dds)

res2 <- data.frame (results(dds, name= "condition_KO_vs_WT"))
head (res2)


res <- merge (res1, res2, by="row.names")
 
colnames (res) <- gsub ("\\.x", ".Parl", colnames (res))
colnames (res) <- gsub ("\\.y", ".Ndufs4", colnames (res))

res <- merge (res, annot, by.x="Row.names", by.y="row.names")
colnames (res)[1] <- "Ensembl_ID"

res <- res[order (res$padj.Parl), ]
head (res)

write.table (res, "muscle_deg_of_the_two_genotypes_for_heatmap.txt", sep="\t", quote=F, row.names=F)



#### heatmap
# see: https://www.reneshbedre.com/blog/heatmap-with-pheatmap-package-r.html?utm_content=cmp-true

library (pheatmap)

a <- read.delim ("muscle_deg_of_the_two_genotypes_for_heatmap.txt", row.names=1)


a <- a[a$padj.Parl < 0.05 | a$padj.Ndufs4 < 0.05, ] 
a <- a[grep ("ENS", row.names (a)), ]
row.names (a) <- a$external_gene_name
a <- unique (a[ ,grep ("log2", colnames (a))])
colnames (a) <- gsub ("log2FoldChange.", "", colnames (a))
a[is.na (a)] <- 0
head (a)

test <- a
paletteLength <- 100

library (RColorBrewer)
               
myColor <- colorRampPalette(rev(brewer.pal(n = 3, name = "RdYlBu"))) (paletteLength)
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
          
pdf ("heatmap muscle.pdf")
pheatmap(a, angle_col = 0, fontsize_row=7, color=myColor, breaks=myBreaks)
dev.off ()




