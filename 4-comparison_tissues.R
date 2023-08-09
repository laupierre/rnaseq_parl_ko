testis <- read.delim ("parl_star_rnaseq_testis_v3.txt")
testis <- testis [ ,c("gene_id", "baseMean", "log2FoldChange", "padj")]
colnames (testis)[2:4] <- paste (colnames (testis)[2:4], "_testis", sep="")

brain <- read.delim ("parl_star_rnaseq_brain_v3.txt")
brain <- brain [ ,c("gene_id", "baseMean", "log2FoldChange", "padj")]
colnames (brain)[2:4] <- paste (colnames (brain)[2:4], "_brain", sep="")

muscle <- read.delim ("parl_star_rnaseq_muscle_v3.txt")
muscle <- muscle [ ,c("gene_id", "baseMean", "log2FoldChange", "padj")]
colnames (muscle)[2:4] <- paste (colnames (muscle)[2:4], "_muscle", sep="")


co <- merge (testis, brain, by="gene_id", all.x=TRUE, all.y=TRUE)
co <- merge (co, muscle, by="gene_id", all.x=TRUE, all.y=TRUE)
co <- co [order (co$log2FoldChange_testis), ]
head (co)
