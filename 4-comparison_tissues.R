annot <- read.delim ("gencode.vM32.annotation.txt")
annot <- annot[ ,grep ("transcript_id", colnames (annot), invert=TRUE)]
annot <- unique (annot)

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

co <- merge (co, annot, by="gene_id")
co <- na.omit (co)
co <- co [order (co$log2FoldChange_testis), ]


co$trend <- "no"

for (i in 1:dim (co)[1]) {
if (
co$log2FoldChange_testis[i] <0 & co$log2FoldChange_brain[i] < 0 & co$log2FoldChange_muscle[i] < 0)
{co$trend[i] <- "down"}
if (
co$log2FoldChange_testis[i] >0 & co$log2FoldChange_brain[i] > 0 & co$log2FoldChange_muscle[i] > 0)
{co$trend[i] <- "up"}
  }



write.table (co, "parl_comparison_tissues_v3.txt", sep="\t", quote=F, row.names=F)








