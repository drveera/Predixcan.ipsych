#!/bin/env Rscript

if ( ! require(qqman)){
install.packages("qqman", repos = "http://cran.us.r-project.org")
}
library(qqman)

#argument 1: association file
#argument 2: ensembl file
#argument 3: output file


args = commandArgs(trailingOnly = TRUE)

assoc = read.table(args[1], header = FALSE)

names(assoc) = c("gene","beta","error","zvalue","p")
assoc = assoc [! is.na(assoc$p), ]
assoc$geneName = gsub("\\..*","",assoc$gene)

ensembl = read.table(args[2], header = TRUE)
assoc.ensembl = merge(assoc, ensembl, by = "geneName")
assoc.ensembl = assoc.ensembl [! duplicated(assoc.ensembl$gene), c("GENE","chrom","txStart","p","zvalue")]
assoc.ensembl$chrom = as.numeric(gsub("chr","",assoc.ensembl$chrom))

assoc.ensembl.sig = assoc.ensembl [ assoc.ensembl$p < (0.05/nrow(assoc.ensembl)),]
write.table(assoc.ensembl.sig, paste0(args[3],"significant"), sep = "\t", row.names = FALSE, quote = FALSE)

pdf(paste0(args[3],".manhattan.plot"))
manhattan(assoc.ensembl, chr = "chrom", bp = "txStart", p = "p", snp = "GENE", ylim=c(0,10) , cex.axis = 0.9,
	col = c("blue4","orange3"), main = basename(args[3]),  suggestiveline = FALSE, genomewideline = -log10(0.05/nrow(assoc.ensembl)))
dev.off()
