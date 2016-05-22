#!usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

#args[1] hapmap file
require(data.table)
hapmap <- fread(args[1])

hapmap.split = split(hapmap, hapmap$V2)

for (i in 1:length(hapmap.split)){
	write.table(hapmap.split[[i]],paste0(names(hapmap.split[i]),".chr.txt"), sep = "\t",
		    quote = FALSE, row.names = FALSE, col.names = FALSE)
}
