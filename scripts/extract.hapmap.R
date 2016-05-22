#Predixcan made on 11-05-2016

#Overall Workflow
#read the chunk info and hapmap file
#merge and subset the hapmap variants
#quality filter
#subset and convert chunks
#once all chunks are converted, combine them chromosome wise
#submit for predictions
#once predictions are done, submit for asscociations


#arguments
args = commandArgs(trailingOnly = TRUE)

#args1 = hapmap file with absolute path
#args2 = chunk file with absolute path
#args3 = chunk info file with full name and full path (since its in seperate folder)
#args4 = outfile name
#args5 = infoscore threshold
#LOAD***********
#packages

library(data.table)

#read the necessary files
#1 read the hapmap file (corresponsing chromosome)
#read only 3 columns(chr,pos,rsid)
hapmap=fread(args[1]) #arg1 -> hapmapfile with absolute location
hapmap = hapmap[,c(2,4,5),with=FALSE]
names(hapmap) = c("chrom","pos","rsid")

#reading chunk map is not required, since we have the variants info in info file itself 

#3 read the info file 
chunk.info = fread(args[3])
#rename the column name: POS to pos
names(chunk.info)[3] <- "pos"

#4 read the chunk
chunk = fread(paste0("zcat ",args[2]), header = TRUE)
#chunk = read.table(paste0(args[2],".gz"), header = TRUE, colClasses = c(rep("character",3), rep("numeric",m.no-3)))


#DEFINE FUNCTIONS *********


#function to convert dosage file to allele dosage format
makeDosage = function(chunk,n){ #n depends on the number columns before the dosage numbers in the resulting merged file
	colNos = seq(n,ncol(chunk),2)

	for(i in colNos){
        print(i)
	chunk[[i]] = chunk[[i]] + 2 * (1- (chunk[[i]] + chunk[[(i-1)]]))
	}
        return(chunk[,c(1:5,colNos), with = FALSE])
}
#function to quality filter 
qfilter = function(info, infoscore = 0.6){  
  info1 = info[info > infoscore,]
  return(info1)
}
  
  
#function to subset dosage file

#PROCESS THE FILES **********

#quality filter chunk info

if (is.null(args[5])){
	infoscore = 0.6
} else {
	infoscore = as.numeric(args[5])
}
info.filtered = qfilter(chunk.info, infoscore = 0.6)

#merge with  hapmap and subset info file
info.hapmap = merge(chunk.info[,c(1:3),with = FALSE],hapmap,by = "pos", sort = FALSE ) #merge1 results in 5 columns: pos CHR SNP chrom rsid

#merge with info.hapmap, and subset chunk
chunk.subset = merge(info.hapmap,chunk,by = "SNP", sort = FALSE) #merge2, SNP pos CHR chrom rsid|A1 A2
chunk.subset = chunk.subset[,c(4,5,2,6:ncol(chunk.subset)),with=FALSE ]

#convert to allele dosage 
chunk.dosage  = makeDosage(chunk.subset,7)

chunk.dosage$maf = apply(chunk.dosage[,c(-1,-2,-3,-4,-5),with = FALSE], 1, function(x) sum(x)/(2*length(x)))

chunk.dosage = chunk.dosage[,c(1:5,ncol(chunk.dosage),6:(ncol(chunk.dosage)-1)), with = FALSE]
write.table(chunk.dosage,args[4], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#OUTPUT**********


#3 merge the map file with hapmap file
#4 quality filter the merged map file
#5 make a list of variants with positions
#6 read the 
