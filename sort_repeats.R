

# this script serves no other purpose other then generating multidimensional tables for genomic repeat contant



rm(list = ls())

setwd("~/Documents/phd/Desktop analyses/new_PCA/script")
# new PCA on chr 16 repeats

require(GenomicRanges)



bin <- read.table("../bins/H_bin.txt", header=T)
bin.gr <- GRanges( seqnames = Rle(bin[,1]),
					ranges = IRanges(start = bin[,2], end = bin[,3])
					)
rep <- read.table("../repeat_files/hg19/hg19_all_chr")

# process the repeat names so they are readable
mamsum <- read.table("../repeat_libraries/human/summary_human2", sep = "\t")
colnames(mamsum) <- "V4"
R <- merge(rep, mamsum)

# set up colnames to make the table subsetable
name <- c("fam", "chrom", "start", "end", "r_start", "r_end", "strand", "score1", "score2", "score3", "zero", "score4", "type", "species")
colnames(R) <- name


# get rid of spaces
for(i in seq(length(R))){
	if(class(R[,i]) == "factor"){
		R[,i] <- gsub(" ", "_", R[,i])
	}
	}

R$type_species <- paste(R$type, R$species,sep ="__")
R$fam_species <- paste(R$fam, R$species,sep ="__")

R.gr <- GRanges(seqnames= Rle( R$chrom),
				ranges = IRanges( start = R$start, end = R$end))


# produce an overlap table that will show all the R entries that will overlap with each bin
bin_Rep_Ol <- as.matrix(findOverlaps(bin.gr, R.gr))

# Sort repeats into their types and their bins

identifier <- unique(R$type_species)
Counts <- data.frame(rep(data.frame(rep(0, dim(bin)[1])), length(identifier)))
colnames(Counts) <- identifier
for(i in 1:dim(bin)[1]){
	b <- R$type_species[bin_Rep_Ol[bin_Rep_Ol[,1] == i,2]]
	for(z in seq(dim(Counts)[2])){
		Counts[i,z] <- length(b[b == colnames(Counts)[z] ])
	}
}
	

# Repeats are in according to an identiffier
# maybe add extra columns for bin chr start stop
# data will then be easier to upload 

write.table(C, file = "../../human-sort_p_gene,species,family", quote = FALSE, sep = "\t", row.names = FALSE, col.names= TRUE)




system.time(R$type_species[bin_Rep_Ol[bin_Rep_Ol[,1] == i,2]])
system.time(Counts[i,z] <- length(b[b == colnames(Counts)[z] ]))


# this thing should take around 14 hours

# there is probably that other table package that can do it much quicker
# on the data sci workshop


b <- R$type_species[bin_Rep_Ol[bin_Rep_Ol[,1] == i,2]]
	for(z in seq(dim(Counts)[2])){
		Counts[i,z] <- length(b[b == colnames(Counts)[z] ])
	}




