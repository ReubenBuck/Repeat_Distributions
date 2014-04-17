
rm(list = ls())
setwd("~/")
setwd("Desktop/new_PCA/script/")
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
R[,14] <- as.factor(R[,14])
R[,15] <- paste(R[,13], R[,14],sep ="_")
R[,16] <- paste(R[,1], R[,15],sep ="_")
R[,15] <- as.factor(R[,15])
R[,16] <- as.factor(R[,16])


R.gr <- GRanges(seqnames= Rle( R[,2]),
				ranges = IRanges( start = R[,3], end = R[,4]),
				rep = R[,1],
				 type = R[,13], 
				 species = R[,14]
)



# find the overlaps

Q <- as.matrix(findOverlaps(bin.gr, R.gr))


v <- length(unique(R[,16]))
C <- NULL
for(i in 1:dim(bin)[1])

{
	
	b <- R[Q[Q[,1] == i,2],]
	
	c <- summary(b[,16], maxsum = v)
	C <- rbind(C,c)
}


write.table(C, file = "../../human-sort_p", quote = FALSE, sep = "\t", row.names = FALSE, col.names= TRUE)

