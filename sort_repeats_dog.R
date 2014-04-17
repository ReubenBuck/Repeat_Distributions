

rm(list = ls())

setwd("Desktop/new_PCA/script/")
# new PCA on chr 16 repeats

require(GenomicRanges)



bin <- read.table("../bins/D_bin.txt", header=T)


bin.gr <- GRanges( seqnames = Rle(bin[,1]),
					ranges = IRanges(start = bin[,2], end = bin[,3])
					)


rep <- read.table("../repeat_files/canFam2/canFam2_all_chr")


# process the repeat names so they are readable




mamsum <- read.table("../repeat_libraries/dog/summary_dog2", sep = "\t")


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



C <- NULL
for(i in 1:dim(bin)[1])

{
	
	b <- R[Q[Q[,1] == i,2],]
	
	c <- summary(b[,16], maxsum = length(unique(R[,16])))
	
	C <- rbind(C,c)
}

write.table(C, file = "../../dog-sort_p", quote = FALSE, sep = "\t", row.names = FALSE, col.names= TRUE)









rownames(C) <- 1:dim(C)[1]


pca <- prcomp(C)
pca.s <- prcomp(C, scale.=T)


# need to select the ones we want

#for now just get the largest 20 families

S <- order(colSums(C), decreasing=TRUE)[c(1:15)]

C1 <- C[,S]

plot(hclust(dist((t(C1)))))

# this program is coming along I,ve got to find better ways at streamlining the process

pca <- prcomp(C1)
pca.s <- prcomp(C1, scale.=T)

C2 <- C1[-823,]

pca <- prcomp(C2)
pca.s <- prcomp(C2, scale.=T)


biplot(pca, xlabs = rep(".", dim(pca$x)[1]), cex = c(2,.5))




# a clsiification tool I may want to use
#M[M[,1] == "L2",15] <- "L2_Eutheria"


# something is wrong with the way repeats are being counted, it could be that the labels are getting mixed up.
# The 

