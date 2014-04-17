

rm(list = ls())

setwd("~/Documents/phd/Desktop analyses/new_PCA/script")
# new PCA on chr 16 repeats

require(GenomicRanges)



bin <- read.table("../bins/H_bin.txt", header=T)


bin.gr <- GRanges( seqnames = Rle(bin[,1]),
					ranges = IRanges(start = bin[,2], end = bin[,3])
					)


rep <- read.table("../repeat_files/hg19/hg19_all_chr")
gene <- read.table("../repeat_files/hg19/refGene.human.txt")

# process the repeat names so they are readable




mamsum <- read.table("../repeat_libraries/human/summary_human2", sep = "\t")


colnames(mamsum) <- "V4"


R <- merge(rep, mamsum)
R[,14] <- as.factor(R[,14])
R[,15] <- paste(R[,13], R[,14],sep ="__")
R[,16] <- paste(R[,1], R[,15],sep ="__")
R[,15] <- as.factor(R[,15])
R[,16] <- as.factor(R[,16])



# add the tissue specific stuff here later
G <- data.frame(gene[,2], gene[,3], gene[,5], gene[,6], NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "gene" , "gene")


colnames(G) <- rep("x", 16)
colnames(R) <- rep("x", 16)
XX <- R
R <- rbind(R,G)



R.gr <- GRanges(seqnames= Rle( R[,2]),
				ranges = IRanges( start = R[,3], end = R[,4]),
				rep = R[,1],
				 type = R[,13], 
				 species = R[,14]
)


# find the overlaps

Q <- as.matrix(findOverlaps(bin.gr, R.gr))


v <- length(unique(R[,15]))
v1 <- data.frame(row.names = unique(R[,15]))
v1[,1] <- NA
s <- t(v1)

C <- data.frame(rep(data.frame(rep(0, dim(bin)[1])), length(s)))
colnames(C) <- colnames(s)

for(i in 1:dim(bin)[1])
{
	b <- R[Q[Q[,1] == i,2],]
	for(z in seq(dim(C)[2])){
		C[i,z] <- dim(b[b[,15] == colnames(C)[z] , ])[1]
	}
}
	
	
	
	c <- summary(b[,15], maxsum = v )
	c2 <- data.frame(c[1:v])
	
	#s[1, match(rownames(c2) , colnames(s))] 
	
	q <- c2[match(colnames(s),rownames(c2) ),1]
	q[is.na(q)] <- 0
		#assign(paste("d", i, sep="") , c2)
	#d <- data.frame(v1, c2[(rownames(c2) %in% as.character(v1[,1])),1] )
	
	C <- rbind(C,q)
}


DDD <- C

colnames(C) <- unique(R[,15])
# summary doesn't really sort the way I want it to

write.table(C, file = "../../human-sort_p_gene,species,family", quote = FALSE, sep = "\t", row.names = FALSE, col.names= TRUE)



rownames(C) <- 1:dim(C)[1]


pca <- prcomp(C)
pca.s <- prcomp(C, scale.=T)


# need to select the ones we want

#for now just get the largest 20 families

S <- order(colSums(C), decreasing=TRUE)[1:20]

C1 <- C[,c(S,1260)]

plot(hclust(dist((t(C1)))))

# this program is coming along I,ve got to find better ways at streamlining the process

pca <- prcomp(C1)
pca.s <- prcomp(C1, scale.=T)

C2 <- C1[-823,]

pca <- prcomp(C2)
pca.s <- prcomp(C2, scale.=T)


biplot(pca.s, xlabs = rep(".", dim(pca$x)[1]), cex = c(2,1))
plot(pca.s$x[,1], pca.s$x[,2])



# a clsiification tool I may want to use
#M[M[,1] == "L2",15] <- "L2_Eutheria"




# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

d <- dist(((C2))) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]

plot(x,y)

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
  main="Metric MDS", type="n")
text(x, y, labels = row.names((C2)), cex=.7)


Y <- colnames(C)
colnames(C) <- gsub(" ", "_", colnames(C))

C2 <- data.frame(C$L1__Eutheria, C$L1__Homo_sapiens, C[,5], C$ERV3__Homo_sapiens, C$ERV3__Eutheria,C[,113], C$CR1__Eutheria) 



biplot(prcomp(C2, scale.=TRUE))







