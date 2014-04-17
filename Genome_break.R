



rm(list = ls())

setwd("Desktop/new_PCA/script/")

x <- 1500000
y <- .1


# select a UCSC genome with ucscGenomes()

library (rtracklayer)
mySession = browserSession("UCSC")
genome(mySession) <- "canFam2"
tbl.gap <- getTable(ucscTableQuery(mySession, track="gap",table="gap"))
ass <- getTable(ucscTableQuery(mySession, track="chainNetEquCab2"))


require("BSgenome.Cfamiliaris.UCSC.canFam2")

sLengths=seqlengths(Cfamiliaris)
L <- sLengths[sLengths > x]

N <- names(L)



B <- NULL

for(i in seq(L)){
	
	s <- as.integer((L[i]%%x)/2)
	
	q <- seq(1:as.integer(L[i]%/%x))
	q2 <- q*x
	q3 <- (q2 + 1)- x
	
	b <- data.frame(N[i], q3,q2)
	#b[,2:3] <- b[,2:3] + s
	
		B <- rbind(B,b)
	
	
}

colnames(B) <- c("chr", "start", "end")


# need to handle gaps
#see if it is possible without GRanges.
#although that does make it easy



gapsgr <- GRanges(seqnames= Rle(tbl.gap[,2]), ranges = IRanges(start= tbl.gap[,3], end = tbl.gap[,4]), type = tbl.gap[,8])

binsgr <- GRanges(seqnames= Rle(B[,1]), ranges = IRanges(start=B[,2], end = B[,3]))


R <- as.matrix(findOverlaps(gapsgr,binsgr, minoverlap = as.integer(y *x)))[,2]



bins <- B[- R,]




write.table(bins, file="../bins/D_bin.txt", quote = FALSE, sep = "\t", row.names=FALSE)





