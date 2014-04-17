# clustering analysis



library("FactoMineR")


rm(list=ls())


C <- read.table("~/Desktop/dog-sort_p", header= TRUE, sep="\t")


#T <- C[,grep("MIR", colnames(C))]

S <- order(colSums(C), decreasing = TRUE)
pie(colSums(C[,S]))
S <- order(colSums(C), decreasing = TRUE)[1:20]
C1 <- C[,S]
pie(colSums(C1), main = "DOG top 20 repeats")

P <- PCA(C1)

Q <- HCPC(P)
plot.HCPC(Q, choice = "map", ind.names=FALSE)






rm(list=ls())



C <- read.table("~/Documents/phd/Desktop analyses/new_PCA/sort results/human-sort_p_gene", header= TRUE, sep="\t")




S <- order(colSums(C), decreasing = TRUE)[1:20]

C1 <- C[,S]
pie(colSums(C1))
P <- PCA(C1)

Q <- HCPC(P)
plot.HCPC(Q, choice = "map", ind.names=FALSE)

P <- PCA(C[,-(grep("SAT", colnames(C)))], scale.unit=FALSE)

Q <- HCPC(P)
plot.HCPC(Q, choice = "map", ind.names=FALSE)


plot.HCPC(T, choice = "map", ind.names=FALSE)

HCPC(T)

T <- C[,grep("L1", colnames(C))]

T <- t(T)

HCPC(T)




head(C)



#####
#lets try MDS

# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

d <- dist(C) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
  main="Metric MDS", type="n")
text(x, y, labels = row.names(mydata), cex=.7)







