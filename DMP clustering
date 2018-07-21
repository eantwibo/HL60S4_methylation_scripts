clus1 <- hclust(dist(data_dmp))
d2 <- data_dmp[clus1$order,]
mycl <- cutree(clus1, h=max(hr$height/2))
myheatcol <- rev(redgreen(75))
mycl  <- mycl[clus1$order]
data_cluster <- cbind(data_dmp, clusterID=mycl)
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
write.table(data_cluster,"DMP_clustered_numbered.csv", quote=F, sep="\t", row.names=T, col.names=F)
