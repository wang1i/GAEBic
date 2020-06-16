miRNA_name = as.matrix(read.table("C:/R-biclust/miRNA nodes.name"))
mRNA_name = as.matrix(read.table("C:/R-biclust/mRNA nodes.name"))

A <- as.matrix(read.table("C:/R-biclust/miRNA-mRNA matrix.txt"))
res <- biclust(A, method=BCSpectral(), normalization="log", numberOfEigenvalues=3, minr=4, minc=6, withinVar=1)

writeBiclusterResults("C:/R-biclust/Spectral results(4X6).txt", res,"Spectral with 4 X 6", miRNA_name, mRNA_name)

