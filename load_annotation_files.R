#' Read Annotation files and read functions.

## Annotation files
xx = read.table("data//hsapiens.ENTREZ.txt", sep="\t", header=F)
Ann.Entrez = as.character(unlist(xx[,2]))
names(Ann.Entrez) = as.vector(unlist(xx[,1]))

xx = read.table("data//hsapiens.SYMBOL.txt",sep="\t", header=F)
Ann.Symbol = as.character(unlist(xx[,2]))
names(Ann.Symbol) = as.vector(unlist(xx[,1]))