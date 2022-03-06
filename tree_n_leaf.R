#' Load a Matrix
#'
#' This function performs a Tree and Leaf network analysis
#'
#' @param infile Path to the input file
#' @export


#read in Mod.Dist.HC
Mod.Dist.HC = as.matrix(read.csv("Mod.Dist.HC.csv", row.names=1))

ITA.TLN = f.tree.leaf.network(Mod.Dist.HC,nodesize=11)
sv = V(ITA.TLN)$size; for(i in 1:length(ITA.modules)){sv[V(ITA.TLN)$name==i] = length(ITA.modules[[i]])}
V(ITA.TLN)$ModuleSize = sv
V(ITA.TLN)$size = sqrt(V(ITA.TLN)$ModuleSize)*0.95; V(ITA.TLN)$size[V(ITA.TLN)$size==min(V(ITA.TLN)$size)] = 0.1
write.table(file="ITA.TLN.Blood_PBMC2.ftM.txt", get.edgelist(ITA.TLN), sep="\t", col.names=F, row.names=F, quote=F)
write.table(file="ITA.TLN.Blood_PBMC2.nodes.txt", get.data.frame(ITA.TLN,"vertices"), sep="\t", col.names=T, row.names=F, quote=F)

A = readLines("ITA.TLN.Blood_PBMC2.ftM.txt.cyjs")
a1 = A[sapply(FUN=function(s){return(grepl(x=s, '\"x\"'))}, X=A)]
a1 = as.numeric(sapply(FUN=function(s){return(substring(s, first=15, last=(nchar(s)-1)))}, X=a1))
a2 = A[sapply(FUN=function(s){return(grepl(x=s, '\"y\"'))}, X=A)]
a2 = as.numeric(sapply(FUN=function(s){return(substring(s, first=15, last=(nchar(s)-1)))}, X=a2))
a3 = A[sapply(FUN=function(s){return(grepl(x=s, '\"name\"'))}, X=A)][2:(length(a1)+1)]
a3 = unname(sapply(FUN=function(s){return(substring(s, first=19, last=(nchar(s)-2)))}, X=a3))
aa = cbind(a1,-a2); rownames(aa) = a3
ITA.TLN$layout = aa[V(ITA.TLN)$name,]

xx = as.matrix(read.table("ModuleNames.txt",sep="\t",header=T))
vx = as.character(xx[,1]); names(vx) =  as.character(xx[,2])
V(ITA.TLN)$label = vx[V(ITA.TLN)$name]

plot(ITA.TLN)
#saveRDS(ITA.TLN, file="/Users/Barrenas/Work_files/Resource_Cadence/ITRs/Blood_PBMC2/TLN/Blood_PBMC2_ITA.TLN.rds")
ITA.TLN = readRDS(file="/Users/Barrenas/Work_files/Resource_Cadence/ITRs/Blood_PBMC2/TLN/Blood_PBMC2_ITA.TLN.rds")

