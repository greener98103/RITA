#####################################
########### GO enrichment ###########
#####################################

xx = read.table("../data/hsapiens.ENTREZ.txt",
                sep="\t", header=F)
Ann.Entrez = as.character(unlist(xx[,2]))
names(Ann.Entrez) = as.vector(unlist(xx[,1]))

f.GOenrich = function(fg, cut=0.01){
  fg.x = Ann.Entrez[fg]; fg.x = fg.x[!is.na(fg.x)]
  bg.x = Ann.Entrez[rownames(TopCorrGene)]; bg.x = bg.x[!is.na(bg.x)]
  A = GOenrichment(genesOfInterest=fg.x, allgenes=bg.x, cutoff=cut)
  a = sort(A$p.values)
  b = A$GOTerms$Term
  names(b) = A$GOTerms$go_id
  return(cbind(a, b[names(a)]))
}

ITA.modules.GOenrich = list()
for(i in 1:length(ITA.modules)){
  ITA.modules.GOenrich[[i]] = read.table(paste("/Users/greener/Desktop/RITA/Module.GOenrich.", i, ".txt", sep=""),sep="\t", quote="", row.names=1)
}
sapply(FUN=function(df){return(unname(df[1:10,2]))}, X=ITA.modules.GOenrich[1:5])
