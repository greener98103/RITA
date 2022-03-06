#' Load a Matrix and Apply User defined filters
#' This function loads pre-compiled correlations scores and values
#' Then a user defined set of genes filters tables and then edges and weights are calculated.
#' Once calculated it will plot the distance between modules, export them to file, and create an
#' R object that can passed on to further functions.
#' #' @param infile Path to the input file
#' @return 
#' @export



load_data <- function(w,x,y,z){
  #w=TopCorrGene.csv
  #x=TopCorrVal.csv
  #y=Blood_PBMC2_TC_DSweight.rds
  #z=User defined list Must have GeneID as a header
  # K is preset at 5
  
  user_list <- read.csv("/data/example_user_list.csv")
  
  if (length(user_list) == 0) {
    stop("User has not provided a appropriate gene list. Using default")
  }
  
  TopCorrGeneFull = read.csv("/data/blood/TopCorrGene.csv", row.names=1)
  TopCorrVal = as.matrix(read.csv("data/blood/TopCorrVal.csv", row.names=1))
  
  
  #TopCorrGene <- y$GeneID %in% row.names(TopCorrGeneFull)
  TopCorrGene <- TopCorrGeneFull
  
  K = 5
  ftM1 = cbind(rep(rownames(TopCorrGene), each=K), as.vector(t(TopCorrGene[,1:K])))
  ftM1val = rep(1:K, times=nrow(TopCorrGene))
  ftM2 = t(apply(FUN=sort, X=ftM1, MARGIN=1))
  ftM3 = unique(cbind(ftM2,ftM1val))
  
  TC.Graph = graph.data.frame(ftM3[,1:2],directed=F)
  E(TC.Graph)$weight = as.numeric(ftM3[,3])
  TC.Graph = simplify(TC.Graph,edge.attr.comb="min")
  
  # return TC.graph 
  # add new function that take the output.  
  
  ## Compiling edge weights from each dataset.
  TC.Graph.DF = get.data.frame(TC.Graph)
  EdgeW.DS = readRDS("Blood_PBMC2_TC_DSweight.edited.rds")
  
  #ignore
  # Ordering filenames
  #nn = sapply(FUN=rownames, X=EdgeW.DS)
  #nn2 = unique(unlist(nn))
  #nx = as.numeric(substring(sapply(FUN=getElement, X=strsplit(nn2,split=".",fixed=T), 2),first=4))
  #nn3 = nn2[order(nx)]
  
  
  # Compiling weights per edge.
  fx = function(ii){
    if(TC.Graph.DF[ii,2]%in%colnames(EdgeW.DS[[TC.Graph.DF[ii,1]]])){
      return(EdgeW.DS[[TC.Graph.DF[ii,1]]][,TC.Graph.DF[ii,2]][nn3])}
    if(TC.Graph.DF[ii,1]%in%colnames(EdgeW.DS[[TC.Graph.DF[ii,2]]])){
      return(EdgeW.DS[[TC.Graph.DF[ii,2]]][,TC.Graph.DF[ii,1]][nn3])}
  };
  
  #ignore
  #TC.Graph.DF.EdgeW.DS = t(sapply(FUN=fx, X=1:nrow(TC.Graph.DF)))
  #Mx = cor(TC.Graph.DF.EdgeW.DS,use="pair")
  #plot(hclust(as.dist(1-Mx),method="ward.D2"),labels=F)
  #plot(hclust(as.dist(1-Mx),method="ward.D2")$height,type="l")
  #lines(diff(hclust(as.dist(1-Mx),method="ward.D2")$height)*4)
  
  #Does this need to be rerun every time? RRG 1/29/22
  
  TC.Gr.Dist = distances(TC.Graph)
  #saveRDS(TC.Gr.Dist, "/Users/Barrenas/Work_files/Resource_Cadence/ITRs/Blood_PBMC2/TopCorr/TopCorrDist.RDS")
  #TC.Gr.Dist = readRDS("/Users/Barrenas/Work_files/Resource_Cadence/ITRs/Blood_PBMC2/TopCorr/TopCorrDist.RDS")
  
  TC.Gr.Dist.HC = hclust(as.dist(TC.Gr.Dist), "average")
  
  #saveRDS(TC.Gr.Dist.HC, "TopCorrHclust.RDS")
  #TC.Gr.Dist.HC = readRDS("TopCorrHclust.RDS")
  
  TC.Gr.Dist.HC.cut = cutreeHybrid(dendro=TC.Gr.Dist.HC, distM=TC.Gr.Dist,
                                   deepSplit=1, minClusterSize = 100)
  table(TC.Gr.Dist.HC.cut$labels)
  ITA.modules = list()
  for(ii in 1:max(TC.Gr.Dist.HC.cut$labels)){
    ITA.modules[[ii]] = TC.Gr.Dist.HC$labels[TC.Gr.Dist.HC.cut$labels==ii]}
  names(ITA.modules) = as.character(1:length(ITA.modules))
  saveRDS(ITA.modules, "ITA.TLN.Modules.RDS")
  
  ## Distances between modules:
  Mod.Dist = matrix(ncol=max(TC.Gr.Dist.HC.cut$labels), nrow=max(TC.Gr.Dist.HC.cut$labels))
  rownames(Mod.Dist) = as.character(1:nrow(Mod.Dist))
  for(ii in 2:nrow(Mod.Dist)){
    for(jj in 1:(ii-1)){
      print(c(ii,jj))
      Mod.Dist[ii,jj] = mean(TC.Gr.Dist[TC.Gr.Dist.HC.cut$labels==ii,TC.Gr.Dist.HC.cut$labels==jj])
      Mod.Dist[jj,ii] = Mod.Dist[ii,jj]
    }
  }
  hist(Mod.Dist[lower.tri(Mod.Dist)],40)
  Mod.Dist.HC = hclust(as.dist(Mod.Dist),"average")
  plot(Mod.Dist.HC)
  return(Mod.Dist.HC)
  write.csv(Mod.Dist.HC, file= "Mod.Dist.HC.csv")
}
