#' TLN construction and plotting
#' 
#' @param infile Path to the input file
#' @export

## TLN construction and plotting
f.tree.leaf.network = function(hc, nodesize=7){
  ftM1 = -hc$merge
  fx = function(i){
    v = ftM1[i,]
    #if(v[1]<0 & v[2]>0){return(cbind(v,c(v[1],-i)))} else {
    return(rbind(v,rep(-i,2)))
  }
  ftM2 = t(matrix(as.vector(sapply(FUN=fx,X=1:nrow(ftM1))),nrow=2))
  g = graph.data.frame(ftM2,directed=F)
  
  aa = hc$height; names(aa) = -(1:length(aa))
  col.x = sapply(FUN=function(s){return(paste(c("grey",s),collapse=""))},X=40:95) 
  E(g)$color = rep(col.x[floor((aa-min(aa))*55/max((aa-min(aa))))+1], each=2)
  E(g)$width = 2
  
  label.index = as.numeric(V(g)$name)
  label = V(g)$name
  label[label.index>0] = hc$label[label.index[label.index>0]]
  label[as.numeric(label.index)<0] = make.names(signif(aa[label[as.numeric(label.index)<0]],3))
  V(g)$label = label
  V(g)$size = nodesize; V(g)$size[igraph::degree(g)>1] = 0
  V(g)$label.family = "sans"
  #g$layout = layout_with_sugiyama(g)$layout
  
  V(g)$size[igraph::degree(g)==1]=nodesize
  V(g)$size[igraph::degree(g)>1]=0.1
  V(g)$label.cex=0.6
  #plot(g)
  return(g)
}
f.TLNplot = function(TLN, ValueM, colorsat=NA, colors=NA, title="NoTitle", markedNodes=NA, valueCap=NA, writeFile=F){
  if(all(is.na(colors))){
    colors=c('darkblue', 'mediumblue', 'dodgerblue', 'white', 'orange', 'red', 'darkred')
  }
  if(is.vector(ValueM)){ValueM = t(as.matrix(ValueM))}
  
  neg.col = rev(colors[1:ceiling(length(colors)/2)])
  pos.col = colors[ceiling(length(colors)/2):length(colors)]
  
  bins = rep(100,(length(neg.col)-1))
  a1 = list()
  a2 = list()
  for(i in 1:length(bins)){
    a1[[i]] = t(colorRamp(neg.col[c(i, i+1)])(seq(from=0, to=1, length.out=bins[i])))
    a2[[i]] = t(colorRamp(pos.col[c(i, i+1)])(seq(from=0, to=1, length.out=bins[i])))
  }
  a1 = matrix(unlist(a1), ncol=3, byrow=T); a2 = matrix(unlist(a2), ncol=3, byrow=T)
  
  b1 = rgb(red=a1[,1], green=a1[,2], blue=a1[,3], max=255)
  b2 = rgb(red=a2[,1], green=a2[,2], blue=a2[,3], max=255)
  color.vector = c(rev(b1), b2, b2[length(b2)])
  
  nn = V(TLN)$name[igraph::degree(TLN)==1]
  ps = lapply(FUN=function(n, x){return(rep(1,n))}, X=1:ncol(ValueM), n=nrow(ValueM))
  
  if(!is.na(valueCap)){ValueM[ValueM<(-valueCap)] = -valueCap; ValueM[ValueM > valueCap] = valueCap}
  if(is.na(colorsat)){mm = max(abs(range(ValueM,na.rm=T)))} else {mm = colorsat}
  if(!is.na(colorsat) & max(abs(range(ValueM,na.rm=T)))>colorsat){print("WARNING: Color saturation too low")}
  
  pc = color.vector[round(ValueM[,nn]*((length(neg.col)-1)*100)/mm)+(((length(neg.col)-1)*100)+1)]
  pc = matrix(pc, nrow=nrow(ValueM))
  pc = lapply(FUN=function(i){return(rev(pc[,i]))}, X=1:ncol(pc))
  TLNx = TLN
  V(TLNx)$label[is.na(suppressWarnings(as.numeric(V(TLNx)$label)))] = ""
  V(TLNx)$shape[igraph::degree(TLNx)==1] = "pie"
  #if(nrow(ValueM)==1){V(TLNx)$shape[igraph::degree(TLNx)==1] = "circle"}
  V(TLNx)$shape[is.na(V(TLNx)$shape)] = "circle"
  V(TLNx)$pie[igraph::degree(TLNx)==1] = ps
  V(TLNx)$pie.color[igraph::degree(TLNx)==1] = pc
  if(nrow(ValueM)==1){V(TLNx)$color[igraph::degree(TLNx)==1] = unlist(pc)}
  
  coord.v = signif(seq(from=-0.5, to=0.5, length.out=8),3)
  cv = color.vector[round(seq(from=1, to=length(color.vector), length.out=7))]
  plotColorScale = function(){
    for(i in 1:length(cv)){rect(xleft=coord.v[i], xright=coord.v[i+1], ytop=-1.2, ybottom=-1.3, col=cv[i], border=NA)}
  }
  
  if(writeFile==F){plot(TLNx, main=title); plotColorScale()}
  if(writeFile==T){
    pdf(paste(c(title, "pdf"),collapse=".")); plot(TLNx); plotColorScale(); title(title, cex.main=2); dev.off()
  }
  #return(TLNx)
}
f.TLNaddSpiderPlot = function(TLN, ValueM, minratio=0.25, sizefactor=2.7, minvalue.override=NULL, maxvalue.override=NULL, spokes=T){
  ## Will add spider plots to nodes, translating node names to ValueM colnames.
  if(is.null(colnames(ValueM))){colnames(ValueM) = as.character(1:ncol(ValueM))}
  ValueM[is.na(ValueM)] = min(ValueM, na.rm=T)
  
  coord.x = layout.norm(TLN$layout)
  angle.x = seq(from=0, to=(2*pi), length.out=nrow(ValueM)+1)[1:nrow(ValueM)]
  angle.x = angle.x+angle.x[2]/2
  
  if(!is.null(minvalue.override)){minv = minvalue.override} else {minv = min(ValueM)}
  if(!is.null(maxvalue.override)){maxv = maxvalue.override} else {maxv = max(ValueM)}
  
  ValueM.norm = (ValueM-minv)/(maxv-minv)*(1-minratio)+minratio
  points(x=coord.x[V(TLN)$size>1,1], y=coord.x[V(TLN)$size>1,2], cex=(V(TLN)$size[V(TLN)$size>1])*minratio/2.35, col="gray60")
  if(spokes==T){
    for(ii in 1:ncol(ValueM)){
      node.n=colnames(ValueM)[ii]
      radius.x = V(TLN)$size[V(TLN)$name==node.n]*sizefactor/diff(range(TLN$layout))
      coord.xx = c()
      for(i in 1:length(angle.x)){coord.xx[i] = coord.x[node.n,1]+sin(angle.x[i])*radius.x}
      coord.xy = c()
      for(i in 1:length(angle.x)){coord.xy[i] = coord.x[node.n,2]+cos(angle.x[i])*radius.x}
      for(i in 1:length(angle.x)){lines(x=c(coord.x[node.n,1],coord.xx[i]),y=c(coord.x[node.n,2],coord.xy[i]),col="gray60")}
    }
  }
  for(ii in 1:ncol(ValueM)){
    node.n=colnames(ValueM)[ii]
    radius.x = V(TLN)$size[V(TLN)$name==node.n]*sizefactor/diff(range(TLN$layout))
    coord.xx = c()
    for(i in 1:length(angle.x)){coord.xx[i] = coord.x[node.n,1]+sin(angle.x[i])*radius.x*ValueM.norm[i,node.n]}
    coord.xy = c()
    for(i in 1:length(angle.x)){coord.xy[i] = coord.x[node.n,2]+cos(angle.x[i])*radius.x*ValueM.norm[i,node.n]}
    lines(x=c(coord.xx, coord.xx[1]), y=c(coord.xy, coord.xy[1]))
  }
  text(paste(rownames(ValueM), collapse="\n"), x=-1,y=1.5, cex=0.7, adj=c(0,1))
}
f.CalculateModRatios = function(DEfiles){
  fa = function(DEfile){
    de.table = as.matrix(read.table(DEfile,row.names=1,header=T,sep="\t"))
    vx = names(Ann.Symbol); names(vx) = Ann.Symbol
    rownames(de.table) = vx[rownames(de.table)]
    de.up = rownames(de.table)[de.table[,2]<0.05 & de.table[,1]>log2(1.5)]
    de.down = rownames(de.table)[de.table[,2]<0.05 & de.table[,1]<(-log2(1.5))]
    fb = function(mod){return((length(intersect(de.up,mod))-length(intersect(de.down,mod)))/length(intersect(rownames(de.table),mod)))}
    return(as.matrix(sapply(FUN=fb, X=ITA.modules)))
  }
  ModRatios = t(sapply(FUN=fa, X=DEfiles))
  colnames(ModRatios) = as.character(1:ncol(ModRatios))
  return(ModRatios)
}