#last updated 10 October 2022

#begin function definitions

## WORKAROUND: https://github.com/rstudio/rstudio/issues/6692
## Revert to 'sequential' setup of PSOCK cluster in RStudio Console on macOS and R 4.0.0
if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}

#function to read in extended newick tree and 
#re-order the data to match node numbers
#takes path to tree file
read.beast.fixed<-function(path, ladderize=F, equal=F){
	tmp<-treeio::read.beast(path)
	tmp@data<-tmp@data[order(as.numeric(tmp@data$node)),]
	
	if(ladderize==T){
	  tmp@phylo<-ladderize(tmp@phylo)
	}
	
	if(equal==T){
	  tmp@phylo$edge.length<-rep(1, length(tmp@phylo$edge.length))
	}
	
	return(tmp)
}

#alt version without resorting data
read.beast.fixed2<-function(path, ladderize=F, equal=F){
  tmp<-treeio::read.beast(path)
  #tmp@data<-tmp@data[order(as.numeric(tmp@data$node)),]
  
  if(ladderize==T){
    tmp@phylo<-ladderize(tmp@phylo)
  }
  
  if(equal==T){
    tmp@phylo$edge.length<-rep(1, length(tmp@phylo$edge.length))
  }
  
  return(tmp)
}

#capture model assignment as tip state data
model2tipdata<-function(output_path){
  gophytree_s4<-read.beast.fixed2(output_path, ladderize=T)
  
  tmp<-as.data.frame(gophytree_s4@data)[as.data.frame(gophytree_s4@data)$node<(length(gophytree_s4@phylo$tip.label)+1),]
  #gophytree_s4@phylo$edge[gophytree_s4@phylo$edge[,2]<(length(gophytree_s4@phylo$tip.label)+1),][,2]
  tmp<-as.data.frame(gophytree_s4@data)
  tmp$node<-as.numeric(tmp$node)
  tmp<-tmp[tmp$node<(length(gophytree_s4@phylo$tip.label)+1),]
  tmp$label<-gophytree_s4@phylo$tip.label
  
  return(tmp)
  
}

#function read in all of the gene trees in a 
#directory in extended newick format
#given directory path
treeload<-function(path, ladderize=F, equal=F, type=NULL){
	files<-list.files(path, pattern=type)
	trees<-list()

	for (i in 1:length(files)){
		print(paste("parsed ", i, "trees out of ", length(files), "or ", round(i/length(files)*100, digits=3), " %"))
		trees[[i]]<-read.beast.fixed(path=paste(path,files[i], sep="/"), ladderize, equal)
		#trees[[i]]@data<-trees[[i]]@data[order(as.numeric(trees[[i]]@data$node)),]
	}
	names(trees)<-unlist(strsplit((files), split=".tre.gophy.results.tre"))
	return(trees)
}

#function to match nodes from consensus 
#to individual gene trees with uneven sampling
#derived from Liam Revell's example-- need to test
match_nodes<-function(t1, t2){
	## step one drop tips
	t1p<-drop.tip(t1,setdiff(t1$tip.label, t2$tip.label))
	t2p<-drop.tip(t2,setdiff(t2 $tip.label, t1$tip.label))

	## step two match nodes "descendants"
	M<-matchNodes(t1p,t2p)

	## step two match nodes "distances"
	M1<-matchNodes(t1,t1p,"distances")
	M2<-matchNodes(t2,t2p,"distances")

	## final step, reconcile
	MM<-matrix(NA,t1$Nnode,2,dimnames=list(NULL,c("left","right")))

	for(i in 1:nrow(MM)){
		MM[i,1]<-M1[i,1]
    	nn<-M[which(M[,1]==M1[i,2]),2]
    if(length(nn)>0){	
    	if(length(which(M2[,2]==nn))>0){
    		MM[i,2]<-M2[which(M2[,2]==nn),1]
    	}
    } else {
    }   
}
return(MM)	
}

#functions from Eliot Miller to generate the matches
{
tableSetup <- function(tree1, tree2){
  output <- data.frame(matrix(nrow=min(c(length(tree1$tip.label),length(tree2$tip.label)))-1, ncol=2))
  names(output) <- c("tree1", "tree2")
  output
}
equivCut <- function(tree1, tree2, results, least.inclusive){
  #identify nodes that appear more than once
  dups <- unique(results[,2][duplicated(results[,2])])
  
  #if there aren't any, return the original results
  if(length(dups)==0)
  {
    return(results)
  }
  
  #loop through, add details on size of clade to data frame
  for(i in 1:length(dups))
  {
    #set aside all the results that pertain to that node
    setAside <- results[results[,2]==dups[i],]
    
    #go through each node and figure out how many taxa descend from that node
    temp <- data.frame(matrix(nrow=dim(setAside)[1], ncol=2))
    names(temp) <- c("node","no.taxa")
    for(j in 1:dim(setAside)[1])
    {
      temp[j,"node"] <- setAside[j,"tree1"]
      temp[j,"no.taxa"] <- length(extract.clade(tree1, setAside[j,"tree1"])$tip.label)
    }
    
    #now you can decide which of these to keep
    if(least.inclusive)
    {
      keep <- temp$node[temp$no.taxa==min(temp$no.taxa)]
      toDrop <- setdiff(temp$node, keep)
      results <- results[!(results$tree1 %in% toDrop),]
    }
    
    else
    {
      keep <- temp$node[temp$no.taxa==max(temp$no.taxa)]
      toDrop <- setdiff(temp$node, keep)
      results <- results[!(results$tree1 %in% toDrop),]
    }
  }
  
  results
}
filter<-function(ref,results){
  output<-matrix(nrow=ref$Nnode, ncol=2)
  output[,1]<-(length(ref$tip.label)+1):(length(ref$tip.label)+length(ref$tip.label)-1)
  rownames(output)<-output[,1]
  rownames(results)<-results[,1]
  
  merged<-merge(output, results, by=0, all=T)
  merged$Row.names<-NULL
  merged$V2<-NULL
  merged$tree1<-NULL
  colnames(merged)<-c("tree1", "tree2")
  merged
}
match_phylo_nodes <- function(tree1, tree2, least.inclusive=TRUE){
  #set the results table up
  results <- tableSetup(tree1, tree2)
  
  #define tree1 nodes here
  nodes1 <- (length(tree1$tip.label)+1):max(tree1$edge)
  
  for(i in 1:length(nodes1))
  {
    #set the correct row in the tree1 column to the node in question
    results[i,1] <- nodes1[i]
    
    #find the descendants of this node
    tips1 <- tips(tree1, nodes1[i])
    
    #drop these to just tips that are in tree 2
    tips1 <- tips1[tips1 %in% tree2$tip.label]
    
    #if there's nothing left, set to no match and move on
    if(length(tips1)==0)
    {
      results[i,2] <- "no.match"
      next()
    }
    
    #find the MRCA of those tips in tree 2
    #if there's only a single taxon left, pull the node it descends from (getMRCA will fail)
    if(length(tips1==1))
    {
      edge2 <- which.edge(tree2, tips1)
      mrca2 <- tree2$edge[edge2,][1]
    }
    else
    {
      mrca2 <- getMRCA(tree2, tips1)
    }
    
    #find the descendants of that node
    tips2 <- tips(tree2, mrca2)
    
    #drop to just taxa that are in tree1
    tips2 <- tips2[tips2 %in% tree1$tip.label]
    
    #if these tips are equivalent, the nodes match
    if(setequal(tips1, tips2))
    {
      results[i,2] <- mrca2
    }
    
    else
    {
      results[i,2] <- "no.match"
    }
  }
  #print(results)
  #consider what you want to do next. do you want to, e.g., add a row for
  #every not-yet-mentioned tree2 node and add "no.match" to the first column,
  #or do you want to return as is, or do you want to only return matches? for
  #now, drop to only matches and re-class both as integer so can demonstrate
  #equality no matter which tree is first or second
  results <- results[results[,2] != "no.match",]
  results[,2] <- as.integer(results[,2])
  
  #run the equivCut function
  results <- equivCut(tree1, tree2, results, least.inclusive)
  
  results <- filter(ref=tree1, results)
  results
}
}

#function to create a list of node
#matches based on the reference and targets
#using match_phylo_nodes
readmatches<-function(ref, targets){
	matches<-list()
	
	for (i in 1:length(targets)){
		matches[[i]]<-match_phylo_nodes(ref@phylo, targets[[i]]@phylo)
		print(paste("parsed ", i, "trees out of ", length(targets), "or ", round(i/length(targets)*100, digits=3), " %"))		
	}
	
	names(matches)<-names(targets)
	return(matches)
}

#parallel version of readmatches
loopmatch<-function(i, targets, ref){
  #matches[[i]]<-
  match<-match_phylo_nodes(ref@phylo, targets[[i]]@phylo)
  #print(paste("parsed ", i, "trees out of ", length(targets), "or ", round(i/length(targets)*100, digits=3), " %"))
  return(match)
}
readmatches.mc<-function(targets, ref, fun){
  
  #set up cluster
  
  #define number of cores to use (max cores - 1)
  no_cores <- detectCores() - 1
  #how long is your data
  n <- length(targets)
  #spool up the cluster
  cl <- makeCluster(no_cores)
  registerDoSNOW(cl)
  
  #for the length of data, run fun_parallel on each in parallel
  result<-foreach(i=1:n, .export=ls(globalenv()), .packages=c('ape', 'phytools', 'geiger'), .verbose=T) %dopar% fun(i, targets, ref)
  
  #kill the cluster
  stopCluster(cl)
  
  #add the names back in
  names(result)<-names(targets)
  
  #return data
  return(result)
}
#parallel version of readmatches


###function to parse the output from readmatches and generate a summary matrix
parsematches<-function(ref, matches){
	matches.matrix<-matrix(nrow=ref@phylo$Nnode, ncol=length(matches)+1)
	matches.matrix[,1]<-matches[[1]][,1]
	colnames(matches.matrix)<-	c("ref", names(matches))

	for (i in 1:length(matches)){
		matches.matrix[,i+1]<-matches[[i]][,2]
	}

	return(matches.matrix)
}

#function to loop across nodes to retreive the uncex values from
#a matching node in a gene
loopone<-function(column, genes, matches.sum){
	
	ref<-matches.sum[,1]
	data<-matches.sum[,2:length(matches.sum[1,])]
	
	col <- data[, column]
	output <- rep(NA, length(data[, 1]))

	for (i in 1:length(ref)){
		#i<-1
		refindex<-ref[i]-(length(ref)+1)
		nodematch<-col[refindex]
	
		if (is.na(nodematch) == T){
			output[i]<-NA
		} else {
			if (colnames(genes[[column]]@data[nodematch,])[3]=="uncex"){
				output[i]<-genes[[column]]@data[nodematch,]$uncex[[1]]
			}
		}
	}
	return(output)
}

#now loop across genes
geneloop<-function(genes, matchsum) {
	summary<-matrix(nrow=length(matchsum[,1]), ncol=length(matchsum[1,])-1)
	colnames(summary)<-names(genes)
	for (i in 1:length(genes)){	
		summary[,i] <- loopone(column=i, genes=genes, matches.sum= matchsum)
	}
	summary<-!is.na(summary)
	return(summary)
}

#function to aggregate the helper functions
match_generator<-function(refpath, targetpath, ladderize=T, kind=""){
  consensus<-read.beast.fixed(refpath, ladderize)
  genetrees<-treeload(targetpath, ladderize, type=kind)
  matches<-readmatches.mc(targets=genetrees, ref=consensus, fun=loopmatch)
  matches.summary<-parsematches(ref=consensus, matches=matches)
  genematch<-geneloop(genes=genetrees, matchsum=matches.summary)
  genematch
}

#calculate concordance as fraction of nodes that exist in reference
#given list of data frame matches
concordance_perc<-function(matches){
  tmp<-do.call("cbind", lapply(matches, "[", 2))
  len<-length(tmp[1,])
  sums<-list()
  for(i in 1:length(tmp[,1])){
    sums[[i]]<-sum(!is.na(tmp[i,]))
    sums<-unlist(sums)
  }
  
  frac<-sums/len
  
}

#function to calculate % of gene trees in which a node is recovered
gene_concordance_perc<-function(refpath, targetpath, ladderize=T, kind=""){
  consensus<-read.beast.fixed(refpath, ladderize)
  genetrees<-treeload(targetpath, ladderize, type=kind)
  matches<-readmatches.mc(targets=genetrees, ref=consensus, fun=loopmatch)
  
  fraction<-concordance_perc(matches)
  
  nodenums<-length(consensus@phylo$tip.label)+1:consensus@phylo$Nnode
 
  output<-data.frame(nodes=nodenums, concordance=fraction, discordance=1-fraction)
  
  output
  
}

#functions to plot nodelabels on a given consensus tree
labelplot1<-function(shifts, minscale, matches.summary){
	for (i in 1:length(shifts[,1])){
		if (rowSums(shifts)[i] > 0){
		par(lwd=0.00001)
		vo=1
		ro=1
		v=rowSums(shifts)[i]
		nodelabels(text=rowSums(shifts)[i], node=matches.summary[,1][i], frame="circle", bg="#000000BF", col="white", cex=(ro*(v/vo)^0.25)*minscale)
		}
	}
}
labelplot<-function(shifts, minscale, color="#000000BF", size=1){
  for (i in 1:length(shifts[,1])){
    if (rowSums(shifts)[i] > 0){
      par(lwd=0.00001)
      vo=1
      ro=1
      v=rowSums(shifts)[i]
      #this line for scaled circles
      #nodelabels(text=rowSums(shifts)[i], node=(length(shifts[,1])+2:(length(shifts[,1])))[i], frame="circle", bg=color, col="white", cex=(ro*(v/vo)^0.25)*minscale)
      #this line for equal sized circles
      nodelabels(text=rowSums(shifts)[i], node=(length(shifts[,1])+2:(length(shifts[,1])))[i], frame="circle", bg=color, col="white", cex=size)
    }
  }
}
labelplot.simple<-function(shifts, minscale, matches.summary){
  for (i in 1:length(shifts[,1])){
      par(lwd=0.00001)
      nodelabels(text=matches.summary[,1][i], node=matches.summary[,1][i], frame="circle", bg="#000000BF", col="white", cex=0.4)
    }
}

#function to plot janustrees
plot_janus<-function(treepath, seed=5, xlim, ladderize=F, equal=F, width=1.0){
	tree<-read.beast.fixed(path = treepath, ladderize, equal)
	
	data<-as.data.frame(tree@data)
	models<-data$model
	
	set.seed(seed=seed)
	cols<-sample(brewer.pal(n=length(unique(models)), name="Paired"))
	
	data$color<-data$model+1
	data$color<-cols[data$color]
	data$node<-as.numeric(data$node)
	noroot<-data#[-395,]
	rownames(noroot)<-(noroot$node)
	noroot<-noroot[order(noroot$node),]
	edge.order<-as.numeric(tree@phylo$edge[,2])
	#plot.new()
	#plotphy<-plot.phylo(tree@phylo, show.tip.label=F, plot=F)
	#print(plotphy$x.lim)
	#print(plotphy$y.lim)
	#plot.window(xlim = c(0.3,0),ylim = plotphy$y.lim)
	par(lwd=1)
	plot.phylo(tree@phylo, edge.color = noroot$color[edge.order], show.tip.label=F, no.margin=F, edge.width=width, x.lim=xlim)
	axisPhylo(backward=F, padj=-0.75)
	abline(v=max(nodeHeights(tree@phylo)), lty=2)
	box()
	
}

#function to plot janustrees from s4 objects already loaded, for time tree
plot_janus_time<-function(tree_s4, seed=5, xlim, width=1.0){
  tree<-tree_s4#read.beast.fixed(path = treepath, ladderize, equal)
  tree@data<-tree@data[order(as.numeric(tree@data$node)),]
  
  data<-as.data.frame(tree@data)
  models<-data$model
  
  set.seed(seed=seed)
  cols<-sample(brewer.pal(n=length(unique(models)), name="Paired"))
  
  data$color<-data$model+1
  data$color<-cols[data$color]
  data$node<-as.numeric(data$node)
  noroot<-data#[-395,]
  rownames(noroot)<-(noroot$node)
  noroot<-noroot[order(noroot$node),]
  edge.order<-as.numeric(tree@phylo$edge[,2])
  #plot.new()
  #plotphy<-plot.phylo(tree@phylo, show.tip.label=F, plot=F)
  #print(plotphy$x.lim)
  #print(plotphy$y.lim)
  #plot.window(xlim = c(0.3,0),ylim = plotphy$y.lim)
  par(lwd=1)
  plot.phylo(tree@phylo, edge.color = noroot$color[edge.order], show.tip.label=F, no.margin=F, edge.width=width, x.lim=xlim)
  axisPhylo(backward=T, padj=-0.75)
  abline(v=max(nodeHeights(tree@phylo)), lty=2)
  box()
  
}

#function to plot janustrees from s4 objects already loaded, for time tree
plot_janus_time_radar<-function(tree_s4, seed=5, xlim, width=1.0){
  tree<-tree_s4#read.beast.fixed(path = treepath, ladderize, equal)
  tree@data<-tree@data[order(as.numeric(tree@data$node)),]
  
  data<-as.data.frame(tree@data)
  models<-data$model
  
  set.seed(seed=seed)
  cols<-sample(brewer.pal(n=length(unique(models)), name="Paired"))
  
  data$color<-data$model+1
  data$color<-cols[data$color]
  data$node<-as.numeric(data$node)
  noroot<-data#[-395,]
  rownames(noroot)<-(noroot$node)
  noroot<-noroot[order(noroot$node),]
  edge.order<-as.numeric(tree@phylo$edge[,2])
  #plot.new()
  #plotphy<-plot.phylo(tree@phylo, show.tip.label=F, plot=F)
  #print(plotphy$x.lim)
  #print(plotphy$y.lim)
  #plot.window(xlim = c(0.3,0),ylim = plotphy$y.lim)
  par(lwd=1, oma=c(0,0,0,0))
  tree@phylo$tip.label<-unlist(lapply(strsplit(tree@phylo$tip.label, split="_"), `[`, 1))
  plot.phylo(tree@phylo, edge.color = noroot$color[edge.order], show.tip.label=T, no.margin=F, edge.width=width, x.lim=xlim, cex=0.25, label.offset = 0.5)
  axisPhylo(backward=T, padj=-0.75)
  #abline(v=max(nodeHeights(tree@phylo)), lty=2)
  box()
  
}

#function to plot janustrees from s4 objects already loaded, for time tree
plot_janus_time_noaxes<-function(tree_s4, seed=5, xlim, width=1.0){
  tree<-tree_s4#read.beast.fixed(path = treepath, ladderize, equal)
  tree@data<-tree@data[order(as.numeric(tree@data$node)),]
  
  data<-as.data.frame(tree@data)
  models<-data$model
  
  set.seed(seed=seed)
  cols<-sample(brewer.pal(n=length(unique(models)), name="Paired"))
  
  data$color<-data$model+1
  data$color<-cols[data$color]
  data$node<-as.numeric(data$node)
  noroot<-data#[-395,]
  rownames(noroot)<-(noroot$node)
  noroot<-noroot[order(noroot$node),]
  edge.order<-as.numeric(tree@phylo$edge[,2])
  #plot.new()
  #plotphy<-plot.phylo(tree@phylo, show.tip.label=F, plot=F)
  #print(plotphy$x.lim)
  #print(plotphy$y.lim)
  #plot.window(xlim = c(0.3,0),ylim = plotphy$y.lim)
  par(lwd=1)
  plot.phylo(tree@phylo, edge.color = noroot$color[edge.order], show.tip.label=F, no.margin=F, edge.width=width, x.lim=xlim)
  #axisPhylo(backward=T, padj=-0.75)
  #abline(v=max(nodeHeights(tree@phylo)), lty=2)
  #box()
  
}

#function to plot janustrees from s4 objects already loaded, for time tree
plot_janus_time_noaxes_radar<-function(tree_s4, seed=5, xlim, width=1.0){
  tree<-tree_s4#read.beast.fixed(path = treepath, ladderize, equal)
  tree@data<-tree@data[order(as.numeric(tree@data$node)),]
  
  data<-as.data.frame(tree@data)
  models<-data$model
  
  set.seed(seed=seed)
  cols<-sample(brewer.pal(n=length(unique(models)), name="Paired"))
  
  data$color<-data$model+1
  data$color<-cols[data$color]
  data$node<-as.numeric(data$node)
  noroot<-data#[-395,]
  rownames(noroot)<-(noroot$node)
  noroot<-noroot[order(noroot$node),]
  edge.order<-as.numeric(tree@phylo$edge[,2])
  #plot.new()
  #plotphy<-plot.phylo(tree@phylo, show.tip.label=F, plot=F)
  #print(plotphy$x.lim)
  #print(plotphy$y.lim)
  #plot.window(xlim = c(0.3,0),ylim = plotphy$y.lim)
  par(lwd=1,oma=c(0,0,0,0))
  tree@phylo$tip.label<-unlist(lapply(strsplit(tree@phylo$tip.label, split="_"), `[`, 1))
  plot.phylo(tree@phylo, edge.color = noroot$color[edge.order], show.tip.label=T, no.margin=F, edge.width=width, x.lim=xlim, cex=0.25, label.offset = 0.5)
  #axisPhylo(backward=T, padj=-0.75)
  #abline(v=max(nodeHeights(tree@phylo)), lty=2)
  #box()
  
}

#function to plot janustrees from s4 objects already loaded, for time tree
plot_janus_time_rect<-function(reference, target){
  tmp<-make_timetree(ref = reference, tar = target@phylo)
  
  tmp2 <- target
  
  tmp2@phylo$edge.length<-tmp$phy$edge.length
  
  plot_janus_time(tree_s4 = tmp2,  xlim=c(0, 86.90747), width=0.000001)
  
  abline(v=max(nodeHeights(tmp$phy))-66, lty=2)
  rect(xleft=max(nodeHeights(tmp$phy))-66-5, ybottom=-100, xright=max(nodeHeights(tmp$phy))-66+5, ytop=250, col="#D2B48C40", border=NA)
  par(new=T)
  plot_janus_time_noaxes(tree_s4=tmp2,  xlim=c(0, 86.90747), width=1.5)
  
}

#function to plot janustrees from s4 objects already loaded, for time tree
plot_janus_time_rect_radar<-function(reference, target){
  tmp<-make_timetree(ref = reference, tar = target@phylo)
  
  tmp2 <- target
  
  tmp2@phylo$edge.length<-tmp$phy$edge.length
  
  plot_janus_time_radar(tree_s4 = tmp2,  xlim=c(-150, 95), width=0.000001)
  
  abline(v=max(nodeHeights(tmp$phy))-66, lty=2)
  rect(xleft=max(nodeHeights(tmp$phy))-66-5, ybottom=-100, xright=max(nodeHeights(tmp$phy))-66+5, ytop=250, col="#D2B48C40", border=NA)
  par(new=T)
  plot_janus_time_noaxes_radar(tree_s4=tmp2,  xlim=c(-150, 95), width=1.5)
  
}

#function to plot janustrees from simmap objects
plot_janus_time_rect_simmap<-function(simmap_object, seed=5, title, width=1, out.line=F, box=T){
  set.seed(seed)
  # Nstates<-length(unique(getStates(simmap_object)))
  # 
  # if(Nstates <=12 ){
  #   cols <- sample(brewer.pal(n=Nstates, name="Paired"))
  # } 
  # 
  # if(Nstates == 13) {
  #   cols <- c(sample(brewer.pal(n=Nstates-1, name="Paired")), sample(brewer.pal(n=Nstates-1, name="Paired")))
  # }
  # 
  # cols<-setNames(cols[1:length(unique(getStates(simmap_object)))], unique(getStates(simmap_object)))
  # 
  cols<-c("#000000", sample(rcartocolor:::carto_pal(n = 12, "Safe")))
  
  cols<-setNames(cols, unique(getStates(simmap_object)))
  
  plot.new()
  par(window)
  par(xpd=F)
  rect(xleft=max(nodeHeights(simmap_object))-66-5, ybottom=0, xright=max(nodeHeights(simmap_object))-66+5, ytop=200, col=make.transparent("#D2B48C",0.1), border=NA) #max(lastPP$yy[-(1:lastPP$Ntip)])
  par(new=T)
  #data(strat2012)
  par(xpd=F)
  data(strat2012)
  axisGeo(GTS = strat2012, unit = c("period","era"), 
          cex = 0.5, gridty = 3, col="#FAFAFA", ages=F)
  
  par(new=T)
  #par(mar=c(3.1, 0.25, 2.1, 0.25))
  plot(simmap_object, ftype="off", colors=make.transparent(cols,1), lwd=width, outline=out.line)
  if(box==T){
    box()
  }
  #title(main=title)
  #axisPhylo()
  
  #abline(v=max(nodeHeights(simmap_object))-66, lty=2)
  #par(xpd=F)
  
  
  #par(new=T, lend=0)
  #this next line uses the "color blind safe palette"
  #comment out if want to
  #col.alt<-c("#000000", rcartocolor:::carto_pal(n = 12, "Safe"))
  #plot(simmap_object, ftype="off", colors=cols, mar=c(3.1, 0.25, 2.1, 0.25), lwd=width, outline=out.line)
  
  #mar=c(3.1, 0.25, 2.1, 0.25)
  #plot(simmap_object, ftype="off", colors=rep("black",13), lwd=width, outline=out.line)

  #data$color<-data$model+1
  #data$color<-cols[data$color]
  
}

#modified axisGeo function from phyloch
axisGeo<-function (GTS, tip.time = 0, unit = c("epoch", "period"), ages = TRUE, 
                   cex = 1, col = "white", texcol = "black", gridty = 0, gridcol = "black") {
  adjustCex <- function(space, string, cex) {
    while (strwidth(string, cex = cex) >= space & cex > 0.001) cex <- cex - 
        0.001
    cex
  }
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  ntips <- lastPP$Ntip
  root <- ntips + 1
  if (lastPP$direction == "rightwards") 
    maxage <- max(lastPP$xx) + tip.time
  if (lastPP$direction == "upwards") 
    maxage <- max(lastPP$yy) + tip.time
  gts <- GTS
  maid <- grep("MA", names(gts))
  gts <- cbind(gts[, 1:maid], c(0, head(gts[, maid], -1)), 
               gts[, ((maid + 1):dim(gts)[2])])
  names(gts)[maid:(maid + 1)] <- c("fromMA", "toMA")
  if (sum(gts[1, maid:(maid + 1)]) == 0) 
    gts <- gts[-1, ]
  ind <- which(gts$fromMA <= maxage & gts$toMA >= tip.time)
  gts <- gts[c(ind, max(ind) + 1), ]
  gts$toMA[1] <- tip.time
  gts$fromMA[dim(gts)[1]] <- maxage
  par(xpd = NA)
  plotGeo <- function(gts, unit, yy) {
    id <- which(names(gts) %in% c(unit, "fromMA", "toMA"))
    gts <- gts[id]
    names(gts) <- c("unit", "from", "to")
    stages <- unique(gts$unit)
    if (length(col) == 1) 
      col <- rep(col, 2)
    col1 <- rep(col, length(stages))
    col1 <- head(col1, length(col1)/2)
    if (length(texcol) == 1) 
      texcol <- rep(texcol, 2)
    col2 <- rep(texcol, length(stages))
    col2 <- head(col2, length(col2)/2)
    xgrid <- NULL
    for (i in seq(along = stages)) {
      cat("\nStage ", i, ": ", stages[i], sep = "")
      from <- maxage - max(gts[gts$unit == stages[i], 2])
      to <- maxage - min(gts[gts$unit == stages[i], 3])
      rect(from, yy[1], to, yy[2], col = col1[i], border = "black", 
           lwd = 0.5)
      xgrid <- c(xgrid, from, to)
      en <- as.character(stages[i])
      if ((to - from) > strwidth(en, cex = cex)) 
        text(mean(c(from, to)), mean(yy), en, cex = cex, 
             col = col2[4])
      else {
        thiscex <- adjustCex(to - from, en, cex) * 0.95
        if (2 * thiscex >= cex) 
          text(mean(c(from, to)), mean(yy), en, cex = thiscex, 
               col = col2[i])
        else {
          while (nchar(en) > 0 & strwidth(en, cex = cex) >= 
                 (to - from)) en <- paste(head(unlist(strsplit(en, 
                                                               "")), -1), collapse = "")
          if (nchar(en) > 1) 
            en <- paste(paste(head(unlist(strsplit(en, 
                                                   "")), -1), collapse = ""), ".", sep = "")
          text(mean(c(from, to)), mean(yy), en, cex = cex, 
               col = col2[i])
        }
      }
    }
    xgrid
  }
  plotGeoToLeft <- function(gts, unit, yy) {
    id <- which(names(gts) %in% c(unit, "fromMA", "toMA"))
    gts <- gts[id]
    names(gts) <- c("unit", "from", "to")
    stages <- unique(gts$unit)
    if (length(col) == 1) 
      col <- rep(col, 2)
    col1 <- rep(col, length(stages))
    col1 <- head(col1, length(col1)/2)
    if (length(texcol) == 1) 
      texcol <- rep(texcol, 2)
    col2 <- rep(texcol, length(stages))
    col2 <- head(col2, length(col2)/2)
    xgrid <- NULL
    for (i in seq(along = stages)) {
      from <- maxage - max(gts[gts$unit == stages[i], 2])
      to <- maxage - min(gts[gts$unit == stages[i], 3])
      rect(yy[2], from, yy[1], to, col = col1[i], border = "black", 
           lwd = 0.5)
      xgrid <- c(xgrid, from, to)
      en <- as.character(stages[i])
      yxr <- (max(lastPP$y.lim) - min(lastPP$y.lim))/(max(lastPP$x.lim) - 
                                                        min(lastPP$x.lim)) * 1.5
      if ((to - from) > strwidth(en, cex = cex * yxr)) 
        text(mean(yy), mean(c(from, to)), en, cex = cex, 
             col = col2[4], srt = 90)
      else {
        asp <- (to - from)/yxr
        thiscex <- adjustCex(asp, en, cex) * 0.95
        if (1.5 * thiscex >= cex) 
          text(mean(yy), mean(c(from, to)), en, cex = thiscex, 
               col = col2[i], srt = 90)
        else {
          while (nchar(en) > 0 & strwidth(en, cex = cex * 
                                          yxr) >= (to - from)) en <- paste(head(unlist(strsplit(en, 
                                                                                                "")), -1), collapse = "")
          if (nchar(en) > 1) 
            en <- paste(paste(head(unlist(strsplit(en, 
                                                   "")), -1), collapse = ""), ".", sep = "")
          text(mean(yy), mean(c(from, to)), en, cex = cex, 
               col = col2[i], srt = 90)
        }
      }
    }
    xgrid
  }
  bh <- -strheight("Ap", cex = cex) * 1.5
  if (lastPP$direction == "rightwards") {
    if (ages) 
      yy <- c(bh, 2 * bh)
    else yy <- c(0, bh)
    for (j in seq_along(unit)) {
      cat("\nPlot unit:", unit[j])
      if (j == 1) 
        xgrid <- plotGeo(gts, unit[j], yy)
      else plotGeo(gts, unit[j], yy)
      yy <- yy + bh
    }
  }
  if (lastPP$direction == "upwards") {
    if (ages) 
      yy <- c(bh, 2 * bh)
    else yy <- c(0, bh)
    for (j in seq(along = unit)) {
      if (j == 1) 
        xgrid <- plotGeoToLeft(gts, unit[j], yy)
      else plotGeoToLeft(gts, unit[j], yy)
      yy <- yy + bh
    }
  }
  xgrid <- unique(sort(xgrid, decreasing = TRUE))
  label <- maxage - xgrid
  id <- TRUE
  for (k in seq(along = xgrid)) {
    if (lastPP$direction == "rightwards") 
      lines(rep(xgrid[k], 2), c(0, ntips + 1), lty = gridty, 
            col = gridcol)
    if (lastPP$direction == "upwards") 
      lines(c(0, ntips + 1), rep(xgrid[k], 2), lty = gridty, 
            col = gridcol)
    if (k < length(xgrid)) {
      spneeded <- strwidth(label[k], cex = cex * 0.8)/2
      spavailable <- xgrid[k] - xgrid[k + 1]
      if (spavailable < spneeded * 1.5) 
        id <- c(id, FALSE)
      else id <- c(id, TRUE)
    }
  }
  id <- c(id, TRUE)
  if (ages) {
    xgrid <- xgrid[id]
    label <- label[id]
    if (lastPP$direction == "rightwards") 
      text(xgrid, -0.2, round(label, digits = 1), cex = cex * 
             0.8)
    if (lastPP$direction == "upwards") 
      text(-0.2, xgrid, round(label, digits = 1), cex = cex * 
             0.8)
  }
  par(xpd = FALSE)
}

#function to plot branch length trees
compare.phylo<-function (t1, t2, limit, width) {
  if (hasArg(colors)) 
    colors <- list(...)$colors
  else colors <- sapply(c("blue", "red"), make.transparent, alpha = 0.4)
  if (hasArg(arr.colors)) 
    arr.colors <- list(...)$arr.colors
  else arr.colors <- sapply(c("blue", "red"), make.transparent, alpha = 0.7)
  h1 <- sapply(1:Ntip(t1), nodeheight, tree = t1)
  h2 <- sapply(1:Ntip(t2), nodeheight, tree = t2)
  plotTree(if (max(h1) > max(h2)) 
    t1
    else t2, plot = FALSE, direction = "rightwards", mar = c(1,0.25,1,0.25))# c(4.1, 1.1, 1.1, 1.1))
  xlim <- limit#get("last_plot.phylo", envir = .PlotPhyloEnv)$x.lim[2:1]
  par(fg = "transparent", new = TRUE)
  plotTree(t1, color = colors[1], xlim = xlim, direction = "rightwards", lwd = width, mar = c(2,0.25,1,1))#c(4.1, 1.1, 1.1, 1.1))
  T1 <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  par(fg = "black")
  axisPhylo(backward=F, padj=-0.75)
  par(fg = "transparent")
  plotTree(t2, color = colors[2], xlim = xlim, add = TRUE, direction = "rightwards", ftype = "off", lwd = width, mar = c(2,0.25,1,1))#c(4.1, 1.1, 1.1, 1.1))
  
  T2 <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  par(fg = "black")
  for (i in 1:t1$Nnode + Ntip(t1)) {
    arrows(T1$xx[i], T1$yy[i], T2$xx[i], T2$yy[i], lwd = width, col = if (T1$xx[i] < T2$xx[i]) arr.colors[2] else arr.colors[1], length = 0.1)
  }
  h <- mapply(function(x, y) if (x < y) 
    x
    else y, x = T1$xx[1:Ntip(t1)], y = T2$xx[1:Ntip(t2)])
  #text(rep(min(h), Ntip(t1)), T1$yy[1:Ntip(t1)], labels = t1$tip.label, 
  #font = 3, pos = 4, offset = 0.1 * max(c(h1, h2)))
  #for (i in 1:Ntip(t1)) lines(c(h[i] + if (h[i] > min(h)) 0.005 * 
  #diff(xlim) else 0, min(h)), rep(T1$yy[i], 2), lty = "dotted")
  abline(v=max(nodeHeights(t1)), lty=2, col="black", lwd=1)
  graphics::box(lty=1, lwd=1, col="black")
  
}

#timetree congruifier
#assumes the taxon names are the same in target+reference
#also that PATHd8 or treePL is installed on user system
make_timetree<-function(ref, tar, method ="treePL"){
  
  #in case there are precision issues
  ultrametric<-force.ultrametric(ref, method="extend")
  
  #print(is.ultrametric(ultrametric))
  
  tax<-as.data.frame(tar$tip.label[tar$tip.label %in% ultrametric$tip.label])
  rownames(tax)<-tax[,1]
  #print(tax)
  
  ultrametric_scale<-congruify.phylo(reference = ultrametric, target = tar, taxonomy = tax, scale = method)
  
  return(ultrametric_scale)
  
}

#function to read in a tree and root it as specified, and then drop the root
readrooter<-function(treepath, outgroup, drop = T){
  tmp<-read.tree(file=treepath)
  tmp<-reroot(tmp, node.number = getMRCA(phy=tmp, tip = outgroup))
  
  if(drop==T){
    tmp<-ape::drop.tip(phy=tmp, tip = outgroup)
  }else{
    tmp
  }
  
  return(tmp)
}

#function to compute substitution distance between allales for the phased trees
allele_dist<-function(tree){
  
  tips<-tree$tip.label
  ## first get the node numbers of the tips
  nodes<-sapply(tips,function(x,y) which(y==x),y=tree$tip.label)
  ## then get the edge lengths for those nodes
  edge.lengths<-setNames(tree$edge.length[sapply(nodes, function(x,y) which(y==x),y=tree$edge[,2])],names(nodes))
  #convert to data frame
  edge.lengths<-as.data.frame(edge.lengths)
  
  #sum each of the two rows to get the approximate nuc_div
  nuc_div<-rowsum(edge.lengths[,1], as.integer(gl(nrow(edge.lengths), 2, nrow(edge.lengths))))
  
  nuc_div<-as.data.frame(nuc_div)
  
  #get the row labels and clean
  rows <- as.data.frame(tips)[seq(1, nrow(as.data.frame(tips)), 2), ]
  
  #make sure all are _1 haplotype labels
  rows<-gsub(pattern="_2", replacement="_1", x=rows)
  
  #add the labels
  nuc_div$label<-rows
  
  #correct column names
  colnames(nuc_div)<-c("nuc_div", "label")
  
  #return the data frame
  return(nuc_div)
  
}

#function to directly compute nuc diversit and variance from the sequences
#intput is list of dnabin objects split into pairs of alleles
nuc_div_pairs<-function(pairs){
  tmp<-list()
  
  for(i in 1:length(pairs)){
    tmp[[i]]<-nuc.div(pairs[[i]], variance=T)
  }
  
  names(tmp)<-names(pairs)
  return(tmp)
}

#function to graft tips for phylogenetic logistic regression
grafter<-function(tree, edge.length=0.0, position=0.0){
  tmp<-tree
  labels<-(seq(from=length(tmp$tip.label), length.out=tmp$Nnode)+1)
  
  for(i in 1:length(seq(from=length(tmp$tip.label), length.out=tmp$Nnode)+1)){
    tmp<-bind.tip(
      tree = tmp,
      tip.label = labels[i],
      edge.length = edge.length,
      where = (seq(from=length(tmp$tip.label), length.out=tmp$Nnode)+1)[i],
      position = position)
  }
  
  tmp
  
}

#convert janus tree to SIMMAP tree
janus_simmap<-function(reference, target){
  #generate a time tree from the reference and target
  tmp<-make_timetree(ref = reference, tar = target@phylo)
  tmp2 <- target
  
  #swap branch lengths from the congruified tree onto the
  #s4 time tree object
  tmp2@phylo$edge.length<-tmp$phy$edge.length
  
  #save the janus data output
  data<-as.data.frame(tmp2@data)
  models<-data$model
  
  #cols represents the values of the model IDs, not colors (old history)
  cols<-sort(unique(data$model))
  #cols<-unique(data$model[c((length(tmp2@phylo$tip.label)+1):(length(tmp2@phylo$tip.label)+length(tmp2@phylo$tip.label)-1), 1:length(tmp2@phylo$tip.label))])
  
  #tmp3<-read.beast.fixed2(path='~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/5d81a027/regular/ALL_MFP_MERGE_MRL3_constraint_RM_UE_UL_M4/ALL_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre', ladderize = F, equal =F)
  #data<-as.data.frame(tmp3@data)
  #data<-data[dim(data)[1]:1,]
  #rownames(data)<-rev(rownames(data))
  
  #reverse order
  #models<-data$model
  #cols<-unique(data$model)
  
  data$color<-data$model+1 #indexing issue because of the zero labeled model
  data$color<-cols[data$color]
  data$node<-as.numeric(data$node)
  noroot<-data#[-395,]
  rownames(noroot)<-(noroot$node)
  noroot<-noroot[order(noroot$node),]
  edge.order<-as.numeric(tmp2@phylo$edge[,2])
  
  #states is a vector of states to assign to all branches in the right order
  #now converted to english words to avoid ambiguity later on
  states<-as.character(as.english(noroot$color[edge.order]))
  
  #get the branch lengths because each state occurs for the whole branch length
  times<-tmp2@phylo$edge.length
  times<-as.list(times)
  
  #assign the model state to the list of branch lengths
  for(i in 1:length(times)){
    names(times[[i]])<-states[i]
  }
  
  sim<-tmp2@phylo
  sim$maps<-times
  
  
  #generate initialized mapped.edge element and assign correct names
  sim$mapped.edge<-matrix(0,nrow(sim$edge),length(unique(states)),dimnames=list(paste(sim$edge[,1],",",sim$edge[,2],sep="")))
  colnames(sim$mapped.edge)<-unique(states)
  
  #copy the maps elements
  mapping<-sim$maps
  
  #duplicate the maps element in the right order
  mapping2<-do.call("cbind", rep(list(unlist(mapping)), length(unique(names(unlist(mapping))))))
  colnames(mapping2)<-unique(names(unlist(mapping)))
  
  #zero out the values of the mapped.edge element that represent
  #states other than the one occupying the entire branch
  for(i in 1:length(unique(names(unlist(mapping))))){
    mapping2[,i][which(!names(mapping2[,i])==unique(names(unlist(mapping)))[i])]<-0
  }
  
  #overwrite the initialized mapped.edge element
  sim$mapped.edge[,c(1:length(unique(names(unlist(mapping)))))]<-mapping2
  
  #rename the attributes of the mapped edge element
  names(dimnames(sim$mapped.edge))<-c("edge","state")
  
  #resort the mapped.edge to follow the map state order
  #sim$mapped.edge<-test$mapped.edge[,names(sort(to_number(colnames(sim$mapped.edge))))]
  
  #add the class definition
  class(sim)<-c("simmap", "phylo")
  return(sim)
  
}

#function for plotting OUwie simulations 
sim_plotter <-
  function(janustree,
           ouwiefit,
           clrs,
           nsim,
           extra,
           width,
           plot = T) {
    H <- nodeHeights(janustree)
    oum.sims <- OUwie.sim(fitted.object = ouwiefit, get.all = T)
    nms <- coalesce(oum.sims$Genus_species, rownames(oum.sims))
    simdat <- setNames(as.numeric(oum.sims$X), nms)
    par(xaxt = "n",
        yaxt = "n",
        mar = c(5.1, 5.1, 2.1, 1.1))
    #plot(simmap.janus.alldata, lwd=0.00001, ftype="off")
    phenogram(
      janustree,
      x = simdat,
      colors = clrs,
      add = F,
      lwd = width,
      spread.labels = F,
      ftype = "off",
      ylim = c(-4, 4)
    )
    par(xaxt = "s",
        yaxt = "s",
        font.lab = 4)
    axis(
      1,
      at = seq(
        from = max(H),
        to = -10,
        by = -20
      ),
      labels = seq(
        from = 0,
        to = max(H) + 10,
        by = 20
      )
    )
    axis(2)
    
    oum.sim <- list()
    for (i in 1:nsim) {
      oum.sim[[i]] <- OUwie.sim(fitted.object = ouwiefit, get.all = T)
      print(paste('simulation', i/nsim, "complete"))
      nms <-
        coalesce(oum.sim[[i]]$Genus_species, rownames(oum.sim[[i]]))
      simdat <- setNames(as.numeric(oum.sim[[i]]$X), nms)
      if (plot == T) {
        phenogram(
          janustree,
          x = simdat,
          colors = clrs,
          add = TRUE,
          lwd = width,
          spread.labels = F
        )
      }
    }
    
    abline(v = max(H) - 66.02, lty = 2)
    
    tmp <- data.frame(Reg = getStates(janustree, type = "tips"))
    for (i in 1:length(oum.sim)) {
      assign(paste("sim", i, sep = ""), oum.sim[[i]][, c(3)][1:198, ])
      tmp <- cbind(tmp,
                   setNames(as.data.frame(eval(
                     parse(text = paste("sim", i, sep = ""))
                   )),
                   paste("sim", i, sep = "")))
      #tmp[[i]]<-sims[[i]][,c(2,3)][1:198,]
    }
    
    #tmp <- lapply(tmp %>% group_by(., Reg) %>%group_split(.), `[`, -c(1))
    #tmp <- lapply(tmp, unlist)
    #tmp <- lapply(tmp, unname)
    
    return(tmp)
  }

#estimate mode from unimodal data 
estimate_mode <- function(x, ...) {
  d <- density(x, bw=0.7)
  mode<-d$x[which.max(d$y)]
  #int<-HDInterval::hdi(d)
  return(c(mode))
}

#checking how many possible nodes are candidates for janus
count_candidates<-function(tree){
  
  x = reorder(tree, "postorder")
  x$edge
  
  ntips = length(x$tip.label)
  m = max(x$edge)
  res = integer(m)
  res[1:ntips] = 1L
  parent = x$edge[,1]
  child = x$edge[,2]
  for(i in 1:nrow(x$edge)){
    res[parent[i]] = res[parent[i]] + res[child[i]]  
  }
  res
  
  
  
}

#"object" is a named charater vector
#"names" is a vector of names from a translation table
#which matches the names of object (not necessarily same order)
#swap is the vector of names from the translation table (a different column)
translation_nameswap<-function(object, names, swap){
  for (i in 1:length(names(object))){
    index<-match(names(object)[i],names)
    names(object)[i]<-swap[index]
  }
  return(object)
}

#same as above but for the readSet function
translation_nameswap_set<-function(object, names, swap){
  for (i in 1:length(object@ranges@NAMES)){
    index<-match(object@ranges@NAMES[i],names)
    object@ranges@NAMES[i]<-swap[index]
  }
  return(object)
}

#get stem ages given MRCAs
#from https://grokbase.com/t/r/r-sig-phylo/142dayv9mk/easiest-way-to-get-the-stem-age-of-a-clade
stem_age<-function(tree, nodenumber=NULL, all=T){
  
  if (all==T){
    tips<-1:length(tree$tip.label)
    nodes<-seq((from=max(tips)+1), length.out=tree$Nnode)
    allnodes<-c(tips,nodes)
    #print(allnodes)
    
    H<-nodeHeights(tree)
    H<-max(H)-H
    ## tips are the species in the clade, or a subset of the species
    ## such that the MRCA of tips is the MRCA of the clade
    
    output<-list()
    for(i in 1:length(allnodes)){
      nn<-allnodes[i]
      output[i]<-H[tree$edge==phytools:::getAncestors(tree,nn,"parent")][1]
    }
    output<-unlist(output)
    names<-tree$tip.label
    
    output<-setNames(output, c(names, nodes))
    
    return(output)
    
  } else{
    H<-nodeHeights(tree)
    H<-max(H)-H
    ## tips are the species in the clade, or a subset of the species
    ## such that the MRCA of tips is the MRCA of the clade
    nn<-nodenumber
    h<-H[tree$edge==phytools:::getAncestors(tree,nn,"parent")][1]
  }
  return(h)
}
stem_node<-function(tree, nodenumber=NULL, all=T){
  
  if (all==T){
    tips<-1:length(tree$tip.label)
    nodes<-seq((from=max(tips)+1), length.out=tree$Nnode)
    allnodes<-c(tips,nodes)
    #print(allnodes)
    
    H<-nodeHeights(tree)
    H<-max(H)-H
    ## tips are the species in the clade, or a subset of the species
    ## such that the MRCA of tips is the MRCA of the clade
    
    output<-list()
    for(i in 1:length(allnodes)){
      nn<-allnodes[i]
      output[i]<-H[tree$edge==phytools:::getAncestors(tree,nn,"parent")][1]
    }
    output<-unlist(output)
    names<-tree$tip.label
    
    output<-setNames(output, c(names, nodes))
    
    return(output)
    
  } else{
    H<-nodeHeights(tree)
    H<-max(H)-H
    ## tips are the species in the clade, or a subset of the species
    ## such that the MRCA of tips is the MRCA of the clade
    nn<-nodenumber
    h<-H[tree$edge==phytools:::getAncestors(tree,nn,"parent")][1]
  }
  return(h)
}

#https://stackoverflow.com/a/25629034
getphylo_x <- function(input, number) {
  if(is.character(number)) {
    number <- which(c(input$tip.label, input$number.label)==number)
  }
  pi <- input$edge[input$edge[,2]==number, 1]
  if (length(pi)) {
    ei<-which(input$edge[,1]==pi & input$edge[,2]==number)
    input$edge.length[ei] + Recall(input, pi)
  } else {
    if(!is.null(input$root.edge)) {
      input$root.edge
    } else {
      0
    }
  }
}
getphylo_y <- function(input, number) {
  if(is.character(number)) {
    number <- which(c(input$tip.label, input$number.label)==number)
  }
  ci <- input$edge[input$edge[,1]==number, 2]
  if (length(ci)==2) {
    mean(c(Recall(input, ci[1]), Recall(input, ci[2])))
  } else if (length(ci)==0) {
    Ntip <- length(input$tip.label)
    which(input$edge[input$edge[, 2] <= Ntip, 2] == number)
  } else {
    stop(paste("error", length(ci)))
  }
}
getcoords<-function(tree, node){

  result<-list()
  result[1]<- getphylo_x(tree, node)
  result[2]<- getphylo_y(tree,node)
  
  return(unlist(result))
}

#quick transformation functions
logitTransform <- function(p) { log(p/(1-p)) }
asinTransform <- function(p) { asin(sqrt(p)) }

#loop ouwie (testing)
OUwie_looper<-function(maps, modlist, df, name, n.cores){
  
  # my.cluster <- parallel::makeCluster(
  #   n.cores, 
  #   type = "PSOCK"
  # )
  # 
  # doParallel::registerDoParallel(cl = my.cluster)
  # 
  output<-list()#rep(NA, length(maps))
  
 # output<-foreach(i = 1:length(maps), .combine='c') %dopar% {OUwie(phy = maps[[i]], data=data.frame(names=rownames(df),regime=getStates(maps[[i]], type = "tips"), trait=get(name,df)), simmap.tree = T, model=mod, diagn=T) }
  
  #parallel::stopCluster(cl = my.cluster)
  
  for(i in 1:length(maps)){
    output[[i]]<-OUwie(phy = maps[[i]], data=data.frame(names=rownames(df),regime=getStates(maps[[i]], type = "tips"), trait=get(name,df)), simmap.tree = T, model=modlist[[i]], diagn=F)
  }
  names(output)<-names(maps)
  return(output)
}

#clean the NaNs
#https://stackoverflow.com/questions/52490552/r-convert-nan-to-na/52490634
is.nan.data.frame <- function(x) {
  do.call(cbind, lapply(x, is.nan))
}

#http://blog.phytools.org/2015/08/removing-node-labels-from-newick-string.html
strip.nodelabels<-function(text){
  obj<-strsplit(text,"")[[1]]
  cp<-grep(")",obj)
  csc<-c(grep(":",obj),length(obj))
  exc<-cbind(cp,sapply(cp,function(x,y) y[which(y>x)[1]],y=csc))
  exc<-exc[(exc[,2]-exc[,1])>1,]
  inc<-rep(TRUE,length(obj))
  if(nrow(exc)>0) for(i in 1:nrow(exc)) 
    inc[(exc[i,1]+1):(exc[i,2]-1)]<-FALSE
  paste(obj[inc],collapse="")
}

#modified pgls plotting
pGLS.plotGrade<- function (Yvar, Xvar, data, tree, model, group, linecol, linelwd, pchlwd, pchcol, pchbg,...) {
  dataTemp <- pruneSample(na.omit(data), tree, group)$data
  treeTemp <- pruneSample(dataTemp, tree, group)$tree
  Y <- dataTemp[, which(colnames(dataTemp) == paste(Yvar))]
  X <- dataTemp[, which(colnames(dataTemp) == paste(Xvar))]
  dataGLS <- as.data.frame(cbind(Y, X))
  rownames(dataGLS) <- rownames(dataTemp)
  switch(model, BM = {
    pGLSTemp <- nlme::gls(Y ~ X, dataGLS, correlation = corBrownian(phy = treeTemp))
  }, lambda = {
    pGLSTemp <- nlme::gls(Y ~ X, dataGLS, correlation = corPagel(1, 
                                                           phy = treeTemp, fixed = FALSE))
  })
  a <- summary(pGLSTemp)$tTable[1, 1]
  b <- summary(pGLSTemp)$tTable[2, 1]
  lines(c(min(X), max(X)), c((a + b * min(X)), (a + b * max(X))), 
        col=linecol, lwd=linelwd,...)
  points(X, Y, col=pchcol, bg=pchbg, lwd=pchlwd,...)
}

#identify edge indices associated with minimum clade size
edge_indices_N<-function(tree, min){
  
  pruned<-tree
  
  theNodes <- length(pruned$tip.label)+1:pruned$Nnode
  results <- c()
  for(i in 1:length(theNodes))
  {
    temp <- extract.clade(pruned, theNodes[i])
    results[i] <- length(temp$tip.label)
  }
  
  names(results) <- theNodes

  edgetable<-pruned$edge
  rownames(edgetable)<-seq(1:length(edgetable[,2]))
  
  as.numeric(rownames(edgetable)[as.numeric(edgetable[,2]) %in% as.numeric(names(results[results >= min]))])
  
}

#identify node indices associated with minimum clade size
node_indices_N<-function(tree, min){
  
  pruned<-tree
  
  theNodes <- length(pruned$tip.label)+1:pruned$Nnode
  results <- c()
  for(i in 1:length(theNodes))
  {
    temp <- extract.clade(pruned, theNodes[i])
    results[i] <- length(temp$tip.label)
  }
  
  names(results) <- theNodes
  
  #results<-results[results>=min]
  #results<-as.numeric(names(results))
  
  return(results)

  #edgetable<-pruned$edge
  #rownames(edgetable)<-seq(1:length(edgetable[,2]))
  
  #as.numeric(rownames(edgetable)[as.numeric(edgetable[,2]) %in% as.numeric(names(results[results >= min]))])
  
}

#identify edge indices associated with particular node numbers
edge_indices_nodes<-function(tree, nodes){
  
  pruned<-tree
  
  edgetable<-pruned$edge
  rownames(edgetable)<-seq(1:length(edgetable[,2]))
  
  as.numeric(rownames(edgetable)[as.numeric(edgetable[,2]) %in% as.numeric(nodes)])
  
}

#identify descendant node number for a given edge number
node_indices_edge<-function(tree, edges){
  edgetable<-tree$edge
  result<-list()
  
  for(i in 1:length(edges)){
    result[i]<- edgetable[edges[i],][2]
  }
  return(unlist(result))
}

#get tip labels descendended from edge number
tips_from_edge<-function(tree, edge){
  node<-node_indices_edge(tree=tree, edges=edge)
  descendants<-getDescendants(tree=tree, node=node)
  tips<-tree$tip.label[descendants]
  tips<-tips[!is.na(tips)]
  return(tips)
  
}

#plot curve
#https://stackoverflow.com/a/66478714
plot_logistic_curve = function(log_mod){
  mod_frame = model.frame(log_mod)
  var_names = names(mod_frame)
  newdat = setNames(data.frame(seq(min(mod_frame[[2]]), max(mod_frame[[2]]), len=100)), var_names[2])
  newdat[var_names[1]] = predict(log_mod, newdata = newdat, type="response")
  plot(mod_frame[[1]] ~ mod_frame[[2]], col = "red4", xlab = var_names[[2]], ylab = var_names[[1]])
  lines(newdat[[var_names[2]]], newdat[[var_names[1]]], col = "green4", lwd = 2)
} 

#function to loop l1ou on n null datasets based on some input vcv
l1ou_nullboot<-function(seed=5, input_tree, boots, input_vcv, ncores, IC, considered_edges, threshold, max.nshifts, timelim){
  set.seed(seed)
  output<-list()
  
 
  for(i in 1:boots){
    input <- adjust_data(input_tree, simRatematrix(tree = input_tree, vcv = input_vcv))
    print("processed input")
    withTimeout(expr={
    tryCatch(output[[i]]<-eModel.par.unconstrained.aic.null <- estimate_shift_configuration(input$tree,
                                                                                   input$Y,
                                                                                   nCores=ncores,
                                                                                   quiet=T,
                                                                                   candid.edges = c(considered_edges),
                                                                                   edge.length.threshold = threshold,
                                                                                   max.nShifts = max.nshifts,
                                                                                   criterion=IC), error=function(e){}) }
    , onTimeout="warning", 
    timeout = timelim,
    elapsed = timelim
    )
    print(paste(i, "out of", boots,"bootstrap reps complete"))
  }
  
  return(output)
  
}
#process nullboot output
l1ou_nullboot_process<-function(boots_output){
  boots.filter<-boots_output[!unlist(lapply(boots_output, is.null))]
  print(paste(length(boots.filter), "samples remain after filtering NAs" ))
  result<-table(unlist(lapply(lapply(boots.filter, function(e){e$shift.configuration}), unname)))/length(boots.filter)
  result<-unlist(as.list(result))
  result$n <- length(boots.filter)
  result<-unlist(result)
  return(result)
}

#function to loop l1ou on n null datasets based on some input vcv
l1ou_nullboot.mvmorph<-function(seed=5, input_tree, boots, model_fit, ncores, IC, considered_edges, threshold, max.nshifts, timelim){
  set.seed(seed)
  output<-list()
  
  
  for(i in 1:boots){
    input <- adjust_data(input_tree, mvSIM(tree=input_tree, nsim=1, model="BM1", param=model_fit ))
    print("processed input")
    withTimeout(expr={
      tryCatch(output[[i]]<-eModel.par.unconstrained.aic.null <- estimate_shift_configuration(input$tree,
                                                                                              input$Y,
                                                                                              nCores=ncores,
                                                                                              quiet=T,
                                                                                              candid.edges = c(considered_edges),
                                                                                              edge.length.threshold = threshold,
                                                                                              max.nShifts = max.nshifts,
                                                                                              criterion=IC), error=function(e){}) }
      , onTimeout="warning", 
      timeout = timelim,
      elapsed = timelim
    )
    print(paste(i, "out of", boots,"bootstrap reps complete"))
  }
  
  return(output)
  
}

#function for looping fisher's test
fisherloop<-function(input_table, alt = "greater"){
  data<-input_table[,-ncol(input_table)]
  names<-colnames(data)
  names<-colnames(data)
  data<-data*input_table[,"n"]
  data<-as.data.frame(data)
  data<-as.data.frame(lapply(data, as.integer))
  
  output<-list()
  #print(data)
  
  for( i in 1:length( data[1,] )){
    #print(rbind(c(data[1,i], input_table[,"n"][1]-data[1,i]), c(data[2,i], input_table[,"n"][2]-data[2,i])))
    output[[i]]<-fisher.test(rbind(c(data[1,i], input_table[,"n"][1]-data[1,i]), c(data[2,i], input_table[,"n"][2]-data[2,i])), 
                             alternative = alt)[c(1,2,3)]
    
  }
  
  output <- unlist(output)
  
  output <- round(output,5)
  
  output <-  as.data.frame(t(data.frame(split(output, ceiling(seq_along(output)/4)))))
  
  rownames(output)<-names#colnames(input_table)[1:length(names(input_table))-1]
  
  #output<-split(output, f=3)
  return(output)
  #return(unlist(setNames(output, names)))
  
}

#https://github.com/mrhelmus/phylogeny_manipulation/blob/master/AIC_func.r
vif<-function(mod, ...){
  UseMethod("vif")
}

#https://github.com/mrhelmus/phylogeny_manipulation/blob/master/AIC_func.r
vif.phyloglm <- function(mod, ...) {
  if (any(is.na(coef(mod))))
    stop ("there are aliased coefficients in the model")
  v <- vcov(mod)
  assign <- attributes(mod$X)$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else {warning("No intercept: vifs may not be sensible.")}
  terms <- colnames(mod$X)[colnames(mod$X)!="(Intercept)"]
  n.terms <- length(terms)
  if (n.terms < 2) stop("model contains fewer than 2 terms")
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) *
      det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {result <- result[, 1]
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  return(result)
}

#moved this to phylolm package temporarily so that modified phylolm function can use it
#moved to phylolm package directly bc of issues with future package
# #https://stat.ethz.ch/pipermail/r-help/2008-August/172323.html
# mode <- function(data) {
#   # Function for mode estimation of a continuous variable
#   # Kernel density estimation by Ted Harding & Douglas Bates (found on
#   #RSiteSearch)	
# 
#     x<-as.numeric(data)
#     lim.inf=min(x)-1; lim.sup=max(x)+1
#     
#     #hist(x,freq=FALSE,breaks=seq(lim.inf,lim.sup,0.2))
#     s<-density(x,from=lim.inf,to=lim.sup,bw=0.2)
#     n<-length(s$y)
#     v1<-s$y[1:(n-2)];
#     v2<-s$y[2:(n-1)];
#     v3<-s$y[3:n]
#     ix<-1+which((v1<v2)&(v2>v3))
#     
#     #lines(s$x,s$y,col="red")
#     #points(s$x[ix],s$y[ix],col="blue")
#     
#     md <- s$x[which(s$y==max(s$y))] 
#     
#     return(md)
# }

#moved to phylolm package directly bc of issues with future package
# #modified phyloglm function to add medians from bootstrapping
# phyloglm.mod<-function (formula, data = list(), phy, method = c("logistic_MPLE", 
#                                                             "logistic_IG10", "poisson_GEE"), btol = 10, log.alpha.bound = 4, 
#                     start.beta = NULL, start.alpha = NULL, boot = 0, full.matrix = TRUE) {
#   if (!inherits(phy, "phylo")) 
#     stop("object \"phy\" is not of class \"phylo\".")
#   if (is.null(phy$edge.length)) 
#     stop("the tree has no branch lengths.")
#   if (is.null(phy$tip.label)) 
#     stop("the tree has no tip labels.")
#   method = match.arg(method)
#   mf = model.frame(formula = formula, data = data)
#   if (is.null(rownames(mf))) {
#     if (nrow(mf) != length(phy$tip.label)) 
#       stop("the number of rows in the data does not match the number of tips in the tree.")
#     warning("the data has no names, order assumed to be the same as tip labels in the tree.\n")
#   }
#   else {
#     taxa_without_data = setdiff(phy$tip.label, rownames(mf))
#     if (length(taxa_without_data) > 0) {
#       warning("will drop from the tree ", length(taxa_without_data), 
#               " taxa with missing data")
#       phy = drop.tip(phy, taxa_without_data)
#     }
#     if (length(phy$tip.label) < 2) 
#       stop("only 0 or 1 leaf with data on all variables: not enough.")
#     taxa_notin_tree = setdiff(rownames(mf), phy$tip.label)
#     if (length(taxa_notin_tree) > 0) {
#       warning(length(taxa_notin_tree), " taxa not in the tree: their data will be ignored")
#       mf = mf[-which(rownames(mf) %in% taxa_notin_tree), 
#               , drop = F]
#     }
#     ordr = match(phy$tip.label, rownames(mf))
#     if (any(is.na(ordr))) 
#       stop("the row names in the data do not match the tip labels in the tree.\n")
#     mf = mf[ordr, , drop = F]
#   }
#   X = model.matrix(attr(mf, "terms"), data = mf)
#   y = model.response(mf)
#   dk = ncol(X)
#   phy = reorder(phy, "pruningwise")
#   original.edge.length = phy$edge.length
#   n <- length(phy$tip.label)
#   N <- dim(phy$edge)[1]
#   ROOT <- n + 1L
#   anc <- phy$edge[, 1]
#   des <- phy$edge[, 2]
#   dis = pruningwise.distFromRoot(phy)
#   if (method %in% c("logistic_MPLE", "logistic_IG10")) {
#     if (sum(!(y %in% c(0, 1)))) 
#       stop("The model by Ives and Garland requires a binary (0 or 1) response (dependent variable).")
#     y = as.numeric(as.factor(y)) - 1
#     if (all(duplicated(y)[-1L])) 
#       stop("the response (dependent variable) is always 0 or always 1.")
#     btouch = 0
#     proposedBetaSD = 0.05
#     D = max(dis[1:n]) - dis[1:n]
#     D = D - mean(D)
#     externalEdge = (des <= n)
#     phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] + 
#       D[des[externalEdge]]
#     times <- pruningwise.branching.times(phy)
#     names(times) <- (n + 1):(n + phy$Nnode)
#     Tmax <- max(times)
#     intern = which(phy$edge[, 2] > n)
#     lok = rep(-1, N)
#     lok[intern] = des[intern] - n
#   }
#   if (method == "poisson_GEE") {
#     if ((!isTRUE(all(y == floor(y))))) 
#       stop("The Poisson regression requires an integer response (dependent variable).")
#     if (sum(y < 0)) 
#       stop("The Poisson regression requires a positive response (dependent variable).")
#   }
#   transf.branch.lengths_poisson_GEE <- function(beta) {
#     if (dk > 1) 
#       g = X %*% beta
#     else g = rep(1, n) * beta
#     mu = as.vector(exp(g))
#     root.edge = 0
#     diag = sqrt(mu/dis[1:n])
#     edge.length = phy$edge.length
#     return(list(edge.length, root.edge, diag))
#   }
#   transf.branch.lengths <- function(B, lL) {
#     if (dk > 1) 
#       g = X %*% B
#     else g = rep(1, n) * B
#     mu = as.vector(1/(1 + exp(-g)))
#     p = mean(mu)
#     alpha = 1/exp(lL)
#     edge.length = numeric(N)
#     distFromRoot <- exp(-2 * alpha * times)
#     tmp = .C("transbranchlengths_IvesGarland2010", as.integer(N), 
#              as.integer(des), as.integer(anc - n), as.integer(lok), 
#              as.double(distFromRoot), as.integer(externalEdge), 
#              as.double(mu), as.double(p), as.double(alpha), as.double(D), 
#              el = as.double(1:N), di = as.double(1:n))
#     edge.length = tmp$el
#     diag = tmp$di
#     root.edge = min(distFromRoot)
#     if (any(is.nan(edge.length))) 
#       stop("edge.length[i] is NaN. Please reduce btol and/or log.alpha.bound.")
#     return(list(edge.length, root.edge, diag))
#   }
#   three.point.compute <- function(trans, y, X) {
#     ole = 4 + 2 * dk + dk * dk
#     tmp = .C("threepoint", as.integer(N), as.integer(n), 
#              as.integer(phy$Nnode), as.integer(1), as.integer(dk), 
#              as.integer(ROOT), as.double(trans[[2]]), as.double(trans[[1]]), 
#              as.integer(des), as.integer(anc), as.double(as.vector(y)), 
#              as.double(as.vector(X)), result = double(ole))$result
#     return(list(vec11 = tmp[2], y1 = tmp[3], yy = tmp[4], 
#                 X1 = tmp[5:(4 + dk)], XX = matrix(tmp[(5 + dk):(ole - 
#                                                                   dk)], dk, dk), Xy = tmp[(ole - dk + 1):ole], 
#                 logd = tmp[1]))
#   }
#   plogregfunct <- function(startB, startlL, y) {
#     convergeflag = 0
#     clL = startlL
#     cB = startB
#     diflL = 100
#     difB = 100
#     counter = 0
#     ttozero = 10^6
#     optss <- list(reltol = .Machine$double.eps^0.5, maxit = 1e+05, 
#                   parscale = 1)
#     while (((diflL > 10^-6) | (difB > 10^-6) | (ttozero > 
#                                                 10^-1)) & (counter < 20)) {
#       counter = counter + 1
#       oldlL = clL
#       oldB = cB
#       olddiflL = diflL
#       olddifB = difB
#       opt <- optim(par = clL, fn = function(par) {
#         plogreglLfunct(cB, par, y)
#       }, method = "L-BFGS-B")
#       clL = as.numeric(opt$par)
#       diflL = (clL - oldlL)^2
#       if (counter >= 10) 
#         clL = (clL + oldlL)/2
#       opt <- optim(par = cB, fn = function(par) {
#         plogregBfunct(par, clL, y)
#       }, method = "L-BFGS-B", control = list(factr = 1e+12))
#       cB = as.vector(opt$par)
#       ttozero = as.numeric(opt$value)
#       if (ttozero > 10^-2) {
#         Btemp = rnorm(dk, startB, proposedBetaSD * pmax(abs(startB), 
#                                                         rep(0.1, dk)))
#         opt <- optim(par = Btemp, fn = function(par) {
#           plogregBfunct(par, clL, y)
#         }, method = "L-BFGS-B", control = list(factr = 1e+12))
#         Btemp = as.vector(opt$par)
#         newttozero = as.numeric(opt$value)
#         if (newttozero < ttozero) {
#           cB = Btemp
#           ttozero = newttozero
#         }
#       }
#       difB = sum((cB - oldB) * (cB - oldB))
#       if (counter >= 10) 
#         cB = (cB + oldB)/2
#     }
#     if (counter >= 19) 
#       if ((max(abs(c(oldlL - clL, oldB - cB))) > 0.1) | 
#           (ttozero > 0.5)) 
#         convergeflag = 1
#     return(list(B = cB, lL = clL, convergeflag = convergeflag))
#   }
#   plogregBfunct <- function(B, lL, y) {
#     if (dk > 1) 
#       g = X %*% B
#     else g = rep(1, n) * B
#     if (any(abs(g) >= btol)) {
#       btouch <<- 1
#       return(1e+06)
#     }
#     mu = as.vector(1/(1 + exp(-g)))
#     temp = transf.branch.lengths(B, lL)
#     dia = temp[[3]]
#     comp = three.point.compute(temp[1:2], (y - mu)/dia, mu * 
#                                  (1 - mu) * X/dia)
#     logdetC = comp$logd + 2 * sum(log(dia)) - sum(log(mu * 
#                                                         (1 - mu)))
#     if (logdetC < -100 * log(10)) 
#       return(1e+06)
#     Z = comp$Xy
#     if (dk == 1) 
#       FirthC = (1 - 2 * mu)/2
#     else {
#       Dx = 0.1
#       infoM = comp$XX
#       invInfoM = solve(infoM)
#       FirthC = rep(NA, dk)
#       for (i in 1:dk) {
#         dB = B
#         dB[i] = dB[i] + Dx
#         g = X %*% dB
#         if (any(abs(g) >= btol)) 
#           return(1e+06)
#         mu = as.vector(1/(1 + exp(-g)))
#         ttemp = transf.branch.lengths(dB, lL)
#         tdiag = ttemp[[3]]
#         tcomp = three.point.compute(ttemp[1:2], (y - 
#                                                    mu)/tdiag, mu * (1 - mu) * X/tdiag)
#         dinfoMp = tcomp$XX
#         dB = B
#         dB[i] = dB[i] - Dx
#         g = X %*% dB
#         if (any(abs(g) >= btol)) 
#           return(1e+06)
#         mu = as.vector(1/(1 + exp(-g)))
#         ttemp = transf.branch.lengths(dB, lL)
#         tdiag = ttemp[[3]]
#         tcomp = three.point.compute(ttemp[1:2], (y - 
#                                                    mu)/tdiag, mu * (1 - mu) * X/tdiag)
#         dinfoMm = tcomp$XX
#         DinfoM = (dinfoMp - dinfoMm)/Dx/2
#         FirthC[i] = sum(diag(invInfoM %*% DinfoM))/2
#       }
#     }
#     tozero = Z + FirthC
#     return(sum(tozero^2))
#   }
#   plogreglLfunct <- function(B, lL, y) {
#     g = X %*% B
#     mu = as.vector(1/(1 + exp(-g)))
#     if (abs(lL - log(Tmax)) >= log.alpha.bound) 
#       return(1e+10)
#     temp = transf.branch.lengths(B, lL)
#     dia = temp[[3]]
#     comp = three.point.compute(temp[1:2], (y - mu)/dia, mu * 
#                                  (1 - mu) * X/dia)
#     LL = (comp$logd + 2 * sum(log(dia)) + comp$yy)/2
#     if (!is.finite(LL)) 
#       LL = 1e+10
#     return(LL)
#   }
#   plogregBSEfunct <- function(B, lL) {
#     g = X %*% B
#     mu = as.vector(1/(1 + exp(-g)))
#     temp = transf.branch.lengths(B, lL)
#     dia = temp[[3]]
#     comp = three.point.compute(temp[1:2], (y - mu)/dia, mu * 
#                                  (1 - mu) * X/dia)
#     infoM = comp$XX
#     covBSE = solve(infoM)
#     BSE = sqrt(diag(covBSE))
#     return(list(BSE = BSE, covBSE = covBSE, info = infoM))
#   }
#   npllh <- function(par, y) {
#     if (abs(par[dk + 1] - log(Tmax)) >= log.alpha.bound) 
#       return(1e+10)
#     g = X %*% par[1:dk]
#     if (any(abs(g) >= btol)) {
#       btouch <<- 1
#       return(1e+10)
#     }
#     mu = as.vector(1/(1 + exp(-g)))
#     temp = transf.branch.lengths(par[1:dk], par[dk + 1])
#     dia = temp[[3]]
#     comp = three.point.compute(temp[1:2], numeric(n), mu * 
#                                  (1 - mu) * X/dia)
#     infoM = comp$XX
#     llk <- .C("logistreglikelihood", as.integer(N), as.integer(n), 
#               as.integer(phy$Nnode), as.integer(ROOT), as.double(original.edge.length), 
#               as.integer(des), as.integer(anc), as.integer(as.vector(y)), 
#               as.double(as.vector(mu)), as.integer(dk), as.double(exp(-par[dk + 
#                                                                              1])), loglik = double(1))$loglik
#     if (dk == 1) 
#       pllik = llk + log(abs(infoM))/2
#     else pllik = llk + log(det(infoM))/2
#     -pllik
#   }
#   llh <- function(mu, alpha) {
#     .C("logistreglikelihood", as.integer(N), as.integer(n), 
#        as.integer(phy$Nnode), as.integer(ROOT), as.double(original.edge.length), 
#        as.integer(des), as.integer(anc), as.integer(as.vector(y)), 
#        as.double(as.vector(mu)), as.integer(dk), as.double(alpha), 
#        loglik = double(1))$loglik
#   }
#   iterate_beta <- function(beta) {
#     difbeta = 1
#     maxint = 10000
#     count = 0
#     curbeta = beta
#     while ((difbeta > 1e-10) && (count < maxint)) {
#       mu = as.vector(exp(X %*% curbeta))
#       temp = transf.branch.lengths_poisson_GEE(curbeta)
#       dia = temp[[3]]
#       if (sum(which(mu == 0)) > 0) 
#         break
#       comp = three.point.compute(temp[1:2], (y - mu)/dia, 
#                                  mu * X/dia)
#       invI = solve(comp$XX)
#       newbeta = curbeta + invI %*% comp$Xy
#       count = count + 1
#       difbeta = sum(abs(newbeta - curbeta))
#       curbeta = newbeta
#     }
#     mu = as.vector(exp(X %*% curbeta))
#     r = (y - mu)/sqrt(mu)
#     phi = sum(r^2)/(n - dk)
#     covBSE = phi * invI
#     BSE = sqrt(diag(covBSE))
#     if (difbeta > 1e-10) 
#       convergeflag = 1
#     else convergeflag = 0
#     return(list(beta = as.vector(curbeta), BSE = BSE, covBSE = covBSE, 
#                 phi = phi, convergeflag = convergeflag))
#   }
#   if (is.null(start.beta)) {
#     if (method %in% c("logistic_MPLE", "logistic_IG10")) {
#       fit = glm(y ~ X - 1, family = binomial)
#       startB = fit$coefficients
#       if (any(abs(X %*% startB) >= btol)) {
#         warning("The estimated coefficients in the absence of phylogenetic signal lead\n  to some linear predictors beyond 'btol'. Increase btol?\n  Starting from beta=0 other than intercept.")
#         startB = numeric(dk)
#         iint = match("(Intercept)", colnames(X))
#         if (!is.na(iint)) 
#           startB[iint] = log(sum(y == 1)/sum(y == 0))
#         if (any(abs(X %*% startB) >= btol)) 
#           startB[iint] = 0
#       }
#     }
#     if (method == "poisson_GEE") {
#       fit = glm(y ~ X - 1, family = poisson)
#       start.beta = fit$coefficients
#     }
#   }
#   else {
#     if (length(start.beta) != dk) 
#       stop(paste("start.beta shoudl be of length", dk))
#     if (method %in% c("logistic_MPLE", "logistic_IG10")) {
#       startB = as.vector(start.beta)
#       if (any(abs(X %*% startB) >= btol)) 
#         stop("With these starting beta values, some linear predictors are beyond 'btol'.\n  Increase btol or choose new starting values for beta.")
#     }
#   }
#   if (method %in% c("logistic_MPLE", "logistic_IG10")) {
#     if (is.null(start.alpha)) 
#       startlL = log(Tmax)
#     else {
#       if (length(start.alpha) != 1) 
#         stop("start.alpha should be a single positive value")
#       if (start.alpha <= 0) 
#         stop("start.alpha should be a positive value")
#       startlL = -log(start.alpha)
#       if (abs(startlL - log(Tmax)) >= log.alpha.bound) {
#         tmp = "start.alpha is outside the bounds, which are\n  exp(+/-log.alpha.bound)/Tmax: "
#         tmp = paste(tmp, signif(exp(-log.alpha.bound)/Tmax, 
#                                 3), ",", signif(exp(log.alpha.bound)/Tmax, 
#                                                 3), " (Tmax=", Tmax, ").", "\n  Change start.alpha or increase log.alpha.bound.", 
#                     sep = "")
#         stop(tmp)
#       }
#     }
#   }
#   if (method %in% c("logistic_MPLE", "logistic_IG10")) {
#     if (method == "logistic_IG10") {
#       plogreg = plogregfunct(startB, startlL, y)
#       lL = plogreg$lL
#       B = plogreg$B
#       convergeflag = plogreg$convergeflag
#     }
#     if (method == "logistic_MPLE") {
#       opt <- optim(par = c(startB, startlL), fn = npllh, 
#                    method = "L-BFGS-B", control = list(factr = 1e+12), 
#                    y = y)
#       B = opt$par[1:dk]
#       lL = opt$par[dk + 1]
#       convergeflag = opt$convergence
#     }
#     alphaWarn = 0
#     if ((lL - log(Tmax) + 0.02) > log.alpha.bound) {
#       warn = paste("the estimate of 'alpha' (", 1/exp(lL), 
#                    ") reached the lower bound (", 1/Tmax/exp(log.alpha.bound), 
#                    ").\n This may reflect a flat likelihood at low alpha values near the lower bound,\n", 
#                    " meaning that the phylogenetic correlation is estimated to be maximal\n", 
#                    " under the model in Ives and Garland (2010).", 
#                    sep = "")
#       warning(warn)
#       alphaWarn = 1
#     }
#     if ((lL - log(Tmax) - 0.02) < -log.alpha.bound) {
#       warn = paste("the estimate of 'alpha' (", 1/exp(lL), 
#                    ") reached the upper bound (", exp(log.alpha.bound)/Tmax, 
#                    ").\n This may simply reflect a flat likelihood at large alpha values,\n", 
#                    " meaning that the phylogenetic correlation is estimated to be negligible.", 
#                    sep = "")
#       warning(warn)
#       alphaWarn = 2
#     }
#     if (btouch == 1) 
#       warning("the boundary of the linear predictor has been reached during the optimization procedure.\nYou can increase this bound by increasing 'btol'.")
#     plogregBSE = plogregBSEfunct(B, lL)
#     results <- list(coefficients = B, alpha = 1/exp(lL), 
#                     sd = plogregBSE$BSE, vcov = plogregBSE$covBSE, convergence = convergeflag)
#   }
#   if (method == "poisson_GEE") {
#     res = iterate_beta(as.vector(start.beta))
#     results <- list(coefficients = res$beta, scale = res$phi, 
#                     sd = res$BSE, vcov = res$covBSE, convergence = res$convergeflag)
#   }
#   if (results$converge) 
#     warning("phyloglm failed to converge.\n")
#   names(results$coefficients) = colnames(X)
#   colnames(results$vcov) = colnames(X)
#   rownames(results$vcov) = colnames(X)
#   results$linear.predictors = as.vector(X %*% results$coefficients)
#   names(results$linear.predictors) = names(y)
#   if (method %in% c("logistic_MPLE", "logistic_IG10")) {
#     if (max(abs(results$linear.predictors)) + 0.01 > btol) 
#       warning("the linear predictor reaches its bound for one (or more) tip.")
#     results$fitted.values = as.vector(1/(1 + exp(-results$linear.predictors)))
#     results$mean.tip.height = Tmax
#     results$logLik = llh(results$fitted.values, results$alpha)
#     results$penlogLik = results$logLik + log(det(as.matrix(plogregBSE$info)))/2
#     results$aic = -2 * results$logLik + 2 * (dk + 1)
#     results$alphaWarn = alphaWarn
#   }
#   if (method == "poisson_GEE") {
#     results$fitted.values = as.vector(exp(-results$linear.predictors))
#     results$logLik = NA
#     results$penlogLik = NA
#     results$aic = NA
#   }
#   names(results$fitted.values) = names(y)
#   results$residuals = y - results$fitted.values
#   results$y = y
#   results$n = n
#   results$d = dk
#   results$formula = formula
#   results$call = match.call()
#   results$method = method
#   results$X = X
#   results$boot = boot
#   if ((boot > 0) && (method %in% c("logistic_MPLE", "logistic_IG10"))) {
#     options(warn = -1)
#     bootobject <- rbinTrait(n = boot, phy = phy, beta = results$coefficients, 
#                             alpha = results$alpha, X = X, model = "LogReg")
#     ncoeff = length(results$coefficients)
#     bootvector <- vector(length = ncoeff + 1)
#     names(bootvector) <- c(names(results$coefficients), "alpha")
#     boot_model <- function(y) {
#       if (method == "logistic_IG10") {
#         bootfit <- try(plogregfunct(startB, startlL, 
#                                     y), silent = TRUE)
#         if (!inherits(bootfit, "try-error")) {
#           bootvector[1:ncoeff] <- bootfit$B
#           bootvector[ncoeff + 1] <- 1/exp(bootfit$lL)
#         }
#       }
#       if (method == "logistic_MPLE") {
#         bootfit <- try(optim(par = c(startB, startlL), 
#                              fn = npllh, method = "L-BFGS-B", control = list(factr = 1e+12), 
#                              y = y), silent = TRUE)
#         if (!inherits(bootfit, "try-error")) {
#           bootvector[1:ncoeff] <- bootfit$par[1:dk]
#           bootvector[ncoeff + 1] <- 1/exp(bootfit$par[dk + 
#                                                         1])
#         }
#       }
#       return(bootvector)
#     }
#     bootmatrix <- future.apply::future_lapply(as.data.frame(bootobject), 
#                                               boot_model)
#     bootmatrix <- do.call(rbind, bootmatrix)
#     ind.na <- which(is.na(bootmatrix[, 1]))
#     if (length(ind.na) > 0) {
#       bootmatrix <- bootmatrix[-ind.na, ]
#       numOnes <- range(apply(bootobject[, ind.na], 2, sum))
#     }
#     bootmean <- apply(bootmatrix, 2, mean)
#     bootmedian <- apply(bootmatrix ,2, median)
#     bootsd <- apply(bootmatrix, 2, sd)
#     bootconfint95 <- apply(bootmatrix, 2, quantile, probs = c(0.025, 
#                                                               0.975))
#     bootmeanAlog <- mean(log(bootmatrix[, ncoeff + 1]))
#     bootsdAlog <- sd(log(bootmatrix[, ncoeff + 1]))
#     bootmode <- apply(bootmatrix, 2, mode)
#     results$bootmean = bootmean
#     results$bootmedian = bootmedian
#     results$bootsd = bootsd
#     results$bootconfint95 = bootconfint95
#     results$bootmeanAlog = bootmeanAlog
#     results$bootsdAlog = bootsdAlog
#     results$bootnumFailed = length(ind.na)
#     results$bootmode = bootmode
#     if (full.matrix) 
#       results$bootstrap = bootmatrix
#     options(warn = 0)
#   }
#   class(results) = "phyloglm"
#   results
# }
# 

#read in directory of fasta files (exons)
exonload<-function(path, treeNames){
  files<-list.files(path)
  loci<-list()
  #print(files)

  for (i in 1:length(files)){
    print(paste("parsed ", i, "locus out of ", length(files), "or ", round(i/length(files)*100, digits=3), " %"))
    
    loci[[i]]<-readSet(file=paste(path,files[i], sep="/"))
  }
  
  
  for( i in 1:length(loci) ) {
    print(paste("resorted ", i, "locus out of ", length(files), "or ", round(i/length(files)*100, digits=3), " %"))
    
    loci[[i]] <- loci[[i]][treeNames[treeNames %in% names(loci[[i]])],]
    
    #loci[[i]]<-loci[[i]][treeNames,]
  }
  
  names(loci)<-unlist(strsplit((files), split="_final_align_NT.fasta"))
  return(loci)
}

Stats <- function(x){
  Mean <- mean(x, na.rm=TRUE)
  Median <- median(x, na.rm=T)
  SD <- sd(x, na.rm=TRUE)
  Min <- min(x, na.rm=TRUE)
  Max <- max(x, na.rm=TRUE)
  return(c(Mean=Mean, Median=Median, SD=SD, Min=Min, Max=Max))
}

#tip.data.molstats$exon_models
# Progress combine function
# f <- function(){
#   pb <- txtProgressBar(min=1, max=n-1,style=3)
#   count <- 0
#   function(...) {
#     count <<- count + length(list(...)) - 1
#     setTxtProgressBar(pb,count)
#     Sys.sleep(0.01)
#     flush.console()
#     c(...)
#   }
# }

multi_anov<-function(exon_stats, tree, models, sims){

  imputed<-list()
  
  for(i in 1:length(exon_stats[1,])){
    print(paste("imputing", colnames(exon_stats)[i]))
    print(paste("completed", (i/length(exon_stats[1,]))*100, "percent"))
    imputed[[i]]<- phylopars(trait_data = data.frame(species=rownames(exon_stats), scuo=exon_stats[,i]), tree = tree)$anc_recon[1:length(rownames(exon_stats))]
  }
  names(imputed)<-colnames(exon_stats)
  #return(imputed)
  
  multi_anov<-list()
  print("completeing dataset")
  
  #create the cluster
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "FORK"
  )
  
  #check cluster definition (optional)
  print(my.cluster)
  
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  
  #run the cluster
  n <- length(imputed)
  multi_anov<-list()
  multi_anov <- foreach(
    i = 1:length(imputed), 
    .combine = 'c', .inorder=F
  ) %dopar% {
    #mod <- 
    #dat <- 
    #print("test")
    phytools::phylANOVA(tree=tree, x = setNames(models, rownames(exon_stats)), y = setNames(imputed[[i]], rownames(exon_stats)), nsim=sims, posthoc = F)$Pf
    
  }
  parallel::stopCluster(cl = my.cluster)

  
  # 
  # for(i in 1:length(imputed)){
  #   print(paste("completed ANOVA", (i/length(imputed))*100, "percent"))
  #   mod <- setNames(models, rownames(exon_stats))
  #   dat <- setNames(imputed[[i]], rownames(exon_stats))
  #   print("test")
  #   multi_anov[[i]]<-phylANOVA(tree=tree, x = mod, y = dat, nsim=sims, posthoc = F)
  # }
  # 
  
  names(multi_anov)<-colnames(exon_stats)
  
  return(multi_anov)
  
}

#base comp for each base
base_comps<-function(DNAString){
  ATGC<-list()
  dnabin<- as.DNAbin(DNAString)
  print("converted to DNAbin")
  #base.freq(dnabin[i])
  
  for(i in 1:length(DNAString)){
    print(i)
    ATGC[[i]]<- base.freq(dnabin[i])
  }
  
  ATGC<- do.call(rbind, ATGC)
  rownames(ATGC)<-names(DNAString)
  return(ATGC)
}

#set up datasets for ggradar for plotting eq base freqs
baseFreq_plotSet<-function(input){
  tmp <- do.call(rbind, input$modelpar)
  tmp <- as.matrix(tmp)
  tmp <- cbind(input[,1], tmp)
  tmp<-tmp[!duplicated(tmp), ]
  colnames(tmp)<-c( "group", "A", "C", "G", "T")
  tmp<-as.data.frame(tmp)
  tmp<-tmp[ order(tmp$group),]
  tmp$group <- paste(tmp$group, "eq", sep='')
  
  return(tmp)
}

#plotting base freq shifts as bars
baseFreq_shiftplot<-function(bf_data, bf_eq, seed=1, title){
  bf.barplot<-t(as.matrix(bf_data[,c(2:5)] - bf_eq[,c(2:5)]))
  colnames(bf.barplot)<- bf_data$group
  set.seed(seed)
  barplot(las=1,
    space=c(0,1.5),
    bf.barplot,
    beside = T,
    col = c("#cacaca", "#00beff", "#ddb310", "#00b25d"),#viridis(4),
    #col = sample(rcartocolor:::carto_pal(n = 4, "Safe"), 4),
    ylim = c(-0.15, 0.1),
    legend.text = T,
    horiz = F,
    args.legend = list(
      bty = 'n',
      horiz = F,
      x = 'bottomleft',
      inset = c(0.0, 0.0)), 
    main=title) #deparse(substitute(bf_data))
  #add h=0
  #abline(h=0, lty=2)
  
}


###quadplot test
quadplot_alt <-
  function (e = NULL,
            f = NULL,
            g = NULL,
            h = NULL,
            angle = 75,
            scale.y = 0.6,
            label = 1:4,
            labelcol = rainbow(4),
            labelpch = 19,
            labelcex = 1.5,
            main = "",
            s3d.control = list(),
            simplex.control = list(),
            legend.control = list(),
            ...)
  {
    corners <- quadtrafo(diag(4) / 1)
    #corners <- corners/
    if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
      message("Package 'scatterplot3d' is required for functionality from the quadplot function")
      return(NULL)
    }
    s3d <- do.call(scatterplot3d::scatterplot3d, c(
      list(
        0.5,
        0.2886751,
        0.2041241,
        type = "n",
        xlim = range(corners[,1]),
        ylim = range(corners[,2]),
        zlim = range(corners[,3]),
        axis = FALSE,
        grid = FALSE,
        angle = angle,
        scale.y = scale.y,
        main = main
      ),
      s3d.control
    ))
    usr <- as.vector(sapply(s3d$xyz.convert(corners[-2,]), range))
    opar <- par(usr = usr, xpd = NA)
    assign("usr", usr, envir = environment(s3d$points3d))
    on.exit(par(opar))
    do.call("quadlines", c(list(e = diag(4)[c(1:4, 1, 3, 2, 4),], sp = s3d), simplex.control))
    do.call("quadpoints", c(
      list(
        e = diag(4),
        sp = s3d,
        pch = labelpch,
        col = labelcol,
        cex = labelcex
      )
    ))
    do.call("legend", c(
      list(
        usr[1],
        usr[4],
        legend = label,
        col = labelcol,
        pch = labelpch,
        cex = labelcex
      ),
      legend.control
    ))
    if (!is.null(e))
      quadpoints(e, f, g, h, sp = s3d, ...)
    return(s3d)
  }



curvemod<-function (expr, from = NULL, to = NULL, n = 101, add = FALSE, 
                    type = "l", xname = "x", xlab = xname, ylab = NULL, log = NULL, 
                    xlim = NULL, ...) {
  sexpr <- substitute(expr)
  if (is.name(sexpr)) {
    expr <- call(as.character(sexpr), as.name(xname))
  }
  else {
    if (!((is.call(sexpr) || is.expression(sexpr)) && xname %in% 
          all.vars(sexpr))) 
      stop(gettextf("'expr' must be a function, or a call or an expression containing '%s'", 
                    xname), domain = NA)
    expr <- sexpr
  }
  if (dev.cur() == 1L && !isFALSE(add)) {
    warning("'add' will be ignored as there is no existing plot")
    add <- FALSE
  }
  addF <- isFALSE(add)
  if (is.null(ylab)) 
    ylab <- deparse(expr)
  if (is.null(from) || is.null(to)) {
    xl <- if (!is.null(xlim)) 
      xlim
    else if (!addF) {
      pu <- par("usr")[1L:2L]
      if (par("xaxs") == "r") 
        pu <- extendrange(pu, f = -1/27)
      if (par("xlog")) 
        10^pu
      else pu
    }
    else c(0, 1)
    if (is.null(from)) 
      from <- xl[1L]
    if (is.null(to)) 
      to <- xl[2L]
  }
  lg <- if (length(log)) 
    log
  else if (!addF && par("xlog")) 
    "x"
  else ""
  if (length(lg) == 0) 
    lg <- ""
  if (grepl("x", lg, fixed = TRUE)) {
    if (from <= 0 || to <= 0) 
      stop("'from' and 'to' must be > 0 with log=\"x\"")
    x <- exp(seq.int(log(from), log(to), length.out = n))
  }
  else x <- seq.int(from, to, length.out = n)
  ll <- list(x = x)
  names(ll) <- xname
  y <- eval(expr, envir = ll, enclos = parent.frame())
  if (length(y) != length(x)) 
    stop("'expr' did not evaluate to an object of length 'n'")
  #if (isTRUE(add)) 
  #  lines(x = x, y = y, type = type, ...)
  #else plot(x = x, y = y, type = type, xlab = xlab, ylab = ylab, 
      #      xlim = xlim, log = lg, ...)
  return(list(x = x, y = y))
}



# #boot confidence intervals
# phylolm_bootstrap_CI <- function(input, lims, interval, breaks){
#   
#   # Create a list of curves
#   tmp <- list()
#   
#   # Loop through the number of bootstrap samples
#   for(i in 1:(input$boot)) {
#     
#     # Extract the bootstrap sample for the current iteration
#     cc <- (input$bootstrap[i,])
#     
#     # Add a transparent curve to the list
#     if(class(input)=='phyloglm'){
#     tmp[[i]] <- curvemod(plogis(cc[1]+cc[2]*x), xlim=lims, n = breaks)
#     }
#   }
#   
#   # Convert the list to a data frame
#   tmp <- as.data.frame(tmp)
#   
#   # Extract the x values and y values from the data frame
#   xvals <- tmp[, 1]
#   columns <- grep("x", names(tmp), invert = T)
#   yvals <- tmp[, columns]
#   
#   # Calculate the 2.5 and 97.5 quantiles for each row
# 
#   quants1 <- apply(yvals, MARGIN = 1, quantile, probs = c((1-interval)/2, (1+interval)/2))
#   quants2 <- apply(yvals, MARGIN = 1, HDInterval::hdi, credMass = interval)
#   
#   # Extract the lower and upper bounds from the quantiles matrix
#   lower_bounds <- quants1[1, ]
#   upper_bounds <- quants1[2, ]
#   
#   lower_dens <- quants2[1, ]
#   upper_dens <- quants2[2, ]
#   
#   # Return a data frame with the x values and lower and upper bounds
#   return(data.frame(xvals, lower_bounds, upper_bounds, lower_dens, upper_dens))
# }

#boot confidence intervals
phylolm_bootstrap_CI <- function(input, lims, interval, breaks, X2=NULL){
  
  # Create a list of curves
  tmp <- list()
  
  # Loop through the number of bootstrap samples
  for(i in 1:(input$boot)) {
    
    # Extract the bootstrap sample for the current iteration
    cc <- (input$bootstrap[i,])
    cc <- cc[-grep("alpha", names(input$bootstrap[1,]))]
    
    # Add a transparent curve to the list
    if(class(input)=='phyloglm'){
      if(length(cc)==3){
        tmp[[i]] <- curvemod(plogis(cc[1]+cc[2]*x+cc[3]*X2), xlim=lims, n = breaks)
      }
      if(length(cc)==2){
        tmp[[i]] <- curvemod(plogis(cc[1]+cc[2]*x), xlim=lims, n = breaks)
      }
    }
  }
  
  # Convert the list to a data frame
  tmp <- as.data.frame(tmp)
  
  # Extract the x values and y values from the data frame
  xvals <- tmp[, 1]
  columns <- grep("x", names(tmp), invert = T)
  yvals <- tmp[, columns]
  
  # Calculate the 2.5 and 97.5 quantiles for each row
  
  quants1 <- apply(yvals, MARGIN = 1, quantile, probs = c((1-interval)/2, (1+interval)/2))
  quants2 <- apply(yvals, MARGIN = 1, HDInterval::hdi, credMass = interval)
  
  # Extract the lower and upper bounds from the quantiles matrix
  lower_bounds <- quants1[1, ]
  upper_bounds <- quants1[2, ]
  
  lower_dens <- quants2[1, ]
  upper_dens <- quants2[2, ]
  
  # Return a data frame with the x values and lower and upper bounds
  return(data.frame(xvals, lower_bounds, upper_bounds, lower_dens, upper_dens))
}

#plot polygon
phylolm_plot_CIs<-function(input, type='quant', n=.nknots.smspl, smooth, color, alpha, outline=T){
  
  # Extract the X and Y coordinates from the data frame
  if(type=='quant'){
    x <- input$xvals
    y1 <- input$lower_bounds
    y2 <- input$upper_bounds
  }
  if(type=='hdi'){
    x <- input$xvals
    y1 <- input$lower_dens
    y2 <- input$upper_dens
  }
  
  # Fit a smooth curve to the X and Y coordinates
  spline1 <- smooth.spline(x, y1, nknots = n, spar=smooth)#, breaks)
  spline2 <- smooth.spline(x, y2, nknots = n, spar=smooth)#, breaks)
  
  # Plot the polygon representing the confidence interval
  #plot(x, y1, type = "n", xlim = range(x), ylim = range(c(y1, y2)), xlab = "X", ylab = "Y")
  #polygon(c(x, rev(x)), c(y1, rev(y2)), col = scales::alpha("gray", 0.25), border = NA)
  #clipping y just a bit
  clip(x1=par('usr')[1], x2=10, y1=0,y2=1)
  polygon(c(spline1$x, rev(spline1$x)), c(spline1$y, rev(spline2$y)), col = scales::alpha(color, alpha), border = NA)
  
  # Add dashed lines to the top and bottom of the polygon
  if(outline==T){
    lines(spline1$x, spline1$y, col = color, lty = "dashed")
    lines(spline2$x, spline2$y, col = color, lty = "dashed")
  }

  
}


#3D surface plotting (Testing)
#https://stackoverflow.com/questions/3979240/r-plotting-a-3d-surface-from-x-y-z
plot_rgl_model_a <- function(fdata, plot_contour = T, plot_points = T, 
                             verbose = F, colour = "rainbow", smoother = F){
  ## takes a model in long form, in the format
  ## 1st column x
  ## 2nd is y,
  ## 3rd is z (height)
  ## and draws an rgl model
  
  ## includes a contour plot below and plots the points in blue
  ## if these are set to TRUE
  
  # note that x has to be ascending, followed by y
  if (verbose) print(head(fdata))
  
  fdata <- fdata[order(fdata[, 1], fdata[, 2]), ]
  if (verbose) print(head(fdata))
  ##
  require(reshape2)
  require(rgl)
  orig_names <- colnames(fdata)
  colnames(fdata) <- c("x", "y", "z")
  fdata <- as.data.frame(fdata)
  
  ## work out the min and max of x,y,z
  xlimits <- c(min(fdata$x, na.rm = T), max(fdata$x, na.rm = T))
  ylimits <- c(min(fdata$y, na.rm = T), max(fdata$y, na.rm = T))
  zlimits <- c(min(fdata$z, na.rm = T), max(fdata$z, na.rm = T))
  l <- list (x = xlimits, y = ylimits, z = zlimits)
  xyz <- do.call(expand.grid, l)
  if (verbose) print(xyz)
  x_boundaries <- xyz$x
  if (verbose) print(class(xyz$x))
  y_boundaries <- xyz$y
  if (verbose) print(class(xyz$y))
  z_boundaries <- xyz$z
  if (verbose) print(class(xyz$z))
  if (verbose) print(paste(x_boundaries, y_boundaries, z_boundaries, sep = ";"))
  
  # now turn fdata into a wide format for use with the rgl.surface
  fdata[, 2] <- as.character(fdata[, 2])
  fdata[, 3] <- as.character(fdata[, 3])
  #if (verbose) print(class(fdata[, 2]))
  wide_form <- dcast(fdata, y ~ x, value.var = "z")
  if (verbose) print(head(wide_form))
  wide_form_values <- as.matrix(wide_form[, 2:ncol(wide_form)])  
  if (verbose) print(wide_form_values)
  x_values <- as.numeric(colnames(wide_form[2:ncol(wide_form)]))
  y_values <- as.numeric(wide_form[, 1])
  if (verbose) print(x_values)
  if (verbose) print(y_values)
  wide_form_values <- wide_form_values[order(y_values), order(x_values)]
  wide_form_values <- as.numeric(wide_form_values)
  x_values <- x_values[order(x_values)]
  y_values <- y_values[order(y_values)]
  if (verbose) print(x_values)
  if (verbose) print(y_values)
  
  if (verbose) print(dim(wide_form_values))
  if (verbose) print(length(x_values))
  if (verbose) print(length(y_values))
  
  zlim <- range(wide_form_values)
  if (verbose) print(zlim)
  zlen <- zlim[2] - zlim[1] + 1
  if (verbose) print(zlen)
  
  if (colour == "rainbow"){
    colourut <- rainbow(zlen, alpha = 0)
    if (verbose) print(colourut)
    col <- colourut[ wide_form_values - zlim[1] + 1]
    # if (verbose) print(col)
  } else {
    col <- "grey"
      if (verbose) print(table(col2))
  }
  
  
  open3d()
  plot3d(x_boundaries, y_boundaries, z_boundaries, 
         box = T, col = "black",  xlab = orig_names[1], 
         ylab = orig_names[2], zlab = orig_names[3])
  
  rgl.surface(z = x_values,  ## these are all different because
              x = y_values,  ## of the confusing way that 
              y = wide_form_values,  ## rgl.surface works! - y is the height!
              coords = c(2,3,1),
              color = col,
              alpha = 1.0,
              lit = F,
              smooth = smoother)
  
  if (plot_points){
    # plot points in red just to be on the safe side!
    points3d(fdata, col = "blue")
  }
  
  if (plot_contour){
    # plot the plane underneath
    flat_matrix <- wide_form_values
    if (verbose) print(flat_matrix)
    y_intercept <- (zlim[2] - zlim[1]) * (-2/3) # put the flat matrix 1/2 the distance below the lower height 
    flat_matrix[which(flat_matrix != y_intercept)] <- y_intercept
    if (verbose) print(flat_matrix)
    
    rgl.surface(z = x_values,  ## these are all different because
                x = y_values,  ## of the confusing way that 
                y = flat_matrix,  ## rgl.surface works! - y is the height!
                coords = c(2,3,1),
                color = col,
                alpha = 1.0,
                smooth = smoother)
  }
}

add_rgl_model <- function(fdata){
  
  ## takes a model in long form, in the format
  ## 1st column x
  ## 2nd is y,
  ## 3rd is z (height)
  ## and draws an rgl model
  
  ##
  # note that x has to be ascending, followed by y
  print(head(fdata))
  
  fdata <- fdata[order(fdata[, 1], fdata[, 2]), ]
  
  print(head(fdata))
  ##
  require(reshape2)
  require(rgl)
  orig_names <- colnames(fdata)
  
  #print(head(fdata))
  colnames(fdata) <- c("x", "y", "z")
  fdata <- as.data.frame(fdata)
  
  ## work out the min and max of x,y,z
  xlimits <- c(min(fdata$x, na.rm = T), max(fdata$x, na.rm = T))
  ylimits <- c(min(fdata$y, na.rm = T), max(fdata$y, na.rm = T))
  zlimits <- c(min(fdata$z, na.rm = T), max(fdata$z, na.rm = T))
  l <- list (x = xlimits, y = ylimits, z = zlimits)
  xyz <- do.call(expand.grid, l)
  #print(xyz)
  x_boundaries <- xyz$x
  #print(class(xyz$x))
  y_boundaries <- xyz$y
  #print(class(xyz$y))
  z_boundaries <- xyz$z
  #print(class(xyz$z))
  
  # now turn fdata into a wide format for use with the rgl.surface
  fdata[, 2] <- as.character(fdata[, 2])
  fdata[, 3] <- as.character(fdata[, 3])
  #print(class(fdata[, 2]))
  wide_form <- dcast(fdata, y ~ x, value.var = "z")
  print(head(wide_form))
  wide_form_values <- as.matrix(wide_form[, 2:ncol(wide_form)])  
  x_values <- as.numeric(colnames(wide_form[2:ncol(wide_form)]))
  y_values <- as.numeric(wide_form[, 1])
  print(x_values)
  print(y_values)
  wide_form_values <- wide_form_values[order(y_values), order(x_values)]
  x_values <- x_values[order(x_values)]
  y_values <- y_values[order(y_values)]
  print(x_values)
  print(y_values)
  
  print(dim(wide_form_values))
  print(length(x_values))
  print(length(y_values))
  
  rgl.surface(z = x_values,  ## these are all different because
              x = y_values,  ## of the confusing way that 
              y = wide_form_values,  ## rgl.surface works!
              coords = c(2,3,1),
              alpha = .8)
  # plot points in red just to be on the safe side!
  points3d(fdata, col = "red")
}

### End function definitions ###
