PrepareBeastMetatree<- function(constraintTree, StratRanges, alignment, myfileName, makeStartTree = TRUE, start.tree.method="mbl", vartime=0.1) {
  library(ape); library(geiger); library(paleotree);

  ## remove extant taxa ##
  StratRanges <- remove.extant(StratRanges)
  
  write.sampling.dates(StratRanges, myfileName)
  
  #td<-treedata(constraintTree, StratRanges)$data
  
  ## get starting values for ages ##
  starting.ages<-apply(StratRanges[,1:2], 1, mean)
  
  write.table(starting.ages, paste(myfileName, "mean.ages.txt", sep="_"),quote = F, sep = "\t")
  
  ## drop tips and make constraint ##
  ## need to remove sequences not  ##
  ##  in tree at this point too    ##
  
  extra_extant<- setdiff(names(alignment), constraintTree$tip.label)
  if(length(extra_extant>0)) {
    print(paste("dropping the following taxa from the alignment as they are not in the tree\n"))
    print(extra_extant)
    alignment <- alignment[-match(setdiff(names(alignment), constraintTree$tip.label), names(alignment))]
  }
 
  make.monophyly.tree(constraintTree, myfileName)

  
  new.alignment <- list()
  nchar <- length(alignment[[1]])
  for(i in 1:length(constraintTree$tip.label)) {
    aa<- match(constraintTree$tip.label[i], names(alignment))
    if(is.na(aa)) {
      new.alignment[[i]] <- rep("?", nchar)
    } else{
      new.alignment[[i]] <- alignment[[aa]]
    }
  }
  names(new.alignment) <- constraintTree$tip.label
  write.nexus.data(new.alignment, paste(myfileName, ".nex", sep=""))
  
  if(makeStartTree==TRUE) {
    extant <- setdiff(constraintTree$tip.label,names(starting.ages))
    timedat <- cbind(c(setNames(rep(0, length(extant)), extant), starting.ages),c(setNames(rep(0, length(extant)), extant), starting.ages))
    
    if(start.tree.method=="cal3") {
      starting.tree <- cal3TimePaleoPhy(constraintTree, timedat, brRate = 0.1, extRate = 0.1,sampRate = 0.1, ntrees=1, anc.wt = 0,root.max=10)
    }
    if(start.tree.method=="mbl")  
      {starting.tree <-timePaleoPhy(constraintTree, timedat, type = "mbl", vartime = vartime, ntrees = 1,
                 timeres = TRUE, add.term = FALSE,
                 inc.term.adj = FALSE, dateTreatment = "firstLast", node.mins = NULL,
                 noisyDrop = TRUE, plot = FALSE)
    }
     write.starting.tree(starting.tree, myfileName)
     write.tree(starting.tree, "startingtree.tre")
  }
  print("done!")
  }


#############################################
remove.extant <- function(d) {
  
  return(d[-which(d[,2]==0), ])
}

#############################################

make.monophyly.tree <- function(tree, outputname="file") {
  
  output1 <- file(paste(outputname, "constraints.txt", sep="_"), "w")  
  output2 <- file(paste(outputname,"logs.txt", sep="_"), "w")  
  ntax <- length(tree$tip.label)
  nodes <- seq(ntax+1, ntax+tree$Nnode)
  all.species <- sapply(nodes, tips, phy=tree)
  names(all.species) <- paste("node", nodes, sep="_")
  
  
  
  for(i in 1:length(all.species)) {
    species <- all.species[[i]]
    
    #id[is.na(match(species, all.species))] <- "idref"
    cat(paste("\t","\t","\t", "<distribution id=\"", names(all.species)[i], ".prior\" spec=\"beast.math.distributions.MRCAPrior\" monophyletic=\"true\" tree=\"x\">", sep=""), file = output1, sep = "\t")
    cat("\n", file = output1)
    cat(paste("\t","\t","\t","\t", "<taxonset id=\"",  names(all.species)[i], "\" spec=\"TaxonSet\">", sep=""), file = output1, sep = "\t")
    cat("\n", file = output1)
    
    if(i==1) {
      for(k in 1:length(species)) {
        cat(paste("\t","\t","\t","\t","\t", "<taxon id =\"", species[k], "\" spec=\"Taxon\"/>", sep=""), file = output1, sep = "\t")
        cat("\n", file = output1)
        
        #all.species <- all.species[-match(species, all.species)]
        
      }
      
      
      
    } else {
      
      for(k in 1:length(species)) {
        cat(paste("\t","\t","\t","\t","\t", "<taxon idref =\"", species[k], "\"/>", sep=""), file = output1, sep = "\t")
        cat("\n", file = output1)
        
        #all.species <- all.species[-match(species, all.species)]
        
      }
    }
      
      cat(paste("\t","\t","\t","\t", "</taxonset>", sep=""), file = output1, sep = "\t")
      cat("\n", file = output1)
      cat(paste("\t","\t","\t", "</distribution>", sep=""), file = output1, sep = "\t")
      cat("\n", "\n", file = output1)
      cat("\n", file = output1)
    
    
    cat(paste("<log idref=\"", names(all.species)[i], ".prior\"/>", sep=""), file=output2, sep="\t")
    cat("\n", file = output2)
    
    
    
  }
  close(output1)
  close(output2)
  
}


#######################################

write.sampling.dates <- function(d, myfileName) {
  
  sd <- file(paste(myfileName, "sample.dates.txt", sep="_"), open="w")
  writeLines("<operator spec=\"SampledNodeDateRandomWalker\" windowSize=\"8\" tree=\"x\" weight=\"25\">", con=sd, sep = "\n")
    
  for(i in 1:nrow(d)) {
    
    writeLines(paste("\t", "\t", "<samplingDates id=\"samplingDate",i,"\" spec=\"beast.evolution.tree.SamplingDate\" taxon=\"", rownames(d)[i] ,"\" upper=\"", d[i,1], "\" lower=\"", d[i,2], "\"/>", sep=""), con=sd, sep = "\n")
    
    
  }
  writeLines("</operator>", con=sd, sep = "\n")
    
  close(sd)
  
}


#######################################

write.starting.tree <- function(tree, myfileName) {
  
  st <- file(paste(myfileName, "starting.tree.txt", sep="_"), open="w")
  writeLines("<init id=\"StartingTree.t:x\" initial=\"@Tree.t:x\" spec=\"beast.util.TreeParser\" IsLabelledNewick=\"true\"> <input name=\"newick\"> ", con=st)
  writeLines(gsub(";", ":0.0 </input> </init> ",write.tree(tree)), con=st, sep = "\n")
  close(st)
  
}
