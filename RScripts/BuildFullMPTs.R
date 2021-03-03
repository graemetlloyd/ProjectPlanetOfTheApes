# SCRIPT TO BUILD METATREE FILES
rm(list=ls())
# Load metatree library:
library(Claddis)
library(metatree)

setwd("ProjectPlanetOfTheApes/Metatree Data/TreeFiles/")

# Safely reinsert exclude taxa and write out to file:
Claddis::safe_taxonomic_reinsertion(input_filename = "STR_MPTs.tre", output_filename  = "MPTs.tre", str_taxa = read.table("PrimateSTRforEXC.txt", header = TRUE, stringsAsFactors = FALSE), multiple_placement_option = "random")

# Get exclude strict consensus and write to file:
ExcludeMPTs <- ape::read.tree("MPTs.tre")
ExcludeSCC <- ape::consensus(ExcludeMPTs)
ExcludeMR50<- ape::consensus(ExcludeMPTs, p=0.5)
ExcludeAdams <- ape:::read.nexus("../ConsensusTrees/Adams_noChunks.nex")

setwd("../Chunks/")
get.tnt.names <- function(file) {
  x <- readLines(file)
  begin <- grep("&", x)+1
  end <- begin+as.numeric(strsplit(x[begin -2], split = " ")[[1]][2])-1
  unlist(lapply(strsplit(x[begin:end], " "), function(x)x[1]))
}


treefiles <- list.files()[grep("strict.nex", list.files())]
tnt.files <- list.files()[grep("tnt", list.files())]


for(i in 1:length(treefiles)) {
  Trees <- read.tree(text = readLines(treefiles[i], warn = FALSE)[grep(",", readLines(treefiles[i], warn = FALSE))])
  tree <- Trees[[length(Trees)]]
  xxx<-get.tnt.names(tnt.files[i])
  yyy <-setdiff(tree$tip.label, xxx)
  if(length(yyy)>0) {
    for(j in yyy) tree$tip.label[match(j, tree$tip.label)] <- xxx[agrep(j, xxx)]
  }
  tree<-drop.tip(tree, "allzero")
  clade<-stringr:::str_to_upper(gsub(pattern = "_mpt_strict.nex",replacement = "",x = treefiles[i]))
  repl <- match(clade, ExcludeMR50$tip.label)
  ExcludeMR50 <- bind.tree(x = ExcludeMR50, y= tree, where=repl, position=1)
  ExcludeMR50<-drop.tip(ExcludeMR50, clade)

  repl <- match(clade, ExcludeSCC$tip.label)
  ExcludeSCC <- bind.tree(x = ExcludeSCC, y= tree, where=repl, position=1)
  ExcludeSCC<-drop.tip(ExcludeSCC, clade)

  repl <- match(clade, ExcludeAdams$tip.label)
  ExcludeAdams <- bind.tree(x = ExcludeAdams, y= tree, where=repl, position=1)
  ExcludeAdams<-drop.tip(ExcludeAdams, clade)

    rm(list=c('Trees', 'tree'))
}



write.tree(ExcludeMR50, "../ConsensusTrees/MRC.tre")
write.tree(ExcludeSCC, "../ConsensusTrees/SCC.tre")
write.tree(ExcludeAdams, "../ConsensusTrees/Adams.tre")

