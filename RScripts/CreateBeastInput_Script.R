rm(list=ls())
library(ape);
library(geiger)
library(paleotree)
#<<<<<<< HEAD
source("ProjectPlanetOfTheApes/BuildBeastXML.R")


### read in constraint tree (MRP consensus, MRL tree etc )
### then read in stratigraphic data in csv format
### strat data should be csv file with name, fad, lad
### finally, read in molecular alignment for extant taxa

constraintTree <- read.tree("MetatreeData/ConsensusTrees/MRC.tre")
StratRanges <- read.csv("final_dates_primates.csv", row.names = 1, stringsAsFactors = F)
alignment <- read.nexus.data("CYTB_add/combined_CYTB.nexus")
setdiff(names(alignment),constraintTree$tip.label)

setdiff(rownames(StratRanges)[which(StratRanges[,2]==0)], names(alignment))
setdiff(names(alignment),rownames(StratRanges)[which(StratRanges[,2]==0)])

# if taxon names in strat range contain spaces, uncomment and run line below
#rownames(StratRanges)<-gsub(" ", "_", rownames(StratRanges))

#if age ranges larger than some amount are not desirable run two lines below
max.age.range <- 10
StratRanges<-StratRanges[-which(StratRanges[,1]-StratRanges[,2]>max.age.range),] ## remove taxa with large strat ranges


name.check(constraintTree, StratRanges)
td <- treedata(constraintTree, StratRanges)



constraintTree <- td$phy
StratRanges <- td$data

fileStem = "Primates_BEAST"

PrepareBeastMetatree(constraintTree = constraintTree, 
                     StratRanges = StratRanges, 
                     alignment = alignment, 
                     myfileName = fileStem,
                     makeStartTree = TRUE, 
                     start.tree.method="mbl", 
                     vartime=0.5)

## output will be a series of files that can be used to create a BEAST xml through beauti
## i.e. nexus file and mean taxon ages
## as well as files that can be pasted into the xml (monophyly constraints, starting tree, age ranges).

mean.ages <-read.table("Primates_BEAST_mean.ages.txt", row.names = 1, stringsAsFactors = F)
846-nrow(mean.ages)

#401 extant taxa
 401/518 

