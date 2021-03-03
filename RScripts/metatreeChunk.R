library(metatree)
PrimatesExclude <- Metatree(MRPDirectory = "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/InputData/MRP", 
                            XMLDirectory = "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/InputData/XML", 
                            TargetClade = "Primates", InclusiveDataList = c(), ExclusiveDataList = c("Andrews_1988a", "Beard_et_MacPhee_1994a", "Begun_et_Kordos_1997a", "Burger_2010aa", "Gebo_etal_2001a", "Kimbel_etal_2004a", "Pattinson_etal_2015a", "Maiolino_etal_2012a", "Seiffert_etal_2015a", "Rose_1997a", "Marivaux_etal_2001a", "Boyer_etal_2017ab", "Boyer_etal_2017ac"), 
                            MissingSpecies = "exclude", VeilLine = TRUE, 
                            IncludeSpecimenLevelOTUs = TRUE, RelativeWeights = c(0, 100, 10, 1), 
                            WeightCombination = "product", ReportContradictionsToScreen = FALSE, 
                            BackboneConstraint = "Springer_etal_2012a")


PrimatesExclude$MonophyleticTaxa


Claddis::WriteMorphNexus(PrimatesExclude$STRMRPMatrix, "~/Dropbox/workingScripts/finalFiles/MRL/PrimateEXCSTR.nex")
Claddis::WriteMorphTNT(PrimatesExclude$STRMRPMatrix, "~/Dropbox/workingScripts/finalFiles/MRL/PrimateEXCSTR.tnt")

write.table(PrimatesExclude$SafelyRemovedTaxa, "~/Dropbox/workingScripts/finalFiles/MRL/PrimateSTRforEXC.txt", row.names = FALSE)


d<-ReadMorphNexus("~/Dropbox/workingScripts/finalFiles/MRL/PrimateEXCSTR.nex")

weights<-(d$Matrix_1$Weights)
round.weights <- round(weights, 0)
write.table(t(round.weights), "~/Dropbox/workingScripts/finalFiles/MRL/PrimateEXCSTR_weights.txt",sep="\n", row.names = F, col.names = F)
write.table((round.weights), "~/Dropbox/workingScripts/finalFiles/MRL/PrimateEXCSTR_weights.csv", row.names = F,  sep="\t", quote=F)


mrp <-d$Matrix_1$Matrix
str(mrp)
mrp[which((is.na(mrp)))] <- "?"
mrp[2,1:10]
d <- list()
for(i in 1:nrow(mrp)) d[[i]] <- mrp[i,]
names(d) <- rownames(mrp)
switch.code <- function(x){
  if(x== "?") return("?")
  if(x== 1) return(0)
  if(x==0) return(1)
}

nchar<-length(d[[1]])
ntax<-length(d)
yy<-sort(sample(x=seq(1, nchar),size = (nchar/2), replace = F))

for(i in 1:ntax) {
  for(j in 1:length(yy)) {
    d[[i]][yy[j]] <- switch.code(d[[i]][yy[j]] )
  }
}
mrp[1,1:10]
unlist(d[[1]][1:10])
unique(unlist(lapply(d, length)))
names(d)

write.nexus.data(d, file="~/Dropbox/workingScripts/finalFiles/MRL/PrimateEXCSTR_mrl.nex", interleaved = F)

tmp <- read.tree("~/Dropbox/workingScripts/finalFiles/MRL/Exclude.tre")
tmp$tip.label[409] <- "Proteropithecia_neuquenensis"
write.tree(multi2di(drop.tip(tmp, setdiff(tmp$tip.label, names(d)))), "~/Dropbox/workingScripts/finalFiles/MRL/StartTree.tre")
## Functions to use with chunking ##
library(geiger)
FindInclusiveMonophyly <- function(taxonomyTree, MonophyleticClades) {
  ## takes the taxonomy tree and vector of monophyletic clades - spits out most inclusive clades (i.e. removes any nested clades)
  mono <- MonophyleticClades
  taxontre<-taxonomyTree

  nodesForTips<-match(mono, taxontre$node.label)

  tipslist <- list()

  for(i in 1:length(mono)) {
    tipslist[[i]] <- tips(phy = taxontre, node = (length(taxontre$tip.label))+ nodesForTips[i])
  }

  names(tipslist) <- mono
  y <- matrix(data=NA, nrow=length(mono), ncol=length(mono))
  rownames(y) <- colnames(y) <- mono
  for(i in 1:length(tipslist)) {
    for(j in 1:length(tipslist)) {
      y[i,j]<- (all(tipslist[[i]]%in%tipslist[[j]]))
    }
  }

  diag(y) <- NA
 # y[lower.tri(y)] <- NA
  z <-which(y==TRUE, arr.ind = TRUE)
  unique.mono <-mono[-match(rownames(z), mono)]
  tipslist<-tipslist[-match(rownames(z),names(tipslist))]
  
  tipslist
  
}


extract_chunk <- function(FullMRPMatrix, otus) {
  # takes full mrp matrix and vector of otus - produces matrices
  FullMRPMatrix$Matrix_1$Matrix <- FullMRPMatrix$Matrix_1$Matrix[c(1,match(otus, rownames(FullMRPMatrix$Matrix_1$Matrix))), ]
  FullMRPMatrix
}


#####################
# load metatree output in RData format

#load("~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/analysis_files/PrimatesExclude.RData")

# first get the most inclusive clades 
tipslist<-FindInclusiveMonophyly(taxonomyTree =PrimatesExclude$TaxonomyTree, MonophyleticClades = PrimatesExclude$MonophyleticTaxa)
## some how is misses nomascus so


# find clade sizes - we will only chunk clades with 4 or more taxa
cladesizes<-unlist(lapply(tipslist, length))
sort(cladesizes)
short.tiplist<-tipslist[-which(cladesizes<4)]
sum(unlist(lapply(short.tiplist, length)))
highertaxaCollapse<-names(short.tiplist)
# now pull out matrices and save 

FullMRPMatrix <- PrimatesExclude$FullMRPMatrix

# note Graham created his own folder here for putting these in called working scripts, so change this if you want to run it
for(i in 1:length(short.tiplist)) {
  otus <- short.tiplist[[i]]
  tmp <- extract_chunk(FullMRPMatrix = FullMRPMatrix, otus =otus)
  Claddis::WriteMorphNexus(tmp, paste("~/Dropbox/workingScripts/sum", names(short.tiplist)[i],".nex", sep=""))
  Claddis::WriteMorphTNT(tmp, paste("~/Dropbox/workingScripts/sum", names(short.tiplist)[i],".tnt", sep=""))
  
}

## now run metatree again with higherTaxaToCollapse set to clades tat have been chunked


# highertaxaCollapse<- c("Apidium","Atelidae","Callithrix","Carpolestinae","Cebus","Cercocebus","Cheirogaleus","Chiropotes",
#                        "Colobus","Crouzeliini","Ekgmowechashalidae","Eulemur","Hylobatidae","Lagothrix","Lepilemur_et_Lepilemuridae",
#                        "Microsyopinae","Microsyops","Nomascus","Notharctus","Omomys","Paromomyidae","Phaner","Piliocolobus",
#                        "Plecturocebus","Plesiadapidae","Presbytis","Propithecus","Purgatorius","Rhinopithecus","Saguinus",
#                        "Sapajus","Soriacebinae","Tetonius","Theropithecus","Trogolemurini","Vectipithex")              

# Build primates metatree excluding missing species:
PrimatesCollapse <- Metatree(MRPDirectory = "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/InputData/MRP", 
                            XMLDirectory = "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/InputData/XML", 
                            TargetClade = "Primates", HigherTaxaToCollapse = highertaxaCollapse, InclusiveDataList = c(), ExclusiveDataList = c("Andrews_1988a", "Beard_et_MacPhee_1994a", "Begun_et_Kordos_1997a", "Burger_2010aa", "Gebo_etal_2001a", "Kimbel_etal_2004a", "Pattinson_etal_2015a", "Maiolino_etal_2012a", "Seiffert_etal_2015a", "Rose_1997a", "Marivaux_etal_2001a", "Boyer_etal_2017ab", "Boyer_etal_2017ac"), 
                            MissingSpecies = "exclude", VeilLine = TRUE, 
                            IncludeSpecimenLevelOTUs = TRUE, RelativeWeights = c(0, 100, 10, 1), 
                            WeightCombination = "product", ReportContradictionsToScreen = FALSE, 
                            BackboneConstraint = "Springer_etal_2012a")

# Write out exclude data:
#Claddis::WriteMorphNexus(PrimatesExclude$FullMRPMatrix, "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/analysis_files/PrimateEXC.nex")
#Claddis::WriteMorphTNT(PrimatesExclude$FullMRPMatrix, "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/analysis_files/PrimateEXC.tnt")
Claddis::WriteMorphNexus(PrimatesCollapse$STRMRPMatrix, "~/Dropbox/workingScripts/product/PrimateEXCSTR.nex")
Claddis::WriteMorphTNT(PrimatesCollapse$STRMRPMatrix, "~/Dropbox/workingScripts/product/PrimateEXCSTR.tnt")
ape::write.tree(PrimatesCollapse$TaxonomyTree, "~/Dropbox/workingScripts/product/PrimatesEXCTaxonomy.tre")
write.table(PrimatesCollapse$SafelyRemovedTaxa, "~/Dropbox/workingScripts/product/PrimateSTRforEXC.txt", row.names = FALSE)
write.table(PrimatesCollapse$MonophyleticTaxa, "~/Dropbox/workingScripts/product/PrimateSTRforEXC_Monophyly.txt", row.names = FALSE)
write.table(PrimatesCollapse$RemovedSourceData, "~/Dropbox/workingScripts/product/PrimateSTRforEXC_excludedData.txt", row.names = FALSE)
