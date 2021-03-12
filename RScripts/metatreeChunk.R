library(metatree)
PrimatesExclude <- Metatree(MRPDirectory = "ProjectPlanetOfTheApes/InputData/MRP", 
                            XMLDirectory = "ProjectPlanetOfTheApes/InputData/XML", 
                            TargetClade = "Primates", InclusiveDataList = c(), ExclusiveDataList = c("Andrews_1988a", "Beard_et_MacPhee_1994a", "Begun_et_Kordos_1997a", "Burger_2010aa", "Gebo_etal_2001a", "Kimbel_etal_2004a", "Pattinson_etal_2015a", "Maiolino_etal_2012a", "Seiffert_etal_2015a", "Rose_1997a", "Marivaux_etal_2001a", "Boyer_etal_2017ab", "Boyer_etal_2017ac"), 
                            MissingSpecies = "exclude", VeilLine = TRUE, 
                            IncludeSpecimenLevelOTUs = TRUE, RelativeWeights = c(0, 100, 10, 1), 
                            WeightCombination = "product", ReportContradictionsToScreen = FALSE, 
                            BackboneConstraint = "Springer_etal_2012a")


PrimatesExclude$MonophyleticTaxa


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

#load("ProjectPlanetOfTheApes/analysis_files/PrimatesExclude.RData")

# first get the most inclusive clades 
tipslist<-FindInclusiveMonophyly(taxonomyTree =PrimatesExclude$TaxonomyTree, MonophyleticClades = PrimatesExclude$MonophyleticTaxa)


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
  Claddis::WriteMorphNexus(tmp, paste("~/Dropbox/workingScripts", names(short.tiplist)[i],".nex", sep=""))
  Claddis::WriteMorphTNT(tmp, paste("~/Dropbox/workingScripts", names(short.tiplist)[i],".tnt", sep=""))
  
}

## now run metatree again with higherTaxaToCollapse set to clades tat have been chunked


highertaxaCollapse<- c("Apidium","Atelidae","Callithrix","Carpolestinae","Cebus","Cercocebus","Cheirogaleus","Chiropotes",
                        "Colobus","Crouzeliini","Ekgmowechashalidae","Eulemur","Hylobatidae","Lagothrix","Lepilemur_et_Lepilemuridae",
                        "Microsyopinae","Microsyops","Nomascus","Notharctus","Omomys","Paromomyidae","Phaner","Piliocolobus",
                        "Plecturocebus","Plesiadapidae","Presbytis","Propithecus","Purgatorius","Rhinopithecus","Saguinus",
                        "Sapajus","Soriacebinae","Tetonius","Theropithecus","Trogolemurini","Vectipithex")              

# Build primates metatree excluding missing species:
PrimatesCollapse <- Metatree(MRPDirectory = "ProjectPlanetOfTheApes/InputData/MRP", 
                            XMLDirectory = "ProjectPlanetOfTheApes/InputData/XML", 
                            TargetClade = "Primates", HigherTaxaToCollapse = highertaxaCollapse, InclusiveDataList = c(), ExclusiveDataList = c("Andrews_1988a", "Beard_et_MacPhee_1994a", "Begun_et_Kordos_1997a", "Burger_2010aa", "Gebo_etal_2001a", "Kimbel_etal_2004a", "Pattinson_etal_2015a", "Maiolino_etal_2012a", "Seiffert_etal_2015a", "Rose_1997a", "Marivaux_etal_2001a", "Boyer_etal_2017ab", "Boyer_etal_2017ac"), 
                            MissingSpecies = "exclude", VeilLine = TRUE, 
                            IncludeSpecimenLevelOTUs = TRUE, RelativeWeights = c(0, 100, 10, 1), 
                            WeightCombination = "product", ReportContradictionsToScreen = FALSE, 
                            BackboneConstraint = "Springer_etal_2012a")

# Write out exclude data:
Claddis::WriteMorphNexus(PrimatesExclude$FullMRPMatrix, "ProjectPlanetOfTheApes/MetatreeData/PrimateEXC.nex")
Claddis::WriteMorphTNT(PrimatesExclude$FullMRPMatrix, "ProjectPlanetOfTheApes/MetatreeData/PrimateEXC.tnt")
Claddis::WriteMorphNexus(PrimatesCollapse$STRMRPMatrix, "ProjectPlanetOfTheApes/MetatreeData/PrimateEXCSTR.nex")
Claddis::WriteMorphTNT(PrimatesCollapse$STRMRPMatrix, "ProjectPlanetOfTheApes/MetatreeData/PrimateEXCSTR.tnt")
ape::write.tree(PrimatesCollapse$TaxonomyTree, "ProjectPlanetOfTheApes/MetatreeData/PrimatesEXCTaxonomy.tre")
write.table(PrimatesCollapse$SafelyRemovedTaxa, "ProjectPlanetOfTheApes/MetatreeData/PrimateSTRforEXC.txt", row.names = FALSE)
write.table(PrimatesCollapse$MonophyleticTaxa, "ProjectPlanetOfTheApes/MetatreeData/PrimateSTRforEXC_Monophyly.txt", row.names = FALSE)
write.table(PrimatesCollapse$RemovedSourceData, "ProjectPlanetOfTheApes/MetatreeData/PrimateSTRforEXC_excludedData.txt", row.names = FALSE)
