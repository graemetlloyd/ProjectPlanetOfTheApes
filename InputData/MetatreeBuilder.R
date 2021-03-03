# Load metatree library:
library(metatree)

# Build primates metatree excluding missing species:
PrimatesExclude <- Metatree(MRPDirectory = "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/InputData/MRP", 
                            XMLDirectory = "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/InputData/XML", 
                            TargetClade = "Primates", InclusiveDataList = c(), ExclusiveDataList = c("Andrews_1988a", "Beard_et_MacPhee_1994a", "Begun_et_Kordos_1997a", "Burger_2010aa", "Gebo_etal_2001a", "Kimbel_etal_2004a", "Pattinson_etal_2015a", "Maiolino_etal_2012a", "Seiffert_etal_2015a", "Rose_1997a", "Marivaux_etal_2001a", "Boyer_etal_2017ab", "Boyer_etal_2017ac"), 
                            MissingSpecies = "exclude", VeilLine = TRUE, 
                            IncludeSpecimenLevelOTUs = TRUE, RelativeWeights = c(0, 100, 10, 1), 
                            WeightCombination = "sum", ReportContradictionsToScreen = FALSE, 
                            BackboneConstraint = "Springer_etal_2012a")

# Write out exclude data:
Claddis::WriteMorphNexus(PrimatesExclude$FullMRPMatrix, "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/analysis_files/PrimateEXC.nex")
Claddis::WriteMorphTNT(PrimatesExclude$FullMRPMatrix, "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/analysis_files/PrimateEXC.tnt")
Claddis::WriteMorphNexus(PrimatesExclude$STRMRPMatrix, "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/analysis_files/PrimateEXCSTR.nex")
Claddis::WriteMorphTNT(PrimatesExclude$STRMRPMatrix, "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/analysis_files/PrimateEXCSTR.tnt")
ape::write.tree(PrimatesExclude$TaxonomyTree, "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/analysis_files/PrimatesEXCTaxonomy.tre")
write.table(PrimatesExclude$SafelyRemovedTaxa, "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/analysis_files/PrimateSTRforEXC.txt", row.names = FALSE)
write.table(PrimatesExclude$MonophyleticTaxa, "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/analysis_files/PrimateSTRforEXC_Monophyly.txt", row.names = FALSE)
write.table(PrimatesExclude$RemovedSourceData, "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/analysis_files/PrimateSTRforEXC_excludedData.txt", row.names = FALSE)
save(PrimatesExclude, file= "~/Dropbox/Mammal_Supertree/ProjectPlanetOfTheApes/analysis_files/PrimatesExclude.RData")
