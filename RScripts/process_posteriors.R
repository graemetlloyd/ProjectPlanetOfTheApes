### Primate Bayes processing ##

setwd("ProjectPlanetOfTheApes/timetree/BeastOutput")
log_1 <- read.table("Primate_BEAST-1602558368884.log", header = T, stringsAsFactors = F)
log_2 <- read.table("Primate_BEAST-1602988762555.log", header = T, stringsAsFactors = F)
log_3 <- read.table("Primate_BEAST-1602988693513.log", header = T, stringsAsFactors = F)
log_1$Sample

library(coda)
log1mcmc<-mcmc(data= log_1, start = 0, end = 500000000, thin = 500000)
log2mcmc<-mcmc(data= log_2, start = 0, end = 500000000, thin = 500000)
log3mcmc<-mcmc(data= log_3, start = 0, end = 500000000, thin = 500000)

comb.mcmc<-mcmc.list(log1mcmc, log2mcmc, log3mcmc)

gelman.diag(comb.mcmc, multivariate = FALSE, autoburnin = F) # Ok!

round(effectiveSize(comb.mcmc), 0)

burnin= ceiling(0.1*nrow(log_1))

comblog<-cbind(log_1[-(1:burnin),], log_2[-(1:burnin),],log_3[-(1:burnin),])
write.table(comblog, file="combined_log.log")
par(mfrow=c(1,3))

### trees 
library(ape)
library(phangorn)
tree_1 <- read.nexus("Primate_BEAST-1602558368884.trees")
tree_2 <- read.nexus("Primate_BEAST-1602988762555.trees")
tree_3 <- read.nexus("Primate_BEAST-1602988693513.trees")

combtrees <- append(append(tree_1[-(1:burnin)], tree_2[-(1:burnin)]), tree_3[-(1:burnin)])
class(combtrees)
write.tree(combtrees, "Primate_Beast_combined.trees")
MAPtree<-combtrees[[which.max(comblog$posterior)]]
write.tree(MAPtree, "Primate_BEAST_MAP.tre")


fossils<-setdiff(MAPtree$tip.label, paleotree:::dropExtinct(MAPtree)$tip.label)
extant <- lapply(combtrees, drop.tip, tip=fossils )
class(extant) <- class(combtrees)
write.tree(extant,"Primate_BEAST_extant.trees")

subsample <-sample(x = seq(1, length(combtrees)), size=100,replace = F)
treesubsample<-combtrees[subsample]
names(treesubsample) <- paste("tree_", subsample, sep = "")

write.tree(treesubsample, "Primate_BEAST_subsample.trees")

#####



phy <- read.tree("Primate_BEAST_median.tree")
extant<-drop.fossil(phy,tol = 0.0005)
write.tree(ladderize(extant), "PrimatesExtantBiogeog.tre")

write.tree(ladderize(drop.tip(phy, extant$tip.label)), "PrimatesFossilBiogeog.tre")
write.tree(ladderize(phy), "PrimatesAllBiogeog.tre")
