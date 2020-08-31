#!/usr/bin/env Rscript
#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO")
library("topGO")
#biocLite("Rgraphviz")
library("Rgraphviz")

#Run from commandline as: Rscript path/topgo_pipe.R path/annotations_topgo.txt path/interestinggenes_file 
args = commandArgs(trailingOnly=TRUE)
geneID2GO <- readMappings(args[1], sep = "\t", IDsep = ",") 
genesOfInterest <- read.table(args[2],header=FALSE)
geneUniverse <- names(geneID2GO)
genesOfInterest <- as.character(genesOfInterest$V1)
head(genesOfInterest)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myBPGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myMFGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myCCGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
#Look at the object
myBPGOdata
myMFGOdata
myCCGOdata
#Access significant genes
#BP
sg_bp <- sigGenes(myBPGOdata)
str(sg_bp)
numSigGenes(myBPGOdata) 
#MF
sg_mf <- sigGenes(myMFGOdata)
str(sg_mf)
numSigGenes(myMFGOdata) 
#CC
sg_cc <- sigGenes(myCCGOdata)
str(sg_cc)
numSigGenes(myCCGOdata)

############ Performing enrichment tests ############ 
#Fishers exact test
#ontology = BP
resultFisherBP <- runTest(myBPGOdata, algorithm="classic", statistic="fisher") 
resultFisherBP
#Selecting 'algorithm=classic' means that the GO hierarchy isn't taken into account, 
#so each GO term is tested independently.
resultFisherBP #The p-values have not been corrected for multiple testing.
#ontology = MF
resultFisherMF <- runTest(myMFGOdata, algorithm="classic", statistic="fisher") 
resultFisherMF #The p-values have not been corrected for multiple testing.
#ontology = CC
resultFisherCC <- runTest(myCCGOdata, algorithm="classic", statistic="fisher") 
resultFisherCC

#The GenTable function returns a data frame containing the top topNodes GO terms identified by the classic algorithm,
#BP
FisherResBP <- GenTable(myBPGOdata, classicFisher = resultFisherBP, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 22)
write.csv(FisherResBP, file = here("annotation/FisherResBP.csv"))
#MF
FisherResMF <- GenTable(myMFGOdata, classicFisher = resultFisherMF, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 22)
write.csv(FisherResMF, file = here("annotation/FisherResMF.csv"))
#CC
FisherResCC <- GenTable(myCCGOdata, classicFisher = resultFisherCC, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 22)
write.csv(FisherResCC, file = here("annotation/FisherResCC.csv"))
