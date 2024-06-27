library(CAGEr)
library(BSgenome.Hkshv.encode.p13)
library(rtracklayer)


iSLK <- readRDS("iSLK_CAGEr_RAMPAGE.rds")

BCBL1 <- readRDS("BCBL1_CAGEr_RAMPAGE.rds")


anno <- import.gff2("gencode.v46.primary_assembly.annotation.gtf")

iSLK <- annotateCTSS(iSLK, anno)
colData(iSLK)[,c("librarySizes", "promoter", "exon", "intron", "unknown")]
plotAnnot(iSLK, "counts")
iSLK <- mergeSamples(iSLK, mergeIndex = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6), 
                   mergedSampleLabels = c("Latent", "h12", "h24", "h48", "h72", "h96"))
iSLK <- annotateCTSS(iSLK, anno)

plotAnnot(iSLK, "counts")
librarySizes(iSLK)

plotReverseCumulatives(iSLK, fitInRange = c(5, 100000), onePlot = TRUE)
iSLK <- normalizeTagCount(iSLK, method = "powerLaw", fitInRange = c(5, 100000), alpha = 1.1, T = 5*10^5)

iSLK[["tagCountMatrix"]]
iSLK <- clusterCTSS( iSLK
                 , threshold = 1
                 , thresholdIsTpm = TRUE
                 , nrPassThreshold = 1
                 , method = "distclu"
                 , maxDist = 20
                 , removeSingletons = TRUE
                 , keepSingletonsAbove = 5)

iSLK <- cumulativeCTSSdistribution(iSLK, clusters = "tagClusters", useMulticore = T)
iSLK <- quantilePositions(iSLK, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
plotInterquantileWidth(iSLK, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)


iSLK <- aggregateTagClusters(iSLK, tpmThreshold = 1, qLow = 0.1, qUp = 0.9, maxDist = 100)
iSLK$outOfClusters / iSLK$librarySizes *100
consensusClustersGR(iSLK)

# the same analysis for BCBL1 cells
BCBL1 <- annotateCTSS(BCBL1, anno)
plotAnnot(BCBL1, "counts")
librarySizes(BCBL1)
BCBL1 <- mergeSamples(BCBL1, mergeIndex = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6), 
		   mergedSampleLabels = c("Latent", "h12", "h24", "h48", "h72", "h96"))
BCBL1 <- annotateCTSS(BCBL1, anno)
plotAnnot(BCBL1, "counts")
plotReverseCumulatives(BCBL1, fitInRange = c(5, 100000), onePlot = TRUE)
BCBL1 <- normalizeTagCount(BCBL1, method = "powerLaw", fitInRange = c(5, 100000), alpha = 1.1, T = 5*10^5)
BCBL1[["tagCountMatrix"]]

BCBL1 <- clusterCTSS( BCBL1
		    , threshold = 1
		    , thresholdIsTpm = TRUE
		    , nrPassThreshold = 1
		    , method = "distclu"
		    , maxDist = 20
		    , removeSingletons = TRUE
		    , keepSingletonsAbove = 5)
BCBL1 <- cumulativeCTSSdistribution(BCBL1, clusters = "tagClusters", useMulticore = T)
BCBL1 <- quantilePositions(BCBL1, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
plotInterquantileWidth(BCBL1, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)
BCBL1 <- aggregateTagClusters(BCBL1, tpmThreshold = 1, qLow = 0.1, qUp = 0.9, maxDist = 100)
BCBL1$outOfClusters / BCBL1$librarySizes *100
consensusClustersGR(BCBL1)
