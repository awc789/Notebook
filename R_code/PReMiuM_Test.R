library(PReMiuM)


inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())

runInfoObj<-profRegr(yModel=inputs$yModel,
                     xModel=inputs$xModel, nSweeps=10, nClusInit=20,
                     nBurn=20, data=inputs$inputData, output="output",
                     covNames = inputs$covNames, outcomeT = inputs$outcomeT,
                     fixedEffectsNames = inputs$fixedEffectNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj) # Distance Matrix
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png", showRelativeRisk = T)


test <- clusSummaryPoissonDiscrete()
