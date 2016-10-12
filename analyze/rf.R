library(randomForest)

rf.mdl <- randomForest(seis.DF$age ~ seis.DF$X_c, ntree=5)
