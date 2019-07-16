library(rpart)
library(gbm)
source("lgd.R")

decompress <- function(data, 
                       actionCntSet, actionFctSet, predictSet,
                       testFlag) {
    res <- list()
    res$train <- subset(data, !testFlag)
    res$test <- subset(data, testFlag)
    res$xTrain <- c()
    res$yTrain <- c()
    res$zTrain <- c()
    res$xTest <- c()
    res$yTest <- c()
    res$zTest <- c()
    
    for (col in actionCntSet) {
        res$zTrain <- cbind(res$zTrain, as.numeric(as.vector(data[[col]]))[!testFlag])
        res$zTest <- cbind(res$zTest, as.numeric(as.vector(data[[col]]))[testFlag])
    }
    for (col in actionFctSet) {
        res$zTrain <- cbind(res$zTrain, as.numeric(as.factor(as.vector(data[[col]])))[!testFlag])
        res$zTest <- cbind(res$zTest, as.numeric(as.factor(as.vector(data[[col]])))[testFlag])
    }
    for (col in predictSet) {
        res$xTrain <- cbind(res$xTrain, as.numeric(as.vector(data[[col]]))[!testFlag])
        res$xTest <- cbind(res$xTest, as.numeric(as.vector(data[[col]]))[testFlag])
    }
    res$yTrain <-  as.numeric(as.vector(data[[target]]))[!testFlag]
    res$yTest <-  as.numeric(as.vector(data[[target]]))[testFlag]
    return(res)
}

c.formula <- function(target, predictSet, actionCntSet = c(), actionFctSet = c()) {
    strPredictSet <- paste(predictSet, collapse="+")
    strActionCntSet <- paste(actionCntSet, collapse="+")
    if (strActionCntSet != "") {
        strActionCntSet <- paste(strActionCntSet, "+", sep="")
    }
    strActionFctSet = ""
    if (!is.null(actionFctSet)) {
        strActionFctSet <- paste(paste("factor(", actionFctSet, ")", sep=""), collapse="+")
        strActionFctSet <- paste(strActionFctSet, "+", sep="")
    } 
    return(as.formula(paste(target, "~", strActionCntSet, strActionFctSet, strPredictSet, sep="")))
}

l2 <- function(a, b) {return(mean((a-b)^2))}

#### DATASET: BEIJINGPM ####
set.seed(42)
id <- "beijingpm"
cat("Working on Dataset: ", id, "...\n")

data <- read.table(paste("Accuracy/", id, ".csv", sep=""), header = TRUE, sep = ",")
data <- subset(data, !is.na(pm2.5))
n <- dim(data)[1]

tenFold <- sample(1:10, size=n, replace=TRUE)
result <- c()
actionCntSet <- c()
actionFctSet <- c("cbwd", "year", "month", "day", "hour")
actionSet <- union(actionCntSet, actionFctSet)
predictSet <- c("DEWP", "TEMP", "PRES", "Iws", "Is", "Ir")
target <- "pm2.5"
gSubsample <- 0.5

modelS <- list()
modelS$diff <- ols.diff
modelS$init <- ols.init.beta
modelS$pred <- ols.pred
modelS$loss <- ols.loss
modelS$dummy <- TRUE
modelS$xscale <- TRUE
modelS$ntree <- 100
modelS$lambda <- 0.05
modelS$method <- "ordinary"
modelS$subsample <- gSubsample

modelS$p <- length(predictSet) + modelS$dummy
modelS$q <- length(actionSet)
modelS$control <- rpart.control(maxdepth = 6, minsplits=100, cp=0.001)
modelS$woods <- list()

modelB <- modelS
modelB$lambda <- 3
modelB$control <- rpart.control(maxdepth = 8, minsplits=5, cp=0.00)
modelB$method <- "boulevard"
modelB$inflate <- 1

for (i in 1:10) {
    gc()
    testFlag <- (tenFold==i)
    simData <- decompress(data, 
                          actionCntSet, actionFctSet, predictSet,
                          testFlag)
    
    #### Linear Model
    tmp <- c()
    lmFit <- lm(c.formula(target, predictSet),
                data = simData$train)
    tmp <- c(tmp, 
             rLmFit = l2(predict(lmFit, newdata = simData$test), simData$yTest))
    
    lmFitSat <- lm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                   data = simData$train)
    tmp <- c(tmp, 
             rLmFitSat = l2(predict(lmFitSat, newdata = simData$test), simData$yTest))
    
    gbmFit <- gbm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                  data = simData$train, 
                  distribution = "gaussian",
                  interaction.depth = modelS$control$maxdepth,
                  bag.fraction = gSubsample,
                  n.trees = modelS$ntree)
    tmp <- c(tmp, 
             rGbmFit = l2(predict(gbmFit, newdata = simData$test, modelS$ntree), simData$yTest))
    
    gc()
    modelS$n = sum(!testFlag)
    svcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model=modelS)
    tmp <- c(tmp,
             rSvcmFit = l2(lgd.predict(x=simData$xTest, z=simData$zTest, model=svcmFit),
                           simData$yTest))
    
    gc()
    modelB$n = sum(!testFlag) 
    bvcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model = modelB)
    tmpPred <- lgd.predict(x=simData$xTest, z=simData$zTest, model=bvcmFit)
    tmp <- c(tmp, rBvcmFitL = l2(tmpPred, simData$yTest),
                  rBvcmFitR = l2(tmpPred * modelB$inflate,  simData$yTest))
    
    result <- rbind(result, tmp)
    cat("SIMULATION ITERATION ", i, "\n", tmp, "\n")
}

save(result, file=paste("Accuracy/", id, ".RData", sep=""))


#### DATASET: BIKEHOUR ####
set.seed(42)
id <- "bikehour"
cat("Working on Dataset: ", id, "...\n")

data <- read.table(paste("Accuracy/", id, ".csv", sep=""), header = TRUE, sep = ",")
data <- subset(data, !is.na(cnt))
n <- dim(data)[1]

tenFold <- sample(1:10, size=n, replace=TRUE)
result <- c()
actionCntSet <- c()
actionFctSet <- c("weathersit", "season", "weekday", "mnth", "hr")
actionSet <- union(actionCntSet, actionFctSet)
predictSet <- c("temp", "atemp", "hum", "windspeed")
target <- "cnt"
gSubsample <- 0.5

modelS <- list()
modelS$diff <- ols.diff
modelS$init <- ols.init.beta
modelS$pred <- ols.pred
modelS$loss <- ols.loss
modelS$dummy <- TRUE
modelS$xscale <- TRUE
modelS$ntree <- 100
modelS$lambda <- 0.05
modelS$method <- "ordinary"
modelS$subsample <- gSubsample

modelS$p <- length(predictSet) + modelS$dummy
modelS$q <- length(actionSet)
modelS$control <- rpart.control(maxdepth = 5, minsplits=100, cp=0.001)
modelS$woods <- list()

modelB <- modelS
modelB$lambda <- 1
modelB$control <- rpart.control(maxdepth = 7, minsplits=5, cp=0.00)
modelB$method <- "boulevard"
modelB$inflate <- 1

for (i in 1:10) {
    gc()
    testFlag <- (tenFold==i)
    simData <- decompress(data, 
                          actionCntSet, actionFctSet, predictSet,
                          testFlag)
    
    #### Linear Model
    tmp <- c()
    lmFit <- lm(c.formula(target, predictSet),
                data = simData$train)
    tmp <- c(tmp, 
             rLmFit = l2(predict(lmFit, newdata = simData$test), simData$yTest))
    
    lmFitSat <- lm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                   data = simData$train)
    tmp <- c(tmp, 
             rLmFitSat = l2(predict(lmFitSat, newdata = simData$test), simData$yTest))
    
    gbmFit <- gbm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                  data = simData$train, 
                  distribution = "gaussian",
                  interaction.depth = modelS$control$maxdepth,
                  bag.fraction = gSubsample,
                  n.trees = modelS$ntree)
    tmp <- c(tmp, 
             rGbmFit = l2(predict(gbmFit, newdata = simData$test, modelS$ntree), simData$yTest))
    
    gc()
    modelS$n = sum(!testFlag)
    svcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model=modelS)
    tmp <- c(tmp,
             rSvcmFit = l2(lgd.predict(x=simData$xTest, z=simData$zTest, model=svcmFit),
                           simData$yTest))
    
    gc()
    modelB$n = sum(!testFlag) 
    bvcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model = modelB)
    tmpPred <- lgd.predict(x=simData$xTest, z=simData$zTest, model=bvcmFit)
    tmp <- c(tmp, rBvcmFitL = l2(tmpPred, simData$yTest),
                  rBvcmFitR = l2(tmpPred * modelB$inflate,  simData$yTest))
    
    result <- rbind(result, tmp)
    cat("SIMULATION ITERATION ", i, "\n", tmp, "\n")
}

save(result, file=paste("Accuracy/", id, ".RData", sep=""))


#### DATASET: STARCRAFT ####
set.seed(42)
id <- "starcraft"
cat("Working on Dataset: ", id, "...\n")

data <- read.table(paste("Accuracy/", id, ".csv", sep=""), header = TRUE, sep = ",")
data <- data[complete.cases(data), ]

result <- c()
actionCntSet <- c("Age", "HoursPerWeek", "UniqueHotkeys", "ComplexUnitsMade")
actionFctSet <- c()
actionSet <- union(actionCntSet, actionFctSet)
predictSet <- c("SelectByHotkeys", "AssignToHotkeys", "MinimapAttacks", "MinimapRightClicks", 
                "GapBetweenPACs", "TotalMapExplored", "WorkersMade", "ComplexAbilitiesUsed")
target <- "LeagueIndex"
for (key in union(actionCntSet, predictSet)) {
    if (is.factor(data[[key]])) {
        data[[key]] = as.numeric(levels(data[[key]]))[data[[key]]]
    }
}
data <- data[complete.cases(data), ]
n <- dim(data)[1]
tenFold <- sample(1:10, size=n, replace=TRUE)

gSubsample <- 0.8

modelS <- list()
modelS$diff <- ols.diff
modelS$init <- ols.init.beta
modelS$pred <- ols.pred
modelS$loss <- ols.loss
modelS$dummy <- TRUE
modelS$xscale <- TRUE
modelS$ntree <- 100
modelS$lambda <- 0.05
modelS$method <- "ordinary"
modelS$subsample <- gSubsample

modelS$p <- length(predictSet) + modelS$dummy
modelS$q <- length(actionSet)
modelS$control <- rpart.control(maxdepth = 4, minsplits=10, cp=0.001)
modelS$woods <- list()

modelB <- modelS
modelB$lambda <- 0.8
modelB$control <- rpart.control(maxdepth = 6, minsplits=5, cp=0.00)
modelB$method <- "boulevard"
modelB$inflate <- 1

for (i in 1:10) {
    gc()
    testFlag <- (tenFold==i)
    simData <- decompress(data, 
                          actionCntSet, actionFctSet, predictSet,
                          testFlag)
    
    #### Linear Model
    tmp <- c()
    lmFit <- lm(c.formula(target, predictSet),
                data = simData$train)
    tmp <- c(tmp, 
             rLmFit = l2(predict(lmFit, newdata = simData$test), simData$yTest))
    
    lmFitSat <- lm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                   data = simData$train)
    tmp <- c(tmp, 
             rLmFitSat = l2(predict(lmFitSat, newdata = simData$test), simData$yTest))
    
    gbmFit <- gbm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                  data = simData$train, 
                  distribution = "gaussian",
                  interaction.depth = modelS$control$maxdepth,
                  bag.fraction = gSubsample,
                  n.trees = modelS$ntree)
    tmp <- c(tmp, 
             rGbmFit = l2(predict(gbmFit, newdata = simData$test, modelS$ntree), simData$yTest))
    
    gc()
    modelS$n = sum(!testFlag)
    svcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model=modelS)
    tmp <- c(tmp,
             rSvcmFit = l2(lgd.predict(x=simData$xTest, z=simData$zTest, model=svcmFit),
                           simData$yTest))
    
    gc()
    modelB$n = sum(!testFlag) 
    bvcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model = modelB)
    tmpPred <- lgd.predict(x=simData$xTest, z=simData$zTest, model=bvcmFit)
    tmp <- c(tmp, rBvcmFitL = l2(tmpPred, simData$yTest),
                  rBvcmFitR = l2(tmpPred * modelB$inflate,  simData$yTest))
    
    result <- rbind(result, tmp)
    cat("SIMULATION ITERATION ", i, "\n", tmp, "\n")
}

save(result, file=paste("Accuracy/", id, ".RData", sep=""))
#### DATASET: ONLINENEWS ####
set.seed(42)
id <- "onlinenews"
cat("Working on Dataset: ", id, "...\n")

data <- read.table(paste("Accuracy/", id, ".csv", sep=""), header = TRUE, sep = ",")
data <- data[complete.cases(data), ]

result <- c()
actionCntSet <- c("n_tokens_title", "n_tokens_content", "num_imgs")
actionFctSet <- c("is_weekend")
actionSet <- union(actionCntSet, actionFctSet)
predictSet <- c("self_reference_avg_sharess", 
                "global_subjectivity", "global_sentiment_polarity", 
                "global_rate_positive_words", "global_rate_negative_words")
target <- "shares"
for (key in union(actionCntSet, predictSet)) {
    if (is.factor(data[[key]])) {
        data[[key]] <- as.numeric(levels(data[[key]]))[data[[key]]]
    }
}
data <- data[complete.cases(data), ]
data[[target]] <- log(data[[target]])
n <- dim(data)[1]
tenFold <- sample(1:10, size=n, replace=TRUE)

gSubsample <- 0.5

modelS <- list()
modelS$diff <- ols.diff
modelS$init <- ols.init.beta
modelS$pred <- ols.pred
modelS$loss <- ols.loss
modelS$dummy <- TRUE
modelS$xscale <- TRUE
modelS$ntree <- 100
modelS$lambda <- 0.05
modelS$method <- "ordinary"
modelS$subsample <- gSubsample

modelS$p <- length(predictSet) + modelS$dummy
modelS$q <- length(actionSet)
modelS$control <- rpart.control(maxdepth = 7, minsplits=10, cp=0.001)
modelS$woods <- list()

modelB <- modelS
modelB$lambda <- 1
modelB$control <- rpart.control(maxdepth = 9, minsplits=10, cp=0.00)
modelB$method <- "boulevard"
modelB$inflate <- 1

for (i in 1:10) {
    gc()
    testFlag <- (tenFold==i)
    simData <- decompress(data, 
                          actionCntSet, actionFctSet, predictSet,
                          testFlag)
    
    #### Linear Model
    tmp <- c()
    lmFit <- lm(c.formula(target, predictSet),
                data = simData$train)
    tmp <- c(tmp, 
             rLmFit = l2(predict(lmFit, newdata = simData$test), simData$yTest))
    
    lmFitSat <- lm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                   data = simData$train)
    tmp <- c(tmp, 
             rLmFitSat = l2(predict(lmFitSat, newdata = simData$test), simData$yTest))
    
    gbmFit <- gbm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                  data = simData$train, 
                  distribution = "gaussian",
                  interaction.depth = modelS$control$maxdepth,
                  bag.fraction = gSubsample,
                  n.trees = modelS$ntree)
    tmp <- c(tmp, 
             rGbmFit = l2(predict(gbmFit, newdata = simData$test, modelS$ntree), simData$yTest))
    
    gc()
    modelS$n = sum(!testFlag)
    svcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model=modelS)
    tmp <- c(tmp,
             rSvcmFit = l2(lgd.predict(x=simData$xTest, z=simData$zTest, model=svcmFit),
                           simData$yTest))
    
    gc()
    modelB$n = sum(!testFlag) 
    bvcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model = modelB)
    tmpPred <- lgd.predict(x=simData$xTest, z=simData$zTest, model=bvcmFit)
    tmp <- c(tmp, rBvcmFitL = l2(tmpPred, simData$yTest),
             rBvcmFitR = l2(tmpPred * modelB$inflate,  simData$yTest))
    
    result <- rbind(result, tmp)
    cat("SIMULATION ITERATION ", i, "\n", tmp, "\n")
}

save(result, file=paste("Accuracy/", id, ".RData", sep=""))

#### DATASET: ENERGYEFFICIENCY ####
set.seed(42)
id <- "energyefficiency"
cat("Working on Dataset: ", id, "...\n")

data <- read.table(paste("Accuracy/", id, ".csv", sep=""), header = TRUE, sep = ",")
data <- subset(data, select = -c(11, 12))
data <- data[complete.cases(data), ]

result <- c()
actionCntSet <- c("X7")
actionFctSet <- c("X6", "X8")
actionSet <- union(actionCntSet, actionFctSet)
predictSet <- c("X1", "X2", "X3")
target <- "Y1"
for (key in union(actionCntSet, predictSet)) {
    if (is.factor(data[[key]])) {
        data[[key]] = as.numeric(levels(data[[key]]))[data[[key]]]
    }
}
data <- data[complete.cases(data), ]
n <- dim(data)[1]
tenFold <- sample(1:10, size=n, replace=TRUE)

gSubsample <- 0.5

modelS <- list()
modelS$diff <- ols.diff
modelS$init <- ols.init.beta
modelS$pred <- ols.pred
modelS$loss <- ols.loss
modelS$dummy <- TRUE
modelS$xscale <- TRUE
modelS$ntree <- 100
modelS$lambda <- 0.05
modelS$method <- "ordinary"
modelS$subsample <- gSubsample

modelS$p <- length(predictSet) + modelS$dummy
modelS$q <- length(actionSet)
modelS$control <- rpart.control(maxdepth = 2, minsplits=10, cp=0.001)
modelS$woods <- list()

modelB <- modelS
modelB$lambda <- 1
modelB$control <- rpart.control(maxdepth = 3, minsplits=10, cp=0.00)
modelB$method <- "boulevard"
modelB$inflate <- 1

for (i in 1:10) {
    gc()
    testFlag <- (tenFold==i)
    simData <- decompress(data, 
                          actionCntSet, actionFctSet, predictSet,
                          testFlag)
    
    #### Linear Model
    tmp <- c()
    lmFit <- lm(c.formula(target, predictSet),
                data = simData$train)
    tmp <- c(tmp, 
             rLmFit = l2(predict(lmFit, newdata = simData$test), simData$yTest))
    
    lmFitSat <- lm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                   data = simData$train)
    tmp <- c(tmp, 
             rLmFitSat = l2(predict(lmFitSat, newdata = simData$test), simData$yTest))
    
    gbmFit <- gbm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                  data = simData$train, 
                  distribution = "gaussian",
                  interaction.depth = modelS$control$maxdepth,
                  bag.fraction = gSubsample,
                  n.trees = modelS$ntree)
    tmp <- c(tmp, 
             rGbmFit = l2(predict(gbmFit, newdata = simData$test, modelS$ntree), simData$yTest))
    
    gc()
    modelS$n = sum(!testFlag)
    svcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model=modelS)
    tmp <- c(tmp,
             rSvcmFit = l2(lgd.predict(x=simData$xTest, z=simData$zTest, model=svcmFit),
                           simData$yTest))
    
    gc()
    modelB$n = sum(!testFlag) 
    bvcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model = modelB)
    tmpPred <- lgd.predict(x=simData$xTest, z=simData$zTest, model=bvcmFit)
    tmp <- c(tmp, rBvcmFitL = l2(tmpPred, simData$yTest),
             rBvcmFitR = l2(tmpPred * modelB$inflate,  simData$yTest))
    
    result <- rbind(result, tmp)
    cat("SIMULATION ITERATION ", i, "\n", tmp, "\n")
}

save(result, file=paste("Accuracy/", id, ".RData", sep=""))

#### DATASET: EGRIDSTAB ####

set.seed(42)
id <- "egridstab"
cat("Working on Dataset: ", id, "...\n")

data <- read.table(paste("Accuracy/", id, ".csv", sep=""), header = TRUE, sep = ",")
data <- data[complete.cases(data), ]

result <- c()
actionCntSet <- c("tau2", "tau3", "tau4")
actionFctSet <- c()
actionSet <- union(actionCntSet, actionFctSet)
predictSet <- c("p2", "p3", "p4", "g2", "g3", "g4")
target <- "stab"
for (key in union(actionCntSet, predictSet)) {
    if (is.factor(data[[key]])) {
        data[[key]] = as.numeric(levels(data[[key]]))[data[[key]]]
    }
}
data <- data[complete.cases(data), ]
n <- dim(data)[1]
tenFold <- sample(1:10, size=n, replace=TRUE)

gSubsample <- 0.5

modelS <- list()
modelS$diff <- ols.diff
modelS$init <- ols.init.beta
modelS$pred <- ols.pred
modelS$loss <- ols.loss
modelS$dummy <- TRUE
modelS$xscale <- TRUE
modelS$ntree <- 100
modelS$lambda <- 0.05
modelS$method <- "ordinary"
modelS$subsample <- gSubsample

modelS$p <- length(predictSet) + modelS$dummy
modelS$q <- length(actionSet)
modelS$control <- rpart.control(maxdepth = 4, minsplits=10, cp=0.001)
modelS$woods <- list()

modelB <- modelS
modelB$lambda <- 1
modelB$control <- rpart.control(maxdepth = 6, minsplits=10, cp=0.00)
modelB$method <- "boulevard"
modelB$inflate <- 1

for (i in 1:10) {
    gc()
    testFlag <- (tenFold==i)
    simData <- decompress(data, 
                          actionCntSet, actionFctSet, predictSet,
                          testFlag)
    
    #### Linear Model
    tmp <- c()
    lmFit <- lm(c.formula(target, predictSet),
                data = simData$train)
    tmp <- c(tmp, 
             rLmFit = l2(predict(lmFit, newdata = simData$test), simData$yTest))
    
    lmFitSat <- lm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                   data = simData$train)
    tmp <- c(tmp, 
             rLmFitSat = l2(predict(lmFitSat, newdata = simData$test), simData$yTest))
    
    gbmFit <- gbm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                  data = simData$train, 
                  distribution = "gaussian",
                  interaction.depth = modelS$control$maxdepth,
                  bag.fraction = gSubsample,
                  n.trees = modelS$ntree)
    tmp <- c(tmp, 
             rGbmFit = l2(predict(gbmFit, newdata = simData$test, modelS$ntree), simData$yTest))
    
    gc()
    modelS$n = sum(!testFlag)
    svcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model=modelS)
    tmp <- c(tmp,
             rSvcmFit = l2(lgd.predict(x=simData$xTest, z=simData$zTest, model=svcmFit),
                           simData$yTest))
    
    gc()
    modelB$n = sum(!testFlag) 
    bvcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model = modelB)
    tmpPred <- lgd.predict(x=simData$xTest, z=simData$zTest, model=bvcmFit)
    tmp <- c(tmp, rBvcmFitL = l2(tmpPred, simData$yTest),
                  rBvcmFitR = l2(tmpPred * modelB$inflate,  simData$yTest))
    
    result <- rbind(result, tmp)
    cat("SIMULATION ITERATION ", i, "\n", tmp, "\n")
}

save(result, file=paste("Accuracy/", id, ".RData", sep=""))

