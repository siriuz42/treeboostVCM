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

ce <- function(b, a) {
    b[b>0.9999] <- 0.9999;
    b[b<0.0001] <- 0.0001;
    return(mean(-a * log(b) - (1-a) * log(1-b)) / (-log(0.5)))
}

cep_DEPRECATED <- function(b, a) {
    b[b>10] <- 10
    b[b<-10] <- -10
    b <- sigmoid(b)
    return(mean(-a * log(b) - (1-a) * log(1-b)) / (-log(0.5)))
}

cep <- function(b, a) {
    b[b>10] <- 10
    b[b<-10] <- -10
    b <- sigmoid(b)
    return(mean( a == (b > 0.5)))
}


#### DATASET: MAGIC ####
set.seed(42)
id <- "magic04"
cat("Working on Dataset: ", id, "...\n")

data <- read.table(paste("Accuracy/", id, ".csv", sep=""), header = FALSE, sep = ",")
data <- data[complete.cases(data), ]
names(data) <- c("fLength", "fWidth", "fSize", "fConc", "fConc1", "fAsym", "fM3Long", "fM3Trans", "fAlpha", "fDist", "class")
data$target <- as.numeric(data$class == "g")

result <- c()
actionCntSet <- c("fSize", "fConc", "fConc1")
actionFctSet <- c()
actionSet <- union(actionCntSet, actionFctSet)
predictSet <- c("fLength", "fWidth", "fAsym", "fM3Long", "fM3Trans", "fAlpha", "fDist")
target <- "target"
gSubsample <- 1

tenFold <- sample(1:10, size=n, replace=TRUE)

modelS <- list()
modelS$diff <- lg.diff
modelS$init <- lg.init.beta
modelS$pred <- lg.pred
modelS$loss <- lg.loss
modelS$dummy <- TRUE
modelS$xscale <- TRUE
modelS$ntree <- 50
modelS$lambda <- 0.01
modelS$method <- "ordinary"
modelS$subsample <- gSubsample
modelS$quantile_truncate <- 0.05

modelS$p <- length(predictSet) + modelS$dummy
modelS$q <- length(actionSet)
modelS$control <- rpart.control(maxdepth = 6, minsplits=100, cp=0.01)
modelS$woods <- list()

modelB <- modelS
modelB$lambda <- 1
modelB$control <- rpart.control(maxdepth = 8, minsplits=100, cp=0.00)
modelB$method <- "boulevard"

for (i in 1:10) {
    gc()
    testFlag <- (tenFold==i)
    simData <- decompress(data, 
                          actionCntSet, actionFctSet, predictSet,
                          testFlag)
    
    #### Linear Model
    tmp <- c()
    lmFit <- glm(c.formula(target, predictSet),
                 family = binomial(link = "logit"),
                 data = simData$train)
    tmp <- c(tmp, 
             rLmFit = cep(predict(lmFit, newdata = simData$test, type="link"), simData$yTest))
    
    lmFitSat <- glm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                    family = binomial(link = "logit"),
                    data = simData$train)
    tmp <- c(tmp, 
             rLmFitSat = cep(predict(lmFitSat, newdata = simData$test, type="link"), simData$yTest))
    
    gbmFit <- gbm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                  data = simData$train, 
                  distribution = "adaboost",
                  interaction.depth = modelS$control$maxdepth,
                  bag.fraction = gSubsample,
                  n.trees = modelS$ntree)
    tmp <- c(tmp, 
             rGbmFit = cep(predict(gbmFit, newdata = simData$test, modelS$ntree, type="link"), simData$yTest))
    
    gc()
    modelS$n = sum(!testFlag)
    svcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model=modelS)
    tmp <- c(tmp,
             rSvcmFit = cep(lgd.predict(x=simData$xTest, z=simData$zTest, model=svcmFit),
                            simData$yTest))
    
    gc()
    modelB$n = sum(!testFlag) 
    bvcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model = modelB)
    tmpPred <- lgd.predict(x=simData$xTest, z=simData$zTest, model=bvcmFit)
    tmp <- c(tmp, rBvcmFitL = cep(tmpPred, simData$yTest),
             rBvcmFitR = cep(tmpPred * (1+modelB$lambda) / modelB$lambda,  simData$yTest))
    
    result <- rbind(result, tmp)
    cat("SIMULATION ITERATION ", i, "\n", tmp, "\n")
}

save(result, file=paste("Accuracy/", id, ".RData", sep=""))

#### DATASET: BANK ####
set.seed(42)
id <- "bank"
cat("Working on Dataset: ", id, "...\n")

data <- read.table(paste("Accuracy/", id, ".csv", sep=""), header = TRUE, sep = ";")
data <- data[complete.cases(data), ]
data$target <- as.numeric(data$y == "yes")
data$pdays <- log(data$pdays + 2)
data$previous <- log(data$previous + 1)

result <- c()
actionCntSet <- c()
actionFctSet <- c("job", "marital", "education", "poutcome", "housing", "loan")
actionSet <- union(actionCntSet, actionFctSet)
predictSet <- c("age", "balance", "duration", "pdays", "previous")
target <- "target"
gSubsample <- 1

tenFold <- sample(1:10, size=n, replace=TRUE)

modelS <- list()
modelS$diff <- lg.diff
modelS$init <- lg.init.beta
modelS$pred <- lg.pred
modelS$loss <- lg.loss
modelS$dummy <- TRUE
modelS$xscale <- TRUE
modelS$ntree <- 100
modelS$lambda <- 0.05
modelS$method <- "ordinary"
modelS$subsample <- gSubsample
modelS$quantile_cut <- 0.05

modelS$p <- length(predictSet) + modelS$dummy
modelS$q <- length(actionSet)
modelS$control <- rpart.control(maxdepth = 7, minsplits=100, cp=0.001)
modelS$woods <- list()

modelB <- modelS
modelB$lambda <- 1
modelB$control <- rpart.control(maxdepth = 8, minsplits=10, cp=0.00)
modelB$method <- "boulevard"

for (i in 1:10) {
    gc()
    testFlag <- (tenFold==i)
    simData <- decompress(data, 
                          actionCntSet, actionFctSet, predictSet,
                          testFlag)
    
    #### Linear Model
    tmp <- c()
    lmFit <- glm(c.formula(target, predictSet),
                 family = binomial(link = "logit"),
                 data = simData$train)
    tmp <- c(tmp, 
             rLmFit = cep(predict(lmFit, newdata = simData$test, type="link"), simData$yTest))
    
    lmFitSat <- glm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                    family = binomial(link = "logit"),
                    data = simData$train)
    tmp <- c(tmp, 
             rLmFitSat = cep(predict(lmFitSat, newdata = simData$test, type="link"), simData$yTest))
    
    gbmFit <- gbm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                  data = simData$train, 
                  distribution = "adaboost",
                  interaction.depth = modelS$control$maxdepth,
                  bag.fraction = gSubsample,
                  n.trees = modelS$ntree)
    tmp <- c(tmp, 
             rGbmFit = cep(predict(gbmFit, newdata = simData$test, modelS$ntree, type="link"), simData$yTest))
    
    gc()
    modelS$n = sum(!testFlag)
    svcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model=modelS)
    tmp <- c(tmp,
             rSvcmFit = cep(lgd.predict(x=simData$xTest, z=simData$zTest, model=svcmFit),
                            simData$yTest))
    
    gc()
    modelB$n = sum(!testFlag) 
    bvcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model = modelB)
    tmpPred <- lgd.predict(x=simData$xTest, z=simData$zTest, model=bvcmFit)
    tmp <- c(tmp, rBvcmFitL = cep(tmpPred, simData$yTest),
             rBvcmFitR = cep(tmpPred * (1+modelB$lambda) / modelB$lambda,  simData$yTest))
    
    result <- rbind(result, tmp)
    cat("SIMULATION ITERATION ", i, "\n", tmp, "\n")
}


save(result, file=paste("Accuracy/", id, ".RData", sep=""))


#### DATASET: OCCUPANCY ####
set.seed(42)
id <- "occupancy"
cat("Working on Dataset: ", id, "...\n")

data <- read.table(paste("Accuracy/", id, ".csv", sep=""), header = TRUE, sep = ",")
data <- data[complete.cases(data), ]
data$month <- strptime(data$date, format = "%Y-%m-%d %H:%M:%S")$mon
data$day <- strptime(data$date, format = "%Y-%m-%d %H:%M:%S")$mday
data$weekday <- strptime(data$date, format = "%Y-%m-%d %H:%M:%S")$wday
data$hour <- strptime(data$date, format = "%Y-%m-%d %H:%M:%S")$hour
data$min <- strptime(data$date, format = "%Y-%m-%d %H:%M:%S")$min
data$target <- as.numeric(data$Occupancy == "1")
data <- data[complete.cases(data), ]

result <- c()
actionCntSet <- c("hour", "min")
actionFctSet <- c("weekday")
actionSet <- union(actionCntSet, actionFctSet)
predictSet <- c("Temperature", "Humidity", "Light", "CO2", "HumidityRatio")
target <- "target"
gSubsample <- 0.5

n <- nrow(data)
tenFold <- sample(1:10, size=n, replace=TRUE)

modelS <- list()
modelS$diff <- lg.diff
modelS$init <- lg.init.beta
modelS$pred <- lg.pred
modelS$loss <- lg.loss
modelS$dummy <- TRUE
modelS$xscale <- TRUE
modelS$ntree <- 100
modelS$lambda <- 0.1
modelS$method <- "ordinary"
modelS$subsample <- gSubsample
modelS$quantile_cut <- 0.05

modelS$p <- length(predictSet) + modelS$dummy
modelS$q <- length(actionSet)
modelS$control <- rpart.control(maxdepth = 7, minsplits=100, cp=0.001)
modelS$woods <- list()

modelB <- modelS
modelB$lambda <- 1
modelB$control <- rpart.control(maxdepth = 8, minsplits=100, cp=0.00)
modelB$method <- "boulevard"

for (i in 1:10) {
    gc()
    testFlag <- (tenFold==i)
    simData <- decompress(data, 
                          actionCntSet, actionFctSet, predictSet,
                          testFlag)
    
    #### Linear Model
    tmp <- c()
    lmFit <- glm(c.formula(target, predictSet),
                 family = binomial(link = "logit"),
                 data = simData$train)
    tmp <- c(tmp, 
             rLmFit = cep(predict(lmFit, newdata = simData$test, type="link"), simData$yTest))
    
    lmFitSat <- glm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                    family = binomial(link = "logit"),
                    data = simData$train)
    tmp <- c(tmp, 
             rLmFitSat = cep(predict(lmFitSat, newdata = simData$test, type="link"), simData$yTest))
    
    gbmFit <- gbm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                  data = simData$train, 
                  distribution = "adaboost",
                  interaction.depth = modelS$control$maxdepth,
                  bag.fraction = gSubsample,
                  n.trees = modelS$ntree)
    tmp <- c(tmp, 
             rGbmFit = cep(predict(gbmFit, newdata = simData$test, modelS$ntree, type="link"), simData$yTest))
    
    gc()
    modelS$n = sum(!testFlag)
    svcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model=modelS)
    tmp <- c(tmp,
             rSvcmFit = cep(lgd.predict(x=simData$xTest, z=simData$zTest, model=svcmFit),
                            simData$yTest))
    
    gc()
    modelB$n = sum(!testFlag) 
    bvcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model = modelB)
    tmpPred <- lgd.predict(x=simData$xTest, z=simData$zTest, model=bvcmFit)
    tmp <- c(tmp, rBvcmFitL = cep(tmpPred, simData$yTest),
             rBvcmFitR = cep(tmpPred * (1+modelB$lambda) / modelB$lambda,  simData$yTest))
    
    result <- rbind(result, tmp)
    cat("SIMULATION ITERATION ", i, "\n", tmp, "\n")
}

save(result, file=paste("Accuracy/", id, ".RData", sep=""))

#### DATASET: SPAMBASE ####
set.seed(42)
id <- "spambase"
cat("Working on Dataset: ", id, "...\n")

data <- read.table(paste("Accuracy/", id, ".csv", sep=""), header = FALSE, sep = ",")
data$target <- as.numeric(data$V58 == "1")
data <- data[complete.cases(data), ]

result <- c()
actionCntSet <- c("V55", "V56", "V57")
actionFctSet <- c()
actionSet <- union(actionCntSet, actionFctSet)
predictSet <- paste("V", 1:54, sep="")
target <- "target"
gSubsample <- 0.5

n <- nrow(data)
tenFold <- sample(1:10, size=n, replace=TRUE)

modelS <- list()
modelS$diff <- lg.diff
modelS$init <- lg.init.beta
modelS$pred <- lg.pred
modelS$loss <- lg.loss
modelS$dummy <- TRUE
modelS$xscale <- TRUE
modelS$ntree <- 100
modelS$lambda <- 0.1
modelS$method <- "ordinary"
modelS$subsample <- gSubsample
modelS$quantile_cut <- 0.05

modelS$p <- length(predictSet) + modelS$dummy
modelS$q <- length(actionSet)
modelS$control <- rpart.control(maxdepth = 5, minsplits=100, cp=0.001)
modelS$woods <- list()

modelB <- modelS
modelB$lambda <- 1
modelB$control <- rpart.control(maxdepth = 7, minsplits=50, cp=0.00)
modelB$method <- "boulevard"

for (i in 1:10) {
    gc()
    testFlag <- (tenFold==i)
    simData <- decompress(data, 
                          actionCntSet, actionFctSet, predictSet,
                          testFlag)
    
    #### Linear Model
    tmp <- c()
    lmFit <- glm(c.formula(target, predictSet),
                 family = binomial(link = "logit"),
                 data = simData$train)
    tmp <- c(tmp, 
             rLmFit = cep(predict(lmFit, newdata = simData$test, type="link"), simData$yTest))
    
    lmFitSat <- glm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                    family = binomial(link = "logit"),
                    data = simData$train)
    tmp <- c(tmp, 
             rLmFitSat = cep(predict(lmFitSat, newdata = simData$test, type="link"), simData$yTest))
    
    gbmFit <- gbm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                  data = simData$train, 
                  distribution = "adaboost",
                  interaction.depth = modelS$control$maxdepth,
                  bag.fraction = gSubsample,
                  n.trees = modelS$ntree)
    tmp <- c(tmp, 
             rGbmFit = cep(predict(gbmFit, newdata = simData$test, modelS$ntree, type="link"), simData$yTest))
    
    gc()
    modelS$n = sum(!testFlag)
    svcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model=modelS)
    tmp <- c(tmp,
             rSvcmFit = cep(lgd.predict(x=simData$xTest, z=simData$zTest, model=svcmFit),
                            simData$yTest))
    
    gc()
    modelB$n = sum(!testFlag) 
    bvcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model = modelB)
    tmpPred <- lgd.predict(x=simData$xTest, z=simData$zTest, model=bvcmFit)
    tmp <- c(tmp, rBvcmFitL = cep(tmpPred, simData$yTest),
             rBvcmFitR = cep(tmpPred * (1+modelB$lambda) / modelB$lambda,  simData$yTest))
    
    result <- rbind(result, tmp)
    cat("SIMULATION ITERATION ", i, "\n", tmp, "\n")
}

save(result, file=paste("Accuracy/", id, ".RData", sep=""))
#### DATASET: ADULT ####
set.seed(42)
id <- "adult"
cat("Working on Dataset: ", id, "...\n")

data <- read.table(paste("Accuracy/", id, ".csv", sep=""), header = FALSE, sep = ",")
data$target <- as.numeric(data$V15 == " <=50K")
data <- data[complete.cases(data), ]

result <- c()
actionCntSet <- c("V1")
actionFctSet <- c("V2", "V8", "V9", "V10")
actionSet <- union(actionCntSet, actionFctSet)
predictSet <- c("V11", "V12", "V13", "V5")
target <- "target"
gSubsample <- 0.5

n <- nrow(data)
tenFold <- sample(1:10, size=n, replace=TRUE)

modelS <- list()
modelS$diff <- lg.diff
modelS$init <- lg.init.beta
modelS$pred <- lg.pred
modelS$loss <- lg.loss
modelS$dummy <- TRUE
modelS$xscale <- TRUE
modelS$ntree <- 100
modelS$lambda <- 0.1
modelS$method <- "ordinary"
modelS$subsample <- gSubsample
modelS$quantile_cut <- 0.05

modelS$p <- length(predictSet) + modelS$dummy
modelS$q <- length(actionSet)
modelS$control <- rpart.control(maxdepth = 6, minsplits=100, cp=0.001)
modelS$woods <- list()

modelB <- modelS
modelB$lambda <- 1
modelB$control <- rpart.control(maxdepth = 7, minsplits=50, cp=0.00)
modelB$method <- "boulevard"

for (i in 1:10) {
    gc()
    testFlag <- (tenFold==i)
    simData <- decompress(data, 
                          actionCntSet, actionFctSet, predictSet,
                          testFlag)
    
    #### Linear Model
    tmp <- c()
    lmFit <- glm(c.formula(target, predictSet),
                 family = binomial(link = "logit"),
                 data = simData$train)
    tmp <- c(tmp, 
             rLmFit = cep(predict(lmFit, newdata = simData$test, type="link"), simData$yTest))
    
    lmFitSat <- glm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                    family = binomial(link = "logit"),
                    data = simData$train)
    tmp <- c(tmp, 
             rLmFitSat = cep(predict(lmFitSat, newdata = simData$test, type="link"), simData$yTest))
    
    gbmFit <- gbm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                  data = simData$train, 
                  distribution = "adaboost",
                  interaction.depth = modelS$control$maxdepth,
                  bag.fraction = gSubsample,
                  n.trees = modelS$ntree)
    tmp <- c(tmp, 
             rGbmFit = cep(predict(gbmFit, newdata = simData$test, modelS$ntree, type="link"), simData$yTest))
    
    gc()
    modelS$n = sum(!testFlag)
    svcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model=modelS)
    tmp <- c(tmp,
             rSvcmFit = cep(lgd.predict(x=simData$xTest, z=simData$zTest, model=svcmFit),
                            simData$yTest))
    
    gc()
    modelB$n = sum(!testFlag) 
    bvcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model = modelB)
    tmpPred <- lgd.predict(x=simData$xTest, z=simData$zTest, model=bvcmFit)
    tmp <- c(tmp, rBvcmFitL = cep(tmpPred, simData$yTest),
             rBvcmFitR = cep(tmpPred * (1+modelB$lambda) / modelB$lambda,  simData$yTest))
    
    result <- rbind(result, tmp)
    cat("SIMULATION ITERATION ", i, "\n", tmp, "\n")
}

save(result, file=paste("Accuracy/", id, ".RData", sep=""))

#### DATASET: EGIRDSTAB_CLS ####

set.seed(42)
id <- "egridstab_cls"
cat("Working on Dataset: ", id, "...\n")

data <- read.table(paste("Accuracy/", id, ".csv", sep=""), header = TRUE, sep = ",")
data <- data[complete.cases(data), ]
data$target <- as.numeric(data$stabf == "stable")

result <- c()
actionCntSet <- c("tau2", "tau3", "tau4")
actionFctSet <- c()
actionSet <- union(actionCntSet, actionFctSet)
predictSet <- c("p2", "p3", "p4", "g2", "g3", "g4")
target <- "target"
for (key in union(actionCntSet, predictSet)) {
    if (is.factor(data[[key]])) {
        data[[key]] = as.numeric(levels(data[[key]]))[data[[key]]]
    }
}
data <- data[complete.cases(data), ]
n <- dim(data)[1]
tenFold <- sample(1:10, size=n, replace=TRUE)

gSubsample <- 0.5

n <- nrow(data)
tenFold <- sample(1:10, size=n, replace=TRUE)

modelS <- list()
modelS$diff <- lg.diff
modelS$init <- lg.init.beta
modelS$pred <- lg.pred
modelS$loss <- lg.loss
modelS$dummy <- TRUE
modelS$xscale <- TRUE
modelS$ntree <- 100
modelS$lambda <- 0.1
modelS$method <- "ordinary"
modelS$subsample <- gSubsample
modelS$quantile_cut <- 0.05

modelS$p <- length(predictSet) + modelS$dummy
modelS$q <- length(actionSet)
modelS$control <- rpart.control(maxdepth = 4, minsplits=100, cp=0.001)
modelS$woods <- list()

modelB <- modelS
modelB$lambda <- 1
modelB$control <- rpart.control(maxdepth = 6, minsplits=50, cp=0.00)
modelB$method <- "boulevard"

for (i in 1:10) {
    gc()
    testFlag <- (tenFold==i)
    simData <- decompress(data, 
                          actionCntSet, actionFctSet, predictSet,
                          testFlag)
    
    #### Linear Model
    tmp <- c()
    lmFit <- glm(c.formula(target, predictSet),
                 family = binomial(link = "logit"),
                 data = simData$train)
    tmp <- c(tmp, 
             rLmFit = cep(predict(lmFit, newdata = simData$test, type="link"), simData$yTest))
    
    lmFitSat <- glm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                    family = binomial(link = "logit"),
                    data = simData$train)
    tmp <- c(tmp, 
             rLmFitSat = cep(predict(lmFitSat, newdata = simData$test, type="link"), simData$yTest))
    
    gbmFit <- gbm(c.formula(target, predictSet, actionCntSet, actionFctSet),
                  data = simData$train, 
                  distribution = "adaboost",
                  interaction.depth = modelS$control$maxdepth,
                  bag.fraction = gSubsample,
                  n.trees = modelS$ntree)
    tmp <- c(tmp, 
             rGbmFit = cep(predict(gbmFit, newdata = simData$test, modelS$ntree, type="link"), simData$yTest))
    
    gc()
    modelS$n = sum(!testFlag)
    svcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model=modelS)
    tmp <- c(tmp,
             rSvcmFit = cep(lgd.predict(x=simData$xTest, z=simData$zTest, model=svcmFit),
                            simData$yTest))
    
    gc()
    modelB$n = sum(!testFlag) 
    bvcmFit <- lgd(x=simData$xTrain, y=simData$yTrain, z=simData$zTrain, model = modelB)
    tmpPred <- lgd.predict(x=simData$xTest, z=simData$zTest, model=bvcmFit)
    tmp <- c(tmp, rBvcmFitL = cep(tmpPred, simData$yTest),
             rBvcmFitR = cep(tmpPred * (1+modelB$lambda) / modelB$lambda,  simData$yTest))
    
    result <- rbind(result, tmp)
    cat("SIMULATION ITERATION ", i, "\n", tmp, "\n")
}

save(result, file=paste("Accuracy/", id, ".RData", sep=""))

