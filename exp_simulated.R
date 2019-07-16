#### simulated data
source("lgd.R")
library(MASS)
library(ggplot2)
#### Simulate f as actual linear span

#### SIM1: DATA ####
n <- 10000
data <- list()
data$z1 <- runif(n)
data$z2 <- runif(n)
data$c1 <- sample(1:10, n, replace=TRUE)
data$c2 <- sample(1:10, n, replace=TRUE)

sigma <- matrix(0, nrow=7, ncol=7)
for (i in 1:7) {
    for (j in 1:7) {
        if (i==j) {
            sigma[i,j] = 1
        } 
        # else if (abs(i-j)==1) {
        #     sigma[i,j] = 0.4
        # }
    }
}

data$x <- matrix(0, nrow=n, ncol=7)
y <- c()
tbeta <- c()
a <- rep(0, 10)
for (i in 1:n) {
    data$x[i, ] <- mvrnorm(1, rep(0, 7), Sigma = sigma)
    if (data$c1[i] < 4) {
        a[1] <- a[1] + 1
        y <- c(y, 1 + 3 * data$x[i, 1] + 7 * data$x[i, 2])
        tbeta <- rbind(tbeta, c(1, 3, 7, 0, 0, 0, 0, 0))
    } else if (data$c1[i] > 8) {
        a[2] <- a[2] + 1
        y <- c(y, -5 + 2 * data$x[i, 1] + 4 * data$x[i, 2] + 6 * data$x[i, 3])
        tbeta <- rbind(tbeta, c(-5, 2, 4, 6, 0, 0, 0, 0))
    } else if (data$c2[i] == 1 | data$c2[i] == 3 | data$c2[i] == 5) {
        a[3] <- a[3] + 1
        y <- c(y, 5 + 5 * data$x[i, 2] + 5 * data$x[i, 3]) 
        tbeta <- rbind(tbeta, c(5, 0, 5, 5, 0, 0, 0, 0))
    } else if (data$z1[i] < 0.5) {
        a[4] <- a[4] + 1
        y <- c(y, 10 + 10 * data$x[i, 4])
        tbeta <- rbind(tbeta, c(10, 0, 0, 0, 10, 0, 0, 0))
    } else if (data$z2[i] < 0.4) {
        a[5] <- a[5] + 1
        y <- c(y, 10 + 10 * data$x[i, 5])
        tbeta <- rbind(tbeta, c(10, 0, 0, 0, 0, 10, 0, 0))
    } else if (data$z1[i] < data$z2[i]) {
        a[6] <- a[6] + 1
        y <- c(y, 5 - 5 * data$x[i, 2] - 10*data$x[i, 3])
        tbeta <- rbind(tbeta, c(5, 0, -5, -10, 0, 0, 0, 0))
    } else {
        a[7] <- a[7] + 1
        y <- c(y, -10 * data$x[i, 1] + 10*data$x[i, 3])
        tbeta <- rbind(tbeta, c(0, -10, 0, 10, 0, 0, 0, 0))
    }
}
y <- y + rnorm(n, 0, 0.5)

#### SIM1: SIM ####
x <- data$x
z <- cbind(data$c1, data$c2, data$z1, data$z2)
model <- list()
model$dummy <- TRUE
model$xscale <- FALSE
model$n <- nrow(x)
model$p <- ncol(x) + model$dummy
model$q <- ncol(z)
model$diff <- ols.diff
model$init <- ols.init.beta
model$pred <- ols.pred
model$loss <- ols.loss
model$lambda <- 0.2
model$subsample <- 0.5
model$ntree <- 100
model$control <- rpart.control(maxdepth = 6, cp=0.0001)
model$woods <- list()

res1 <- lgd(y=y, x=x, z=z, model=model)

model$method <- "boulevard"
model$lambda <- 1
model$ntree <- 100
model$control <- rpart.control(maxdepth = 6, cp=0.0001)
model$inflate <- 1

res2 <- lgd(y=y, x=x, z=z, model=model)

for (i in 0:7) {
    png(paste("nbeta_", i, ".png", sep=""), width=240, height=200)
    dat <- rbind(data.frame(x = tbeta[, i+1], id = "gt"),
                  data.frame(x = res1$beta[, i+1], id = "vcm1"),
                  data.frame(x = res2$beta[, i+1] * res2$inflate, id = "vcm2"))
    print(
        ggplot(dat, aes(x = x)) + 
        geom_histogram(data=subset(dat, id=="gt"), fill="black", alpha=0.3, bins=50) +
        geom_histogram(data=subset(dat, id=="vcm1"), fill="red", alpha=0.5, bins=50) +
        # geom_histogram(data=subset(dat, id=="vcm2"), fill="blue", alpha=0.5, bins=50) +
        labs(x="Coefficient", y="Count", title=bquote(beta[.(i)]))
    )
    dev.off()
}


