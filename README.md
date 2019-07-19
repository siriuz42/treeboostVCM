# Readme: Manual

This repo holds the supplementary code for the arXiv submission *Tree Boosted Varying Coefficient Models*

https://arxiv.org/abs/1904.01058


## List of files in this repo:

+ lgd.R : The main tree boosted VCM (local gradient descent) library that implements all supported functions.

+ honestRpart.R : A dependent lib. 

+ Exp_BeijingHousing.R : Code we used to get the results on Beijing Housing data. A good starting point.

+ exp_simulated.R : Code we used to get synthetic data and evaluate our method.

+ accuracy.R : Code for evaluating the predictive accuracy of the method on several UCI datasets.

## How to train Tree Boosted VCM

1. Load lgb.R
2. Create a lgb model. We take Beijing Housing data as an example:
```r
model <- list()
model$dummy <- TRUE
model$xscale <- TRUE
model$n <- nrow(x)
model$p <- ncol(x) + model$dummy
model$q <- ncol(z)
model$diff <- ols.diff
model$init <- ols.init.beta
model$pred <- ols.pred
model$loss <- ols.loss
model$lambda <- 0.05
model$ntree <- 200
model$subsample <- 0.5
model$control <- rpart.control(maxdepth = 5, cp = 0)
model$woods <- list()
```
