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
model <- list()           # Create an empty model.
model$dummy <- TRUE       # [TRUE] if we want an intercept term in the linear part. [FALSE] otherwise. 
model$xscale <- TRUE      # [TRUE] to scale the covariates for better boosting result. Always set to [TRUE].
model$n <- nrow(x)        # [n] is the number of training points.
model$p <- ncol(x) + model$dummy    # [p] is the number of all predictive covariates.
model$q <- ncol(z)        # [q] is the number of action covariates (effect modifiers).
model$diff <- ols.diff    # [diff] is the derivative of the loss function. [ols.diff] for ordinary least square (regression).
model$init <- ols.init.beta         # [init] is the initial values of the parameters. [ols.init.beta] for regression.
model$pred <- ols.pred    # [pred] is the prediction function. [ols.pred] for regression. 
model$loss <- ols.loss    # [loss] is the loss function. [ols.loss] for regression.
model$lambda <- 0.05      # [lambda] is the learning rate for boosting.
model$ntree <- 200        # [ntree] is the number of trees used for boosting.
model$subsample <- 0.5    # [subsample] is the subsampling rate for boosting.
model$control <- rpart.control(maxdepth = 5, cp = 0)  # [control] is an [rpart.control] object that tells how to trees.
model$woods <- list()     # [woods] is initialized as emplty for remembering the boosting history. 
```
3. 
