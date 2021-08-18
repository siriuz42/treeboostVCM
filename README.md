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

1. Load `lgd.R`.

2. Create an empty `lgd` model. We take Beijing Housing data as an example:
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
Another set of arguments can be used for logistic regression, including `lg.diff`, `lg.init.beta`, `lg.pred`, and `lg.loss`. 

3. There are some other arguments that can be set:
```r
model$savetree            # [savetree] is whether the ensemble is saved. Set to [FALSE] when no need to make predictions.
model$method              # [method] is the method to boost. Now only [ordinary] works. Please leave empty.
```

4. Train by calling `lgd()` to get a trained lgd model. For instance
```r
res <- lgd(y=y, x=x, z=z, model=model)
```
Here `x` is the predictive covariates, arranged in an n by p matrix. `z` is the action covariates in an n by q matrix, `y` is a n-vector representing the response.

5. Inspect the result. For instance, the structure of this trained `res` is 
```r
res$yhat    # [yhat] is the fitted y values of the training data. 
res$lc      # [lc] is the training loss as a function of the boosting iterations.
res$beta    # [beta] is an n by p matrix showing the varying coefficients for every training point.
res$woods   # [woods] is the boosted trees.
```

## How to make new predictions with Tree Boosted VCM
Use the function `lgd.predict`. There is no such example in Beijing Housing file, but it would be something as the following
```r
newPred <- lgd.predict(x=predX, z=predZ, model=res)
```
Here `x` is the new predictive covariates, `z` is the new action covaraites, and `model` is the trained model. There are a couple of examples in `accuracy.R`.

Current there is no output for the predicted coefficients. Please contact the authors (yichenzhou% at google% dot com% with % removed) and they will patch the code asap.
