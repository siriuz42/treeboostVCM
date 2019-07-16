source("lgd.R")

data <- read.csv("BeijingHousing/new.csv",sep=',',stringsAsFactors=FALSE, fileEncoding="latin1")
colnames(data)
data$"floor" <- as.numeric(gsub("[^0-9]", "", data$"floor"))
y <- as.numeric(as.vector(data$"totalPrice")) / as.numeric(as.vector(data$"square"))
x <- cbind(as.numeric(as.vector(data$"floor")),
           as.numeric(as.vector(data$"bathRoom")),
           as.numeric(as.vector(data$"livingRoom")),
           as.numeric(as.vector(data$"elevator")),
           as.numeric(as.numeric(as.vector(data$"buildingType")) == 1),
           as.numeric(as.vector(data$"renovationCondition")) - 2
           )
x[x[, 6] < 0, 6] <- 0

z <- cbind(as.numeric(as.vector(data$"Lng")),
           as.numeric(as.vector(data$"Lat")))

na.list <- which(is.na(y))
for (i in 1:ncol(x)) {
    na.list <- c(na.list, which(is.na(x[, i])))
}
for (i in 1:ncol(z)) {
    na.list <- c(na.list, which(is.na(z[, i])))
}

x <- x[-na.list, ]
z <- z[-na.list, ]
y <- y[-na.list]

xsample <- sample(nrow(x), 30000)
x <- x[xsample, ]
z <- z[xsample, ]
y <- y[xsample]

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

res <- lgd(y=y, x=x, z=z, model=model)
ssample <- sample(res$n, 10000)

summary(res$beta[, 1])
plot_ly(x=z[model$subsample, 1], y=z[model$subsample, 2],
        color=res$beta[model$subsample, 1],
        type="scatter", mode="markers", colors=c("red", "lightgrey", "blue"),
        marker = list(size = 3.5)) %>%
    layout(title = "Intercept",
           xaxis = list(title="Longtitude"),
           yaxis = list(title="Latitude"))

qplot(x=z[ssample, 1], y=z[ssample, 2], geom="point", 
      color = res$beta[ssample, 1],
      main = "Intercept",
      xlab = "Longitude",
      ylab = "Latitude") + 
    scale_color_distiller(palette="RdBu", direction=1) +
    labs(colour = 'Val') +
    geom_point(shape=".")
dev.off()

summary(res$beta[, 1])
plot_ly(x=z[model$subsample, 1], y=z[model$subsample, 2], 
        color=res$beta[model$subsample, 2], 
        type="scatter", mode="markers", colors=c("red", "lightgrey", "blue"),
        marker = list(size = 3.5)) %>%
    layout(title = "Floor",
           xaxis = list(title="Longtitude"),
           yaxis = list(title="Latitude"))

plot_ly(x=z[model$subsample, 1], y=z[model$subsample, 2], 
        color=res$beta[model$subsample, 3], 
        type="scatter", mode="markers", colors=c("red", "lightgrey", "blue"),
        marker = list(size = 3.5)) %>%
    layout(title = "Bathroom",
           xaxis = list(title="Longtitude"),
           yaxis = list(title="Latitude"))

plot_ly(x=z[model$subsample, 1], y=z[model$subsample, 2], 
        color=res$beta[model$subsample, 4], 
        type="scatter", mode="markers", colors=c("red", "lightgrey", "blue"),
        marker = list(size = 3.5)) %>%
    layout(title = "Living Room",
           xaxis = list(title="Longtitude"),
           yaxis = list(title="Latitude"))

plot_ly(x=z[model$subsample, 1], y=z[model$subsample, 2], 
        color=res$beta[model$subsample, 5], 
        type="scatter", mode="markers", colors=c("red", "lightgrey", "blue"),
        marker = list(size = 3.5)) %>%
    layout(title = "Elevator",
           xaxis = list(title="Longtitude"),
           yaxis = list(title="Latitude"))

plot_ly(x=z[model$subsample, 1], y=z[model$subsample, 2], 
        color=res$beta[model$subsample, 6], 
        type="scatter", mode="markers", colors=c("red", "lightgrey", "blue"),
        marker = list(size = 3.5)) %>%
    layout(title = "Tower",
           xaxis = list(title="Longtitude"),
           yaxis = list(title="Latitude"))

plot_ly(x=z[model$subsample, 1], y=z[model$subsample, 2], 
        color=res$beta[model$subsample, 7], 
        type="scatter", mode="markers", colors=c("red", "lightgrey", "blue"),
        marker = list(size = 3.5)) %>%
    layout(title = "Refurbishment",
           xaxis = list(title="Longtitude"),
           yaxis = list(title="Latitude"))

#### GGPLOT ####
png(filename="BeijingHouse_Plot/plot_bj_intercept.png", width=300, height=240)
dat <- data.frame(Longitude = z[ssample, 1], 
                  Latitude = z[ssample, 2],
                  Val = res$beta[ssample, 1])
g <- ggplot(dat, aes(Longitude, Latitude)) + 
    geom_point(aes(color = Val), alpha=1, cex=0.2) + 
    scale_color_distiller(palette="RdBu", direction=1) + 
    ggtitle("Intercept")
print(g)
dev.off()

png(filename="BeijingHouse_Plot/plot_bj_floor.png", width=300, height=240)
dat <- data.frame(Longitude = z[ssample, 1], 
                  Latitude = z[ssample, 2],
                  Val = res$beta[ssample, 2])
g <- ggplot(dat, aes(Longitude, Latitude)) + 
    geom_point(aes(color = Val), alpha=1, cex=0.2) + 
    scale_color_distiller(palette="RdBu", direction=1) + 
    ggtitle("Floor")
print(g)
dev.off()

png(filename="BeijingHouse_Plot/plot_bj_bathroom.png", width=300, height=240)
dat <- data.frame(Longitude = z[ssample, 1], 
                  Latitude = z[ssample, 2],
                  Val = res$beta[ssample, 3])
g <- ggplot(dat, aes(Longitude, Latitude)) + 
    geom_point(aes(color = Val), alpha=1, cex=0.2) + 
    scale_color_distiller(palette="RdBu", direction=1) + 
    ggtitle("Bathroom")
print(g)
dev.off()

png(filename="BeijingHouse_Plot/plot_bj_living.png", width=300, height=240)
dat <- data.frame(Longitude = z[ssample, 1], 
                  Latitude = z[ssample, 2],
                  Val = res$beta[ssample, 4])
g <- ggplot(dat, aes(Longitude, Latitude)) + 
    geom_point(aes(color = Val), alpha=1, cex=0.2) + 
    scale_color_distiller(palette="RdBu", direction=1) + 
    ggtitle("Living Room")
print(g)
dev.off()

png(filename="BeijingHouse_Plot/plot_bj_elevator.png", width=300, height=240)
dat <- data.frame(Longitude = z[ssample, 1], 
                  Latitude = z[ssample, 2],
                  Val = res$beta[ssample, 5])
g <- ggplot(dat, aes(Longitude, Latitude)) + 
    geom_point(aes(color = Val), alpha=1, cex=0.2) + 
    scale_color_distiller(palette="RdBu", direction=1) + 
    ggtitle("Elevator")
print(g)
dev.off()

png(filename="BeijingHouse_Plot/plot_bj_refurbish.png", width=300, height=240)
dat <- data.frame(Longitude = z[ssample, 1], 
                  Latitude = z[ssample, 2],
                  Val = res$beta[ssample, 7])
g <- ggplot(dat, aes(Longitude, Latitude)) + 
    geom_point(aes(color = Val), alpha=1, cex=0.2) + 
    scale_color_distiller(palette="RdBu", direction=1) + 
    ggtitle("Refurbishment")
print(g)
dev.off()
