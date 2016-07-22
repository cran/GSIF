\name{fit.regModel-methods}
\alias{fit.regModel}
\alias{fit.regModel,formula,data.frame,SpatialPixelsDataFrame,character-method}
\title{Fits a regression model to spatial data}
\description{Fits a regression or a trend model (e.g. a GLM) and, if not available, a variogram for the response residuals using the default settings.}
\usage{
\S4method{fit.regModel}{formula,data.frame,SpatialPixelsDataFrame,character}(formulaString, rmatrix, predictionDomain, 
      method = list("GLM", "rpart", "randomForest", "quantregForest", "lme",
       "xgboost", "ranger"), 
      dimensions = NULL, fit.family = gaussian(), stepwise = TRUE, rvgm, 
      GLS = FALSE, steps, subsample, subsample.reg, \dots)
}
\arguments{
  \item{formulaString}{object of class \code{"formula"} --- regression model} 
  \item{rmatrix}{object of class \code{"data.frame"}; regression matrix produced as a result of spatial overlay}
  \item{predictionDomain}{object of class \code{"SpatialPixelsDataFrame"}; spatial domain of interest}
  \item{method}{character; family of methods considered e.g. \code{"GLM"}, \code{"rpart"} (regression trees), \code{"randomForest"} (random forest)}
  \item{dimensions}{character; \code{"2D"}, \code{"3D"}, \code{"2D+T"}, or \code{"3D+T"}}
  \item{fit.family}{family to be passed to the \code{glm} (see examples below)}
  \item{stepwise}{specifies whether to run step-wise regression on top of GLM to get an optimal subset of predictors}
  \item{rvgm}{residual variogram (to avoid fitting the variogram set as NULL)}
  \item{GLS}{fit trend model using Generalized Least Squares implemented in the nlme package}
  \item{steps}{integer; the maximum number of steps to be considered for step-wise regression; see \code{stats::step} for more details}
  \item{subsample}{integer; maximum number of observations to be taken for variogram model fitting (to speed up variogram fitting)}
  \item{subsample.reg}{integer; maximum number of observations to be taken for regression model fitting (especially important for randomForest modelling)}
  \item{\dots}{other optional arguments that can be passed to \code{gstat::fit.variogram}}
}
\details{Produces an object of class \code{"gstatModel"} that contains: (1) fitted regression model (e.g. a GLM, cubist model, or randomForest model), (2) fitted variogram, and (c) object of class \code{"SpatialPoints"} with observation locations. To combine overlay and model fitting operations, consider using \code{\link{fit.gstatModel}}.}
\author{ Tomislav Hengl, Mario Antonio Guevara Santamaria and Bas Kempen}
\seealso{ \code{\link{fit.gstatModel}}, \code{stats::glm}, \code{gstat::fit.variogram}, \code{randomForest::randomForest} }
\examples{
library(sp)
library(rpart)
library(nlme)
library(gstat)
library(randomForest)
library(quantregForest)
library(xgboost)
library(caret)
library(scales)
library(AICcmodavg)

## load the Meuse data set:
demo(meuse, echo=FALSE)

## prepare the regression matrix:
ov <- over(meuse, meuse.grid)
ov <- cbind(data.frame(meuse["om"]), ov)
## skip variogram fitting:
m <- fit.regModel(om~dist+ffreq, rmatrix=ov, meuse.grid, 
    fit.family=gaussian(log), method="GLM", rvgm=NULL)
m@regModel
m@vgmModel
plot(m)
## fit a GLM with variogram:
m1 <- fit.regModel(om~dist+ffreq, rmatrix=ov, meuse.grid, 
    fit.family=gaussian(log), method="GLM")
m1@vgmModel
plot(m1)
rk1 <- predict(m1, meuse.grid)
plot(rk1)

## fit a regression tree with variogram:
m2 <- fit.regModel(log1p(om)~dist+ffreq, rmatrix=ov, meuse.grid,
     method="rpart")
plot(m2)
rk2 <- predict(m2, meuse.grid)
plot(rk2)

## fit a lme model with variogram:
m3 <- fit.regModel(log1p(om)~dist, rmatrix=ov, meuse.grid,
     method="lme", random=~1|ffreq)
plot(m3)
#rk3 <- predict(m3, meuse.grid)
#plot(rk3)

## fit a randomForest model with variogram
## NOTE: no transformation required
m4 <- fit.regModel(om~dist+ffreq, rmatrix=ov, meuse.grid,
     method="randomForest")
plot(m4)
rk4 <- predict(m4, meuse.grid)
plot(rk4)

## RF is very sensitive to the 'mtry' argument:
m4b <- fit.regModel(om~dist+ffreq, rmatrix=ov, meuse.grid,
     method="randomForest", mtry=2)
plot(m4b)

\dontrun{
## RF with uncertainty (quantregForest package)
m5 <- fit.regModel(om~dist+ffreq, rmatrix=ov, meuse.grid,
     method="quantregForest")
plot(m5)
rk5 <- predict(m5, meuse.grid)
plot(rk5)
dev.off()
## plot prediction error map:
spplot(rk5@predicted["var1.var"])

## ranger
m6 <- fit.regModel(om~dist+ffreq, rmatrix=ov, meuse.grid,
     method="ranger", rvgm=NULL)
plot(m6)
rk6 <- predict(m6, meuse.grid)
plot(rk6)

## XGBoost
m7 <- fit.regModel(om~dist+ffreq, rmatrix=ov, meuse.grid,
     method="xgboost", rvgm=NULL)
plot(m7)
rk7 <- predict(m7, meuse.grid)
plot(rk7)
}
}
