\name{autopredict-methods}
\alias{autopredict}
\alias{autopredict,SpatialPointsDataFrame,SpatialPixelsDataFrame-method}
\title{Auto predict numeric or factor type variables}
\description{Fits either geostatistical model via the \code{\link{fit.gstatModel}} function (in the case of numeric variable) or a multinomial logistic regression model via the \code{\link{spmultinom}} function (factor-type variable) and generates predictions.}
\usage{
\S4method{autopredict}{SpatialPointsDataFrame,SpatialPixelsDataFrame}(target,covariates, 
     auto.plot=TRUE, \dots)
}
\arguments{
  \item{target}{object of class \code{"SpatialPointsDataFrame"} containing observations of the target variable} 
  \item{covariates}{object of class \code{"SpatialPixelsDataFrame"}; spatial covariates}
  \item{auto.plot}{logical; specifies whether to immediately plot the data via the plotKML function}
  \item{\dots}{other optional arguments that can be passed to \code{fit.gstatModel} or \code{spmultinom}}
}
\author{ Tomislav Hengl}
\seealso{ \code{\link{fit.gstatModel}}, \code{\link{spmultinom}} }
\examples{
## Ebergotzen data:
library(sp)
library(gstat)
library(randomForest)
library(plotKML)

## load input data:
data(eberg)
eberg <- eberg[runif(nrow(eberg))<.1,]
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")

## predict sand content:
SNDMHT_A <- autopredict(eberg["SNDMHT_A"], eberg_grid,
   auto.plot=FALSE)
plot(SNDMHT_A)

## predict soil types:
soiltype <- autopredict(eberg["soiltype"], eberg_grid, 
  auto.plot=FALSE)
spplot(soiltype@predicted)
}
