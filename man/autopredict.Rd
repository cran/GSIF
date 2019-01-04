\name{autopredict-methods}
\alias{autopredict}
\alias{makePixels}
\alias{autopredict,SpatialPointsDataFrame,SpatialPixelsDataFrame-method}
\title{Autopredict numeric or factor type variables}
\description{Fits a spatial prediction model via the \code{\link{fit.gstatModel}} function (in the case of numeric variable), or some classification algorithm (factor-type variable), and generates predictions. By defaults uses machine learning method (random forest) as implemented in the ranger package.}
\usage{
\S4method{autopredict}{SpatialPointsDataFrame,SpatialPixelsDataFrame}(target, covariates, 
     auto.plot=TRUE, spc=TRUE, buffer.dist=TRUE, \dots)
}
\arguments{
  \item{target}{object of class \code{"SpatialPointsDataFrame"} containing observations of the target variable} 
  \item{covariates}{object of class \code{"SpatialPixelsDataFrame"}; spatial covariates}
  \item{auto.plot}{logical; specifies whether to immediately plot the data via the plotKML function}
  \item{spc}{logical; specifies whether to generate spatial predictive components}
  \item{buffer.dist}{logical; specifies whether to use buffer distances as covariates}
  \item{\dots}{other optional arguments that can be passed to \code{fit.gstatModel}}
}
\details{For factor-type variables, classes with <5 observations are automatically removed from analysis. Generation of principal components via the \code{spc} function and buffer distances can be time-consuming for large rasters.}
\author{ Tomislav Hengl}
\seealso{ \code{\link{fit.gstatModel}} }
\examples{
\dontrun{
library(sp)
library(gstat)
library(ranger)
library(plotKML)
library(raster)

## Ebergotzen data set:
data(eberg)
## subset to 8%
eberg <- eberg[runif(nrow(eberg))<.08,]
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")

## predict sand content:
SNDMHT_A <- autopredict(eberg["SNDMHT_A"], eberg_grid,
   auto.plot=FALSE, rvgm=NULL)
plot(raster(SNDMHT_A$predicted["SNDMHT_A"]), col=SAGA_pal[[1]])

## predict soil types:
soiltype <- autopredict(eberg["soiltype"], eberg_grid, 
  auto.plot=FALSE)
spplot(soiltype$predicted, col.regions=R_pal[[2]])
## most probable class:
eberg_grid$soiltype <- as.factor(apply(soiltype$predicted@data, 1, which.max))
levels(eberg_grid$soiltype) = names(soiltype$predicted@data)
spplot(eberg_grid["soiltype"])

## Meuse data set:
demo(meuse, echo=FALSE)
zinc <- autopredict(meuse["zinc"], meuse.grid[c("dist","ffreq")], 
  auto.plot=FALSE, rvgm=NULL)
spplot(zinc$predicted["zinc"])
}
}