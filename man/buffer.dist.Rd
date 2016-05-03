\name{buffer.dist-methods}
\alias{buffer.dist}
\alias{buffer.dist,SpatialPointsDataFrame,SpatialPixelsDataFrame-method}
\title{Derive buffer distances to a set of points}
\description{Derive buffer distances using the \code{raster::distance} function, so that these can be used as predictors for spatial prediction i.e. to account for spatial proximity to low, medium and high values.}
\usage{
\S4method{buffer.dist}{SpatialPointsDataFrame,SpatialPixelsDataFrame}(observations, predictionDomain, classes, width, \dots)
}
\arguments{
  \item{observations}{object of class \code{"SpatialPointsDataFrame"} containing observations of the target variable} 
  \item{predictionDomain}{object of class \code{"SpatialPixelsDataFrame"}; prediction domain for which distances are estimated}
  \item{classes}{factor; split of the points}
  \item{width}{numeric; maximum search radius}
  \item{\dots}{other optional arguments that can be passed to \code{raster::distance}}
}
\note{Number of breaks (numeric resolution) should be higher than the number of bins, for example, estimated for the histogram display. Machine learning techniques can be quite sensitive to blunders / artifacts in the input point data. 
}
\author{ Tomislav Hengl}
\seealso{ \code{\link{fit.gstatModel}} }
\examples{
library(sp)
library(raster)
library(gstat)
library(randomForest)
library(quantregForest)
library(plotKML)
library(scales)

## Load the Meuse data set:
demo(meuse, echo=FALSE)
classes <- cut(meuse$om, breaks=seq(0, 17, length=8))
## are these optimal splits?
grid.dist <- buffer.dist(meuse["om"], meuse.grid[1], classes)
plot(stack(grid.dist))
## quantregForest as a 'replacement' for kriging:
dn <- paste(names(grid.dist), collapse="+")
fm <- as.formula(paste("om ~", dn))
m <- fit.gstatModel(meuse, fm, grid.dist, 
    method="quantregForest")
plot(m)
dev.off()
## Residual variogram shows no spatial structure
rk.m <- predict(m, grid.dist, vgmmodel=NULL)
plot(rk.m)
dev.off()
## prediction error:
plot(sqrt(raster(rk.m@predicted[2])))
points(meuse, pch="+")
\dontrun{
plotKML(rk.m@predicted["om"], colour_scale = SAGA_pal[[1]])
kml(meuse, file.name="om_points.kml", colour=om, labels=meuse$om)
kml_View("om_points.kml")
meuse$classes <- classes
plotKML(meuse["classes"])
}

## Example with zinc:
meuse.gridT <- meuse.grid
meuse.gridT@data <- cbind(meuse.grid@data, grid.dist@data)
fm1 <- as.formula(paste("om ~", dn, "+soil+dist+ffreq"))
m1 <- fit.gstatModel(meuse, fm1, meuse.gridT, 
     method="quantregForest", rvgm=NULL)
## no need to fit variogram in this case
plot(m1)
dev.off()
rk.m1 <- predict(m1, meuse.gridT, vgmmodel=NULL)
plot(rk.m1)
varImpPlot(m1@regModel)
dev.off()
\dontrun{
plotKML(rk.m1@predicted["om"], 
   file.name="rk_combined.kml", 
   colour_scale = SAGA_pal[[1]])
}

## Zinc:
classes2 <- cut(meuse$zinc, 
   breaks=seq(min(meuse$zinc), max(meuse$zinc), length=10))
grid.dist2 <- buffer.dist(meuse["zinc"], meuse.grid[1], classes2)
dn2 <- paste(names(grid.dist2), collapse="+")
meuse.gridT2 <- meuse.grid
meuse.gridT2@data <- cbind(meuse.grid@data, grid.dist2@data)
fm2 <- as.formula(paste("zinc ~", dn2, "+soil+dist+ffreq"))
m2 <- fit.gstatModel(meuse, fm2, meuse.gridT2, 
      method="quantregForest")
varImpPlot(m2@regModel)
rk.m2 <- predict(m2, meuse.gridT2, vgmmodel=NULL)
plot(rk.m2)
dev.off()
## prediction error:
plot(raster(rk.m2@predicted[2]))
\dontrun{
plotKML(rk.m2@predicted["zinc"], 
    file.name="rk_combined_zinc.kml", 
    colour_scale = SAGA_pal[[1]])
kml(meuse, colour=zinc, 
    file.name="zinc_points.kml", labels=meuse$zinc)
kml_View("zinc_points.kml")
}
}