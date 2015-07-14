\name{sample.grid}
\alias{sample.grid}
\alias{sample.grid.SpatialPoints}
\alias{sample.grid.SpatialPointsDataFrame}
\alias{sample.grid,SpatialPoints-method}
\alias{sample.grid,SpatialPointsDataFrame-method}
\title{sample spatial points by grids}
\description{Get a subset of a object of class \code{"SpatialPoints"} or \code{"SpatialPointsDataFrame"} avoiding spatial clustering.
}
\usage{
\S4method{sample.grid}{SpatialPoints}(obj, cell.size, n, bbox, \dots)
\S4method{sample.grid}{SpatialPointsDataFrame}(obj, cell.size, n, bbox, \dots)
}
\arguments{
  \item{obj}{\code{"SpatialPoints*"} object}
  \item{cell.size}{numeric; the cell size of the overlayed \code{"SpatialGridDataFrame"} in the form of \code{c(x,y)}}
  \item{n}{integer; specifies maximum number points in each grid}
  \item{bbox}{matrix; the bounding box of output \code{"SpatialPoints"} or \code{"SpatialPointsDataFrame"}; it is set the same as the \code{obj} if missing}  
  \item{\dots}{other optional arguments that can be passed to \code{over}}
}
\note{Spatial points are overlayed with spatial grids with a specified cell size and then get a subset from each grid with a specified number at most. If one grid has less points than the specified number, all the points are taken. If one grid has more points than the specified number, only this number of points are taken by \code{\link{sample}}. This function can be used when there are too much point observations to be handled, especially for spatially clustered observations. The total number of sampled points are determined by \code{cell.size} and \code{n} together. You will get fewer the sampled points when \code{cell.size} is larger, or/and when \code{n} is smaller. Similar sample sizes can be achieved by differen combination of \code{cell.size} and \code{n}.}
\value{
Returns a list of two objects: (1) an object of type \code{"SpatialPoints"} or \code{"SpatialPointsDataFrame"} that contains a subset of the {obj}, and (2) resulting grid.
}
\author{Wei Shangguan}
\examples{
library(sp)
data(isis)
profs <- isis[["sites"]]
coordinates(profs) <- ~  LONWGS84 + LATWGS84
proj4string(profs) <- CRS("+proj=longlat +datum=WGS84")
## sample SpatialPointsDataFrame: 
#bbox <- matrix(c(-180, -90, 180, 90), nrow=2)
prof1 <- sample.grid(profs, cell.size = c(5,5), n = 1)
l0 <- list("sp.points", profs, pch=1, col="red")
l1 <- list("sp.points", prof1$subset, pch="+", col="black", cex=1.2)
spplot(prof1$grid, scales=list(draw=TRUE), 
   col.regions="grey", sp.layout=list(l0, l1))
## Subsampling ratio:
round(length(prof1$subset)/length(profs)*100, 1)
}
\keyword{methods}