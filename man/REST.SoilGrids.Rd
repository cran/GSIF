\name{REST.SoilGrids-class}
\docType{class}
\alias{REST.SoilGrids-class}
\alias{REST.SoilGrids}
\alias{over,REST.SoilGrids,SpatialPoints-method}
\title{A class for SoilGrids REST API}
\description{A class for \href{https://rest.soilgrids.org/}{SoilGrids REST API} Service. Can be used to overlay points or fetch grid values from SoilGrids Soil Information System.}
\section{Slots}{
  \describe{
    \item{\code{server}:}{object of class \code{"character"}; contains the location of the server that executes REST.SoilGrids calls}
    \item{\code{query}:}{object of class \code{"list"}; contains parameters or REST.SoilGrids query}
    \item{\code{stream}:}{object of class \code{"character"}; contains parameters or REST.SoilGrids stream operation}
  }
}
\section{Methods}{
  \describe{
  \item{over}{\code{signature(x = "REST.SoilGrids", y = "SpatialPoints")}: overlays spatial points and the target grids defined via the REST.SoilGrids-class (point-by-point) and returns list of objects of \code{"SpatialPixelsDataFrame"}-class}
  }
}
\note{More examples of overlay and download functions are available via \url{http://rest.soilgrids.org/}. \code{over} method is not recommended for large point data sets.}
\references{
\itemize{
\item SoilGrids --- a system for automated soil mapping (\url{https://soilgrids.org})
\item REST API SoilGrids (\url{http://rest.soilgrids.org})
}
}
\author{ Tomislav Hengl & Jorge S. Mendes de Jesus}
\seealso{ \code{\link{SoilGrids-class}}, \code{\link{WPS-class}}
}
\examples{
\dontrun{
library(rjson)
library(sp)
## 2 points:
pnts <- data.frame(lon=c(10.65,5.36), lat=c(51.81,51.48), id=c("p1","p2"))
coordinates(pnts) <- ~lon+lat
proj4string(pnts) <- CRS("+proj=longlat +datum=WGS84")
pnts
## REST example:
soilgrids.r <- REST.SoilGrids(c("ORCDRC","PHIHOX"))
ov <- over(soilgrids.r, pnts)
str(ov)
}
}
\keyword{classes}