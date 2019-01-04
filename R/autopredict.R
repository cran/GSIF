setMethod("autopredict", signature(target = "SpatialPointsDataFrame", covariates = "SpatialPixelsDataFrame"), function(target, covariates, auto.plot=TRUE, spc=TRUE, buffer.dist=TRUE, ...){
  
  ## parent call:
  parent_call <- as.list(substitute(list(...)))[-1]
  ## TH: TO-DO estimate processing time
  
  if(buffer.dist==TRUE){
    message("Generating buffer distances...")
    if(is.factor(target@data[,1])){
      ## Drop levels with <5 observations
      xg = summary(target@data[,1], maxsum=length(levels(target@data[,1])))
      selg.levs = attr(xg, "names")[xg > 5]
      if(length(selg.levs)<length(levels(target@data[,1]))){
        flev <- as.factor(target@data[,1])
        flev[which(!target@data[,1] %in% selg.levs)] <- NA
        target@data[,1] <- droplevels(flev)
      }
      classes <- target@data[,1]
    } else {
      classes <- cut(target@data[,1], breaks=seq(min(target@data[,1], na.rm=TRUE), max(target@data[,1], na.rm=TRUE), length=2*round(log1p(length(target)))))
    }
    grid.dist <- buffer.dist(target, covariates[1], classes)
    covariates@data <- cbind(covariates@data, grid.dist@data)
  }
  
  ## (optional) predictive components:
  if(spc==TRUE){
    spc.fm <- as.formula(paste("~", paste(names(covariates), collapse = "+")))
    covariates <- spc(covariates, spc.fm)@predicted
  }
  
  ## Model:
  fm <- as.formula(paste(names(target)[1], '~', paste(names(covariates), collapse="+")))
  if(is.factor(target@data[,1])){
    ov <- over(x=target, y=covariates)
    ov <- cbind(as.data.frame(target), ov)
    ov <- ov[complete.cases(ov[,all.vars(fm)]),]
    message("Fitting a random forest model using 'ranger'...")
    m <- ranger::ranger(fm, ov, importance="impurity", write.forest=TRUE, probability=TRUE)
    p <- covariates[1]
    message("Generating predictions...")
    p@data <- data.frame(round(predict(m, covariates@data, probability=TRUE, na.action = na.pass)$predictions*100))
    if(auto.plot==TRUE){
      kml_open(file.name=paste0(names(target)[1], "_predicted.kml"))
      for(j in names(p)){
        g = p[paste(j)]
        names(g) = "colour"
        kml_layer(g, subfolder.name=paste(j), png.width=g@grid@cells.dim[1]*3, png.height=g@grid@cells.dim[1]*3, colour="colour", z.lim=c(0,50), raster_name=paste0(names(target)[1], "_predicted_", j, ".png"), colour_scale=R_pal[["heat_colors"]])
      }
      kml_layer(target, subfolder.name="observed", points_names=paste(target@data[,1]))
      kml_close(paste0(names(target)[1], "_predicted.kml"))
      kml_View(file.name=paste0(names(target)[1], "_predicted.kml"))
    }
    p <- list(m, p, covariates)
    names(p) = c("model","predicted","covariates")
    return(p)
  }
  
  if(is.numeric(target@data[,1])){
    if(!any(names(parent_call) %in% "method")){
      m <- fit.gstatModel(target, fm, covariates, method="ranger", ...)
    } else {
      m <- fit.gstatModel(target, fm, covariates, ...) 
    }
    ## predict:
    p <- predict(m, covariates)
    if(auto.plot==TRUE){ 
      plotKML(p, folder.name=names(target)[1], png.width=p@predicted@grid@cells.dim[1]*3, png.height=p@predicted@grid@cells.dim[1]*3, colour_scale=SAGA_pal[[1]], file.name=paste0(names(target)[1], "_predicted.kml"))
    }
    p <- list(m, p@predicted, covariates)
    names(p) = c("model","predicted","covariates")
    return(p)
  }
})

makePixels <- function(x, y, factors, pixel.size=1e2, t.dens=.2, min.dim=50, max.dim=2000, gdalwarp, sigma=1e3, show.progress=TRUE, area.poly=NULL, remove.files=TRUE, ncores=1){
  
  if(requireNamespace("spatstat", quietly = TRUE)&requireNamespace("maptools", quietly = TRUE)&requireNamespace("snowfall", quietly = TRUE)){
    if(!is.null(x)){
      if(!class(x)=="SpatialPoints"){ stop("object of class 'SpatialPoints' expected") }
      if(is.na(proj4string(x))){ stop("Object 'x' missing SRS / proj4string") }
    }
    if(missing(gdalwarp)){
      gdalwarp <- .programPath(utility="gdalwarp")
    }
    if(any(!file.exists(y))){ stop("Missing files in 'y'") }
    if(missing(factors)){ factors = rep(FALSE, length(y)) }
    
    if(class(area.poly)=="SpatialPixelsDataFrame"){
      if(!is.null(x)){
        if(!proj4string(area.poly)==proj4string(x)){ stop("proj4string for object 'area.poly' not consistent with object 'x'") }
      }
      dens <- area.poly[1]
    } else {
      if(!is.null(x)){
        ## determine spatial domain:
        ncol = round(abs(diff(x@bbox[1,]))/pixel.size); nrow = round(abs(diff(x@bbox[2,]))/pixel.size)
        if(nrow>max.dim|ncol>max.dim){ stop("Grid too large. Consider increasing 'pixel.size'") }
        if(nrow<min.dim|ncol<min.dim){ stop("Grid too small. Consider reducing 'pixel.size'") }
        grid.owin <- spatstat::owin(x@bbox[1,], x@bbox[2,], mask=matrix(TRUE,nrow,ncol))
        x.ppp <- spatstat::ppp(x=x@coords[,1], y=x@coords[,2], window=grid.owin)
        dens <- spatstat::density.ppp(x.ppp, sigma=sigma)
        dens <- maptools::as.SpatialGridDataFrame.im(dens)
        if(is.null(area.poly)){
          dx <- quantile(dens@data[,1], t.dens, na.rm=TRUE)
          dens@data[,1] <- ifelse(dens@data[,1]<dx, NA, dens@data[,1])
        } else {
          if(!class(area.poly)=="SpatialPolygons"){ stop("object of class 'SpatialPolygons' expected") }
          area.poly <- as(rasterize(area.poly, raster(dens), field=1), "SpatialGridDataFrame")
          dens@data[,1] <- ifelse(is.na(area.poly@data[,1]), NA, dens@data[,1])
        }
        ## generate mask map:
        dens <- as(dens[1], "SpatialPixelsDataFrame")
        if(length(dens)>1e6){ warning("Grid contains more than 1e6 elements. Consider reducing 'pixel.size'.") }
        proj4string(dens) = proj4string(x)
        ## copy grid parameters
        attr(dens@bbox, "dimnames")[[1]] = attr(x@bbox, "dimnames")[[1]]
        attr(dens@coords, "dimnames")[[2]] = attr(x@coords, "dimnames")[[2]]
      } else {
        stop("Either 'SpatialPoints' object or mask map required")
      }
    }
    ## resample all layers to target grid:
    if(ncores==1){
      message(paste("Resampling values from", length(y), "rasters..."))
      out.fname <- .resampleGDAL(y, factors, dens, pixel.size, show.progress, gdalwarp)
    } else {
      snowfall::sfInit(parallel=TRUE, cpus=ncores)
      snowfall::sfExport("y", "factors", "dens", "pixel.size", "gdalwarp", ".resampleGDAL")
      snowfall::sfLibrary("rgdal", character.only=TRUE)
      snowfall::sfLibrary("RSAGA", character.only=TRUE)
      out.fname <- unlist(snowfall::sfClusterApplyLB(1:length(y), fun=function(i){ try( .resampleGDAL(y[i], factors[i], dens, pixel.size, show.progress=FALSE, gdalwarp) )}))
      snowfall::sfStop()
    }
    
    ## Import the resampled maps
    ## round up numbers for density map
    dens@data[,1] <- (dens@data[,1]/max(dens@data[,1], na.rm=TRUE))*1e2
    for(i in 1:length(y)){
      dens@data[,i+1] <- rgdal::readGDAL(out.fname[[i]])$band1[dens@grid.index]
    }
    names(dens) = c("dens", unlist(out.fname))
    ## Assign factors where necessary
    for(i in 1:length(y)){
      if(factors[i]){ dens@data[,i+1] = as.factor(dens@data[,i+1]) }
    }
    ## Remove all layers without any variation:
    rm.v <- !sapply(dens@data, function(x){var(as.numeric(x), na.rm=TRUE)})>0
    if(any(rm.v)){
      dens <- dens[!rm.v]
    }
    if(remove.files==TRUE){
      unlink(unlist(out.fname))
    }
    return(dens)
  } else {
    stop("Missing packages 'maptools' and/or 'spatstat'")
  }
}

.resampleGDAL <- function(y, factors, dens, pixel.size, show.progress, gdalwarp){
  if (show.progress) { pb <- txtProgressBar(min=0, max=length(y), style=3) }
  out.fname <- list(NULL)
  for(i in 1:length(y)){
    if(.Platform$OS.type == "windows"){
      if(!file.exists(gdalwarp)){ stop("'gdalwarp' program not found") }
      fname <- normalizePath(y[i], winslash="/")
    } else {
      fname <- normalizePath(y[i])
    }
    ## decide about whether to aggregate or downscale
    gi <- GDALinfo(fname)
    if(gi[["res.x"]]<pixel.size & !factors[i]){ 
      resample.method = "average"
    } else {
      if(gi[["res.x"]]>pixel.size & !factors[i]){ 
        resample.method = "cubicspline"
      } else {
        resample.method = "near"
      }
    } 
    out.fname[[i]] = set.file.extension(paste0("r_", basename(fname)), ".tif")
    if(!file.exists(out.fname[[i]])){ system(paste0(gdalwarp, ' ', fname, ' ', out.fname[[i]], ' -t_srs \"', proj4string(dens), '\" -co \"COMPRESS=DEFLATE\" -r \"', resample.method, '\" -tr ', pixel.size, ' ', pixel.size,' -te ', paste(as.vector(dens@bbox), collapse=" "))) }
    if (show.progress) { setTxtProgressBar(pb, i) }
  }
  if (show.progress) { close(pb) }
  return(out.fname)
}
