# This will run on starup, detect your device,
# and then change different function default values.
.onLoad <- function(libname, pkgname){
    systemType <- Sys.info()[1]
    if( systemType != "Windows" ){
        formals(tcd)$info <<- F
        formals(tcd)$bcex <<- 0.5
        
        windows <<- cairoDevice::Cairo
        formals(windows)$pointsize <<- 7
        
        formals(PeakFunc7)$bcex <<- 1.5
        formals(RDView)$wh <<- 11
        formals(RDView)$hh <<- 6
    }
}