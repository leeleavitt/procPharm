#' Function to track the various users data during analysis
#' The goal is to append a line to the log with 
#' 
#' Needs to happen before the function goes
#' time1 <- proc.time()
#' Fill up your additional info...
#' additionalInfo <- choices
#' timeInFunction <- (proc.time() - time1)[3]
#' functionName <- as.character(match.call())[1]
#' additionalInfo <- choices
#' logger(functionName, timeInFunction, additionalInfo)
#' @export
logger <- function(functionName, timeInFunction, additionalInfo = NA, rate = F){
    pathSplit <- strsplit(getwd(), "/")[[1]]
    driveLoc <- pathSplit[1]

    # 1 
    dateTime <- as.character(Sys.time())
    # 2
    user <- pathSplit[2]
    # 3
    experiment <- ls(pattern = "^RD[.]", envir = .GlobalEnv)
    # 4
    experimentLocation <- getwd()
    # 5
    functionName <- functionName
    # 6
    timeInFunction <- timeInFunction

    # 7
    if(rate){
        alarm()
        flush.console()
        cat("\nRate quality of your scoring 1 to 10\n")
        quality <- scan(n=1, quiet=TRUE)
    }else{
        quality <- NA
    }

    # 8 
    additionalInfo <- paste0(additionalInfo, collapse = '_')

    logToAdd <- paste0(c(dateTime, user, experiment, experimentLocation, functionName, timeInFunction, quality, additionalInfo), collapse =',')

    # Save to the specified location
    logLocation <- "/Lee Leavitt/spyLog.txt"
    logLocation <- paste0(driveLoc, logLocation)

    write(logToAdd, logLocation, append = TRUE)
}
