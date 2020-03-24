#Funciton to save the work along with create a unique savehistory

saveRD <- function(dat){
    cat("\nDO NOT CLOSE UNTIL I SAY YOU CAN!\nWait for the sound...")
    bringToTop(-1)
    Sys.sleep(1)

    #History Saver
    experimentorsName <- strsplit(getwd(),'/')[[1]][2]
    historyName <- paste(experimentorsName, Sys.time(), 'History.r')
    historyName <- gsub(":", '_',historyName)
    savehistory(historyName)

    #Exp Saver
    expName <- deparse(substitute(dat))
    expToSave <- get(expName, envir = .GlobalEnv)
    save(expToSave, file=paste0(expName,".Rdata") )
    alarm()
    cat('\nYou can now close. Please consider cleaning up the file,\n',historyName,'\n')
}

