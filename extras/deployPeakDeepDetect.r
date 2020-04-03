require(procPharm)
require(reticulate)
pyPharm <- import('python_pharmer')

tmpRD <- get(load("./extras/RD.200309.30.m.m3.p1.Rdata"))

# Find where the AITC is located and go out 120 points

fancyBin<- function(dat){
    pulsesWithNN <- c('^[bB]ob.*', "^AITC.*", "^[cC]aps.*", "^[mM]enth.*", "[kK][.]40.*")
    nnNames <- c('blob','aitc', 'menthol', 'capsaicin', 'k40')

    for( i in 1:length(pulsesWithNN)){
        print(i)
        
        # Make sure the pulse exists
        pulse <- grep(pulsesWithNN[i], dat$w.dat$wr1)
        
        if(length(pulse) > 0){
            minWin <- min( pulse )
            maxWin <- minWin + 119

            # Snag the pulse for all cells
            pulseToScore <- as.data.frame(t(dat$blc[minWin:maxWin,-1]))

            # Now use the python score all the responses of interest
            featureFrame <- pyPharm$featureMaker(pulseToScore, 10)
            featureScores <- pyPharm$modelRunner(featureFrame, nnNames[i])

            # Transfer these scoring to the binary dataframe
            binName <- grep(pulsesWithNN[i], names(dat$bin), value=T)
            dat$bin[binName] <- featureScores
        }

    }
    return(dat)
}

tmpRD <- fancyBin(tmpRD)
