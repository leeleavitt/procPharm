#' Function to perform the basline correction
#' @export
ProcConstPharm <- function(dat,shws=2,phws=20,bl.meth="SNIP"){
    if(class(dat)=="data.frame"){
        (dat1<-dat)
    }else{
        dat1 <- dat$t.dat
    }
    t.names <- names(dat1)[-1]#Time in first column
    dat1.snr <- dat1 #peak calls stored as SNR
    dat1.snr[,t.names] <- 0
    dat1.bc <- dat1.snr #baseline corrected data

    for(i in t.names){
        p1 <- PeakFunc2(dat1, i, shws=shws, phws=phws, Plotit=F, bl.meth=bl.meth)
        dat1.snr[match(MALDIquant::mass(p1$peaks),dat1[,1]),i] <- MALDIquant::snr(p1$peaks)
        dat1.bc[i] <- MALDIquant::intensity(p1$dat)
    }
    return(list(snr=dat1.snr,blc=dat1.bc))

}

#' @export 
modelMaker <- function(dat){
    datum <- read.csv("../python_packages/python_pharmer/python_pharmer/peakDeepDetect/models/AITC.100uM.h5")
    return(datum)
}

#' binary score for all cells for the regions of interest bScore
#' returns the scoring for all cells subject to the above parameters.
#' as well as the sum for the snr scores and the sd for the snr scores.
#' @param blc is the baseline corrected data
#' @param snr is the signal to noise peak data
#' @param snr.lim is the threshold for significance on the peaks
#' @param blc.lim is the intensity above baseline theshold
#' @param levs indicates the regions of interest. (e.g. the response windows for which the cells will be scored)
#' @param wr indicates the response windows. 
#' @param cnames indicates the cells to score (if null all cells will be scored)
#' @export
bScore <- function(blc, snr, snr.lim, blc.lim, levs, wr, cnames=NULL){
    
    notzero <- function(x){
        as.integer(sum(x) > 0)
    }

    if(is.null(cnames)){
        cnames <- names(blc)[-1]
    }
    wr2 <- wr[is.element(wr,levs)]
    b.snr <- snr[is.element(wr,levs),cnames]
    b.blc <- blc[is.element(wr,levs),cnames]
    b.call <- b.blc
    b.call[,] <- 0
    b.call[b.snr > snr.lim & b.blc > blc.lim] <- 1
    b.score <- data.frame(tot=apply(b.snr,2,sum))
    b.score["sd"] <- apply(b.snr,2,sd)
    for(i in levs){
        b.score[i] <- apply(b.call[wr2==i,],2,notzero)
    }
    return(b.score)
}

#'  Binary scoring dependent upon score const pharm talbe values
#'  Best way to determine parameters is to look through trace click before hand
#' @param snr.min = minimun signal to noise value
#' @param max.min= minimun above baseline threshold
#' @param tot.min= area minimun to consider
#' @param wm.min= which max, Where within the window region does the maximun value occur
#' @param wm.max= where to stop looking for the maximun value
#' @export
bscore2 <- function(dat, levs.1=NULL, snr.lim = 5 , blc.lim = 0.05){
    if(is.null(dat$bin)){
        dat$bin <- data.frame(matrix(nrow = dim(dat$c.dat)[1], ncol = 0))
        row.names(dat$bin) <- row.names(dat$c.dat)
    }

    levs <- setdiff(unique(as.character(dat$w.dat$wr1)),"")
    
    if(is.null(levs.1)){
        levs.1<-levs
    }else{
        levs.1<-levs.1
    }

    for(i in 1:length(levs.1)){
        snr.name <- paste(levs.1[i],".snr", sep="")
        max.name <- paste(levs.1[i],".max", sep="")

        logic <-    dat$scp[,snr.name] >= snr.lim &
                    dat$scp[,max.name] >= blc.lim 
            
        dat$bin[logic, levs.1[i] ] <- 1
        dat$bin[!logic, levs.1[i] ] <- 0

    }
    return(dat)
}

#' Function to use the neural networks within the python package
#' @param dat is the RD file input
#'
fancyBin <- function(dat){
    pyPharm <- reticulate::import('python_pharmer')
    print(pyPharm)
    pulsesWithNN <- c('^[bB]ob.*', "^AITC.*", "^[cC]aps.*", "^[mM]enth.*", "[kK][.]40.*")
    nnNames <- c('blob','aitc','capsaicin', 'menthol', 'k40')

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
            #tryCatch({
                featureScores <- pyPharm$modelRunner(featureFrame, nnNames[i])
                # Transfer these scoring to the binary dataframe
                binName <- grep(pulsesWithNN[i], names(dat$bin), value=T)
                dat$bin[binName] <- featureScores
                #}
            #    , error=function(e) print(paste("Could not score", pulsesWithNN[i]))
            #)
        }
    }
    return(dat)
}

#' Function to return probability scores from the models
#' @param dat is the RD file input
#' @export
probMaker <- function(dat){
    pyPharm <- reticulate::import('python_pharmer')
    pulsesWithNN <- c("^AITC.*", "^[cC]aps.*", "^[mM]enth.*", "[kK][.]40.*")
    nnNames <- c('aitc','capsaicin', 'menthol', 'k40')

    for( i in 1:length(pulsesWithNN)){
        print(i)
        
        # Make sure the pulse exists
        pulse <- grep(pulsesWithNN[i], dat$w.dat$wr1)
        
        if(length(pulse) > 0){
            # Grab the model
            model <- pyPharm$modelLoader(nnNames[i])
            
            # Snag the pulse for all cells
            minWin <- min( pulse )
            maxWin <- minWin + 119

            pulseToScore <- as.data.frame(t(dat$blc[minWin:maxWin,-1]))


            # Now use the python score all the responses of interest
            tryCatch({
                featureFrame <- pyPharm$featureMaker(pulseToScore, 10)
                probs <- model$predict(featureFrame)
                colnames(probs) <- c(0,1)
                # Transfer these scoring to the binary dataframe
                pulseName <- grep(pulsesWithNN[i], names(dat$bin), value=T)
                dat[['probs']][[pulseName]] <- probs
                }
                , error=function(e) print(paste("Could not score", pulsesWithNN[i]))
            )
        }
    }
    return(dat)
}

#' calculate a table of cell characteristics globally and 
#' within specific windows
#' these specifics should include
#' mean and sd, sum of in window peaks, sum of out of window peaks
#' 1)  some measure of dead cell
#' 2)  yes/no peak response for each window
#' 3) peak height
#' 4) max peak SNR
#' 5) peak timing in window
#' 6)
#' variance of smoothed - raw in window
#' define and number blank windows.
#' @export
ScoreConstPharm <- function(dat, blc=NULL, snr=NULL, der=NULL, snr.lim=5, blc.lim=.05, shws=2, t.type = NA){
    t.dat<-dat$t.dat
    if(is.null(blc)){
        if(!is.na(t.type)){
            blc <- dat[[ t.type ]]
        }else{
            stop("There is no blc input or a t.type specified\nTry:\nblc = dat$blc\nor t.type = 'blc'")
        }
        blc <- dat$blc
    }else{
        blc <- blc
    }
    if(is.null(snr)){
        snr<-dat$snr
    }else{
        snr<-snr
        }

    wr<-dat$w.dat$wr1

    gtfunc <- function(x,alph){sum(x > alph,na.rm=T)}
    
    lt5func <- function(x,y){
        ltfunc <- function(i){summary(lm(y[i:(i+5)] ~ x[i:(i+5)]))$coefficients[2,3]}
        iseq <- 1:(length(x)-5)
        res <- sapply(iseq,ltfunc)
        return(range(res))
    }

    levs <- setdiff(unique(wr),"")
    cnames <- names(t.dat)[-1]
    res.tab <- data.frame(mean=apply(blc[,cnames],2,mean))
    res.tab["sd"] <- apply(blc[,cnames],2,sd)
    res.tab["snr.iws"] <- apply(snr[is.element(wr,levs),cnames],2,sum)
    res.tab["snr.ows"] <- apply(snr[!is.element(wr,levs),cnames],2,sum)
    res.tab["snr.iwc"] <- apply(snr[is.element(wr,levs),cnames],2,gtfunc,alph=snr.lim)
    res.tab["snr.owc"] <- apply(snr[!is.element(wr,levs),cnames],2,gtfunc,alph=snr.lim)
    
    tryCatch({
        for(i in cnames){

            s1 <- MALDIquant::createMassSpectrum(t.dat[,"Time"],t.dat[,i])
            s3 <- MALDIquant::smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=shws)
            bl.th <- MALDIquant::estimateBaseline(s3, method="TopHat")[,"intensity"]
            bl.snp <- MALDIquant::estimateBaseline(s3, method="SNIP")[,"intensity"]
            eseq <- 1:ceiling((nrow(t.dat)/2))
            lseq <- max(eseq):nrow(t.dat)
            res.tab[i,"bl.diff"] <- mean(bl.th-bl.snp)
            res.tab[i,"earl.bl.diff"] <- mean(bl.th[eseq]-bl.snp[eseq])
            res.tab[i,"late.bl.diff"] <- mean(bl.th[lseq]-bl.snp[lseq])        
        }
    }, error = function(e)cat("\nCannot Make Stats from MALDIquant\n"))
    
    for(i in levs){
        res.tab[paste(i,".snr",sep="")] <- apply(snr[wr==i,cnames],2,max)
        res.tab[paste(i,".tot",sep="")] <- apply(blc[wr==i,cnames],2,sum)
        res.tab[paste(i,".max",sep="")] <- apply(blc[wr==i,cnames],2,max)
        res.tab[paste(i,".wm",sep="")] <- apply(blc[wr==i,cnames],2,which.max)
    }
    return(res.tab)
}

### These functions create automatic stats but only if the window regions are completely fill out and corred 
# Indirect effect stat maker
#' Function taht creates a ide stat automatically if you have all control and test window regions in the data.
#' To make this work you can't half ass your window regions you must include all control window regions
#' @param statType "max", "snr", "tot" are entries into this
#' @param testPulseNames this is the test pulse names that contains the experiments, examples include "^[kK].*40" or "^[aA][cC][hH].*300[uU][mM]" or "ATP"
#' @param controlToUse this is an integer input to specify the testPulse to use for control is lefft NA, the pulse prior to the pulse is control pulse
#' @export
ideStatMaker <- function(dat, statType = "max", testPulseExp = "^[kK][.]20", controlToUse = c(1,2,3)){
    # statType <- ".max"
    # testPulseNames <- "^[kK][.]20"
    # controlToView <- c(1,2)

    levs <- setdiff(unique(as.character(dat$w.dat$wr1)), "")

    pulses <- grep(testPulseExp, 
        levs, 
        value = F, 
        ignore.case = T
    )
    
    pulsesNames <- grep(testPulseExp, 
        levs, 
        value = T, 
        ignore.case = T
    )

    # Within this region we will compute a specific minmax norm stat
    testPulses <- pulsesNames[2:length(pulsesNames)]
    testPulsesScp <- paste0(testPulses, ".", statType)

    
    if(!is.na(controlToUse)){
        controlPulses <- pulsesNames[controlToUse]
        controlPulsesScp <- paste0(controlPulses, ".", statType)
    }else{
        controlPulses <- pulsesNames
        controlPulsesScp <- paste0(controlPulses, ".", statType)
    }

    for(i in 1:length(testPulses) ){
        # First answer is there a response in the test or control. 
        
        if(!is.na(controlToUse)){
            test <- dat$scp[testPulsesScp[i]]
            control <- apply(dat$scp[controlPulsesScp], 1, mean)
            
            respLogic <- 
            !(
                # Test pulses
                dat$bin[ testPulses[i] ] == 1 |
                apply(dat$bin[ controlPulses ] == 1, 1, any)
            )
        }else{
            test <- dat$scp[ testPulsesScp[i] ]
            control <- dat$scp[ controlPulsesScp[i] ]
            
            respLogic <- 
            !(
                # Test pulses
                dat$bin[ testPulses[i] ] == 1 |
                dat$bin[ controlPulses[i] ] == 1
            )
        }

        # Now compute the F2 Stat
        stat <- (test - control)/(test + control)

        newName <- paste0(levs[pulses[i] + 1], "_", statType, "_", "ide", "_mmnorm")
        dat$scp[newName] <- stat
        dat$scp[respLogic, newName] <- NA
    }
    return(dat)
}

#' Function to create a direct effect stat. 
#' this stat is automatically tested against the control to use.
#' Due to nature of these experiments there are ususally only one
#' window region really compare to, to observe the direct effect.
#' a major weakness of this statistic is taht direct effects cannot be automatically scored
#' if they are automatically score this would be significantly better.
#' since we are unable to compute this we compare window regions of anything. This confuses the outpu
#' of this analysis
deStatMaker <- function(dat, statType = "tot", testPulseNames = "^[kK][.]30", controlToUse = 1){
    # Within this region we will compute a specific minmax norm stat
    testPulses <- grep(testPulseNames, names(dat$bin), value = F, ignore.case = T)
    windows <- names(dat$bin)[testPulses + 1]
    comparisons <- windows[-1]
    controlWindows <- grep('control', windows, value = T)[controlToUse]

    for(i in 2:length(comparisons)){
        # If the cell doesn't respond to the potassium pulse or the potassium pulse following
        # the comparison, then dis regard it
        # respLogic <- !(dat$bin[testPulses[i]] == 1 |
        #     dat$bin[testPulses[0]] == 1)

        statsToComp <- dat$scp[paste0(c(controlWindows, windows[i]),'.', statType)]
        
        #statsToComp[respLogic,] <- NA

        collumnstat <- (statsToComp[2] - statsToComp[1] ) / (statsToComp[2] + statsToComp[1])
        newName <- paste0(comparisons[i], ".",statType,".","de", ".mmnorm")

        dat$scp[newName] <- collumnstat
    }
    return(dat)
}

#' @export
drStatMaker <- function(dat, statType = "tot", testPulseNames=paste0("^[kK].", c(10,15,20,30,40)) ){
    
    testPulses <- paste0(testPulseNames, ".*", statType)
    
    # Conduct the peak height stat
    for(i in 1:length(testPulseNames) ){
        testNames <- rev(rev(grep(testPulseNames[i], names(dat$bin), value = T, ignore.case = T))[1:2])
        testNamesLoc <- rev(rev(grep(testPulseNames[i], names(dat$bin), value = F, ignore.case = T))[1:2])

        respLogic <- 
            !(
                # Test pulses
                dat$bin[ testNames[1] ] == 1 |
                dat$bin[ testNames[2] ] == 1
            )

        # Find the statistic
        statNamesToGrab <- rev(rev(grep(testPulses[i], names(dat$scp), value = T, ignore.case = T))[1:2])
        statsToComp <- dat$scp[statNamesToGrab]
        statsToComp[respLogic,] <- NA

        collumnstat <- (statsToComp[2] - statsToComp[1] ) / (statsToComp[2] + statsToComp[1])[1]
        newName <- paste0(testNames[1], "_", names(dat$bin)[ testNamesLoc[2] - 1 ], "_", statType, "_", "dr", "_mmnorm")
        dat$scp[newName] <- collumnstat
    }

    return(dat)
}

