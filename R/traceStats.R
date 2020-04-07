#' Function to perform the basline correction
#' @export
ProcConstPharm <- function(dat,shws=2,phws=20,bl.meth="SNIP"){
    if(class(dat)=="data.frame"){(dat1<-dat)}else{dat1 <- dat$t.dat}
    t.names <- names(dat1)[-1]#Time in first column
    dat1.snr <- dat1 #peak calls stored as SNR
    dat1.snr[,t.names] <- 0
    dat1.bc <- dat1.snr #baseline corrected data

    for(i in t.names)
    {
        p1 <- PeakFunc2(dat1,i,shws=shws,phws=phws,Plotit=F,bl.meth=bl.meth)
        dat1.snr[match(mass(p1$peaks),dat1[,1]),i] <- snr(p1$peaks)
        dat1.bc[i] <- intensity(p1$dat)
    }
    dat1.der<-dat1.bc[-1,]-dat1.bc[-nrow(dat1.bc),]
    dat1.der <- sweep(dat1.der[,-1],1,dat1.der[,1],'/')

    #    dat1.crr <- allCRR(dat1,t.names,Plotit=F) #leave off advanced processing for now
    return(list(snr=dat1.snr,blc=dat1.bc, der=dat1.der))
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
    
    notzero <- function(x){as.integer(sum(x) > 0)}
    if(is.null(cnames)){cnames <- names(blc)[-1]}
    wr2 <- wr[is.element(wr,levs)]
    b.snr <- snr[is.element(wr,levs),cnames]
    b.blc <- blc[is.element(wr,levs),cnames]
    b.call <- b.blc
    b.call[,] <- 0
    b.call[b.snr > snr.lim & b.blc > blc.lim] <- 1
    b.score <- data.frame(tot=apply(b.snr,2,sum))
    b.score["sd"] <- apply(b.snr,2,sd)
    for(i in levs)
    {
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
bscore2<-function(dat, levs.1=NULL, snr.min=2.8, max.min=.03, wm.min=0, wm.max=600){
    scp<-dat$scp
    levs<-setdiff(unique(as.character(dat$w.dat[,2])),"")
    if(is.null(levs.1)){levs.1<-levs}
    else{levs.1<-levs.1}
    #dat2<-matrix(0, nrow=length(dat$c.dat[,1]), ncol=length(levs))
    dat2<-dat$bin[levs]
    #row.names(dat2)<-dat$c.dat[,1]
    #colnames(dat2)<-levs
    x.names<-dat$c.dat[,1]
    for(j in x.names){	
        for(i in levs.1){
            snr.name<-grep(paste(i,".snr", sep=""), names(dat$scp), value=T)
            tot.name<-grep(paste(i,".tot", sep=""), names(dat$scp), value=T)
            max.name<-grep(paste(i,".max", sep=""), names(dat$scp), value=T)
            wm.name<-grep(paste(i,".wm", sep=""), names(dat$scp), value=T)
            
            if(dat$scp[j,snr.name]>=snr.min &
                dat$scp[j,max.name]>=max.min &
                dat$scp[j,wm.name]>=wm.min &
                dat$scp[j,wm.name]<=wm.max)
            {dat2[j,i]<-1}
            else{dat2[j,i]<-0}
            }
            }
            return(dat2)
}

#' Function to use the neural networks within the python package
#' @param dat is the RD file input
#' @export
fancyBin<- function(dat){
    pyPharm <- reticulate::import('python_pharmer')

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
ScoreConstPharm <- function(dat, blc=NULL, snr=NULL, der=NULL, snr.lim=3, blc.lim=.03, shws=2){
    t.dat<-dat$t.dat
    if(is.null(blc)){
        blc<-dat$blc
    }else{
        blc<-blc
    }
    if(is.null(snr)){
        snr<-dat$snr
    }else{
        snr<-snr
    }
    if(is.null(der)){
        der<-dat$der
    }else{
        der<-der
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

    dat.der<-der
    
    for(i in cnames){
        s1 <- createMassSpectrum(t.dat[,"Time"],t.dat[,i])
        s3 <- smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=shws)
        bl.th <- estimateBaseline(s3, method="TopHat")[,"intensity"]
        bl.snp <- estimateBaseline(s3, method="SNIP")[,"intensity"]
        eseq <- 1:ceiling((nrow(t.dat)/2))
        lseq <- max(eseq):nrow(t.dat)
        res.tab[i,"bl.diff"] <- mean(bl.th-bl.snp)
        res.tab[i,"earl.bl.diff"] <- mean(bl.th[eseq]-bl.snp[eseq])
        res.tab[i,"late.bl.diff"] <- mean(bl.th[lseq]-bl.snp[lseq])        
    }
    
    for(i in levs){
        res.tab[paste(i,".snr",sep="")] <- apply(snr[wr==i,cnames],2,max)
        res.tab[paste(i,".tot",sep="")] <- apply(blc[wr==i,cnames],2,sum)
        res.tab[paste(i,".max",sep="")] <- apply(blc[wr==i,cnames],2,max)
        res.tab[paste(i,".ph.a.r",sep="")] <-res.tab[paste(i,".tot",sep="")]/res.tab[paste(i,".max",sep="")]

        res.tab[paste(i,".wm",sep="")] <- apply(blc[wr==i,cnames],2,which.max)
        
        ## Derviative measures
        #res.tab[paste(i,".der.tot",sep="")] <- apply(dat.der[wr==i,cnames],2,sum)
        res.tab[paste(i,".der.tot",sep="")] <- apply(dat.der[wr==i,cnames],2,sum)
        #res.tab[paste(i,".der.tot",sep="")] <- apply(na.omit(dat.der[wr==i,cnames]),2,function(x){sum(x[x>0])})
        res.tab[paste(i,".der.max",sep="")] <- apply(na.omit(dat.der[wr==i,cnames]),2,max)
        res.tab[paste(i,".der.min",sep="")] <- apply(na.omit(dat.der[wr==i,cnames]),2,min)
        res.tab[paste(i,".der.wmax",sep="")] <- apply(na.omit(dat.der[wr==i,cnames]),2,which.max)#function(x){which.max(x[5:length(row.names(x))])})
        res.tab[paste(i,".der.wmin",sep="")] <- apply(na.omit(dat.der[wr==i,cnames]),2,which.min)
        #res.tab[c(paste(i,".dn5",sep=""),paste(i,".up5",sep=""))] <- t(apply(t.dat[wr==i,cnames],2,lt5func,x=t.dat[wr==i,1]))
        #res.tab[paste(i,".dn5",sep="")] <- apply(blc[wr==i,cnames],2,dn5func)                
    }
    return(res.tab)
}
#' calculate a table of cell characteristics globally and 
#' within specific windows
#' these specifics should include
#' mean and sd, sum of in window peaks, sum of out of window peaks
#' @export
ScoreConstPharm.2 <- function(dat,t.type=NULL, snr=NULL, der=NULL, snr.lim=3,blc.lim=.03,shws=2){
    require(MALDIquant)
    t.dat<-dat$t.dat
    if(is.null(t.type)){t.type<-'blc'
    }else{t.type<-t.type}
    if(is.null(snr)){snr<-dat$snr
    }else{snr<-snr}
    if(is.null(der)){der<-dat$der
    }else{der<-der}


    wr<-dat$w.dat$wr1

        gtfunc <- function(x,alph){sum(x > alph,na.rm=T)}
        
    lt5func <- function(x,y)
    {
        ltfunc <- function(i){summary(lm(y[i:(i+5)] ~ x[i:(i+5)]))$coefficients[2,3]}
        iseq <- 1:(length(x)-5)
        res <- sapply(iseq,ltfunc)
        return(range(res))
    }

    levs <- setdiff(unique(wr),"")
    cnames <- names(t.dat)[-1]
    res.tab <- data.frame(mean=apply(dat[['blc']][,cnames],2,mean))
    res.tab["sd"] <- apply(dat$blc[,cnames],2,sd)
    res.tab["snr.iws"] <- apply(snr[is.element(wr,levs),cnames],2,sum)
    res.tab["snr.ows"] <- apply(snr[!is.element(wr,levs),cnames],2,sum)
    res.tab["snr.iwc"] <- apply(snr[is.element(wr,levs),cnames],2,gtfunc,alph=snr.lim)
    res.tab["snr.owc"] <- apply(snr[!is.element(wr,levs),cnames],2,gtfunc,alph=snr.lim)

    dat.der<-der
    
    for(i in cnames)
    {
        s1 <- createMassSpectrum(t.dat[,"Time"],t.dat[,i])
        s3 <- smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=shws)
        bl.th <- estimateBaseline(s3, method="TopHat")[,"intensity"]
        bl.snp <- estimateBaseline(s3, method="SNIP")[,"intensity"]
        eseq <- 1:ceiling((nrow(t.dat)/2))
        lseq <- max(eseq):nrow(t.dat)
        res.tab[i,"bl.diff"] <- mean(bl.th-bl.snp)
        res.tab[i,"earl.bl.diff"] <- mean(bl.th[eseq]-bl.snp[eseq])
        res.tab[i,"late.bl.diff"] <- mean(bl.th[lseq]-bl.snp[lseq])        
    }
    for(i in levs)
    {
        res.tab[paste(i,".snr",sep="")] <- apply(snr[wr==i,cnames],2,max)
        res.tab[paste(i,".tot",sep="")] <- apply(dat[[t.type]][wr==i,cnames],2,sum)
        res.tab[paste(i,".max",sep="")] <- apply(dat[[t.type]][wr==i,cnames],2,max)
        res.tab[paste(i,".ph.a.r",sep="")] <-res.tab[paste(i,".tot",sep="")]/res.tab[paste(i,".max",sep="")]

        res.tab[paste(i,".wm",sep="")] <- apply(dat[[t.type]][wr==i,cnames],2,which.max)
    }
    return(res.tab)
}
