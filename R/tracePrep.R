smoothfunc <- function(y,ts,pts=min(13,sum(!is.na(y)))){
    yloe <- loess(y ~ ts,span=pts/length(y))
    ymp <- predict(yloe,newdata=data.frame(ts=ts))
    return(ymp)
}

DespikeSmooth <- function(tmp,ulim=NULL,dlim= NULL){
    #print("Starting Despike-smooth timing")
    start_time <- Sys.time()
    wt <- tmp$t.dat
    ts <- tmp$t.dat[,1]
    wtd <- wt[-1,] - wt[-nrow(wt),]
    #print(paste("step 1 time",Sys.time()-start_time))
        wtd <- sweep(wtd[,-1],1,wtd[,1],'/')
    #print(paste("step 2 time", Sys.time()-start_time))
        wtm <- wtd[-1,]*wtd[-nrow(wtd),]
    #print(paste("step 3 time", Sys.time()-start_time))
        wrb <- NumBlanks(tmp$w.dat[,"wr1"])
    #print(paste("step 4 time", Sys.time()-start_time)	)
        if(is.null(ulim) | is.null(dlim))
        {
            qvals <- quantile(
                        as.vector(
                            as.matrix(
                                wtm[grep("blank",wrb[1:nrow(wtm)]),]
                            )
                        )
                    ,probs=c(0,.001,.5,.999,1))
        }
    #print(paste("step 5 time",Sys.time()-start_time))
        if(is.null(dlim)){dlim <- qvals[2]}
        if(is.null(ulim)){ulim <- qvals[4]}	
        wtm.z <- wtm[1,]
        wtm.z[,] <- 0
        wtm <- rbind(wtm.z,wtm,wtm.z)
        wtrm <- wtm < dlim
        x <- wt[,1]
        wt <- wt[,-1]
        wt[wtrm] <- NA
        mp <- sapply(wt, smoothfunc, ts=x)
    #print(paste("step 6 time",Sys.time()-start_time))
        tmp$mp <- tmp$t.dat
        tmp$mp[,-1] <- mp
        return(tmp)
}

LinearBLextend <- function(x,y,intvl=3,plotit=F){
    require(MASS)
    #break into intervals of 3 minute.
    time.tot <- max(x)-min(x)
    gaps <- ceiling(time.tot/intvl)
    gap.fac <- sort(rep(seq(1,gaps),length.out=length(x)))
    gap.i <- tapply(y,gap.fac,which.min)
    gap.i <- gap.i + match(unique(gap.fac),gap.fac)-1
    mindat <- data.frame(x=x[gap.i],y=y[gap.i])
    rlt <- rlm(y ~ x,data=mindat)	
    if(plotit)
    {
        plot(x,y)
        points(x[gap.i],y[gap.i],pch=16,col="red")
        abline(rlt,lty=2,col="blue",lwd=2)
        #lines(x[gap.i],predict(xloe))
        #points(x1,predict(rlt,newdata=data.frame(x=x1)),pch=16,col="blue",cex=2)
    }
    return(coefficients(rlt))
}

#add points to the end of the t.dat
PadTdat <- function(tmp,n=5){
    xdat <- tmp$t.dat
    r1 <- sapply(tmp$t.dat[,-1],LinearBLextend,x=tmp$t.dat[,1],plotit=F)
    r1 <- data.frame(t(r1))
    tmp$scp[row.names(r1),"rlm.a"] <- r1[,1]
    tmp$scp[row.names(r1),"rlm.b"] <- r1[,2]
    x <- xdat[,1]
    dx <- x[-1]-x[-length(x)]
    
    tmax <- max(x)
    tseq <- seq(1,n)*median(dx)+tmax
    
    r2 <- data.frame(t(r1[,"x"] %*% t(tseq) + r1[,1]))
    names(r2) <- row.names(r1)
    r2[,"Time"] <- tseq
    
    xdat <- rbind(xdat,r2[,names(xdat)])
    w1 <- tmp$w.dat[1:n,]
    w1[,"Time"] <- tseq	
    
    for(i in names(w1)[sapply(w1,is.character)]){
        w1[,i] <- "epad"
    }
    tmp$w.dat <- rbind(tmp$w.dat,w1)
    row.names(tmp$w.dat) <- tmp$w.dat$Time #Lee's additions
    tmp$t.dat <- xdat
    row.names(tmp$t.dat) <- tmp$t.dat$Time #Lee's additions
    return(tmp)
}

TraceNormal<-function(dat, t.type='blc'){
    tmp.dat<-dat[[t.type]]

    for(k in 1:length(colnames(dat$mp))){

        tmp.dat[,k]<-dat$mp[,k]-min(dat$mp[,k])
        tmp.dat[,k]<-tmp.dat[,k]/max(tmp.dat[,k])

    }
    tmp.dat[,1]<-dat$t.dat[,1]

    dat$t.norm<-tmp.dat
    return(dat)
}

#' take a defined window vector and
#' number of the contiguos blank regions ("")
NumBlanks <- function(x){
    nw <-  as.character(x)
    xlen <- length(x)
    bl.cnt <- 1
    mi <- match("",nw)
    while(!is.na(mi) & bl.cnt < 20)
    {
        mi2 <- mi+1
        while((x[mi2] == "") & mi2 <= xlen){mi2 <- mi2+1}
        nw[mi:(mi2-1)] <- paste("blank",bl.cnt,sep="")
        bl.cnt <- bl.cnt+1
        mi <- match("",nw)
    }
    return(as.factor(nw))
}

#Display the analysis of a single trace 
#dat is the trace dataframe with "Time" in the first column and cell trace intensities in subsequent columns
#i is the index column to be analyzed and displayed.
#shws is the smoothing half window size
#Plotit is a flag indicating that the results should be ploted or not.
#wr is the response window factor 
#SNR.lim is the signal to noise ratio limit for peak detection
#bl.meth is the method for baseline correction.
PeakFunc2 <- function(dat, i, shws=2, phws=20, Plotit=F, wr=NULL, SNR.lim=2, bl.meth="TopHat", lmain=NULL){
    s1 <- MALDIquant::createMassSpectrum(dat[,"Time"],dat[,i])
    
    if(shws > 1){
        s3 <- MALDIquant::smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=shws)
    }else{
        s3 <- s1
    }

    s4 <- MALDIquant::removeBaseline(s3, method=bl.meth)
    Baseline <- MALDIquant::estimateBaseline(s3, method=bl.meth)
    p <- MALDIquant::detectPeaks(s4, method="MAD", halfWindowSize=phws, SNR=SNR.lim)
    if(Plotit){
        bSnip <- MALDIquant::estimateBaseline(s3, method="SNIP")
        bTopHat <- MALDIquant::estimateBaseline(s3, method="TopHat")

        xlim <- range(mass(s1)) # use same xlim on all plots for better comparison
        ylim <- c(-.1,1.4)
        #ylim <- range(intensity(s1))
        plot(s1, main=paste(lmain,i),xlim=xlim,ylim=ylim,xlab="Time (min)", xaxt="n")
        axis(1, at=seq(0, length(dat[,1]), 5))  
        if(length(wr) > 0)
        {
            levs <- setdiff(unique(wr),"")
            levs <- setdiff(levs,grep("blank",levs,value=T))
            x1s <- tapply(dat[,"Time"],as.factor(wr),min)[levs]
            x2s <- tapply(dat[,"Time"],as.factor(wr),max)[levs]
            y1s <- rep(min(ylim)-.2,length(x1s))
            y2s <- rep(max(ylim)+.2,length(x1s))
            #cols <- rainbow(length(x1s))
            rect(x1s,y1s,x2s,y2s,col="lightgrey")
            #points(dat[,"Time"],as.integer(wr=="")*-1,pch=15,cex=.6)
            ## for(j in levs)
            ## {
            ##     x1 <- mass(s3)[min(grep(j,wr))]
            ##     x2 <- mass(s3)[max(grep(j,wr))]
            ##     y1 <- min(ylim)-.2
            ##     y2 <- max(ylim)+.2
            ##     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col="lightgrey",lwd=.1)
            ## }
            text(dat[match(levs,wr),"Time"],rep(-.1,length(levs)),levs,pos=4,offset=0,cex=.5)
        }
        
        lines(s3,lwd=3,col="cyan")
        lines(s1)
        lines(bSnip, lwd=2, col="red")
        lines(bTopHat, lwd=2, col="blue")
        lines(s4,lwd=2)
        
        if((length(p) > 0)){
            points(p)
            ## label top 40 peaks
            top40 <- intensity(p) %in% sort(intensity(p), decreasing=TRUE)[1:40]
            labelPeaks(p, index=top40, underline=TRUE,labels=round(snr(p)[top40],2))
        }
    }

    return(list(peaks=p,baseline=Baseline,dat=s4))
}

#' This is our trace cleaning protocol
#' @export 
#' @param dat is the RD.experiment to input
#' @param bscore logical, if true all scores will be replaced by bscore2 function
#' @param blcPrep logical, if this is true, an expensive baseline correction schema occurs. If false only some SNR occurs
TraceBrewer <- function(dat, bscore = F, blcPrep = T, verbose = F){

    tmp.rd <- dat
    start.time <- proc.time()
    if(blcPrep){
        if(verbose){
        cat('
#The current flow of our trace cleaning protocol is as follows, and this is\n
#what the function automatically fills in for the RD list

#t.dat>                 Raw data                   t.dat
#t.dat.pad>             3 End points added at end  NA
#t.dat.pad.ds.s>        despike and smooth         mp
#t.dat.pad.ds.s.n>      normalize 0 to 1           t.norm
#t.dat.pad.ds.s.n.blc>  Baseline Corrected         blc
        ')
        }
        # Kevin has created a new way of creating cleaned up traces.
        # Add a 3 point padding to the end of the experiment which return the trace back to baseline
        # This helps preserve the response shape
        tmp.rd <- PadTdat(tmp.rd)
        #print(paste("Completed Padding at:",(proc.time()-start.time))[3])
        # Kevin now uses this to despike and smooth the data
        tmp.rd <- DespikeSmooth(tmp.rd)
        #print(paste("Completed Despiking at:",(proc.time()-start.time)[3]))
        # Now what we need to do is provide the analysis with some type of normalized trace
        tmp.rd <- TraceNormal(tmp.rd, 'mp')
        #print(paste("Completed Normalizing at:",(proc.time()-start.time)[3]))
        # Now do a baseline correction based on the trace padding, and the despike smooth function, and
        # the normalized trace 
        pcp.tmp <- ProcConstPharm(tmp.rd$t.norm)
        #print(paste("Completed Baseline Correction at:",(proc.time()-start.time)[3]))
        tmp.rd$blc <- pcp.tmp$blc
        tmp.rd$snr <- pcp.tmp$snr
    }else{
        # We need to grab the snr data to proceed
        dat <- tmp.rd$blc
        t.names <- names(dat)[-1]#Time in first column
        dat1.snr <- dat #peak calls stored as SNR
        dat1.snr[,t.names] <- 0
        for(i in 1:length(t.names)){
            s1 <- MALDIquant::createMassSpectrum(dat[,"Time"],dat[,t.names[i]])
            peaks <- MALDIquant::detectPeaks(s1, method="MAD", halfWindowSize=20, SNR=2)
            dat1.snr[match(MALDIquant::mass(peaks),dat[,1]),t.names[i]] <- MALDIquant::snr(peaks)
        }
        tmp.rd$snr <- dat1.snr
    }
    
    # Now perform trace statistics on the specified trace
    tmp.rd$scp <- ScoreConstPharm(tmp.rd, t.type = 'blc')

    # Skip scoring if specified. This helps to not overwrite hardwork done.
    if(bscore){
        tmp.rd <- bscore2(tmp.rd)
    }

    cat("\nCompleted Brew. CHEERS!:",round((proc.time()-start.time)[3], digits = 3),' seconds\n')

    # This is to get rid of the padding
    epadRM <- tmp.rd$w.dat$wr1 != "epad"
        
    tmp.rd$blc <- tmp.rd$blc[epadRM,]
    tmp.rd$w.dat <- tmp.rd$w.dat[epadRM,]
    tmp.rd$t.dat <- tmp.rd$t.dat[epadRM,]
    return(tmp.rd)
}

#' This fucntion stiches two videos together
#' 
#' Make sure both Videos are in seperate Folders
#' ALSO HAVE SEPERATE WR1's that match  exactly what happened
#' 1 Now tht we know the videos are ok we will simply do the cell profiler pipeline
#' 2 we copy everything from the cell profiler pipeline into the 2 folders that contain the videos
#' 3 Now do the pharming harvest and select both folders that contain the seperate videos
#' @param dat1Loc is the flocation of the frist RD fil ex. "./1 video that yeeted itself/RD.1.Rdata"
#' @param dat2Loc is the second rd file taht needs to stich to "./2 R3J, TTAA, agonist vid/RD.2.Rdata"
#' @param timeBuffer is the time separation in minutes between the frist trace and the second trace
#' @param newName is the name of the experiment to save 
#' @export
traceSticher <- function(dat1Loc, dat2Loc, timeBuffer = 3, newName = NULL){
    cat("
    # MAKE SURE BOTH VIDEOS ARE IN SEPERATE FOLDERS, 
    # ALSO HAVE SEPERATE WR1's that match  exactly what happened
    # 1 Now tht we know the videos are ok we will simply do the cell profiler pipeline
    # 2 we copy everything from the cell profiler pipeline into the 2 folders that contain the videos
    # 3 Now do the pharming harvest and select both folders that contain the seperate videos
    
    traceSticher(
        './1 video that yeeted itself/RD.1.Rdata', 
        './2 R3J, TTAA, agonist vid/RD.2.Rdata', 
        3, 
        'RD.bob')
    ")
    dat1 <- get(load(dat1Loc))
    dat2 <- get(load(dat2Loc))

    # fName <- paste0(newName, ".Rdata")
    # assign(newName, dat1)
    # assign(newName, dat1, envir = .GlobalEnv)  
    # save(list = newName ,file=fName )

    #How long were your experiments seperated?
    #I think at least 7 min is mandatory
    #time_buffer <- 3
    #look at the max time value in the t.dat on the first experiemnt
    #add it to the time_buffer
    #look at the max time value in the t.dat on the first experiemnt
    #add it to the time_buffer
    time_to_change<-max(dat1$t.dat[1])+timeBuffer
    #now increase every time value in the second experiment by the 
    #value above
    changed_time <- dat2$t.dat[,1] + time_to_change

    #Now we need make unique names of rhte windew regions
    names_to_change<-setdiff(unique(dat2$w.dat[,"wr1"]),"")
    new_names<-paste(names_to_change,"_2", sep="")
    #dat2$w.dat["ignore"]<-0
    for(i in 1:length(names_to_change)){
        dat2$w.dat[ dat2$w.dat["wr1"]==names_to_change[i] ,2]<-new_names[i]
    }

    #Select the t.dat, blc and, w.dat
    t_to_view<-c('t.dat', 't.340', 't.380', 'w.dat')
    ##now merg ethe datasets
    for(i in 1:length(t_to_view)){
        #we need to get rid of the epad
        notEpadLogic <- dat1$w.dat$wr1 != "epad"
        dat1[[ t_to_view[i] ]] <- dat1[[ t_to_view[i] ]][notEpadLogic, ]
        #change the row names first
        
        notEpadLogic <- dat2$w.dat$wr1 != "epad"
        dat2[[ t_to_view[i] ]] <- dat2[[ t_to_view[i] ]][notEpadLogic, ]
        row.names( dat2[[ t_to_view[i] ]] ) <- as.character(changed_time)
        #change the first collumn value
        dat2[[ t_to_view[i] ]][1]<-changed_time
        #combine the experiments trace dataframes together
        dat1[[ t_to_view[i] ]]<-rbind(dat1[[ t_to_view[i] ]], dat2[[ t_to_view[i] ]])
    }

    #And reprocess the data
    levs<-setdiff(unique(as.character(dat1$w.dat[,2])),"")
        snr.lim=5;hab.lim=.05;sm=2;ws=3;blc="SNIP"
    pcp <- ProcConstPharm(dat1,sm,ws,blc)
    scp <- ScoreConstPharm(dat1,pcp$blc,pcp$snr,pcp$der,snr.lim,hab.lim,sm)
    bin <- bScore(pcp$blc,pcp$snr,snr.lim,hab.lim,levs,dat1$w.dat[,"wr1"])

    dat1$scp<-scp
    dat1$snr<-pcp$snr
    dat1$blc<-pcp$blc
    dat1$bin<-bin
    
    #dat1 <- TraceBrewer(dat1)

    fName <- paste0(newName, '.Rdata')
    assign(newName, dat1)
    assign(newName, dat1, envir = .GlobalEnv)  
    save(list = newName ,file=fName )
}


