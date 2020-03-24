smoothfunc <- function(y,ts,pts=min(13,sum(!is.na(y)))){
    yloe <- loess(y ~ ts,span=pts/length(y))
    ymp <- predict(yloe,newdata=data.frame(ts=ts))
    return(ymp)
}

DespikeSmooth <- function(tmp,ulim=NULL,dlim= NULL){
    print("Starting Despike-smooth timing")
    start_time <- Sys.time()
    wt <- tmp$t.dat
    ts <- tmp$t.dat[,1]
    wtd <- wt[-1,] - wt[-nrow(wt),]
    print(paste("step 1 time",Sys.time()-start_time))
        wtd <- sweep(wtd[,-1],1,wtd[,1],'/')
    print(paste("step 2 time", Sys.time()-start_time))
        wtm <- wtd[-1,]*wtd[-nrow(wtd),]
    print(paste("step 3 time", Sys.time()-start_time))
        wrb <- NumBlanks(tmp$w.dat[,"wr1"])
    print(paste("step 4 time", Sys.time()-start_time)	)
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
    print(paste("step 5 time",Sys.time()-start_time))
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
    print(paste("step 6 time",Sys.time()-start_time))
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
    for(i in names(w1)[sapply(w1,is.character)]){w1[,i] <- "epad"}
    tmp$w.dat <- rbind(tmp$w.dat,w1)
    row.names(tmp$w.dat)<-tmp$w.dat$Time #Lee's additions
    tmp$t.dat <- xdat
    row.names(tmp$t.dat)<-tmp$t.dat$Time #Lee's additions
    tmp$bin['epad']<-0 #Lee's additions
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

TraceBrewer<-function(dat){
    cat('
    #The current flow of our trace cleaning protocol is as follows, and this is
    #what the function automatically fills in for the RD list

    #t.dat>                 Raw data                   t.dat
    #t.dat.pad>             3 End points added at end  NA
    #t.dat.pad.ds.s>        despike and smooth         mp
    #t.dat.pad.ds.s.n>      normalize 0 to 1           t.norm
    #t.dat.pad.ds.s.n.blc>  Baseline Corrected         blc
    ')
    tmp.rd<-dat
    start.time<-proc.time()
    # Kevin has created a new way of creating cleaned up traces.
    # Add a 3 point padding to the end of the experiment which return the trace back to baseline
    # This helps preserve the response shape
    tmp.rd<-PadTdat(tmp.rd)
    print(paste("Completed Padding at:",(proc.time()-start.time))[3])
    # Kevin now uses this to despike and smooth the data
    tmp.rd<-DespikeSmooth(tmp.rd)
    print(paste("Completed Despiking at:",(proc.time()-start.time)[3]))
    # Now what we need to do is provide the analysis with some type of normalized trace
    tmp.rd<-TraceNormal(tmp.rd,'mp')
    print(paste("Completed Normalizing at:",(proc.time()-start.time)[3]))

    # Now do a baseline correction based on the trace padding, and the despike smooth function, and
    # the normalized trace 
    pcp.tmp<-ProcConstPharm(tmp.rd$t.norm)
    print(paste("Completed Baseline Correction at:",(proc.time()-start.time)[3]))
    tmp.rd$blc<-pcp.tmp$blc
    tmp.rd$snr<-pcp.tmp$snr
    # Now perform trace statistics on the specified trace
    tmp.rd$scp<-ScoreConstPharm.2(tmp.rd,'blc')
    print(paste("Completed Window Statistics at:",(proc.time()-start.time)[3]))
    return(tmp.rd)
}

# MAKE SURE BOTH VIDEOS ARE IN SEPERATE FOLDERS,
# ALSO HAVE SEPERATE WR1's that match  exactly what happened
# 1 Now tht we know the videos are ok we will simply do the cell profiler pipeline
# 2 we copy everything from the cell profiler pipeline into the 2 folders that contain the videos
# 3 Now do the pharming harvest and select both folders that contain the seperate videos

# dat1 is the flocation of the frist RD fil ex. "./1 video that yeeted itself/RD.1.Rdata"
# dat2 is the second rd file taht needs to stich to "./2 R3J, TTAA, agonist vid/RD.2.Rdata"
# timeBuffer is the time separation in minutes between the frist trace and the second trace
# newName is the name of the experiment to save 
# Example of how to use
#traceSticher(
#    './1 video that yeeted itself/RD.1.Rdata',
#   './2 R3J, TTAA, agonist vid/RD.2.Rdata',
#   3,
#    'RD.bob')
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

    fName <- paste0(newName, ".Rdata")
    assign(newName, dat1)
    assign(newName, dat1, envir = .GlobalEnv)  
    save(list = newName ,file=fName )

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


