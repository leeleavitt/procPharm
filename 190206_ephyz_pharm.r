#################################################################
#Package to look at the data collected during ca img and ephyz
#################################################################

###########################################
#function to import data
#if the data is properly structured then this should work
#make sure you are in the right working directory
if( !library(abf2, logical.return=T) ){
    install.packages('abf2')
}else{}

if( !library(readABF, logical.return=T) ){
    install.packages('readABF')
}else{}


ephyz_import<-function(dat, abName= NULL, plotz=F){
    if(is.null(abName)){
		ab_fname<-select.list(list.files(pattern=".abf"))
	}else{ab_fname <- abName}
    dat$ab<-abf2::abfload(ab_fname)
    if(plotz){
        cat("\nI'm going to plot the trace for you\n")
        ephyz_total_plotter(dat)
    }else{}
    return(dat)
}
    
ephyz_total_plotter<-function(dat){
    bringToTop(-1)
    cat("\nThis function plots both the Calcium imaging trace and the \n electrophysiology trace it may take some time\n")

    ###########################################
    #Plotting Entire Trace
    dat.name<-deparse(substitute(dat))
    dev.new(width=14,height=8)
    par(mfrow=c(2,1))
    PeakFunc7(dat, "X.1",t.type="t.dat", bcex=1.2, zf=80, dat.n=dat.name )
    par(xpd=F)
    abline(h=0)
    cat("\n The Data has been imported, i am working on plotting it\n")
    start_time<-Sys.time()
    v_region_to_view<-min(which(dat$ab$s > min(dat$t.dat$Time)*60, arr.ind=T)):max(which(dat$ab$s<max(dat$t.dat$Time)*60, arr.ind=T))
    plot(dat$ab$s[v_region_to_view]/60,dat$ab$traces[1,v_region_to_view], ylim=c(-7.9943e+01,10), type='l', ylab="Vm (mV)", xlab="Time (s)")

    
    #plot(dat$ab$s/60, dat$ab$traces[1,], ylim=c(-7.9943e+01,10), type='l', ylab="Vm (mV)", xlab="Time (s)")
    end_time<-Sys.time()
    print(paste("Time to plot=",round(end_time-start_time, digits=3) ) )
    bringToTop(-1)
    #Adding window region to the ephyz plot
    wr<-dat$w.dat$wr1
    levs<-setdiff(unique(as.character(wr)),"")
    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs]
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs]
    y1s <- rep(par("usr")[4],length(x1s))
    y2s <- rep(par("usr")[3],length(x1s))
    rect(x1s,y1s,x2s,y2s)
    par(xpd=T)
    levs.loc<-tapply(dat$w.dat[,"Time"],as.factor(wr),mean)[levs]
    text(levs.loc,par("usr")[3],levs,pos=3,offset=-4.2,cex=1.2, srt=90)	
}

#Trace zoom
ephyz_zoomer<-function(dat, plot.new=F, t_type="t.dat", ephyz_ymin=-70){
    cat("\nThis function allows you to zoom in on the calcium imaging Trace\n")
    if(plot.new){
        graphics.off()
        dev.new(width=14,height=8)

    }else{}
    main_win<-dev.cur()
    layout(rbind(c(1,1),c(2,2)))
    #plot the main trace
    PeakFunc7(dat, "X.1",t.type=t_type, bcex=1.2, zf=80 )
    #identify region to zoom
    bringToTop(-1)
    cat("\nPlease identify the region you would like to more closely observe\nSelect the right side first followed by the left side\nPRESS ENTER TO CONTINUE\n")
    scan(n=1)
    
    continue<-"T"
    plot_counter<-1
    while(continue){
        click1<-locator(n=1,type='p')
        click2<-locator(n=1,type='p')
        par(xpd=F)
        abline(v=click2$x, col="blue", lwd=3)
        abline(v=click1$x, col="blue", lwd=3)
        if(plot_counter>1){dev.new(width=12, height=4)} 
        
        region_to_view<-min(which(dat$t.dat$Time>click1$x, arr.ind=T)):max(which(dat$t.dat$Time<click2$x, arr.ind=T))
        
        v_region_to_view<-min(which(dat$ab$s>click1$x*60, arr.ind=T)):max(which(dat$ab$s<click2$x*60, arr.ind=T))
        
        plot(dat$ab$s[v_region_to_view]/60,dat$ab$traces[1,v_region_to_view],type="l", ylab="", col="red", xlab="",yaxt="n",xaxt='n',ylim=c(ephyz_ymin,10) )
        axis(side = 4)
        mtext("Vm (mV)", side = 4, line = 3)

        par(new=T)

        plot(dat$t.dat$Time[ region_to_view ],dat$t.dat[region_to_view,"X.1"], type="l", xlab="Time", ylab=expression(paste(Delta,"F/F")))

        points(dat$t.dat$Time[region_to_view],dat$t.dat[region_to_view,"X.1"], pch=15)

        continue<-scan(n=1, what='character')
        if(length(continue)==0){continue<-"T"}
        plot_counter<-plot_counter+1
        dev.set(main_win)
    }
}

#TRACE ALIGNMENT
ephyz_ca_aligner<-function(dat, delay_time=3.2){
    cat("\nThis function requires you to align the traces\n\nIt is best if you use ephyz_zoomer to get an idea on how much to shift the ephyz trace by\nWe start by aligning the traces of calcium and ephyz. The delay time is then used to shift the ephyz trace to the proper location.")

    ###########################################################################
    #TRACE ALIGNMENT

    #Which has the max value for starting poit
    ( start_time<-max( min(dat$t.dat$Time)*60, min(dat$ab$s) ) )

    #find minimun of endtimes
    ( end_time<-min( max(dat$t.dat$Time) *60, max(dat$ab$s) ) )

    #get start sample for ephyz and ca imaging
    ( sampling_rate_ca_img<- 1/round( ( dat$t.dat$Time[3]*60-dat$t.dat$Time[2]*60 ), digits=0 ) )
    ( sampling_rate_ephyz<- 1/ ( dat$ab$s[2] - dat$ab$s[1] ) )

    ( start_sample_ca_img<- floor( start_time * sampling_rate_ca_img ) )
    ( start_sample_ephyz<- floor(start_time * sampling_rate_ephyz) )

    ( end_sample_ca_img<- floor( end_time * sampling_rate_ca_img ) )
    ( end_sample_ephyz<- floor( end_time * sampling_rate_ephyz) )

    #Adjusted calcium imaging Trace
    
    ( ca_img_time_seq <- dat$t.dat$Time[ start_sample_ca_img : end_sample_ca_img ] )
    ( dat$t.dat <- dat$t.dat[ start_sample_ca_img : end_sample_ca_img,  ] )

    #Adjusted ephyz trace
    ( dat$ab$s <- dat$ab$s[start_sample_ephyz : end_sample_ephyz] )
    ( dat$ab$traces <- dat$ab$traces[ ,start_sample_ephyz : end_sample_ephyz] )


    ###########################################################################
    #EPHYZ TRACE ADVANCE
    #########################################################################
    #Now that the traces have been aligned we need to create a trace delay function

    ( delay_samples <- floor( delay_time * sampling_rate_ephyz ) )

    #hold the time in a diff location
    #time_tmp<-dat$ab$s

    #first duplicate the first row and append it to the top of the array
    start_time<-Sys.time()

    dat$ab$traces <- cbind( dat$ab$traces, dat$ab$traces[,1:delay_samples] )
    dat$ab$traces <- dat$ab$traces[ ,-(1:delay_samples) ]

    end.time<-Sys.time()
    print(paste("It took:", end.time-start_time ) )
    return(dat)
}

######################################################################################
# Series for analyzing current steps
# See "Y:/Halen/Ephyz/191210.40.m.m3.p1 L2 TTA and RIIIJ/c.860" For example
#Steps Plotter
stepPlotter<-function(){
	#install.packages('readABF')
	#require(readABF)

	abfFiles <- list.files(pattern='[.]abf$')
	experiments <- select.list(abfFiles, multiple=T)
	for(j in 1:length(experiments)){
		tmpabf<- readABF(experiments[j])
		totalSteps <-length(tmpabf$data)

		#png(paste0(experiments[j],'.png'), 2048, 4096, res=300)
		pdf(paste0(experiments[j],'.pdf'), 8, 12)
		par(mfrow=c(totalSteps,2))

		for(i in totalSteps:1){
			par(mai=c(0,.5,0,0))
			plot(tmpabf$data[i][[1]][,2], type='l', ylim=c(-200,200))
			plot(tmpabf$data[i][[1]][,1], type='l', ylim=c(-200,50),xaxt='n')
		}

		dev.off()
	}
}

#Load in the package to read in the abf files
require(readABF)

#Function To count action Potentials
# apInterval: Shortest interval inbetween action potentials 
# apThresh: When action potential Fires
apCounter <- function(trace, apInterval = 100, apThresh = -10){
    #Trace is a step for firing.
    #Increas if ap has fired
    apCounter <- 0
    # Increase  Until ap, Then Reset.
    timeCounter <- 0

    for(i in 1:length(trace)){
        timeCounter <- timeCounter + 1
        # If thye value rises above 0 an action potential has been generated
        # If it hasn't been more than 1000 
        # Now in crease our counter
        if(trace[i] > apThresh && timeCounter > apInterval){
            apCounter <- apCounter + 1
            timeCounter <- 0 
        }
    }

    return(apCounter)
}

# Function to Find Steps
stepFinder <- function(abData){
    steps <- c()
    for(i in 1:length(abData$data)){
        #Max value
        maxVal <- max(abData$data[[i]][,2])
        minVal <- min(abData$data[[i]][,2])

        if(abs(maxVal) > abs(minVal)){
            steps[i] <- round(maxVal, digits=-1)
        }else{
            steps[i] <- round(minVal, digits=-1)
        }
    }
    return(steps)
}

# Function To count action potentials in each sweep
# Takes in an ab object loaded from the readABF function
stepCounter <- function(abDf, apInterval = 100, apThresh = -10){
    steps <-  stepFinder(abDf)
    apPerStep <- c()
    for(i in 1:length(abDf$data)){
        apPerStep[i] <- apCounter(abDf$data[[i]][,1],apInterval, apThresh) 
    }
    names(apPerStep)<-steps

    return(apPerStep)
}

#Function to plot differential comparison.
#Input needs to be a differenc
diffAbsPlotter <- function(expDiff){
    expDiffMat <- as.matrix(t(expDiff))

    #COLOR
    colorbins <- 10
    colorFunction <- colorRampPalette(c('red','white','blue'))
    cols <- colorFunction(colorbins)[as.numeric(cut(expDiffMat, seq(-1.1,1.1,length.out=colorbins)))]

    censusBarPlot <- barplot(as.vector(expDiffMat), 
        col=cols, 
        xlim = c(-1.1,1.1), 
        space = 1, 
        xlab=c('Loss of AP\'s'),
        ylab='Current Step nA',
        horiz=T,
        main='Change in Action Potential frequency'
        )

    #Label it
    par(xpd=T)
    x_pos <- expDiffMat
    x_pos_inch <- x_pos + ( sign(x_pos) * xinch(0.25) )

    y_pos <- censusBarPlot
    barLabels <- names(expDiff)
    text(x_pos_inch, y_pos, barLabels, font=1)
}

diffPlotter <- function(expComps){
    par(xpd=FALSE)
    lims <- c(-1*max(expComps), max(expComps))
    cols <- c('red', 'blue')

    mainName <- paste(colnames(expComps)[1], 'vs', colnames(expComps)[2])
    barplot(-1*expComps[,1], col = cols[1], xlim = lims, horiz=T, yaxt = 'n', xaxt='n', main = mainName, xlab = '# of Action Potentials')
    abline(v=1)
    barplot(expComps[,2], col = cols[2], horiz=T, add=T, las=1,  xaxt='n')
    axis(side = 1, seq(-max(expComps),max(expComps),1), las=2)
    abline(v= seq(-max(expComps),max(expComps),1))
}

#This Function onlt works when you are located in the directory
#apThresh is the threshhold for the action potentials
stepsSpikeComparer <- function(apInterval=50, apThresh=-20){
    abFiles <- list.files(pattern='abf$')

    abToCompare <- c()
    abToCompare[1] <- select.list(abFiles, multiple=F, title="SELECT First")
    abToCompare[2] <- select.list(abFiles, multiple=F, title="SELECT Second")

    abList <- list()
    expCompDf <- list()
    for(i in 1:length(abToCompare)){
        abName <- sub('[.]abf','',abToCompare[i])
        abList[[abName]] <- readABF(abToCompare[i])
        expCompDf[[i]] <- stepCounter(abList[[abName]],apInterval,apThresh)
    }
    expComps <- Reduce(cbind,expCompDf)
    colnames(expComps) <- names(abList)

    print(expComps)

    ##################################
    expDiff <- (expComps[,1] - expComps[,2]) / (expComps[,1] + expComps[,2])

    dev.new(width = 20, height = 10)
    par(mfrow = c(1,2))
    diffAbsPlotter(expDiff)
    diffPlotter(expComps)
}