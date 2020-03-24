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

ephyz_import<-function(dat, plotz=T){
    ab_fname<-list.files(pattern=".abf")
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
    plot(dat$ab$s[v_region_to_view]/60,dat$ab$traces[1,v_region_to_view], ylim=c(-7.9943e+01,10), type='l', ylab="Vm (mV)", xlab="Time (s)", col='red')

    
    #plot(dat$ab$s/60, dat$ab$traces[1,], ylim=c(-7.9943e+01,10), type='l', ylab="Vm (mV)", xlab="Time (s)")
    end_time<-Sys.time()
    print(paste("Time to plot=",round(end_time-start_time, digits=3) ) )
    bringToTop(-1)
    #Adding window region to the ephyz plot
    wr<-dat$w.dat$wr1
    levs<-setdiff(unique(as.character(wr)),"")
    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs]
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs]
    y1s <- rep(par("usr")[3]+yinch(.11),length(x1s))
    y2s <- rep(par("usr")[3]+yinch(.01),length(x1s))
    rect(x1s,y1s,x2s,y2s, col='black', border='white')
    par(xpd=T)
    levs.loc<-tapply(dat$w.dat[,"Time"],as.factor(wr),mean)[levs]
	
	levs_cex <- nchar(levs)
	levs_cex[ levs_cex<=12 ] <- 1
	levs_cex[ levs_cex > 12 ] <- 12/levs_cex[ levs_cex>12 ]*1.3
    text(levs.loc, par("usr")[3]-yinch(1), levs,pos=3, cex=levs_cex, srt=90)	
}

###################################################################
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

###########################################################################
#TRACE ALIGNMENT
ephyz_ca_aligner<-function(dat, delay_time=3.2, deelay=T){
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
    if(deelay){
        dat$ab$traces <- cbind( dat$ab$traces, dat$ab$traces[,1:delay_samples] )
    }else{
        dat$ab$traces <- cbind(dat$ab$traces[,1:delay_samples], dat$ab$traces)
    }

    dat$ab$traces <- dat$ab$traces[ ,-(1:delay_samples) ]

    end.time<-Sys.time()
    print(paste("It took:", end.time-start_time ) )
    return(dat)
}

