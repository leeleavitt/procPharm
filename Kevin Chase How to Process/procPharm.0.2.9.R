#Display the analysis of a single trace 
#dat is the trace dataframe with "Time" in the first column and cell trace intensities in subsequent columns
#i is the index column to be analyzed and displayed.
#shws is the smoothing half window size
#Plotit is a flag indicating that the results should be ploted or not.
#wr is the response window factor 
#SNR.lim is the signal to noise ratio limit for peak detection
#bl.meth is the method for baseline correction.
PeakFunc2 <- function(dat,i,shws=2,phws=20,Plotit=F,wr=NULL,SNR.lim=2,bl.meth="TopHat",lmain=NULL)
{
    library("MALDIquant")
    s1 <- createMassSpectrum(dat[,"Time"],(dat[,i]-(min(dat[,i])-.1)))
    if(shws > 1)
        s3 <- smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=shws)
    else
        s3 <- s1
    if(Plotit)
    {
        bSnip <- estimateBaseline(s3, method="SNIP")
        bTopHat <- estimateBaseline(s3, method="TopHat")
    }
    s4 <- removeBaseline(s3, method=bl.meth)
    Baseline <- estimateBaseline(s3, method=bl.meth)
    p <- detectPeaks(s4, method="MAD", halfWindowSize=phws, SNR=SNR.lim)
    if(Plotit)
    {
        xlim <- range(mass(s1)) 
        ylim <- c(-.1,max(intensity(s1)))
#        ylim <- range(intensity(s1))
        plot(s1, main=paste(lmain,i),xlim=xlim,ylim=ylim,xlab="Time (min)")
#		axis(1, at=seq(0, length(dat[,1]), 5))  
        if(length(wr) > 0)
        {
            levs <- setdiff(unique(wr),"")
            levs <- setdiff(levs,grep("blank",levs,value=T))
            x1s <- tapply(dat[,"Time"],as.factor(wr),min)[levs]
            x2s <- tapply(dat[,"Time"],as.factor(wr),max)[levs]
            y1s <- rep(min(ylim)-.2,length(x1s))
            y2s <- rep(max(ylim)+.2,length(x1s))
#            cols <- rainbow(length(x1s))
            rect(x1s,y1s,x2s,y2s,col="lightgrey",lwd=.1)
#            points(dat[,"Time"],as.integer(wr=="")*-1,pch=15,cex=.6)
            ## for(j in levs)
            ## {
            ##     x1 <- mass(s3)[min(grep(j,wr))]
            ##     x2 <- mass(s3)[max(grep(j,wr))]
            ##     y1 <- min(ylim)-.2
            ##     y2 <- max(ylim)+.2
            ##     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col="lightgrey",lwd=.1)
            ## }
            text(dat[match(levs,wr),"Time"],rep(c(-.05,-.1),length=length(levs)),levs,pos=4,offset=0,cex=.5)
        }
        
        lines(s3,lwd=1,col="cyan")
        lines(s1)
        lines(bSnip, lwd=.5, col="red")
        lines(bTopHat, lwd=.5, col="blue")
        points(s4,pch=16,cex=.75)
    }
    if((length(p) > 0)&Plotit)
    {
        points(p)
        ## label top 40 peaks
        top40 <- intensity(p) %in% sort(intensity(p), decreasing=TRUE)[1:40]
        labelPeaks(p, index=top40, underline=TRUE,labels=round(snr(p)[top40],2))
    }
    return(list(peaks=p,baseline=Baseline,dat=s4))
}

#read three files.
#t.name is the name of the csv file with the trace data
#w.name is the name of the csv file with response windows
#c.name is the name of the file with cell details.

#this function reads csv formated trace files exported from the NIS program
#rows are time intervals columns are individual cell data
#in particular there should be a "Time" column and the name of each column
#should be something like: "X.178..Ratio.340.nm.380.nm"  where the number after the "X"
#is the unique identifier for that cell in this experiment.
ReadTraceFile <- function(fname)
{
    dat <- read.csv(fname)
    names(dat) <- sub("\\.\\..*","",names(dat))#Time column may be duplicated, remove all but the first
    g.names <- grep("X",names(dat),value=T)
    if(length(g.names) == 0)
    {stop("Incorrect cell names in trace input file")}
    ti <- grep("Time",names(dat),value=T)
    if(length(ti) == 0)
    {stop("No Time column in the trace input file")}
    dat <- dat[c(ti[1],g.names)]#this will fail with no Time column
    
    dat["Time"] <- ConvertTime(as.single(dat[,"Time"]))
#    g.names <- sub("X",fname,names(dat))
#    names(dat) <- g.names
    return(dat)
}

#convert time to minutes
ConvertTime <- function(x)
{
	if(is.numeric(x))
	{
		if(max(x) > 1000000)#in ms
		{
			retval <- x/60000
		}
		else if(max(x) > 1500) #in seconds
		{
			retval <- x/60	
		}
	}
	else
	{			
	vals <- strsplit(x,":")[[1]]
	retval <- NA
	if(length(vals)==3)
	{
		retval <- as.integer(vals[1])*60+as.integer(vals[2])+as.single(vals[3])/60
	}
	if(length(vals)==2)
	{
		retval <- as.integer(vals[1])+as.single(vals[2])/60
	}
	}
	return(retval)
}
#Read whole data dump file 
#fname is the name of the dump file
#wrdef is the name of the csv files that defines the window regions.
#if wrdef is left out a fake wr1 will be made.
#rd.name is the name resulting data list (THIS SHOULD START WITH RD)
#if this is left out RD"Date" will be used.
ReadDataDump <- function(fname,wrdef=NULL,rd.name=NULL,sep=";",encode="UCS-2LE")
{
	tmp <- read.delim(fname,fileEncoding=encode,sep=sep,header=T)
	all.names <- names(tmp)
	time.name <- grep("Time",all.names,value=T,ignore=T)[1]
	if(time.name != "Time..ms."){warning(paste(time.name,"assumed to be in ms"))}
	id.name <- grep("ROI.ID",all.names,value=T,ignore=T)[1]
	if(id.name != "ROI.ID"){warning(paste(id.name,"assumed to be it ROI.ID"))}
	area.name <- grep("Area",all.names,value=T,ignore=T)[1]
	if(is.na(area.name)){stop("no ROI.Area data")}
	else{if(area.name != "ROI.Area"){warning(paste(area.name,"assumed to be ROI.Area"))}}
	ratio.name <- grep("Ratio",all.names,value=T,ignore=T)
	if(is.na(ratio.name)){stop("no ratio data")}
	else{if(ratio.name != "Ratio.340.380"){warning(ratio.name,"assumed to be Ratio data")}}
	cx.name <- grep("Center.X",all.names,value=T,ignore=T)
	if(is.na(cx.name)){stop("no Center X data")}
	else{if(cx.name != "Center.X"){warning(cx.name,"assumed to be Center X data")}}
	cy.name <- grep("Center.Y",all.names,value=T,ignore=T)
	if(is.na(cy.name)){stop("no Center Y data")}
	else{if(cy.name != "Center.Y"){warning(cy.name,"assumed to be Center Y data")}}
	
	x.names <- unique(tmp[,id.name])
	x.tab <- table(tmp[,id.name])
	if(max(x.tab) != min(x.tab)){error("all ids do not have the same number of data points")}
	x.row <- max(x.tab)
	t.dat <- matrix(tmp[,ratio.name],byrow=FALSE,nrow=x.row)
	time.val <- as.character(tmp[tmp[,id.name]==x.names[1],time.name])
	
	if(length(grep(":",time.val[1]))==0)
	{
		x <- as.single(time.val)
		if(max(x) > 1000000)#in ms
		{
			x <- x/60000
		}
		else if(max(x) > 1500) #in seconds
		{
			x <- x/60	
		}		
		time.val <- x #in minutes
	}
	else{time.val <- sapply(as.character(time.val),ConvertTime)}
	t.dat <- cbind(time.val,t.dat) #note assumption of ms
	t.dat <- as.data.frame(t.dat)
	names(t.dat) <- c("Time",paste("X.",x.names,sep=""))
	c.names <- c(area.name,cx.name,cy.name)
	o.names <- setdiff(all.names,c(time.name,id.name,area.name,ratio.name,cx.name,cy.name))
	if(length(o.names) > 0){warning(paste(o.names,"added to c.dat"));c.names <- c(c.names,o.names)}

	c.dat <- tmp[match(x.names,tmp[,id.name]),c.names]
	c.dat <- cbind(paste("X.",x.names,sep=""),c.dat)
	c.dat <- data.frame(c.dat)
	names(c.dat)[1:4] <- c("id","ROI.Area","Center.X","Center.Y") 
	row.names(c.dat) <- c.dat[,"id"]
	if(!is.null(wrdef))
	{
		wr <- ReadResponseWindowFile(wrdef)  #complete and revise this section
		w.dat <- MakeWr(t.dat,wr)
	}
	else
	{
		w.dat <- t.dat[,1:2]
		w.dat[2] <- ""
		names(w.dat)[2] <- "wr1"
		w.dat[1:10,"wr1"] <- rep("noDef",n=10)
	}
	if(is.null(rd.name)){rd.name <- paste("RD",make.names(date()),sep="")}
	tmp.rd <- list(t.dat=t.dat,w.dat=w.dat,c.dat=c.dat)
	f.name <- paste(rd.name,".Rdata",sep="")
	assign(rd.name,tmp.rd)
	save(list=rd.name,file=f.name)
	return(paste(nrow(tmp.rd$c.dat),"traces read saved to ",f.name))
	#save as RD file
}

#read only the ratio data from datadump
ReadRatioDump <- function(fname,wrdef=NULL,cd.name=NULL,rd.name=NULL,sep="\t",encode="UCS-2LE")
{
	tmp <- read.delim(fname,fileEncoding=encode,sep=sep,header=T)
	all.names <- names(tmp)
	time.name <- grep("Time",all.names,value=T,ignore=T)[1]
	if(time.name != "Time..ms."){warning(paste(time.name," assumed to be in ms"))}
	id.name <- grep("ROI.ID",all.names,value=T,ignore=T)[1]
	if(is.na(id.name)){stop("no id column detected")}
	if(id.name != "ROI.ID"){warning(paste(id.name,"assumed to be it ROI.ID"))}
	ratio.name <- grep("Ratio",all.names,value=T,ignore=T)
	if(is.na(ratio.name)){stop("no ratio data detected")}
	else{if(ratio.name != "Ratio.340.380"){warning(ratio.name," assumed to be Ratio data")}}

	x.names <- unique(tmp[,id.name])
	x.tab <- table(tmp[,id.name])
	if(max(x.tab) != min(x.tab)){stop("all ids do not have the same number of data points")}
	x.row <- max(x.tab)
	t.dat <- matrix(tmp[,ratio.name],byrow=FALSE,nrow=x.row)
	time.val <- as.character(tmp[tmp[,id.name]==x.names[1],time.name])	
	if(length(grep(":",time.val[1]))==0)
	{
		x <- as.single(time.val)
		if(max(x) > 1000000)#in ms
		{
			x <- x/60000
		}
		else if(max(x) > 1500) #in seconds
		{
			x <- x/60	
		}		
		time.val <- x #in minutes
	}
	else{time.val <- sapply(as.character(time.val),ConvertTime)}
	t.dat <- cbind(time.val,t.dat) #note assumption of ms
	t.dat <- as.data.frame(t.dat)
	t.dat <- t.dat[unique(row.names(t.dat)),]
	names(t.dat) <- c("Time",paste("X.",x.names,sep=""))
	#c.dat
	if(!is.null(cd.name))
	{
		c.dat <- ReadCellDetailsFile(cd.name)
		if(!setequal(row.names(c.dat),names(t.dat)[-1])){stop("cell details ids don't match trace data ids")}
	}
	else
	{
		c.dat <- data.frame(id=names(t.dat)[-1],ROI.Area=NA,Center.X=NA,Center.Y=NA)
		row.names(c.dat) <- c.dat[,"id"]
	}		
	
	if(!is.null(wrdef))
	{
		wr <- ReadResponseWindowFile(wrdef)  #complete and revise this section
		w.dat <- MakeWr(t.dat,wr)
	}
	else
	{
		w.dat <- t.dat[,1:2]
		w.dat[2] <- ""
		names(w.dat)[2] <- "wr1"
		w.dat[1:10,"wr1"] <- rep("noDef",n=10)
	}
	if(is.null(rd.name)){rd.name <- paste("RD",make.names(date()),sep="")}
	tmp.rd <- list(t.dat=t.dat,w.dat=w.dat,c.dat=c.dat)
	f.name <- paste(rd.name,".Rdata",sep="")
	assign(rd.name,tmp.rd)
	save(list=rd.name,file=f.name)
	return(paste(nrow(tmp.rd$c.dat),"traces read saved to ",f.name))
	#save as RD file
}


ReadResponseWindowFile <- function(fname)
{
    dat <- read.csv(fname)
    return(dat)
}

ReadCellDetailsFile <- function(fname)
{
	dat <- read.delim(fname,sep="\t",header=T,fileEncoding="UCS-2LE")
    #dat <- read.csv(fname,header=T)
    id.name <- grep("id",names(dat),ignore.case=T,value=F)
    if(is.na(id.name)){stop(paste("no id column in",fname))}
    id.name <- id.name[1]
    names(dat)[id.name] <- "id"
    dat["id"] <- paste("X.",dat[,"id"],sep="")
    dat$id <- as.character(dat$id)
    row.names(dat) <- dat$id
    centx <- grep("cent.*xpx$",names(dat),ignore.case=T,value=F)
    if(length(centx)==0)
    	{dat["Center.X"] <- NA}
    else
    	{names(dat)[centx[1]] <- "Center.X"}
    centy <- grep("cent.*ypx$",names(dat),ignore.case=T,value=F)
    if(length(centy)==0)
    	{dat["Center.Y"] <- NA}
    else
    	{names(dat)[centy[1]] <- "Center.Y"}
    roiarea <- grep("area",names(dat),ignore.case=T,value=F)
    if(length(roiarea)==0)
    	{dat["ROI.Area"] <- NA}
    else
    	{names(dat)[roiarea[1]] <- "ROI.Area"}
    #remove some extraneous columns
	c.names <- setdiff(names(dat),c("X","Item"))
	dat <- dat[,c.names]
    
    return(dat)
}


#check for id alignment
#make provisions for null windows or details
ReadExperiment <- function(t.name=NULL,w.name=NULL,c.name=NULL)
{
    t.dat <- NULL
    w.dat <- NULL
    c.dat <- NULL
    if(is.null(t.name))
    {
        stop("Must have a trace file")
    }
    else
    {
        t.dat <- ReadTraceFile(t.name)
    }
    if(!is.null(w.name))
    {
        w.dat <- ReadResponseWindowFile(w.name)
        if(nrow(w.dat) != nrow(t.dat))
        {stop("Number of rows for the window response file must match the number of rows from the trace file")}
    }
    if(!is.null(c.name))
    {
        c.dat <- ReadCellDetailsFile(c.name)
        g.names <- grep("X",names(t.dat),value=T)
        if(length(setdiff(g.names,c.dat$id)) > 0)
        {warning(c("All cells not characterized in Cell Details File ",setdiff(g.names,c.dat$id)))}
        if(length(setdiff(c.dat$id,g.names)) > 0)
        {warning(c("Cell details given for cells with no trace data ",setdiff(c.dat$id,g.names)))}
 
    }
    return(list(t.dat = t.dat, w.dat = w.dat, c.dat=c.dat))
}

WindowMax <- function(mval,hws=20,levs=c("1","2"))
{
    tmp.window <- rep("",length(mval))
    wc <- rep(seq(1,length(levs)),ceiling(length(mval)/length(levs)))
    wc <- sort(wc[1:length(mval)])

    for(i in 1:length(levs))
    {
        j1 <- match(i,wc)
        jmax <- 0
        me <- 0
        for(j in (j1+1+hws):(j1+sum(wc==i)-hws))
        {
            me <- mean(mval[(j-hws):(j+hws)],na.rm=T)
            if(me >= jmax)
            {
                jmax <- me
                jmaxj <- j
            }
        }
        lowi <- jmaxj-hws
        highi <- jmaxj+(hws-1)
        tmp.window[lowi:highi] <- levs[i]
    }
    return(tmp.window)
}


#wr file should be
#NEW format for wr file (three column table, treatment, at, duration)
GetWr <- function(fname)
{
	wr1 <- read.csv(fname)
	return(wr1)
}

MakeWr <- function(t.dat,wr1,padL=0,padR=0)
{
	w.dat <- t.dat[,1:2]
	names(w.dat)[2] <- "wr1"
	w.dat["wr1"] <- ""
	wr1["treatment"] <- make.names(wr1[,"treatment"],unique=T)
	for(i in 1:nrow(wr1))
	{
		x1 <- which.min(abs(wr1[i,"at"]-t.dat[,"Time"]))
		x2 <- which.min(abs((wr1[i,"at"]+wr1[i,"duration"])-t.dat[,"Time"]))
		w.dat[max((x1-padL),1):min((x2+padR),nrow(t.dat)),"wr1"] <- as.character(wr1[i,"treatment"])
	}
	return(w.dat)
}

#fill forward for flevs in the window region.
FillWR <- function(wr1,flevs)
{
	u.names <- unique(wr1)
	wr2 <- NumBlanks(wr1)

	u2.names <- unique(wr2)
	b.names <- grep("blank",u2.names,value=T)
	for(i in flevs)
	{
		for(j in 1:(length(u2.names)-1))
		{
			if(u2.names[j]==i & is.element(u2.names[j+1],b.names) )
			{
				wr1[wr2==u2.names[j+1]] <- i
			}
		}
	}
	return(wr1)
}

#adjust the windows to maximize shift regions and peak regions
#try to minimize the false positive rates but growing/shinking windows
#works reasonably well, but is only counting peaks.  It is not accountinf
#for shape aspects of the trace.
WrAdjust <- function(dat,pcp=NULL,wr=NULL,wr.levs=NULL,snrT=4,minT=10)
{
	gtrfunc <- function(x,a){sum(x>a)}
	if(is.null(wr)){wr <- dat$w.dat[,"wr1"]}
	wr.new <- wr
	wrb <- NumBlanks(wr)
	wi <- 1:(length(wrb))
	x.names <- names(dat$t.dat[,-1])
	if(is.element("bin",names(dat)))
		if(is.element("drop",names(dat$bin)))
		{
			x.names <- row.names(dat$bin[dat$bin[,"drop"]==0,])
		}
	if(is.null(wr.levs))
	{
		wr.levs <- unique(wr)
		wr.levs <- wr.levs[wr.levs != ""]
	}
	if(is.null(pcp))
	{
		pcp <- ProcConstPharm(dat)
	}
	#OK expand/contract each window to give best false positive ratio.
	#keep a min width.
	hits <- apply(pcp$snr[,x.names],1,gtrfunc,a=snrT)
	wrb.levs <- unique(wrb)
	b.levs <- grep("blank",wrb.levs,value=T)
	for(i in wr.levs[wr.levs != wrb.levs[length(wrb.levs)]])
	{
		i1 <- match(i,wrb.levs)
		if(is.element(wrb.levs[i1+1],b.levs))
		{
			targs <- hits[wrb==i | wrb==wrb.levs[i1+1]]
			tval <- NULL
			endT <- length(targs)
			lp <- 0
			for(j in minT:(endT-1))
			{
				lp <- lp+1
				#tval[lp] <- mean(targs[1:j])/((sum(targs[(j+1):endT])+1)/length(targs[(j+1):endT]))
				tval[lp] <- 1/((sum(targs[(j+1):endT])+1)/length(targs[(j+1):endT]))				
			}			
			iopt <- match(i,wr)+which.max(tval)+(minT-1)
		}
		else
		{iopt <- max(wi[wr==i])}
		wr.new[wr==i] <- ""
		wr.new[match(i,wr):iopt] <- i
	}	
	return(wr.new)
}

#given a directory read the files within and process all the experiments
ReadMultipleExperiments <- function(dir.name=NULL,padL=2,padR=20)
{
	if(is.null(dir.name)){stop("not a directory")}
	setwd(dir.name)
	f.names <- list.files(pattern="\\.csv$")
	td.names <- grep("^td",f.names,value=T)
	cd.names <- grep("^cd",f.names,value=T)
	if(!setequal(sub("^td","",td.names),sub("^cd","",cd.names)))
	{
		print(td.names)
		print(cd.names)
		stop("File Sets not equal")
	}
	#ah now what to name the datafile.
	RD.names <- make.names(sub("^td","RD",sub("\\.csv","",td.names)))
	RD.f.names <- paste(RD.names,".Rdata",sep="")
	wr.name <- grep("^wr",f.names,value=T)
	if(length(wr.name)==0){stop("no wr.file")}
	wr1 <- GetWr(wr.name) #NEW format for wr file (three column table, treatment, at, duration) time in minutes
	for(i in 1:length(td.names))
	{
		tmp <- ReadExperiment(td.names[i],,cd.names[i])
		if(max(tmp$t.dat[,1]) > 5000)
			{tmp$t.dat[,1] <- tmp$t.dat[,1]/60000}#probably in ms
		else if(max(tmp$t.dat[,1] > 500))
			{tmp$t.dat[,1] <- tmp$t.dat[,1]/60}#probably in s

		tmp$w.dat <- MakeWr(tmp$t.dat,wr1,padL,padR) #place the treatments
		assign(RD.names[i],tmp)
		save(list=RD.names[i],file=RD.f.names[i])	
	}
	return(RD.f.names)
}

#should probably break this into ScoreMulti and ReviewMulti
#Score all RD...Rdata files in a given directory with review
#check for an existing bin file and just review that.
#add a "drop" column to the bin file
ScoreMulti <- function(dir.name=NULL,snr.lim=4,hab.lim=.05,sm=3,ws=30,review=T,wr.i="wr1",bl.meth="TopHat",t.name="t.dat")
{
	if(is.null(dir.name)){dir.name <- getwd()}
	setwd(dir.name)
	f.names <- list.files(pattern="RD.*\\.Rdata$")
	if(length(f.names) == 0){stop("no RD...Rdata files in given directory")}
	rd.list <- sub("\\.Rdata*","",f.names)
	RD.names <- rd.list #paste(rd.list,".b",sep="")
	RD.f.names <- paste(RD.names,".Rdata",sep="")
	sel.i <- menu(rd.list,title="Select Data to review")			
	while(sel.i != 0)
	{

		j <- sel.i
		load(f.names[j])
		i <- rd.list[j]
		tmp <- get(i)
		tlevs <- c(as.character(unique(tmp$w.dat[,wr.i])),"drop")
		tlevs <- tlevs[!is.element(tlevs,c("tot","sd",""))]
		tmp$t.dat <- tmp[[t.name]]
		if(is.null(tmp$bin))
		{
		tmp.pcp <- ProcConstPharm(tmp,sm,ws,bl.meth)
		tmp.scp <- ScoreConstPharm(tmp$t.dat,tmp.pcp$blc,tmp.pcp$snr,snr.lim,hab.lim,tmp$w.dat[,wr.i],sm)
		tmp.bin <- bScore(tmp.pcp$blc,tmp.pcp$snr,snr.lim,hab.lim,tlevs,tmp$w.dat[,wr.i])
		tmp.bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
		}
		else
		{
			tmp.bin <- tmp$bin
			tmp.scp <- tmp$scp
			#tmp.blc <- tmp$blc
		}
		if(review)
		{
		tmp.bin <- ScoreReview1(tmp$t.dat,tmp.bin[,tlevs],tmp$w.dat[,wr.i])
		tmp.bin <- ScoreReview0(tmp$t.dat,tmp.bin[,tlevs],tmp$w.dat[,wr.i])
		}
		tmp$bin <- tmp.bin[,tlevs]
		tmp$scp <- tmp.scp
		#tmp$blc <- tmp.blc
		
		assign(RD.names[j],tmp)
		save(list=RD.names[j],file=RD.f.names[j])
		print(paste("DONE REVIEWING ",RD.names[j]," CHANGES SAVED TO FILE.",sep=""))
		sel.i <- menu(rd.list,title="Select Data to review")			
	}
	return(RD.f.names)		
}

#this is not complete
#condi is the indicator for the conditional frequency table
#this is bad
#####add selection section of selection of experiments to include/exclude
#####conditional expresion tables.
SummarizeMulti <- function(dir.name=NULL,condi=1,recur=F,f.names=NULL,rd.list=NULL)
{
	if(is.null(dir.name)){stop("not a directory")}
	setwd(dir.name)
	if(is.null(f.names) & is.null(rd.list))
	{
	f.names <- list.files(pattern=".*RD.*\\.Rdata$",recursive=recur,full.names=T)
	f.names <- select.list(f.names,multiple=T,title="Select Experiments For Analysis")
	if(length(f.names) == 0){stop("no RD...Rdata files in given directory")}
	}
	if(is.null(rd.list))
	{
	for(i in f.names){load(i)}
	rd.list <- sub("\\.Rdata*","",basename(f.names))
	RD.names <- ls(pat="^RD")
	RD.names <- intersect(rd.list,RD.names)
	if(!setequal(RD.names,rd.list)){stop("dataframes loaded do not match files listed in directory")}
	RD.f.names <- paste(RD.names,".Rdata",sep="")
	}
	i <- rd.list[1]
	tmp <- get(i)
	if(sum(is.element(c("bin","scp"),names(tmp))) < 2){stop("Data frame has not been scored")}

	if(names(tmp$bin)[c(1,2)]==c("tot","sd"))
	{tmp$bin <- tmp$bin[,-c(1,2)]}
	freq.tab <- data.frame(mean=apply(tmp$bin[tmp$bin[,"drop"]==0,],2,mean))
	kfreq.tab <- data.frame(mean=apply(tmp$bin[tmp$bin[,"drop"]==0 & tmp$bin[,condi]==1,],2,mean))

	b.names <- row.names(freq.tab)[row.names(freq.tab) != "drop"]
	q.names <- paste(b.names,".max",sep="")
	resp.tab <- data.frame(mean=apply(tmp$scp[tmp$bin[,"drop"]==0,q.names],2,mean))
	for(rn in row.names(resp.tab)){resp.tab[rn,"mean"] <- mean(tmp$scp[tmp$bin[,"drop"]==0 & tmp$bin[,sub("\\.max$","",rn)]==1,rn],na.rm=T)}
	pf.tot <- data.frame(str = apply(tmp$bin[tmp$bin[,"drop"]==0,names(tmp$bin)!="drop"],1,paste,collapse=""))
	pf.tot["exp"] <- i
	for(j in 2:length(rd.list))
	{
		i <- rd.list[j]
		tmp <- get(i)
		if(names(tmp$bin)[c(1,2)]==c("tot","sd"))
		{tmp$bin <- tmp$bin[,-c(1,2)]}
		
		m1 <- apply(tmp$bin[tmp$bin[,"drop"]==0,],2,mean)
		freq.tab[i] <- m1[row.names(freq.tab)]
		m2 <- apply(tmp$bin[tmp$bin[,"drop"]==0 & tmp$bin[,condi]==1,],2,mean)
		kfreq.tab[i] <- m2[row.names(kfreq.tab)]
		resp.tab[i] <- NA
		for(rn in intersect(row.names(resp.tab),names(tmp$scp))){resp.tab[rn,i] <- mean(tmp$scp[tmp$bin[,"drop"]==0 & tmp$bin[,sub("\\.max$","",rn)]==1,rn],na.rm=T)}
		pf.tmp <- data.frame(str = apply(tmp$bin[tmp$bin[,"drop"]==0,names(tmp$bin)!="drop"],1,paste,collapse=""))		
		pf.tmp["exp"] <- i
		pf.tot <- rbind(pf.tot,pf.tmp)
}
	names(freq.tab)[1] <- rd.list[1]
	names(kfreq.tab)[1] <- rd.list[1]
	names(resp.tab)[1] <- rd.list[1]
	pf.tab <- table(pf.tot[,1],pf.tot[,2])
	return(list(freq.tab=freq.tab,kfreq.tab=kfreq.tab,resp.tab=resp.tab,pf.tab=pf.tab))
}


#RoughReview
#the first argument is an experiment object with t.dat, w.dat and c.dat
#note that only c.dat is not really necessary at this point.
#the second argument is the smoothing half window size (defaults to 2)
#the third argument is the peak detection half window size (defaults to 20)
#the fourth argument is the column from the w.dat to use as the window response regions.  This can be the name of the column or the number
#the fifth argument is the number of rough groups to divide the data into. (defaults to 20).  Generall weird traces will end up in groups by themselves.

#this function opens two windows: a bulk group window that shows all raw traces for cells in that group.
#and a single cell window that shows the raw data and the processed data.
#this function allows you to step through the groups, the cells in the group and to change the peak detection
#parameters.  The current values for the parameters is shown below the single cell graph.

RoughReview <- function(dat,shws=2,phws=20,wr.i=2,pam.k=20)
{
    library(cluster)
    graphics.off()
    dev.new(width=14,height=4)    
    dev.new(width=14,height=4)  
    devs <- dev.list()
    if(length(devs) < 2){stop("windows will not open")}
	devs <- devs[1:2]
    linesfunc <- function(n.names)
    {
        dev.set(devs[2])
        ylim <- c(-.1,1.4)        
        plot(dat$t.dat[,"Time"],dat$t.dat[,n.names[1]],type="n",ylim=ylim,xlab="Time (min)",ylab="Ratio",xaxt="n")
		axis(1, at=seq(0, length(dat$t.dat[,1]), 5))  
        for(i in n.names)
        {lines(dat$t.dat[,"Time"],dat$t.dat[,i])}
        dev.set(devs[1])
    }
    lines.flag <- 1
#    pam.k <- min(as.integer(ncol(dat$t.dat)/10),30)
    g.names <- names(dat$t.dat[,-1])
    dat.pam <- pam(t(dat$t.dat[,g.names]),k=pam.k)
    group.i <- 1
    cell.i <- 1
    click.i <- 1
    basel <- "TopHat"

    while(click.i != 10)
    {
        g.num <- sum(dat.pam$clustering==group.i)
        group.names <- g.names[dat.pam$clustering==group.i]
        cell.pick <- g.names[dat.pam$clustering==group.i][cell.i]
        dev.set(devs[1])
        p1 <- PeakFunc2(dat$t.dat,cell.pick,shws=shws,phws=phws,Plotit=T,wr=dat$w.dat[,wr.i],SNR.lim=2,bl.meth=basel)
        
        stext <- c("Group+","Group-","Cell+","Cell-","SHWS+","SHWS-","PHWS+","PHWS-","BaseL","Done")
		xs <- rep(dat$t.dat[1,"Time"],length(stext))
		ys <- seq(max(dat$t.dat[,cell.pick]-(min(dat$t.dat[,cell.pick])-.1)),-.1,length=(length(stext)))

        if(lines.flag==1){linesfunc(group.names);lines.flag <- 0}
        title(sub=paste("Group ",group.i," n=",g.num," Cell ",cell.i," SHWS ",shws," PHWS ",phws," BaseL ",basel,sep=""))
        text(x=xs,y=ys,labels=stext,pos=2,cex=.5)
        points(xs,ys,pch=16)
        click.i <- identify(x=xs,y=ys,n=1,plot=F)
        if(click.i==1)
        {group.i <- group.i + 1;if(group.i > pam.k){group.i <- 1};cell.i<-1;lines.flag <- 1}
        if(click.i==2)
        {group.i <- group.i - 1;if(group.i < 1){group.i <- pam.k};cell.i<-1;lines.flag <- 1}
        if(click.i==3)
        {cell.i <- cell.i + 1;if(cell.i > g.num){cell.i <- 1}}
        if(click.i==4)
        {cell.i <- cell.i - 1;if(cell.i < 1){cell.i <- g.num}}
        if(click.i==6)
        {shws <- max(shws-1,0)}
        if(click.i==5)
        {shws <- shws+1}
        if(click.i==8)
        {phws <- max(phws-1,2)}
        if(click.i==7)
        {phws <- phws+1}
        if(click.i==9)
        {basel <- setdiff(c("TopHat","SNIP"),basel)}

    }
}

#the first argument is the raw data
#the second argument is the halfwindow size for smoothing (shws)
#the third argument is the peak detection halfwindow size (phws)
#the last argument is the baseline correction method (TopHat = blue line SNIP = red line)
#Note that you should use the RoughReview function to determine the best values for
#arguments 2,3 and 4.

#returns a list with two dataframes: snr and blc.
#snr has the peaks detected for all cells, blc has the baseline corrected data for all cells. 

ProcConstPharm <- function(dat,shws=2,phws=20,bl.meth="TopHat")
{
    dat1 <- dat$t.dat
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
#    dat1.crr <- allCRR(dat1,t.names,Plotit=F) #leave off advanced processing for now
    return(list(snr=dat1.snr,blc=dat1.bc))
}

#binary score for all cells for the regions of interest bScore
#argument 1 is the baseline corrected data
#argument 2 is the snr peak data
#argument 3 is the threshold for significance on the peaks
#argument 4 is the intensity above baseline theshold
#argument 5 indicates the regions of interest. (e.g. the response windows for which the cells will be scored)
#argument 6 indicates the response windows. 
#argument 7 indicates the cells to score (if null all cells will be scored)
#returns the scoring for all cells subject to the above parameters.
#as well as the sum for the snr scores and the sd for the snr scores.
bScore <- function(blc,snr,snr.lim,blc.lim,levs,wr,c.names=NULL,delta.lim=5)
{
    notzero <- function(x){as.integer(sum(x) > 0)}
    if(is.null(c.names)){c.names <- names(blc)[-1]}
    wt <- blc
	wtd <- wt[-1,]-wt[-nrow(wt),]
	wtd <- sweep(wtd[,-1],1,wtd[,1],'/')
	wtd <- rbind(wtd,rep(0,ncol(wtd)))
    wr2 <- wr[is.element(wr,levs)]
    b.snr <- snr[is.element(wr,levs),c.names]
    b.blc <- blc[is.element(wr,levs),c.names]
    wtd <- wtd[is.element(wr,levs),c.names]
    b.call <- b.blc
    b.call[,] <- 0
    b.call[(b.snr > snr.lim & b.blc > blc.lim)] <- 1
    b.score <- data.frame(tot=apply(b.snr,2,sum))
    b.score["sd"] <- apply(b.snr,2,sd)
    for(i in levs)
    {
        b.score[i] <- apply(b.call[wr2==i,],2,notzero)
		x.names <- row.names(b.score)[b.score[i]==1]
		if(length(x.names) > 5)
		{
			p.pos <- median(apply(wtd[wr2==i,x.names],2,which.max))
			y.names <- setdiff(row.names(b.score),x.names)
			bplus <- apply(wtd[wr2==i,y.names],2,max)
			bpos <- apply(wtd[wr2==i,y.names],2,which.max)
			bdif <- abs(bpos-p.pos)
			z.names <- y.names[bdif < 2 & bplus > delta.lim]
			if(length(z.names) > 0){b.score[z.names,i] <- 1}
		}
    }
    return(b.score)
}

# calculate a table of cell characteristics globally and 
# within specific windows
# these specifics should include
# mean and sd, sum of in window peaks, sum of out of window peaks
# 1)  some measure of dead cell
# 2)  yes/no peak response for each window
# 3) peak height
# 4) max peak SNR
# 5) peak timing in window
# 6)
# variance of smoothed - raw in window
# define and number blank windows.
ScoreConstPharm <- function(raw,blc,snr,snr.lim,blc.lim,wr,shws=2)
{
    gtfunc <- function(x,alph){sum(x > alph,na.rm=T)}
    
lt5func <- function(x,y)
{
    ltfunc <- function(i){summary(lm(y[i:(i+5)] ~ x[i:(i+5)]))$coefficients[2,3]}
    iseq <- 1:(length(x)-5)
    res <- sapply(iseq,ltfunc)
    return(range(res))
}

    levs <- setdiff(unique(wr),"")
    c.names <- names(raw)[-1]
    res.tab <- data.frame(mean=apply(blc[,c.names],2,mean))
    res.tab["sd"] <- apply(blc[,c.names],2,sd)
    res.tab["snr.iws"] <- apply(snr[is.element(wr,levs),c.names],2,sum)
    res.tab["snr.ows"] <- apply(snr[!is.element(wr,levs),c.names],2,sum)
    res.tab["snr.iwc"] <- apply(snr[is.element(wr,levs),c.names],2,gtfunc,alph=snr.lim)
    res.tab["snr.owc"] <- apply(snr[!is.element(wr,levs),c.names],2,gtfunc,alph=snr.lim)

    for(i in c.names)
    {
        s1 <- createMassSpectrum(raw[,"Time"],raw[,i])
        s3 <- smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=shws)
        bl.th <- estimateBaseline(s3, method="TopHat")[,"intensity"]
        bl.snp <- estimateBaseline(s3, method="SNIP")[,"intensity"]
        eseq <- 1:ceiling((nrow(raw)/2))
        lseq <- max(eseq):nrow(raw)
        res.tab[i,"bl.diff"] <- mean(bl.th-bl.snp)
        res.tab[i,"earl.bl.diff"] <- mean(bl.th[eseq]-bl.snp[eseq])
        res.tab[i,"late.bl.diff"] <- mean(bl.th[lseq]-bl.snp[lseq])        
    }
    for(i in levs)
    {
        res.tab[paste(i,".snr",sep="")] <- apply(snr[wr==i,c.names],2,max)
#peak intensity        res.tab[paste(i,".snr",sep="")] <- apply(snr[wr==i,c.names],2,max)        
        res.tab[paste(i,".tot",sep="")] <- apply(blc[wr==i,c.names],2,sum)
        res.tab[paste(i,".max",sep="")] <- apply(blc[wr==i,c.names],2,max)
        res.tab[paste(i,".wm",sep="")] <- apply(blc[wr==i,c.names],2,which.max)
#        res.tab[c(paste(i,".dn5",sep=""),paste(i,".up5",sep=""))] <- t(apply(raw[wr==i,c.names],2,lt5func,x=raw[wr==i,1]))
#        res.tab[paste(i,".dn5",sep="")] <- apply(blc[wr==i,c.names],2,dn5func)                
    }
    return(res.tab)
}

#topdown parsing of all traces
TraceChase <- function(dat,blc=NULL,levs=NULL,x.names=NULL,scale=T)
{
	library(cluster)
	if(is.null(blc)){
	if(is.element("blc",names(dat))){blc <- dat$blc}
	else
	{tmp.pcp <- ProcConstPharm(dat);blc <- tmp.pcp$blc}}
	if(is.null(levs))
	{
		levs <- unique(dat$w.dat[,"wr1"])
		levs <- select.list(levs,multiple=T,title="Select Regions for clustering")
	}
	dmat <- t(scale(blc[is.element(dat$w.dat[,"wr1"],levs),-1],scale=scale,center=scale))
	a.names <- names(blc)[-1]
	if(!is.null(x.names)){a.names <- intersect(x.names,names(blc))}
	done=FALSE
	while(!done)
	{
		if(length(a.names) < 21)
		{
			x.names <- TraceSelect(dat$t.dat,,a.names,dat$w.dat[,"wr1"],levs,"Final Select")
			done=TRUE
		}
		else
		{
			#pam20 <- pam(dmat[a.names,],k=20)
			clmb20 <- ClimbTree(dmat[a.names,],k=20)
			lmain <- paste("Select Traces (all or none to end) n=",length(a.names))
			#x.names <- SmashSelect(blc[c("Time",a.names)],pam20$clustering,row.names(pam20$medoids),dat$w.dat[,"wr1"],levs,lmain=lmain)				
			x.names <- SmashSelect(blc[c("Time",a.names)],clmb20,names(clmb20)[match(1:length(unique(clmb20)),clmb20)],dat$w.dat[,"wr1"],levs,lmain=lmain)							
			if(length(a.names)==length(x.names)){done = TRUE}
			if(length(x.names)==0){done= TRUE}
			a.names <- x.names
		}
	
	}
	return(x.names)	
}

#given a set of traces (or trace seqments)
#calculate the distances and group into K groups
#by height of tree cutting. One of the K groups will
#be a catch-all for all small groups
ClimbTree <- function(x,k=20)
{
	tabstat <- function(x){return(list(mean=mean(x),length=length(x),median=median(x),sd=sd(x),gt5c=sum(x>5)))}
	library(cluster)
	d1 <- dist(x)
	h1 <- hclust(d1)
	q1 <- quantile(h1$height,probs=1:10/10)
	clust <- cutree(h1,h=q1[5])
	clust.tab <- table(clust)
	clust.tab <- clust.tab[order(clust.tab,decreasing=T)]	
	new.num <- clust.tab
	new.num[] <- seq(1,length(new.num))
	clust[] <- new.num[as.character(clust)]
	clust.tab <- table(clust)
	if(length(clust.tab) > k)
	{
		in.grp <- names(clust.tab[1:(k-1)])
		out.grp <- setdiff(names(clust.tab),in.grp)
		clust[is.element(clust,out.grp)] <- k	
	}
	return(clust)
	# clust.stat <- data.frame(tabstat(clust.tab))
	# for(i in 2:length(q1))
	# {
		# clust <- cutree(h1,h=q1[i])
		# clust.tab <- table(clust)
		# clust.stat[i,] <- tabstat(clust.tab)
	# }
	# return(clust.stat)
	
}

GroupReview <- function(dat,shws=2,phws=20,wr.i=4,clust,bl.meth="TopHat")
{
    library(cluster)
    graphics.off()
    dev.new(width=14,height=4)    
    dev.new(width=14,height=4)    
    linesfunc <- function(n.names)
    {
        dev.set(dev.list()[2])
        ylim <- c(-.1,1.4)
        plot(dat$t.dat[,"Time"],dat$t.dat[,n.names[1]],type="n",ylim=ylim,xlab="Time (min)",ylab="Ratio", xaxt="n")
        axis(1, at=seq(0, length(dat$t.dat[,1]), 5))        
        for(i in n.names)
        {lines(dat$t.dat[,"Time"],dat$t.dat[,i])}
        dev.set(dev.list()[1])
    }
    clust.name <- unique(clust)
    lines.flag <- 1
    g.names <- names(dat$t.dat[,-1])
    pam.k <- max(clust)
    group.i <- 1
    cell.i <- 1
    click.i <- 1
    while(click.i != 5)
    {
        g.num <- sum(clust==group.i)
        if(g.num > 0)
        {
        group.names <- g.names[clust==group.i]
        cell.pick <- g.names[clust==group.i][cell.i]
        dev.set(2)
        p1 <- PeakFunc2(dat$t.dat,cell.pick,shws=shws,phws=phws,Plotit=T,wr=dat$w.dat[,wr.i],SNR.lim=2,bl.meth=bl.meth)
        if(lines.flag==1){linesfunc(group.names);lines.flag <- 0}
        title(sub=paste("Group ",group.i," n=",g.num," Cell ",cell.i,sep=""))
        xs <- rep(dat$t.dat[50,"Time"],5)
        points(x=xs,y=c(1.2,1.1,1.0,.9,.8),pch=16)
        text(x=xs,y=c(1.2,1.1,1.0,.9,.8),labels=c("Group +","Group -","Cell +","Cell -","Done"),pos=2,cex=.5)
        click.i <- identify(x=xs,y=c(1.2,1.1,1.0,.9,.8),n=1,plot=F)
        }
        if(click.i==1)
        {group.i <- group.i + 1;if(group.i > pam.k){group.i <- 1};cell.i<-1;lines.flag <- 1}
        if(click.i==2)
        {group.i <- group.i - 1;if(group.i < 1){group.i <- pam.k};cell.i<-1;lines.flag <- 1}
        if(click.i==3)
        {cell.i <- cell.i + 1;if(cell.i > g.num){cell.i <- 1}}
        if(click.i==4)
        {cell.i <- cell.i - 1;if(cell.i < 1){cell.i <- g.num}}
    }
}

#TraceSelectLarge takes a large list of traces
#subsets it and passes each on to Trace select
TraceSelectLarge <- function(t.dat,snr=NULL,m.names,wr,levs=NULL,lmain="",subset.n=10,rtag=NULL)
{
	sel.names <- NULL
	s <- ceiling(length(m.names)/subset.n)
	for(i in 1:s)
	{
		x1 <- (i-1)*subset.n+1
		x2 <- min(length(m.names),x1+subset.n)
		x.sel <- TraceSelect(t.dat,snr,m.names[x1:x2],wr,levs,lmain,rtag[x1:x2])
		sel.names <- union(sel.names,x.sel)
	}		
	return(sel.names)	
}

TraceSelect <- function(t.dat,snr=NULL,m.names,wr,levs=NULL,lmain="",rtag=NULL)
{
	if(!is.null(rtag)){names(rtag) <- m.names}
	sf <- max(t.dat[,m.names])*.8
    library(RColorBrewer)
    m.names <- intersect(m.names,names(t.dat))
    lwds <- 3
    if(length(m.names) == 0)
    {stop("no named traces exist in trace dataframe.")}
    
    xseq <- t.dat[,1]
    cols <-brewer.pal(8,"Dark2")
    cols <- rep(cols,ceiling(length(m.names)/length(cols)))
    cols <- cols[1:length(m.names)]
    dev.new(width=14,height=8)
    m.pca <- prcomp(t(t.dat[,m.names]),scale=F,center=T)
    m.names <- m.names[order(m.pca$x[,1],decreasing=sum(m.pca$rot[,1]) < 0)]
    hbc <- length(m.names)*sf+max(t.dat[,m.names])
    hb <- ceiling(hbc)
    
    plot(xseq,t.dat[,m.names[1]],ylim=c(-sf,hbc),xlab="Time (min)",ylab="Ratio with shift",main=lmain,type="n", xaxt="n")
	axis(1, at=seq(0, length(t.dat[,1]), 5))

    if(length(wr) > 0)
    {
    	if(is.null(levs)){levs <- setdiff(unique(wr),"")}
        x1s <- tapply(xseq,as.factor(wr),min)[levs]
        x2s <- tapply(xseq,as.factor(wr),max)[levs]
        y1s <- rep(-.3,length(x1s))
        y2s <- rep(hbc+.2,length(x1s))
        rect(x1s,y1s,x2s,y2s,col="lightgrey",border=NA)
        text(xseq[match(levs,wr)],rep(c(sf/2,-sf/2),length(levs)),levs,pos=4,offset=0,cex=1)
    }
    x.sel <- NULL
    xs <-c(rep(0,length(m.names)),c(.1,.1,.1))
    ys <- seq(1,length(m.names))*sf+t.dat[1,m.names]
    ys <- as.vector(c(ys,c(sf,0,-sf)))
#    xs[(length(xs)-2):length(xs)] <- c(0,5,10)
    p.names <- c(m.names,"ALL","NONE","FINISH")
    done.n <- length(p.names)
    none.i <- done.n-1
    all.i <- none.i-1
    p.cols <- c(cols,c("black","black","black"))
    for(i in 1:length(m.names))
    {
        lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lwds)
        if(!is.null(snr))
        {
        pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
        pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
        points(xseq[pp1],t.dat[pp1,m.names[i]]+i*sf,pch=1,col=cols[i])
        points(xseq[pp2],t.dat[pp2,m.names[i]]+i*sf,pch=0,col=cols[i])
        }
    }
	text(x=xs,y=ys,labels=p.names,pos=2,cex=.7,col=p.cols)
    points(x=xs,y=ys,pch=16,col=p.cols)
    if(!is.null(rtag))
        {
        	rtag <- rtag[m.names]
	        text(rep(xseq[length(xseq)-10],length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag,cex=.8,col=cols,pos=4)
        }
        
    click.i <- 1    
    while(click.i != done.n)
    {
        click.i <- identify(xs,ys,n=1,plot=F)
        if(click.i < (length(m.names)+1) & click.i > 0)
        {
            i <- click.i
            if(is.element(i,x.sel))
            {
                lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lwds)
                x.sel <- setdiff(x.sel,i)
            }
                else
                {
	    	    lines(xseq,t.dat[,m.names[i]]+i*sf,col="black",lwd=lwds)
                #lines(xseq,t.dat[,m.names[i]]+i*sf,col="white",lwd=2,lty=2)
                x.sel <- union(x.sel,i)
            }
        }
        if(click.i == none.i)
        {
        	x.sel <- NULL
	    	for(i in 1:length(m.names))
		    {
    		    lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lwds)
	    	}
	    }
        if(click.i == all.i)	
        {
        	x.sel <- seq(1,length(m.names))
	    	for(i in 1:length(m.names))
		    {
    		    lines(xseq,t.dat[,m.names[i]]+i*sf,col="black",lwd=lwds)
	    	}
        	
        }
    }
    dev.off()
    return(m.names[x.sel])
}

#smash select plot the smashes and return the selected.
#all data in t.dat is ploted (1st col must be time)
#m.names are taken to be the medoids of the clusters
SmashSelect <- function(t.dat,clust,m.names,wr,levs=NULL,lmain="")
{	
	rtag <- table(clust)
	names(rtag) <- m.names[order(clust[m.names])]
	sf <- 1
	gcol <- rgb(10,10,10,alpha=120,max=255)
	#gcol <- "grey"
	x <- t.dat[,-1]
	xm <- apply(x,2,max)
	xn <- scale(x,center=F,scale=xm)
	for(i in 1:nrow(xn)){xn[i,] <- xn[i,]+clust}
	
    library(RColorBrewer)
    lwds <- 2
    
    xseq <- t.dat[,1]
    cols <-brewer.pal(8,"Dark2")
    cols <- rep(cols,ceiling(length(m.names)/length(cols)))
    cols <- cols[1:length(m.names)]
    dev.new(width=14,height=9)
    op <- par(yaxt="n",bty="n",mar=c(4,0,2,1),cex=1)
    plot(xseq,xn[,m.names[1]],ylim=c((min(xn)-2),max(xn)),xlab="Time (min)",ylab="Ratio with shift",main=lmain,type="n", xaxt="n")
	axis(1, at=seq(0, length(t.dat[,1]), 5))
	apply(xn,2,lines,x=xseq,col=gcol,lwd=2)
	hbc <- 1
    if(length(wr) > 0)
    {
    	if(is.null(levs)){levs <- setdiff(unique(wr),"")}
        x1s <- tapply(xseq,as.factor(wr),min)[levs]
        x2s <- tapply(xseq,as.factor(wr),max)[levs]
        y1s <- rep(-.3,length(x1s))
        y2s <- rep(hbc+.2,length(x1s))
        #rect(x1s,y1s,x2s,y2s,col="lightgrey")
        text(xseq[match(levs,wr)],rep(c(.2,-.2),length.out=length(levs)),levs,pos=4,offset=0,cex=1)
    }
    x.sel <- NULL
    xs <-c(rep(0,length(m.names)),c(.1,.1,.1))
    ys <- xn[1,m.names]
    ys <- as.vector(c(ys,c(sf*.9,0,-sf*.9)))
#    xs[(length(xs)-2):length(xs)] <- c(0,5,10)
    p.names <- c(rep(" ",length(m.names)),"ALL","NONE","FINISH")
    done.n <- length(p.names)
    none.i <- done.n-1
    all.i <- none.i-1
    p.cols <- c(cols,c("black","black","black"))
    for(i in 1:length(m.names))
    {
  	    #lines(xseq,xn[,m.names[i]],col="black",lwd=lwds*.5)
        lines(xseq,xn[,m.names[i]],col=cols[i],lwd=lwds)
    }
    text(x=rep(max(xseq),length(m.names)),y=xn[nrow(xn),m.names],cex=.9,rtag,pos=4,col=p.cols)
	text(x=xs,y=ys,labels=p.names,pos=2,cex=.7,col=p.cols)
    points(x=xs,y=ys,pch=16,col=p.cols,cex=1.5)
    click.i <- 1    
    while(click.i != done.n)
    {
        click.i <- identify(xs,ys,n=1,plot=F)
        if(click.i < (length(m.names)+1) & click.i > 0)
        {
            i <- click.i
            if(is.element(i,x.sel))
            {
                lines(xseq,xn[,m.names[i]],col=cols[i],lwd=lwds)
                x.sel <- setdiff(x.sel,i)
            }
                else
                {
	    	    lines(xseq,xn[,m.names[i]],col="black",lwd=lwds)
                x.sel <- union(x.sel,i)
            }
        }
        if(click.i == none.i)
        {
        	x.sel <- NULL
	    	for(i in 1:length(m.names))
		    {
    		    lines(xseq,xn[,m.names[i]],col=cols[i],lwd=lwds)
	    	}
	    }
        if(click.i == all.i)	
        {
        	x.sel <- seq(1,length(m.names))
	    	for(i in 1:length(m.names))
		    {
    		    lines(xseq,xn[,m.names[i]],col="black",lwd=lwds)
	    	}
        	
        }
    }
    c.sel <- clust[m.names[x.sel]]
    x.ret <- names(clust[is.element(clust,c.sel)])
    dev.off()
    return(x.ret)
}

#this simply finds the traces in t.dat that are similar to targs
#note this is "complete" similarity other options may be
#"average" and "best"
GetCloser <- function(t.dat,targs,k=20)
{
	x.names <- setdiff(names(t.dat),targs)
	ct <- cor(t.dat[,x.names],t.dat[,targs])
	x.max <- apply(ct,1,min)
	y.names <- x.names[order(x.max,decreasing=T)[1:k]]
	return(y.names)	
}

#this is a bit raw still
#Given a set of traces (t.dat) and a list of targets (targs)
#identify the 20 most similar traces using wr and the select levs.
#allow the user to select from those to add to the master list.
SimilarSelect <- function(t.dat,targs,wr,levs=NULL,sel.cnt=20,lmain="Similar Select")
{
	
	plot(t.dat[,1],t.dat[,targs[1]],type="n",ylim=c(min(t.dat[-1]),(length(targs)+50)*.2),xlab="Time",ylab="Ratio")
	sf <- 0
	for(i in targs){lines(t.dat[,1],t.dat[,i]+sf);sf<-sf+.2}
	a.names <- setdiff(names(t.dat)[-1],targs)
	rjct <- rep(0,length(a.names))
	names(rjct) <- a.names
	done=FALSE
	tps <- seq(1:nrow(t.dat))
	if(!is.null(levs)){tps <- tps[is.element(wr,levs)]}
	while(!done)
	{	
		if(sum(rjct==0) < 21)
		{done=TRUE}
		else
		{
			x.names <- GetCloser(t.dat[tps,c(a.names[rjct==0],targs)],targs,k=sel.cnt)
			rjct[x.names] <- 1
			y.names <- TraceSelect(t.dat,,x.names,wr,lmain=lmain)
	
			if(length(y.names)==0){done=TRUE}
			else
			{
				targs <- c(targs,y.names)
				for(i in y.names){lines(t.dat[,1],t.dat[,i]+sf);sf<-sf+.2}
			}
		}
	}
	return(targs)
	#plot targs and allow user to
	#paint region of interest if you can do this it makes a very good window adjust function.
	#find matches within t.dat
	#show matches in trace select allow user to choose.
	#merge all selected and return that list.
	
	
}

ScoreSelect <- function(t.dat,snr=NULL,m.names,wr,levs=NULL,lmain="")
{
	sf <- .2
    library(RColorBrewer)
    m.names <- intersect(m.names,names(t.dat))
    lwds <- 3
    if(length(m.names) == 0)
    {stop("no named traces exist in trace dataframe.")}
    
    xseq <- t.dat[,1]
    cols <-brewer.pal(8,"Dark2")
    cols <- rep(cols,ceiling(length(m.names)/length(cols)))
    cols <- cols[1:length(m.names)]
    dev.new(width=14,height=8)
    m.pca <- prcomp(t(t.dat[,m.names]),scale=F,center=T)
    m.names <- m.names[order(m.pca$x[,1],decreasing=sum(m.pca$rot[,1]) < 0)]
    hbc <- length(m.names)*sf+max(t.dat[,m.names])
    hb <- ceiling(hbc)
    
    plot(xseq,t.dat[,m.names[1]],ylim=c(-sf,hbc),xlab="Time (min)",ylab="Ratio with shift",main=lmain,type="n", xaxt="n")
	axis(1, at=seq(0, length(t.dat[,1]), 5))

    if(length(wr) > 0)
    {
    	if(is.null(levs)){levs <- setdiff(unique(wr),"")}
        x1s <- tapply(xseq,as.factor(wr),min)[levs]
        x2s <- tapply(xseq,as.factor(wr),max)[levs]
        y1s <- rep(-.3,length(x1s))
        y2s <- rep(hbc+.2,length(x1s))
        rect(x1s,y1s,x2s,y2s,col="lightgrey")
        text(xseq[match(levs,wr)],rep(-.1,length(levs)),levs,pos=4,offset=0,cex=1)
    }
    x.sel <- NULL
    xs <-rep(0,(length(m.names)+4))
    ys <- seq(1,length(m.names))*sf+t.dat[1,m.names]
    ys <- as.vector(c(ys,c(2*sf,sf,0,-sf)))
#    xs[(length(xs)-2):length(xs)] <- c(0,5,10)
    p.names <- c(m.names,"ALL","NONE","FINISH","DROP")
	drop.i <- length(p.names)
    done.n <- drop.i-1
    none.i <- drop.i-2
    all.i <- drop.i-3

    p.cols <- c(cols,c("black","black","black","black"))
    for(i in 1:length(m.names))
    {
        lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lwds)
        if(!is.null(snr))
        {
        pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
        pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
        points(xseq[pp1],t.dat[pp1,m.names[i]]+i*sf,pch=1,col=cols[i])
        points(xseq[pp2],t.dat[pp2,m.names[i]]+i*sf,pch=0,col=cols[i])
        }
    }
	text(x=xs,y=ys,labels=p.names,pos=2,cex=.7,col=p.cols)
    points(x=xs,y=ys,pch=16,col=p.cols)
    click.i <- 1    
    while(click.i < done.n)
    {
        click.i <- identify(xs,ys,n=1,plot=F)
        if(click.i < (length(m.names)+1) & click.i > 0)
        {
            i <- click.i
            if(is.element(i,x.sel))
            {
                lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lwds)
                x.sel <- setdiff(x.sel,i)
            }
                else
                {
	    	    lines(xseq,t.dat[,m.names[i]]+i*sf,col="black",lwd=lwds)
                #lines(xseq,t.dat[,m.names[i]]+i*sf,col="white",lwd=2,lty=2)
                x.sel <- union(x.sel,i)
            }
        }
        if(click.i == none.i)
        {
        	x.sel <- NULL
	    	for(i in 1:length(m.names))
		    {
    		    lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lwds)
	    	}
	    }
        if(click.i == all.i)	
        {
        	x.sel <- seq(1,length(m.names))
	    	for(i in 1:length(m.names))
		    {
    		    lines(xseq,t.dat[,m.names[i]]+i*sf,col="black",lwd=lwds)
	    	}
        	
        }
    }
    return(list(cells=m.names[x.sel],click = p.names[click.i]))
}




LinesSome <- function(t.dat,snr=NULL,m.names,wr,levs,lmain="",pdf.name=NULL,morder=NULL,subset.n=5,sf=.25,lw=2,bcex=.6)
{
	library(cluster)
	if(length(m.names) < subset.n)
	{stop("group size lower than subset size")}
	pam5 <- pam(t(t.dat[,m.names]),k=subset.n)
	s.names <- row.names(pam5$medoids)
	if(!is.null(morder))
	{
		names(morder) <- m.names
		morder <- morder[s.names]
	}
	pam5.tab <- table(pam5$clustering)
	tags <- paste(paste("#",names(pam5.tab),sep=""),as.vector(pam5.tab),sep=":")
	LinesEvery(t.dat=t.dat,snr=snr,m.names=s.names,wr=wr,levs=levs,lmain=lmain,pdf.name=pdf.name,morder=morder,rtag=tags,sf=sf,lw=lw,bcex=bcex)
	return(pam5$clustering)
}

LinesEvery <- function(t.dat,snr=NULL,m.names,wr,levs=NULL,lmain="",pdf.name=NULL,morder=NULL,rtag=NULL,sf=.25,lw=3,pw=.3,bcex=1,p.ht=7,p.wd=10,sig.line=NULL,scy=.05,scx=15/60,bw.flag=FALSE)
{
    m.names <- intersect(m.names,names(t.dat))
    xseq <- t.dat[,1]
    library(RColorBrewer)

    if(length(m.names) > 0)
    {
        if(is.null(pdf.name))
        {dev.new(width=14,height=8)}
        else
        {if(length(grep("\\.pdf",pdf.name))>0){pdf(pdf.name,width=p.wd,height=p.ht)}else{png(pdf.name,width=1200,height=600)}}#pdf(pdf.name,width=28,height=16)}
        if(is.null(morder))
        {
            m.pca <- prcomp(t(t.dat[,m.names]),scale=F,center=T)
            morder <- m.pca$x[,1] * c(1,-1)[(sum(m.pca$rot[,1]) < 0)+1]
            #m.names <- m.names[order(m.pca$x[,1],decreasing=sum(m.pca$rot[,1]) < 0)]
        }
        m.names <- m.names[order(morder)]
        
        hbc <- length(m.names)*sf+max(t.dat[,m.names])
        hb <- ceiling(hbc)
        #cols <- rainbow(length(m.names),start=.55)
		cols <-brewer.pal(8,"Dark2")
        cols <- rep(cols,ceiling(length(m.names)/length(cols)))
        cols <- cols[1:length(m.names)]
        if(bw.flag==T){cols[] <- "black"}
        par(mar=c(4,1,4,1))
        plot(xseq,t.dat[,m.names[1]],ylim=c(0,hbc),xlab="Time (min)",main=lmain,type="n", xaxt="n",yaxt="n",xlim=c(min(xseq)-1.5,max(xseq)+1.5))#-sf
        axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
		if(length(wr) > 0)
        {
        	if(!is.null(levs))
        	{
	           	#levs <- setdiff(unique(wr),"")
	            x1s <- tapply(xseq,as.factor(wr),min)[levs]
	            x2s <- tapply(xseq,as.factor(wr),max)[levs]
	            y1s <- rep(-.3,length(x1s))
	            y2s <- rep(hbc+.2,length(x1s))
	            rect(x1s,y1s,x2s,y2s,col=NA,border="darkgrey")
	            cpx <- xseq[match(levs,wr)+round(table(wr)[levs]/2,0)]
	            offs <- nchar(levs)*.5
	            text(cpx,rep(c(sf/2,sf),length=length(levs)),levs,pos=1,cex=bcex)#,offset=-offs
            }
        }
        for(i in 1:length(m.names))
        {
            lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lw)
            points(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],pch=15,cex=pw)
            if(!is.null(snr))
            {
	            pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
	            pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
	                                        #                pp3 <- dat$crr[,m.names[i]] > 0
	            points(xseq[pp1],t.dat[pp1,m.names[i]]+i/10,pch=1,col=cols[i])
	            points(xseq[pp2],t.dat[pp2,m.names[i]]+i/10,pch=0,col=cols[i])
	                                        #                points(xseq[pp3],t.dat[pp3,m.names[i]]+i/10,pch=2,col=cols[i],cex=.5)
            }    
			if(!is.null(sig.line))
			{
				
				lines(c(xseq[1],xseq[1]),c(t.dat[1,m.names[i]]+i*sf,t.dat[1,m.names[i]]+i*sf+sig.line),lwd=6,col="red")
			}
            
        }
        text(rep(0,length(m.names)),seq(1,length(m.names))*sf+t.dat[1,m.names],m.names,cex=.8*bcex,col=cols,pos=2)
        if(!is.null(rtag))
        {
        	rtag <- rtag[order(morder)]
	        text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag,cex=.8*bcex,col=cols,pos=4)
        }
      if(!is.null(pdf.name))
        {dev.off()}
    }
    
}

#Group summarry
#generate pdfs with line graphs
#table of means and frequencies for all c.dat
#THIS MUST BE CLEANED UP 040314
GroupSummary <- function(dat,snr,c.dat,wr,levs,groups,pref="Group")
{
    g.levs <- unique(groups)
    for(i in g.levs)
    {
        c.names <- names(groups[groups==i])
        pdf.name <- paste(pref,i,".pdf",sep="")
        lmain <- paste(pref,i,sep="")
        LinesEvery(dat,snr,c.names,wr,levs,lmain,pdf.name)
        dev.off()
    }
    res.tab <- data.frame(mean=apply(c.dat[names(groups),],2,mean))
    res.tab["sd"] <- apply(c.dat[names(groups),],2,sd)
    for(i in g.levs)
    {
        c.names <- names(groups[groups==i])
        res.tab[paste(pref,i,".mean",sep="")] <- apply(c.dat[c.names,],2,mean)
        res.tab[paste(pref,i,".sd",sep="")] <- apply(c.dat[c.names,],2,sd)
    }
    tab.name <- paste(pref,".table.csv",sep="")
    write.csv(res.tab,file=tab.name)
    #lines figure similar to boxplot
    ## tmp <- scale(c.dat[names(groups),],center=T,scale=T)
    ## tmp.mn <- data.frame(t(apply(tmp,2,function(x){tapply(x,as.factor(groups),mean)})))
    ## tmp.sd <- data.frame(t(apply(tmp,2,function(x){tapply(x,as.factor(groups),sd)})))
    ## tmp.se <- t(t(tmp.sd)/sqrt(summary(as.factor(groups))))
    ## ylim <- c(min(tmp.mn)-2,max(tmp.mn))
    ## miny <- min(ylim)+1
    ## dev.new()
    ## par(xaxt="n",mar=c(2,4,4,2))
    ## plot(seq(1,nrow(tmp.mn)),tmp.mn[,1],ylim=ylim,xlim=c(0,(nrow(tmp.mn)+1)),type="n",ylab="Normalized Mean +- SE",xaxt="n")
    ## cols <- rainbow(ncol(tmp.mn),start=.3)
    ## names(cols) <- names(tmp.mn)
    ## nudge <- 0
    ## ## for(i in names(tmp.mn))
    ## {
    ##     xseq <- seq(1,nrow(tmp.mn))
    ##     rect(nudge+seq(1,nrow(tmp.mn))-.05,tmp.mn[,i]-tmp.se[,i],nudge+seq(1,nrow(tmp.mn))+.05,tmp.mn[,i]+tmp.se[,i],col=cols[i],border=NA)
    ##     points(nudge+seq(1,nrow(tmp.mn)),tmp.mn[,i],pch=16,col=cols[i],lwd=2,type="b")        
    ##     nudge <- nudge+.1
    ## }
    ## text(rep(nrow(tmp.mn),ncol(tmp.mn)),tmp.mn[nrow(tmp.mn),],paste(pref,names(tmp.mn),sep=""),cex=.8,col=cols,pos=4)
    ## text(seq(1,nrow(tmp.mn))+.25,miny,names(c.dat),srt=90,pos=3)
    c.mn <- data.frame(t(apply(c.dat,2,function(x){tapply(x,as.factor(groups),mean)})))
    c.sd <- data.frame(t(apply(c.dat,2,function(x){tapply(x,as.factor(groups),sd)})))
    c.se <- t(t(c.sd)/sqrt(summary(as.factor(groups))))
    return(list(mean=c.mn,sd=c.sd,se=c.se))   
}
        

#take a defined window vector and
#number of the contiguos blank regions ("")
NumBlanks <- function(x)
{
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

NumAny <- function(x,targ)
{
    nw <-  as.character(x)
    xlen <- length(x)
    bl.cnt <- 1
    mi <- match(targ,nw)
    while(!is.na(mi) & bl.cnt < 20)
    {
        mi2 <- mi+1
        while((x[mi2] == targ) & mi2 <= xlen){mi2 <- mi2+1}
        nw[mi:(mi2-1)] <- paste(targ,bl.cnt,sep="")
        bl.cnt <- bl.cnt+1
        mi <- match(targ,nw)
    }
    return(as.factor(nw))	
}

#measure the variance about the smooth
SmoothVar <- function(x,shws=2)
{
    s1 <- createMassSpectrum(1:length(x),x)
    s3 <- smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=shws)
    return(sd(x-intensity(s3)))
}

#Given two experiments dat1 and dat2 merge into one experiment
#aligning the treatment window variables specified by name in wr1 and wr2
#treatments in wr1 and wr2 MUST be the same in the same order.
#THIS FUNCTION SHOULD BE USED WITH CAUTION
#ALIGNMENT WILL NOT BE EXACT UNLESS THE EXPERIMENTS WERE IDENTICAL
#align the two start points. 
#identify the points of lowest variation within an experiment in the interval to the next
#window.  duplicate points until the first windows align.  repeat for the next interval
#PROBABLY SHOULD IDENTIFY NEW WINDOW POINTS AT THE MIDPOINT OF BLANK REGIONS.
#THIS WILL BE BEST FOR SLIPPING POINTS.
Merge2 <- function(dat1,dat2,wr1,wr2)
{
	w1 <- dat1$w.dat[,wr1]
	w2 <- dat2$w.dat[,wr2]
	w1.lev <- unique(w1)[-1]
	w2.lev <- unique(w2)[-1]
	if(sum(w1.lev != w2.lev)>0){stop("windows not identical")}
	#Define a midpoints between w1 and w2
	p1 <- as.integer((c(match(w1.lev,w1),length(w1))+c(1,match(w1.lev,w1)))/2)
	p2 <- as.integer((c(match(w2.lev,w2),length(w2))+c(1,match(w2.lev,w2)))/2)

	s1 <- seq(1,length(w1))
	s2 <- seq(1,length(w2))
	#first gap from start to first window
#	i <- w1.lev[1]
#	gap1 <- seq(1,match(i,w1))
#	gap2 <- seq(1,match(i,w2))
	gap1 <- seq(1,p1[1])
	gap2 <- seq(1,p2[1])
	v1 <- apply(dat1$t.dat[gap1,-1],1,sd,na.rm=T)
	v2 <- apply(dat2$t.dat[gap2,-1],1,sd,na.rm=T)
	if(length(gap1) > length(gap2))
	{
		#v <- apply(dat2$t.dat[gap2,-1],1,sd,na.rm=T)
		#gap2 <- SmearFunc(gap2,length(gap1),v)
		gap2 <- gap2[AlignFunc(v2,v1)]
	}
	if(length(gap2) > length(gap1))
	{
		#v <- apply(dat1$t.dat[gap1,-1],1,sd,na.rm=T)
		#gap1 <- SmearFunc(gap1,length(gap2),v)
		gap1 <- gap1[AlignFunc(v1,v2)]
	}
	seq1 <- gap1
	seq2 <- gap2
	for(i in 1:(length(p1)-1))
	{
		gap1 <- seq(p1[i]+1,p1[(i+1)])
		gap2 <- seq(p2[i]+1,p2[(i+1)])
		v1 <- apply(dat1$t.dat[gap1,-1],1,sd,na.rm=T)
		v2 <- apply(dat2$t.dat[gap2,-1],1,sd,na.rm=T)		
		if(length(gap1) > length(gap2))
		{
			#v <- apply(dat2$t.dat[gap2,-1],1,sd,na.rm=T)
			#gap2 <- SmearFunc(gap2,length(gap1),v)
			gap2 <- gap2[AlignFunc(v2,v1)]
		}
		if(length(gap2) > length(gap1))
		{
			#v <- apply(dat1$t.dat[gap1,-1],1,sd,na.rm=T)
			#gap1 <- SmearFunc(gap1,length(gap2),v)
			gap1 <- gap1[AlignFunc(v1,v2)]
		}
		seq1 <- c(seq1,gap1)
		seq2 <- c(seq2,gap2)		
	}	
	#now to the end.
	gap1 <- seq(max(p1)+1,length(w1))
	gap2 <- seq(max(p2)+1,length(w2))
	v1 <- apply(dat1$t.dat[gap1,-1],1,sd,na.rm=T)
	v2 <- apply(dat2$t.dat[gap2,-1],1,sd,na.rm=T)		
	
	if(length(gap1) > length(gap2))
	{
#		v <- apply(dat2$t.dat[gap2,-1],1,sd,na.rm=T)
#		gap2 <- SmearFunc(gap2,length(gap1),v)
			gap2 <- gap2[AlignFunc(v2,v1)]

	}
	if(length(gap2) > length(gap1))
	{
#		v <- apply(dat1$t.dat[gap1,-1],1,sd,na.rm=T)
#		gap1 <- SmearFunc(gap1,length(gap2),v)
			gap1 <- gap1[AlignFunc(v1,v2)]

	}
	seq1 <- c(seq1,gap1)
	seq2 <- c(seq2,gap2)		
	tot.dat <- dat1
	tot.dat$w.dat <- cbind(dat1$w.dat[seq1,],dat2$w.dat[seq2,])
	tot.dat$t.dat <- cbind(dat1$t.dat[seq1,],dat2$t.dat[seq2,-1])
	names(tot.dat$t.dat) <- c("Time",paste(names(dat1$t.dat)[-1],".1",sep=""),paste(names(dat2$t.dat)[-1],".2",sep=""))
	tot.dat$t.dat[,"Time"] <- seq(1,length(seq1))
	c.names <- union(names(dat1$c.dat),names(dat2$c.dat))
	for(i in setdiff(c.names,names(dat1$c.dat))){dat1$c.dat[i] <- NA}
	for(i in setdiff(c.names,names(dat2$c.dat))){dat2$c.dat[i] <- NA}
	tot.dat$c.dat <- rbind(dat1$c.dat[,c.names],dat2$c.dat[,c.names])
	c.names <- c(paste(dat1$c.dat[,"id"],".1",sep=""),paste(dat2$c.dat[,"id"],".2",sep=""))
	tot.dat$c.dat[,"id"] <- c.names
	row.names(tot.dat$c.dat) <- c.names
	tot.dat$c.dat["exp"] <- c(rep(1,nrow(dat1$c.dat)),rep(2,nrow(dat2$c.dat)))
	return(tot.dat)
}
#smear vector x out to length n
#using the lowest variances indicated in v
#NO NAs
SmearFunc <- function(x,n,v)
{
	n1 <- length(x)
	if(n1 >= n)
	return(x)
	difn <- n-n1
	min.i <- which.min(v)
	n.seq <- c(seq(1,min.i),rep(min.i,(difn-1)),seq(min.i,n1))
	return(x[n.seq])	
}

#shift compare and align 
#align the two sequences
#sv must me the shorter of the two.
#not allowed to shift outside the boundaries.
AlignFunc <- function(sv,lv)
{	
	corfunc <- function(x){cor(c(rep(NA,x),sv,rep(NA,dif.n-x)),lv,use="pair")}
	dif.n <- length(lv)-length(sv)
	ct <- sapply(0:dif.n,corfunc)
	best.i <- which.max(ct)-1
	return(c(rep(1,best.i),seq(1,length(sv)),rep(length(sv),dif.n-best.i)))
	
}

#given a traces region estimate the complexity (e.g. baseline deviation)
#currently this is the amount of variation explained by a polynomial regression 
#of degree 3.
#consider using the slope of the regression line for increasing degree on the
#polynomial.
#names of the tdat are the individual ids.
#note no time column
EstimateComplexity <- function(tdat,np=3)
{
	res <- data.frame(mean=apply(tdat,2,mean))
	res["sd"] <- apply(tdat,2,sd)
	nr <- nrow(tdat)
	nc <- ncol(tdat)
	i.names <- names(tdat)
	x <- seq(1,nr)
	for(p in 1:np)
	{
		p.name <- paste("p",p,"r2",sep=".")
	for(i in i.names)
	{
		lt <- lm(tdat[,i] ~ poly(x,p))
		lt.sum <- summary(lt)
		res[i,p.name] <- lt.sum$r.squared
	}
	}
	res[,"pr2.slope"] <- NA
	if(ncol(res) < 4){return(res)}
	p.seq <- 3:(ncol(res))
	for(i in i.names)
	{
	  lt <- lm(t(res[i,p.seq]) ~ seq(1,length(p.seq)))
	  res[i,"pr2.slope"] <- lt$coefficients[2]		
	}
	return(res)
}
#estimate complexity as the sum of log10(p-values)
#except for excluded coefficients ex.coef
EstimateComplexity2 <- function(tdat,np=10,ex.coef=c(1,2,3))
{
	res <- data.frame(mean=apply(tdat,2,mean))
	res["sd"] <- apply(tdat,2,sd)
	nr <- nrow(tdat)
	nc <- ncol(tdat)
	i.names <- names(tdat)
	x <- seq(1,nr)
	for(i in i.names)
	{
		lt <- lm(tdat[,i] ~ poly(x,np))
		lt.sum <- summary(lt)
		res[i,"sum.logp"] <- sum(-log10(lt.sum$coefficients[,4][-ex.coef]))
		lt <- lm(tdat[,i] ~ poly(x,(length(ex.coef))))
		res[i,"res.sd"] <- sd(residuals(lt))
	}

	return(res)
}

#report the start of all monotonic increases or decreases
#for the sequence of points x.
TrendFunc <- function(x,tlen=5)
{
	res <- NULL
	lp <- 1
	for(i in 1:(length(x)-tlen))
	{
		if(sum(x[i:(i+tlen-1)] < x[(i+1):(i+tlen)]) == tlen)
		{
			res[lp] <- i
			lp <- lp+1
		}
	}
	return(res)	

}

#all blanks estimate the ide for all levels indicated in blevs
All.DE <- function(t.dat,wr,blevs,tleft=0,tright=0)
{
	blevs <- as.character(blevs)
	bdat <- t.dat[wr==blevs[1],-1]
	if(tleft > 0 & tleft < nrow(bdat))
	{bdat <- bdat[tleft:nrow(bdat),]}
	if(tright > 0 & tright < nrow(bdat))
	{bdat <- bdat[1:(nrow(bdat)-tright),]}
	ec2 <- EstimateComplexity2(bdat,np=3,ex.coef=c(1))
	if(mean(ec2[,"sum.logp"]) > 3 || max(ec2[,"sum.logp"]) > 10){warnings("significant higher-order effects")}
	res <- data.frame(ec2[,"res.sd"])	
	for(i in blevs[-1])
	{
		bdat <- t.dat[wr==i,-1]
		if(tleft > 0 & tleft < nrow(bdat))
		{bdat <- bdat[tleft:nrow(bdat),]}
		if(tright > 0 & tright < nrow(bdat))
		{bdat <- bdat[1:(nrow(bdat)-tright),]}
		
		ec2 <- EstimateComplexity2(bdat,np=3,ex.coef=c(1))
		if(mean(ec2[,"sum.logp"]) > 3 || max(ec2[,"sum.logp"]) > 10){warnings("significant higher-order effects")}
        	res[i] <- ec2[,"res.sd"]    			
	}
	names(res) <- paste(blevs,".DE",sep="")
	return(res)	
}

#give a set of pulses estimate the deviation of each pulse from the average.
#for each individual estimate mean and sd using all other pulses then compute
#the #standard deviation distance for each cell.
All.IDE <- function(scp.dat)
{
	res.dat <- scp.dat
	names(res.dat) <- paste(names(res.dat),".IDE",sep="")
	t.names <- names(scp.dat)
	for(i in 1:length(t.names))
	{
		mn.cnt <- apply(scp.dat[,-i],1,mean)
		sd.cnt <- apply(scp.dat[,-i],1,sd)
		res.dat[i] <- (scp.dat[,i]-mn.cnt)/mn.cnt
	}
	return(res.dat)
	
}

#no NAs for now please
Size.Bin <- function(x)
{
	x.bin <- rep("M",length(x))
	x.bin[x < 300] <- "S"
	x.bin[x > 600] <- "L"
	return(x.bin) 
}

#give a set of traces divide them into clusters of 5 or more
#graph the clusters with transparency and color
PeakShapeView <- function(t.dat,min.n=4,cut.h=.2,cut.k=NULL,aval=100,wd.val=4,main.lab="",pdf.name=NULL)
{
	library(RColorBrewer)
	
	t.dat[1] <- round(t.dat[,1]*60,0) #convert time to seconds
	cr <- cor(t.dat[,-1])
	cr1 <- cor(t.dat[seq(1,nrow(t.dat),by=2),-1])
	cr2 <- cor(t.dat[seq(2,nrow(t.dat),by=2),-1])
	cval <- apply(cbind(cr[lower.tri(cr)],cr1[lower.tri(cr1)],cr2[lower.tri(cr)]),1,median)
	#cval <- cr[lower.tri(cr)]
	cdist <- dist(t(t.dat[,-1]))
	cdist[] <- 1-cval
	chclust <- hclust(cdist)
	if(is.null(pdf.name))
	{
	dev.new()
	plot(chclust)
	}
	if(!is.null(cut.k))
		{cutc <- cutree(chclust,k=cut.k)}
	else
		{cutc <- cutree(chclust,h=cut.h)}
	
	cut.tab <- table(cutc)	
	ng <- sum(cut.tab > min.n)
	#cols <- rainbow(ng,alpha=.1)
	cols <-brewer.pal(8,"Dark2")
	col.mat <- col2rgb(cols)
	for(i in 1:length(cols)){cols[i] <- rgb(col.mat[1,i],col.mat[2,i],col.mat[3,i],alpha=aval,max=255)}
    cols <- rep(cols,ceiling(ng/length(cols)))
    cols <- cols[1:ng]
	gr <- (1:length(cut.tab))[cut.tab > min.n]
	
	med.dat <- t.dat[,1:ng]
	lp<-1
	for(i in gr)
	{
		x.names <- names(cutc[cutc==i])
		med <- row.names(pam(t(t.dat[,x.names]),k=1)$medoids)
		med.dat[lp] <- t.dat[,med]
		lp <- lp+1
	}
	cols <- cols[order(apply(med.dat,2,mean))]
	op <- par(xaxp=c(range(t.dat[,1]),1))

	if(is.null(pdf.name)){dev.new()}else{pdf(pdf.name)}
	
	plot(t.dat[,1],t.dat[,2],type="n",ylim=range(med.dat)*c(1,1.1),xlab="seconds",ylab="fitted ratio",main=main.lab)
	lp <- 1
	for(i in gr)
	{
		x.names <- names(cutc[cutc==i])
		med <- row.names(pam(t(t.dat[,x.names]),k=1)$medoids)
		y.names <- setdiff(x.names,med)
		lines(t.dat[,1],t.dat[,med],type="l",lwd=wd.val,col=cols[lp])
		for(i in y.names)
		{
			lt <- lm(t.dat[,med] ~ t.dat[,i])
			lines(t.dat[,1],predict(lt),type="l",lwd=wd.val,col=cols[lp])
		}
		lp <- lp+1	
	}
	if(!is.null(pdf.name)){dev.off()}
	return(cutc)
}

##review binary scoring file and toggle 1/0
##names of binary scoring bin must be in wr
##NO NAs
ScoreReview1 <- function(tdat,bin,wr,maxt=20)
{
	subD <- function(xdat)#trace dat with names NO TIME COL
	{
		s.x <- apply(xdat,2,sum)
		s.names <- names(xdat)[order(s.x)]
		sub.list <- list()
		sub.i <- seq(1,ncol(xdat),by=(maxt+1))
		if(length(sub.i) > 1)
		{
		for(i in 1:(length(sub.i)-1))
		{
			sub.list[[i]] <- s.names[sub.i[i]:(sub.i[i]+maxt)]
		}
		}
		i <- length(sub.i)
		sub.list[[i]] <- s.names[sub.i[i]:(ncol(xdat))]
		return(sub.list)
	}
	
	b.levs <- names(bin)[names(bin) != "drop"]
	drop <- rep(0,nrow(bin))
	#if(is.element("drop",names(bin))){drop <- bin[,"drop"]}
	names(drop) <- row.names(bin)
	for(i in b.levs)
	{
		lmain <- paste("Scored as 1 for ",i,sep="")
		b.1 <- row.names(bin)[bin[,i]==1 & drop==0]
		if(length(b.1) > 0)
		{
		if(length(b.1) < maxt){sub1 <- list(b.1)}else{sub1 <- subD(tdat[wr==i,b.1])}
		for(x.names in sub1)
		{
			no.names <- NULL
			dropit <- TRUE
			while(dropit==TRUE & (length(x.names) > 0))
			{
				inp <- ScoreSelect(tdat,,x.names,wr,i,lmain)
				no.names <- inp[["cells"]]
				dropit <- (inp[["click"]]=="DROP")
				if(dropit){drop[no.names] <- 1;x.names <- setdiff(x.names,no.names)}
				dev.off()
			}
			if(length(no.names) > 0)
			{
				bin[no.names,i] <- 0
			}
		}
		}
	}
	bin["drop"] <- drop
	return(bin)
		
}


ScoreReview0 <- function(tdat,bin,wr,maxt=20)
{
	subD <- function(xdat)#trace dat with names NO TIME COL
	{
		s.x <- apply(xdat,2,sum)
		s.names <- names(xdat)[order(s.x)]
		sub.list <- list()
		sub.i <- seq(1,ncol(xdat),by=(maxt+1))
		if(length(sub.i) > 1)
		{
		for(i in 1:(length(sub.i)-1))
		{
			sub.list[[i]] <- s.names[sub.i[i]:(sub.i[i]+maxt)]
		}
		}
		i <- length(sub.i)
		sub.list[[i]] <- s.names[sub.i[i]:(ncol(xdat))]
		return(sub.list)
	}
	
	b.levs <- names(bin)[names(bin) != "drop"]
	drop <- rep(0,nrow(bin))
	if(is.element("drop",names(bin))){drop <- bin[,"drop"]}
	names(drop) <- row.names(bin)
	for(i in b.levs)
	{
		lmain <- paste("Scored as 0 for ",i,sep="")
		b.1 <- row.names(bin)[bin[,i]==0 & drop==0]
		if(length(b.1) > 0)
		{
		if(length(b.1) < maxt){sub1 <- list(b.1)}else{sub1 <- subD(tdat[wr==i,b.1])}			
		for(x.names in sub1)
		{
			no.names <- NULL
			dropit <- TRUE
			while(dropit==TRUE & (length(x.names)>0))
			{
				inp <- ScoreSelect(tdat,,x.names,wr,i,lmain)
				no.names <- inp[["cells"]]
				dropit <- (inp[["click"]]=="DROP")
				if(dropit){drop[no.names] <- 1;x.names <- setdiff(x.names,no.names)}
				dev.off()
			}
			if(length(no.names) > 0)
			{
				bin[no.names,i] <- 1
			}
		}
		}
	}
	bin["drop"] <- drop
	return(bin)
}

#Given a list of traces allow for partitioning into 4 groups
#presumably these will be drop, unk, class1, class2
#drops are removed
#class1 traces are shown in separate panel (updated each shift change)
#class2 traces are shown in separate panel (updated each shift change)
ClassReview <- function(tdat,bin,wr,maxt=10)
{
	subD <- function(xdat)#trace dat with names NO TIME COL
	{
		s.x <- apply(xdat,2,sum)
		s.names <- names(xdat)[order(s.x)]
		sub.list <- list()
		sub.i <- seq(1,ncol(xdat),by=(maxt+1))
		if(length(sub.i) > 1)
		{
		for(i in 1:(length(sub.i)-1))
		{
			sub.list[[i]] <- s.names[sub.i[i]:(sub.i[i]+maxt)]
		}
		}
		i <- length(sub.i)
		sub.list[[i]] <- s.names[sub.i[i]:(ncol(xdat))]
		return(sub.list)
	}
	
	b.levs <- names(bin)[names(bin) != "drop"]
	drop <- rep(0,nrow(bin))
	if(is.element("drop",names(bin))){drop <- bin[,"drop"]}
	names(drop) <- row.names(bin)
	for(i in b.levs)
	{
		lmain <- paste("Scored as 0 for ",i,sep="")
		b.1 <- row.names(bin)[bin[,i]==0 & drop==0]
		if(length(b.1) > 0)
		{
		if(length(b.1) < maxt){sub1 <- list(b.1)}else{sub1 <- subD(tdat[wr==i,b.1])}			
		for(x.names in sub1)
		{
			no.names <- NULL
			dropit <- TRUE
			while(dropit==TRUE & (length(x.names)>0))
			{
				inp <- ScoreSelect(tdat,,x.names,wr,i,lmain)
				no.names <- inp[["cells"]]
				dropit <- (inp[["click"]]=="DROP")
				if(dropit){drop[no.names] <- 1;x.names <- setdiff(x.names,no.names)}
				dev.off()
			}
			if(length(no.names) > 0)
			{
				bin[no.names,i] <- 1
			}
		}
		}
	}
	bin["drop"] <- drop
	return(bin)
}


SpikyNorm <- function(xdat)
{
		shapfunc <- function(x){shapiro.test(x)$p.value}
		i1 <- seq(1,nrow(xdat))
		s1 <- xdat[c(1,i1[-length(i1)]),] #shift 1 time interval forward
		s2 <- xdat[c(i1[-1],i1[length(i1)]),] #shift 1 time interval back
		s3 <- xdat-((s1+s2)/2)
		s.x <- apply(abs(s3),2,shapfunc)
	return(s.x)	
}

DropPick <- function(tdat,bin,wr,maxt=10,s.x=NULL,lmain="Select Cells to drop")
{
	#order traces by spikey trait.
	#allow drop selection until 0 selected.
	#spikes are defined as single point deviations from previous and next.
	subD <- function(s.x)#trace dat with names NO TIME COL
	{
		s.names <- names(s.x)[order(s.x)]
		sub.list <- list()
		sub.i <- seq(1,length(s.x),by=(maxt+1))
		if(length(sub.i) > 1)
		{
		for(i in 1:(length(sub.i)-1))
		{
			sub.list[[i]] <- s.names[sub.i[i]:(sub.i[i]+maxt)]
		}
		}
		i <- length(sub.i)
		sub.list[[i]] <- s.names[sub.i[i]:(length(s.x))]
		return(sub.list)
	}
	
	b.levs <- c("drop") #names(bin)[names(bin) != "drop"]
	drop <- rep(0,nrow(bin))
	if(is.element("drop",names(bin))){drop <- bin[,"drop"]}
	names(drop) <- row.names(bin)
	for(i in b.levs)
	{

		b.1 <- row.names(bin)[bin[,i]==0 & drop==0]
		if(is.null(s.x)){s.x <- SpikyNorm(tdat[,-1])}

		if(length(b.1) > 0)
		{
			s.x <-s.x[b.1]
		if(length(b.1) < maxt){sub1 <- list(b.1)}else{sub1 <- subD(s.x)}

		for(x.names in sub1)
		{
			no.names <- NULL
			dropit <- TRUE
			nd <- 0			
			while(dropit==TRUE & (length(x.names)>0))
			{

				inp <- ScoreSelect(tdat,,x.names,wr,,lmain)
				no.names <- inp[["cells"]]
				dropit <- (inp[["click"]]=="DROP")
				if(dropit){drop[no.names] <- 1;x.names <- setdiff(x.names,no.names);nd=1}
				dev.off()
			}
			if(length(no.names) > 0)
			{
				drop[no.names] <- 1
			}
			if(length(no.names)==0 & nd==0)
			{break}
		}
		}
	}
	return(drop)
}

DropTestList <- function(tmp)
{
		#tmp <- get(rd.name)
		#x1 <- DropPick(tmp$t.dat,tmp$bin,tmp$w.dat[,"wr1"],lmain="Select spikey traces to Drop") #defaults to spiky test
		#tmp$bin[,"drop"] <- x1
		#x1 <- DropPick(tmp$t.dat,tmp$bin,tmp$w.dat[,"wr1"],s.x= -apply(tmp$scp[,"snr.owc",drop=F],1,mean),lmain="Select out of window peaks to Drop")
		#tmp$bin[,"drop"] <- x1		
		x1 <- DropPick(tmp$t.dat,tmp$bin,tmp$w.dat[,"wr1"],s.x= -apply(tmp$scp[,"bl.diff",drop=F],1,mean),lmain="Select Baseline Drops")
		tmp$bin[,"drop"] <- x1		
		# if(sum(x1 > 2)) #check highest correlations with dropped cells.
		# {
			# d.names <- names(x1[x1==1])
			# ct <- cor(tmp$t.dat[,-1])
			# mn <- -apply(ct[,d.names],1,max)
			# x1 <- DropPick(tmp$t.dat,tmp$bin,tmp$w.dat[,"wr1"],s.x= mn,lmain="Correlated with other drops")
			# tmp$bin[,"drop"] <- x1		
		# }
		return(tmp)
}

DropTestMulti <- function(dir.name=NULL,snr.lim=4,hab.lim=.05,sm=3,ws=30,review=F)
{
	if(is.null(dir.name)){dir.name <- getwd()}
	setwd(dir.name)
	f.names <- list.files(pattern="RD.*\\.Rdata$")
	if(length(f.names) == 0){stop("no RD...Rdata files in given directory")}
	rd.list <- sub("\\.Rdata*","",f.names)
	RD.names <- rd.list #paste(rd.list,".b",sep="")
	RD.f.names <- paste(RD.names,".Rdata",sep="")
	sel.i <- menu(rd.list,title="Select Data to review")
	while(sel.i != 0)
	{
		j <- sel.i
		load(f.names[j])
		i <- rd.list[j]
		tmp <- get(i)
		tlevs <- c(as.character(unique(tmp$w.dat[,"wr1"])[-1]),"drop")
		if(is.null(tmp$bin))
		{
		tmp.pcp <- ProcConstPharm(tmp,sm,ws,"TopHat")
		tmp.scp <- ScoreConstPharm(tmp$t.dat,tmp.pcp$blc,tmp.pcp$snr,snr.lim,hab.lim,tmp$w.dat[,"wr1"],sm)
		tmp.bin <- bScore(tmp.pcp$blc,tmp.pcp$snr,snr.lim,hab.lim,tlevs,tmp$w.dat[,"wr1"])
		tmp.bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
		}
		else
		{
			tmp.bin <- tmp$bin
			tmp.scp <- tmp$scp
			#tmp.blc <- tmp$blc
		}

		tmp$bin <- tmp.bin[,tlevs]
		tmp$scp <- tmp.scp
		#tmp$blc <- tmp.blc
		
		tmp <- DropTestList(tmp)
		if(review)
		{
		  	tmp.bin <- ScoreReview1(tmp$t.dat,tmp.bin[,tlevs],tmp$w.dat[,"wr1"])
		  	tmp.bin <- ScoreReview0(tmp$t.dat,tmp.bin[,tlevs],tmp$w.dat[,"wr1"])
	    	tmp$bin <- tmp.bin[,tlevs]
		}
				
		assign(RD.names[j],tmp)		
		save(list=RD.names[j],file=RD.f.names[j])
		print(paste("DONE REVIEWING ",RD.names[j]," CHANGES SAVED TO FILE.",sep=""))
		sel.i <- menu(rd.list,title="Select Data to review")			
	}
	return(RD.f.names)		
}
	
			

#return the counts for each pairwise combination of columns
#below the diagnal and significance above the diagnal 0/1 values only
CountTab <- function(x)
{
	ct <- cor(x)
	ct[,] <- NA
	for(i in 1:(ncol(x)-1))
	{
		for(j in (i+1):ncol(x))
		{
			f.tab <- table(x[,i],x[,j])
			if(ncol(f.tab)==2 & nrow(f.tab)==2)
			{
			f.fish <- fisher.test(f.tab)
			ct[i,j] <- -log10(f.fish$p.value)
			ct[j,i] <- f.tab[2,2]
			}
		}
	}
	return(ct)	
}

semfunc <- function(x)
{
	n <- sum(!is.na(x))
	if(n < 3){return(NA)}
	return(sd(x,na.rm=T)/sqrt(n))
}

#given an rd data list and a list of traces
#output a summary of the traces.
#plain plot
#fancy plot
#binary score
#peak heights
#raw trace
#blc trace
#roi

TraceSummary <- function(dat,blc=NULL,x.names)
{
	
}


#1 dat, trace data the first column must be Time 
#2 x.names, which traces to plot.
#3 png.name, the file name
#4 f.ht, the image height
#5 f.wt, image width
#6 sfp, the amount to shift each trace plot. (spacing between traces)
#7 lwd, The line width.
#8 xlim, range of x-values (defualts to all)
#9 rtag, the label to show on to the right of the traces
PlainPlot <- function(dat,x.names,png.name="tmp.png",f.ht=800,f.wt=1200,sfp=0,lwd=4,xlim=NULL,rtag=NULL)
{
	png(png.name,height=f.ht,width=f.wt,bg = "transparent")
	op <- par(yaxt="n",bty="n",mar=c(4,0,0,1),cex=2)
	xval <- dat[,1]
	if(is.null(xlim)){xlim <- range(xval)}
	x <- dat[,x.names]
	if(sfp==0){sfp <- max(x)} 	
	ylim <- c(0,sfp*length(x.names))
	plot(dat[,1],dat[,x.names[1]],xlim=xlim,ylim=ylim,type="n",ylab="",xlab="Minutes")
	#mtext("Minutes",side=1)
	sf <- 0
	for(i in 1:ncol(x))
	{
		lines(xval,x[,i]+sf,lwd=4)
		sf <- sf+sfp	
	}
	if(!is.null(rtag))
	{
		xs <- rep(max(xval),length(x.names))
		ys <- (0:(length(x.names)-1))*sfp
		text(xs,ys,format(rtag,width=4,justify="right"),pos=4,cex=.8,offset=0.1)
	}
	dev.off()				
}


#given a raw data dataframe name
#bin and scp output a summary of the experiment
ExperimentSummary <- function(dat.name,wr.i="wr1")
{
	
	tmp <- get(dat.name)
	if(!is.element("bin",names(tmp))){stop("no binary scoring (bin)")}
	if(!is.element("scp",names(tmp))){stop("no quantitative scoring (scp)")}
	#output pdf overview

	x.names <- row.names(tmp$bin)[tmp$bin[,"drop"]==0]
	p.levs <- as.character(unique(tmp$w.dat[,wr.i]))
	p.levs <- p.levs[p.levs != ""]
	pmain <- paste(dat.name," Cells n=",length(x.names),sep="")
	pdfname <- paste(dat.name,".ex.pdf",sep="")
	rlab <- NULL
	if(is.element("ROI.Area",names(tmp$c.dat))){rlab <- round(tmp$c.dat[x.names,"ROI.Area"],0)}    
	LinesSome(tmp$t.dat,,x.names,tmp$w.dat[,wr.i],p.levs,sf=.5,pdf=pdfname,lw=3,bcex=.6,subset.n=20,lmain=pmain)

	x.names <- row.names(tmp$bin)[tmp$bin[,"drop"]==1]
	if(length(x.names) > 0)
	pmain <- paste(dat.name," Dropped Cells n=",length(x.names),sep="")
	pdfname <- paste(dat.name,".ex.D.pdf",sep="")
	rlab <- NULL
	if(is.element("ROI.Area",names(tmp$c.dat))){rlab <- round(tmp$c.dat[x.names,"ROI.Area"],0)}    
	if(length(x.names) > 20)
	{LinesSome(tmp$t.dat,,x.names,tmp$w.dat[,wr.i],p.levs,sf=.5,pdf=pdfname,lw=3,bcex=.8,lmain=pmain,subset.n=20)}
	else	
	{LinesEvery(tmp$t.dat,,x.names,tmp$w.dat[,wr.i],p.levs,sf=.5,pdf=pdfname,lw=3,bcex=.8,lmain=pmain,rtag=rlab)}	

	#output frequency binary scoring and magnitude of response means
	tab.name <- paste(dat.name,".FreqMag.csv",sep="")
	b.names <- names(tmp$bin)
	b.names <- b.names[!is.element(b.names,c("drop","sd","tot"))]
	t.names <- paste(b.names,".max",sep="")
	d <- tmp$bin[,"drop"]==0
	dk <- d
	ki <- NULL
	if(length(grep("KCl",names(tmp$bin))) >0){ki <- grep("KCl",names(tmp$bin),ignore.case=T)[1];dk <- tmp$bin[,"drop"]==0 & tmp$bin[,ki]==1}
	freq.tab <- data.frame(mean=apply(tmp$bin[d,b.names],2,mean))
	freq.tab["sem"] <- apply(tmp$bin[d,b.names],2,semfunc)
	mag.tab <- tmp$scp[,t.names]
	mag.tab[tmp$bin[,b.names]==0] <- NA
	freq.tab["mag"] <- apply(mag.tab[d,t.names],2,mean,na.rm=T)		
	freq.tab["mag.sem"] <- apply(mag.tab[d,t.names],2,semfunc)				
	
	if(!is.null(ki))
	{
		freq.tab["KCl.mean"] <- apply(tmp$bin[dk,b.names],2,mean)
		freq.tab["KCl.sem"] <- apply(tmp$bin[dk,b.names],2,semfunc)		
		freq.tab["KCl.mag"] <- apply(mag.tab[dk,t.names],2,mean,na.rm=T)		
		freq.tab["KCl.mag.sem"] <- apply(mag.tab[dk,t.names],2,semfunc)				
	}
	write.csv(freq.tab,file=tab.name)
	#what to do about cell details file


}

#given a raw data data frame
#print all non-dropped cells to a pdf ameniable to printing.
PrintAll <- function(dat.name,mlab=NULL)
{
	pdf.name <- paste(dat.name,".all.so.pdf",sep="")
	tmp <- get(dat.name)
	x.names <- row.names(tmp$bin[tmp$bin[,"drop"]==0,])
	size <- round(tmp$c.dat[x.names,"ROI.Area"],0)
#	ib4 <- tmp$c.dat[x.names,"IB4.manual"]
#	ib4[is.na(ib4)] <- 0
#	ib4[ib4 > 0] <- 1
#	rlab <- paste(size,ib4,sep=":")
	rlab <- size
	levs <- unique(tmp$w.dat[,"wr1"])[-1]
	p.ht <- round(length(x.names)/6,0)
	mlab <- paste(dat.name,mlab,sep=":")
	LinesEvery(tmp$t.dat,,x.names,tmp$w.dat[,"wr1"],levs,sf=.5,lw=2,pdf.name = pdf.name,lmain=mlab,morder=size,rtag=size,p.ht=p.ht,bcex=.8)
}

PrintAll2 <- function(dat.name,mlab=NULL,x.names=NULL)
{
	tmp <- get(dat.name)
	if(is.null(x.names)){x.names <- row.names(tmp$bin[tmp$bin[,"drop"]==0,])}
	size <- tmp$c.dat[x.names,"ROI.Area"]
	x.names <- x.names[order(size)]
	size <- round(tmp$c.dat[x.names,"ROI.Area"],0)
	x1s <- seq(1,length(x.names),by=20)
	x2s <- x1s+19
	x2s[length(x2s)] <- length(x.names)
	mlab <- paste(dat.name,mlab,sep=":")
	levs <- unique(tmp$w.dat[,"wr1"])[-1]	
	for(i in 1:length(x1s))
	{
		pdf.name <- paste(dat.name,".all.p",i,".pdf",sep="")		
		xn <- x.names[x1s[i]:x2s[i]]
		rlab <- size[x1s[i]:x2s[i]]
		LinesEvery(tmp$t.dat,,xn,tmp$w.dat[,"wr1"],levs,sf=.5,lw=2,pdf.name = pdf.name,lmain=mlab,rtag=rlab,bcex=.8)		
	}

}

#XYtrace positions
XYtrace <- function(rdat)
{
	levs <- unique(rdat$w.dat[,"wr1"])
	col="black"
	pch=16
	cex=1
	lmain="X Y ROI"
	dev.new(height=4,width=14)
	rr.dev <- dev.cur()
	dev.new()
	plot(rdat$c.dat[,"Center.X"],rdat$c.dat[,"Center.Y"],col=col,pch=pch,cex=cex,main=lmain)
	i <- identify(rdat$c.dat[,"Center.X"],rdat$c.dat[,"Center.Y"],n=1,plot=F)
	my.dev <- dev.cur()
	while(length(i) > 0)
	{
		x.names <- row.names(tmp$c.dat)[i]		
		#points(lookat[i,"x"],lookat[i,"y"],pch=8,cex=.5)
		lmain <- i
		dev.set(which=rr.dev)
		PeakFunc2(rdat$t.dat,x.names,3,30,TRUE,rdat$w.dat[,"wr1"],lmain=lmain)
		dev.set(which=my.dev)
		i <- identify(rdat$c.dat[,"Center.X"],rdat$c.dat[,"Center.Y"],n=1,plot=F)
	}

}
#show traces from multiple collected.
PointTrace <- function(lookat,png=F,col=rep("black",nrow(lookat)),pch=16,cex=1,lmain="PointTrace",x.trt=NULL,y.trt=NULL,wr="wr1",t.names=NULL)
{
	if(!is.null(x.trt)){lookat["x"] <- lookat[,x.trt]}
	if(!is.null(y.trt)){lookat["y"] <- lookat[,y.trt]}	
	dev.new(height=4,width=14)
	rr.dev <- dev.cur()
	dev.new(height=4,width=4)
	plot(lookat[,"x"],lookat[,"y"],col=col,pch=pch,cex=cex,main=lmain,xlab=x.trt,ylab=y.trt)
	ret.list <- NULL
	i <- identify(lookat[,"x"],lookat[,"y"],n=1,plot=F)
	my.dev <- dev.cur()
	while(length(i) > 0)
	{
		x.names <- lookat[i,"trace.id"]
		#points(lookat[i,"x"],lookat[i,"y"],pch=8,cex=.5)
		rn.i <- row.names(lookat)[i]
		tmp <- get(lookat[i,"rd.name"])
		levs <- unique(tmp$w.dat[,"wr1"])
		lmain <- paste(i,lookat[i,"rd.name"])
		#LinesEvery(tmp$t.dat,,x.names,tmp$w.dat[,"wr1"],levs,lmain=lmain)
		dev.set(which=rr.dev)
		PeakFunc2(tmp$t.dat,x.names,3,30,TRUE,tmp$w.dat[,wr],lmain=lookat[i,"rd.name"])
		if(!is.null(t.names)){mtext(paste(t.names,tmp$c.dat[x.names,t.names],collapse=":"))}
		if(png==TRUE)
		{
			f.name <- paste(lookat[i,"rd.name"],lookat[i,"trace.id"],"png",sep=".")
			png(f.name,heigh=600,width=1200)
			PeakFunc2(tmp$t.dat,x.names,3,30,TRUE,tmp$w.dat[,wr],lmain=lookat[i,"rd.name"])
			dev.off()
		}

		dev.set(which=my.dev)
		if(is.element(rn.i,ret.list))
			{points(lookat[i,"x"],lookat[i,"y"],col=col[i],pch=pch,cex=cex);ret.list <- setdiff(ret.list,rn.i)}		
		else
			{points(lookat[i,"x"],lookat[i,"y"],col="red",pch=pch,cex=cex);ret.list <- union(rn.i,ret.list)}		
		i <- identify(lookat[,"x"],lookat[,"y"],n=1,plot=F)
	}
	return(ret.list)

}


TrimClust <- function(clust,min.n=5)
{
	clust.tab <- table(clust)
	x.names <- names(clust.tab[clust.tab < min.n])
	clust[is.element(clust,x.names)] <- 0
	return(clust)
}

ClusterAnalysis <- function(dat.name)
{
	library(cluster)
	tmp <- get(dat.name)

	if(!is.element("bin",names(tmp))){stop("no binary scoring (bin)")}
	if(!is.element("scp",names(tmp))){stop("no quantitative scoring (scp)")}
	if(!is.element("ROI.Area",names(tmp$c.dat))){stop("no ROI.Area (c.dat)")}
	x.names <- row.names(tmp$bin[tmp$bin[,"drop"]==0,])
	#three similarity scores (binary, trace correlation, scp)
	b.names <- names(tmp$bin)[names(tmp$bin)!="drop"]
	b.names <- select.list(b.names,multiple=T,title="Select responses for clustering")
	if(length(b.names) == 0){stop("no variables selected")}
	dist.b <- dist(tmp$bin[x.names,b.names],method="manhattan")
	q.names <- paste(b.names,".max",sep="")#,paste(b.names,".tot",sep=""))
	# dist.q <- dist(tmp$scp[x.names,q.names]) 
	# dist.c <- dist.b
	# cr.tot <- rep(1,length(dist.c))
	# for(i in b.names)
	# {
		# cr <- cor(tmp$t.dat[tmp$w.dat[,"wr1"]==i,x.names])
		# cr.tot <- cr.tot+(1-cr[lower.tri(cr)])
	# }
	# dist.c[] <- cr.tot
	# glt.b <- tmp$bin[x.names,b.names]
	# for(i in b.names)
	# {
		# g.dat <- as.matrix(t(tmp$t.dat[tmp$w.dat[,"wr1"]==i,x.names]))
		# glt <- glm(tmp$bin[x.names,i] ~ g.dat)
		# glt.b[i] <- predict(glt)
	# }
	# dist.d <- dist(glt.b,method="manhattan")
# #	return(cbind(dist.b,dist.q,dist.c,dist.d))
	#binary cluster
	h1 <- hclust(dist.b)
	clust1 <- cutree(h1,h=.1)
	clust1 <- TrimClust(clust1,5)
	res.dat <- cbind(tmp$scp[x.names,q.names],tmp$c.dat[x.names,-1])
	res.tab <- data.frame(n=tapply(res.dat[,1],clust1,length))
	for(i in names(res.dat))
	{
		glt <- glm(scale(res.dat[,i]) ~ as.factor(clust1))
		res.tab[i] <- summary(glt)$coefficients[,3]
	}
	return(res.tab)
	#plot all clusters > 5
	#table of means (t.values,z.values from glm) scale variable or ignore 0 group.
	#heat graph of group characterization.

}

#given a list of file names collect and merge all bin scp and c.dat data
CollectMulti <- function(f.names,rd.names=NULL)
{
	if(is.null(rd.names))
	{
		rd.names <- sub("\\.rdata$","",sub(".*\\/","",f.names),ignore.case=T)
		for(i in f.names){load(i)}
	}
	
	b.names <- NULL
	s.names <- NULL
	c.names <- NULL
	for(i in rd.names)
	{
		tmp <- get(i)
		names(tmp$bin) <- make.names(names(tmp$bin))
		names(tmp$scp) <- make.names(names(tmp$scp))		
		names(tmp$c.dat) <- make.names(names(tmp$c.dat))		
		b.names <- union(b.names,names(tmp$bin))
		s.names <- union(s.names,names(tmp$scp))
		c.names <- union(c.names,names(tmp$c.dat))
	}
	c.names <- setdiff(c.names,b.names)
	s.names <- setdiff(s.names,b.names)
	c.names <- setdiff(c.names,s.names)
	
	tot.names <- c(b.names,s.names,c.names,"rd.name","trace.id")
	ret.dat <- data.frame(matrix(rep(1,length(tot.names)),ncol=length(tot.names)))
	names(ret.dat) <- tot.names
	for(i in rd.names)
	{
		tmp <- get(i)
		names(tmp$bin) <- make.names(names(tmp$bin))
		names(tmp$scp) <- make.names(names(tmp$scp))		
		names(tmp$c.dat) <- make.names(names(tmp$c.dat))			
		ret.tmp <- data.frame(cbind(tmp$bin,tmp$scp,tmp$c.dat))
		ret.tmp["rd.name"] <- i
		ret.tmp["trace.id"] <- row.names(tmp$bin)
#		ret.dat <- merge(ret.dat,ret.tmp)
		i.names <- setdiff(tot.names,names(ret.tmp))

		for(j in i.names)
		{
			ret.tmp[j] <- NA
		}

		ret.add <- ret.tmp[,tot.names]
		ret.dat <- rbind(ret.dat,ret.add)
	}
	ret.dat <- ret.dat[-1,]
	return(ret.dat)	
}

#given a collected multi (a.tot) and a list of names to plot
#plot all in sets of subset.n
LinesEveryMulti <- function(a.tot,x.names,subset.n=10,rtag.name=NULL,sf=.6,name.tag="")
{
	rd.names <- unique(a.tot[x.names,"rd.name"])
	for(i in rd.names)
	{
		tmp <- get(i)
		y.names <- x.names[a.tot[x.names,"rd.name"]==i]
		yseq <- sort(rep(seq(1,ceiling(length(y.names)/subset.n)),length.out=length(y.names)))
		y.lst <- split(y.names,yseq)
		
		for(j in 1:length(y.lst))
		{
			lmain <- paste(i,name.tag,j,sep=".")
			z.names <- a.tot[y.lst[[j]],"trace.id"]
			pdf.name <- paste(i,name.tag,j,"pdf",sep=".")
			rtag <- NULL
			p.ht <- (length(z.names)/subset.n)*7+3
			
			if(!is.null(rtag.name)){rtag <- a.tot[y.lst[[j]],rtag.name]}
			LinesEvery(tmp$t.dat,,z.names,tmp$w.dat[,"wr1"],levs=unique(tmp$w.dat[,"wr1"])[-1],lmain=lmain,pdf.name=pdf.name,rtag=rtag,sf=sf,sig.line=.05,bcex=.65,p.ht=p.ht)
			##LinesEvery
		}
	}
	
}
	
	
	
	
