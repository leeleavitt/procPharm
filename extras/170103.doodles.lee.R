
#t.dat is the trace data with Time in the first column.
#should realy be baseline corrected data.(at least minimum 0)
#review one trace at a time sending up down drop or stay
#smash the selected up and down.
#color code regions of correlation to the up and down group.
#forget the layout stuff and just use one plot region.
#overprint in white to erase the previous trace.
TracePartition <- function(t.dat,m.names,wr,levs=NULL,ctype=factor(levels=c("drop","unk","class1","class2"),ordered=T))
{
	RndLine <- function(i,ly,hy,lcol){lines(xseq,t.dat[,i]+runif(1,ly,hy),col=lcol)}
		
	StepLine1 <- function(i){lines(xseq,(t.dat[,i]-min(t.dat[,i]))+4+(sum(ctype==c1)/(length(ctype)/2)),col="blue")}
	StepLine2 <- function(i){lines(xseq,(t.dat[,i]-min(t.dat[,i]))+2-(sum(ctype==c2)/(length(ctype)/2)),col="green")}
	
	sf <- 1 
    m.names <- intersect(m.names,names(t.dat))
    lwds <- 2
    if(length(m.names) == 0)
    {stop("no named traces exist in trace dataframe.")}
    init.plot=T
	if(length(ctype) < length(m.names)){ctype=factor(rep("unk",length(m.names)),levels=c("drop","unk","class1","class2"),ordered=T)}
	names(ctype) <- m.names
    xseq <- t.dat[,1]

	#levels of the ctype factor must be "drop < unk < class1 < class2" 
	#the last two can be user defined.
	if(length(levels(ctype))!=4){stop("wrong number of levels for ctype")}
	if(levels(ctype)[1] != "drop"){stop("wrong levels for ctype")}
	if(levels(ctype)[2] != "unk"){stop("wrong levels for ctype")}

	c1 <- levels(ctype)[3]
	c2 <- levels(ctype)[4]


	c.cnt <- summary(ctype)
    dev.new(width=14,height=8)
    par(mar=c(0,0,0,0),xaxt="n",yaxt="n")
	plot(xseq,rep(1,length(xseq)),type="n",ylim=c(-.5,7))
	
    if(length(wr) > 0)
    {
    	if(is.null(levs)){levs <- setdiff(unique(wr),"")}
        x1s <- tapply(xseq,as.factor(wr),min)[levs]
        x2s <- tapply(xseq,as.factor(wr),max)[levs]
        y1s <- rep(-.3,length(x1s))
        y2s <- rep(.1,length(x1s))
        rect(x1s,y1s,x2s,y2s,col="lightgrey",border=T)
        text(xseq[match(levs,wr)],rep(-.1,length(levs)),levs,pos=4,offset=0,cex=1)
    }
    text(0,6.5,c1,col="blue")
    text(0,.5,c2,col="green")    

    p.names <- c("Drop",c1,c2,"NEXT","DONE")
    xs <- rep(0,length(p.names))+c(0,2,2,4,0)/2
    ys <- c(3.5,3.7,3.3,3.5,3)-.5
    p.cols <- c("red","blue","green","black","black")
    for(i in sample(names(ctype[ctype!="drop"]),10,replace=T)){RndLine(i,2.5,3,"lightgrey")}
	for(i in names(ctype[ctype==c1])) {RndLine(i,4,6,"blue")}
	for(i in names(ctype[ctype==c2])) {RndLine(i,0,2,"green")}	
	
	for(i in names(ctype[ctype!="drop"]))
	{
		yseq <- t.dat[,i]
		lines(xseq,yseq+3,lwd=lwds,col="black")
		points(xs,ys,col=p.cols,cex=2,pch=16)
		click.i <- identify(xs,ys,n=1,plot=F)
		if(click.i == 1){ctype[i] <- "drop"}
		if(click.i == 2){ctype[i] <- c1;StepLine1(i)}
		if(click.i == 3){ctype[i] <- c2;StepLine2(i)}
		if(click.i == 4){ctype[i] <- "unk"}
		if(click.i == 5){dev.off();return(ctype)}
		lines(xseq,yseq+3,lwd=lwds*1.5,col="white")		
	}
	dev.off()
	return(ctype)
	
}


#a.tot has rd.name, trace.id, ctype.
#all RDs must be loaded. ctype must be a factor compatilble with TracePartition.
#usage a.tot["ctype"] <- factor(rep("unk",nrow(a.tot)),levels=c("drop","unk","class1","class2"),ordered=T)
#rename class one and class two as you see fit.
MultiType <- function(a.tot)
{
	for(i in unique(a.tot[,"rd.name"]))
	{
		tmp <- get(i)
		levs <- unique(tmp$w.dat[,"wr1"])[-1]
		ct <- a.tot[a.tot[,"rd.name"]==i,"ctype"]
		x.names <- a.tot[a.tot[,"rd.name"]==i,"trace.id"]
		names(ct) <- x.names
		pt1 <- TracePartition(tmp$t.dat,x.names,tmp$w.dat[,"wr1"],levs,ctype=ct)
		a.tot[a.tot[,"rd.name"]==i,"ctype"] <- pt1
	}
	return(a.tot)
}

#remove single point spikes that go up then back down
#or down then back up.
#wt is $t.dat from RD 
#zthresh is the z transformed threshold to be exceeded in both directions
#gthresh is the absolute change in ratio.
SpikeTrim <- function(wt,zthresh=NULL,gthresh=1)
{
	if(is.null(zthresh)){zthresh <- abs(qnorm(1/(nrow(wt)*(ncol(wt)-1))))/2}
	
	wtd <- wt[-1,]-wt[-nrow(wt),]
	wtd <- sweep(wtd[,-1],1,wtd[,1],'/')
	#last delta value is missing or unkown (log transform?)
	#edges of blanks may have response
	mn.cnt <- apply(wtd,2,mean)
	sd.cnt <- apply(wtd,2,sd)
	#z-scores
	wtz <- sweep(wtd,2,mn.cnt,'-')
	wtz <- sweep(wtz,2,sd.cnt,'/')
	wtp <- wtz > zthresh | wtd > gthresh #significant increase	
	wtn <- wtz < -zthresh | wtd < -gthresh #significant decrease
	wtrm <- wtp[-1,]+wtn[-nrow(wtn),] #remove all 2s
	wtrm2 <- wtn[-1,]+wtp[-nrow(wtp),] #remove all 2
	wtrm <- wtrm==2 | wtrm2==2
	#now impute missing.
	rm.cnt <- apply(wtrm,2,sum)
	mp.names <- names(rm.cnt)[rm.cnt>0]
	
	for(i in mp.names)
	{
		na.x <- row.names(wtrm)[-nrow(wtrm)][wtrm[-1,i]]
		wt[na.x,i] <- NA
		xspan <- min(20/nrow(wt),.5)
		xloe <- loess(substitute(i ~ Time,list(i=as.name(i))),span=xspan,data=wt)
		loepred <- predict(xloe,newdata=wt[,1])
		names(loepred) <- row.names(wt)
		wt[na.x,i] <- loepred[na.x]
	}
	return(wt)
}	

SpikeTrim2 <- function(wt,ulim=NULL,dlim= NULL)
{
	
	wtd <- wt[-1,]-wt[-nrow(wt),]
	wtd <- sweep(wtd[,-1],1,wtd[,1],'/')
	if(is.null(ulim) | is.null(dlim))
	{
		qvals <- quantile(as.vector(as.matrix(wtd)),probs=c(0,.01,.5,.99,1))
	}
	if(is.null(dlim)){dlim <- qvals[2]}
	if(is.null(ulim)){ulim <- qvals[4]}	
	wt.up <- wtd > ulim
	wt.dn <- wtd < dlim
	wt.ud <- wt.up[-nrow(wt.up),] + wt.dn[-1,]
	wt.du <- wt.up[-1,] + wt.dn[-nrow(wt.dn),]
	wt.na <- wt[2:(nrow(wt)-1),-1]
	wt.na[wt.ud==2] <- NA
	wt.na[wt.du==2] <- NA	
	sum(is.na(wt.na))
	wt[2:(nrow(wt)-1),-1] <- wt.na

	#impute missing using mean of flanking.
	#consider replicating first and last columns and doing this all as a vector
	
	return(wt)
}

#each point is replaced with the mean of the two neighboring points
Mean3 <- function(wt)
{
	wt.mn <- (wt[-c(1,2),]+wt[-c(nrow(wt),(nrow(wt)-1)),])/2
	wt[2:(nrow(wt)-1),] <- wt.mn
	return(wt)
}
#tmp is an RD file.  levs are the levels to refine (defaults to all)
#maxws is the maximum window size in minutes.
#padL and padR are padding to the left and right.
WindowRefine <- function(tmp,levs=NULL,maxws=1,wr="wr1",zthresh=2,cthresh=NULL,padL=0,padR=0)
{
	if(is.null(levs)){levs <- unique(tmp$w.dat[,wr]);levs <- levs[levs!=""]}
	if(is.null(cthresh)){cthresh <- max(nrow(tmp$c.dat)*.03,2)}
	if(is.element("bin",names(tmp))){x.names <- row.names(tmp$bin[tmp$bin[,"drop"]==0,])}else{x.names <- names(tmp$t.dat)[-1]}
	#define null response
	wr.trim <- tmp$w.dat[,wr]
	if(padL != 0 | padR != 0)
	{
		wr.out <- wr.trim
		wr.out[] <- ""
		for(i in levs)
		{
			wL <-  max(min((1:length(wr.trim))[wr.trim==i])-padL,1)
			wR <- min(max((1:length(wr.trim))[wr.trim==i])+padR,length(wr.trim))
			wr.out[wL:wR] <- i
		}
		wr.trim <- wr.out
	}
	wr2 <- NumBlanks(wr.trim)
	wr2 <- wr2[-length(wr2)]
	b.levs <- grep("blank",unique(wr2),value=T)
	if(length(b.levs) == 0){stop("No blank regions to establish null threshold")}
	#each cell will have an delta increase that signals a non random response
	#This is a CA increase only seen in 2.5% of the time in blank regions.
	#the first timeframe in the window for which > 4 cells exceed this is the 
	#putative start of response to the challenge.
	#each cell has a mean and sd within each window.
	wt <- tmp$t.dat
	wtd <- wt[-1,]-wt[-nrow(wt),]
	wtd <- sweep(wtd[,-1],1,wtd[,1],'/')
	#last delta value is missing or unkown (log transform?)
	#edges of blanks may have response
	mn.cnt <- apply(wtd[is.element(wr2,b.levs),],2,mean)
	sd.cnt <- apply(wtd[is.element(wr2,b.levs),],2,sd)
	#z-scores
	wtz <- sweep(wtd,2,mn.cnt,'-')
	wtz <- sweep(wtz,2,sd.cnt,'/')
	wtt <- wtz > zthresh	
	for(i in levs)
	{
		t.cnt <- apply(wtt[wr2==i,x.names],1,sum)
		if(sum(t.cnt > cthresh) ==0)
		{warning(paste("threhsold for leading edge not met for",i))}
		else
		{
			wr.start <- wt[wr.trim==i,1][t.cnt > cthresh][1] #this is the first time point for which 5 cells exceed the 97.5 threshold
			wr.end <- which.min(abs(wt[,1]-(wr.start + maxws)))
			wr.start <- match(wr.start,wt[,1])
			wr.trim[wr.trim==i] <- ""
			wr.trim[wr.start:wr.end] <- i
		}
	}
	return(wr.trim)	
}


TypeReview <- function(a.tot)
{
rs.levs <- names(a.tot)[c(1:4,15)]	
for(my.type in c("CIMpAm","CIMpAp","CSMmAm","HTMpAm","HTMpAp","LTMpAm","MTMpAp"))	
{
for(i in rd.list)
{
		tmp <- get(i)
		tmp$w.dat["wr2"] <- tmp$w.dat[,"wr1"]
		tmp.wr <- WindowRefine(tmp,rs.levs,.8,"wr2",zthresh=2)		
		tmp$w.dat["wr1"] <- tmp.wr
		x.names <- a.tot[a.tot[,"rd.name"]==i & a.tot[,"ctype3"]==my.type,"trace.id"]
		if(length(x.names) > 0)
		{
		rtag <- round(a.tot[a.tot[,"rd.name"]==i & a.tot[,"ctype3"]==my.type,"ROI.Area"],0)

	z.names <- TraceSelectLarge(tmp$t.dat,,x.names,tmp$w.dat[,"wr1"],unique(tmp$w.dat[,"wr1"])[-1],lmain=paste("select to remove from",my.type,i),rtag=rtag)
	if(length(z.names) > 0)
	{
		a.tot[a.tot[,"rd.name"]==i & is.element(a.tot[,"trace.id"],z.names),"ctype3"] <- "unk"	
		x.names <- setdiff(x.names,z.names)
	}
	
	}
}
}	
	return(a.tot)
	
}


getCommonBin <-function(rd.list)
{
	tmp <- get(rd.list[1])
	bn.tot <- names(tmp$bin)
	for(i in rd.list){tmp <- get(i);bn.tot <- intersect(bn.tot,names(tmp$bin))}
	return(bn.tot)
}

#maybe use loess
#step 1 loess estimates of all cells on a common time frame
#step 2 shift ex1 with respect to ex2 until max correlation
#step 3 repeat with all other rds.
#Time must be the first column
fitPoly <- function(tdat,dgre = 3)
{
	lmpolyfunc <- function(x,y){lt <- lm(y ~ poly(x,dgre));lt.sum <- summary(lt);return(lt.sum$coefficients[,1])}
	res.tab <- data.frame(mean=apply(tdat[,-1],2,mean))
	res.tab["sd"] <- apply(tdat[,-1],2,sd)

	
}	
#create the best sync'd polynomial fit across
#all named levels align the wr across all rds
PolySync <- function(a.tot,levs=NULL,wr="wr1")
{
	if(!is.element("rd.name",names(a.tot))){stop("rd.name required for PolySync")}
	if(is.null(levs)){levs <- getCommonBin(unique(a.tot[,"rd.name"]))}	
	ply.stat <- fitPoly(tmp$t.dat[tmp$w.dat[,wr]==lv,])
		
		
}

#find the best alignment of two trace sets x and y
#and loess fit into targ time slots
#targ must a subset of x[,1],y[,1] #time var for traces
#returns may be smaller due to shift and align


TraceAlign <- function(x,y,targ=NULL,x.names=NULL,y.names=NULL)
{
	matdist <- function(rx,ry)
	{
		xydist <- dist(rbind(t(xt[rx,-1]),t(yt[ry,-1])))
		return(mean(as.vector(as.matrix(xydist)[xi,yi]))/length(rx))
	}		
	if(is.null(targ))
	{
		x[,1] <- x[,1]-min(x[,1])
		y[,1] <- y[,1]-min(y[,1])
		mint <- min(max(x[,1]),max(y[,1]))
		targ <- data.frame(Time=seq(0,mint,by=1/120))
	}
	#impute into the targ time frame
	#both data sets will end up on the same time scale
	#with the same number of points.
	if(sum(x[,"Time"] != targ[1:nrow(x),"Time"]) > 0)
	{
	xt <- targ
	for(i in names(x[,-1]))
	{
		xspan <- 5/nrow(x)
		xloe <- loess(substitute(i ~ Time,list(i=as.name(i))),span=xspan,data=x)
		xt[i] <- predict(xloe,newdata=targ)
		xt[i] <- xt[,i]-min(xt[,i])
	}
	}
	else
	{
		xt <- x
	}
	yt <- targ
	for(i in names(y[,-1]))
	{
		xspan <- 5/nrow(y)
		xloe <- loess(substitute(i ~ Time,list(i=as.name(i))),span=xspan,data=y)
		yt[i] <- predict(xloe,newdata=targ)
		yt[i] <- yt[,i]-min(yt[,i])
	}
	#minimize distance between two sets
	#let's try using ccf on the delta
	xtd <- xt[-1,]-xt[-nrow(xt),]
	xtd <- sweep(xtd[,-1],1,xtd[,1],'/')	
	xtd[xtd<0] <- 0
	sdx.cnt <- apply(xtd[,-1],2,sd)
	if(is.null(x.names)){x.names <- names(sdx.cnt[sdx.cnt > quantile(sdx.cnt,.9)[[1]]])}
	
	ytd <- yt[-1,]-yt[-nrow(yt),]
	ytd <- sweep(ytd[,-1],1,ytd[,1],'/')
	ytd[ytd<0] <- 0
	sdy.cnt <- apply(ytd[,-1],2,sd)
	if(is.null(y.names)){y.names <- names(sdy.cnt[sdy.cnt > quantile(sdy.cnt,.9)[[1]]])}
	ccf.mat <- ccf(xtd[,x.names[1]],ytd[,y.names[1]],plot=F,lagmax=20)$acf
	for(i in x.names){for(j in y.names){ccf.mat <- cbind(ccf.mat,ccf(xtd[,i],ytd[,j],plot=F)$acf)}}

	si <- 20-which.max(apply(ccf.mat,1,mean))
	xi <- seq(1,nrow(xt))
	yi <- seq(1,nrow(yt))
	if(si > 0)
	{
		yi <- seq((si+1),nrow(yt))
		xi <- seq(1,(nrow(xt)-si))
	}
	if(si < 0)
	{
		xi <- seq((-si+1),nrow(xt))
		yi <- seq(1,(nrow(yt)+si))		
	}	
	xt <- xt[xi,]
	yt <- yt[yi,]
	xt[,1] <- targ[1:nrow(xt),"Time"]
	yt[,1] <- targ[1:nrow(yt),"Time"]	
	xtmean <- apply(xt[,x.names],1,mean)
	ytmean <- apply(yt[,y.names],1,mean)
	return(list(ccf.mat = ccf.mat,xt = xt,yt=yt,x.names = x.names,y.names=y.names))
	#now how to get all ccfs between xtd and ytd
		
# # 	xi <- seq(1,(ncol(xt)-1))
	# yi <- max(xi) + seq(1,(ncol(yt)-1))	
	# rx <- seq(1,nrow(xt))
	# ry <- seq(1,nrow(xt))
	# wiggle.dist <- data.frame(xi=1,yi=1,xydist=matdist(rx,ry))
	# lp <- 1
	# for(i in 1:wiggle)
	# {
		# rx <- seq(1,nrow(xt)-i)
		# ry <- seq(1+i,nrow(xt))
		# lp <- lp+1
		# wiggle.dist[lp,] <- c(1,1+i,matdist(rx,ry))
		# ry <- seq(1,nrow(xt)-i)
		# rx <- seq(1+i,nrow(xt))
		# lp <- lp+1
		# wiggle.dist[lp,] <- c(1+i,1,matdist(rx,ry))		
	# }	
	# wi <- which.min(wiggle.dist[,3])
	# i <- wiggle.dist[wi,1]
	# j <- wiggle.dist[wi,2]
	# xt <- xt[i:(nrow(xt)-(j-1)),]
	# yt <- yt[j:(nrow(yt)-(i-1)),]
	# xt[,1] <- xt[,1]-min(xt[,1])
	# yt[,1] <- yt[,1]-min(yt[,1])
	# return(list(wiggle=wiggle.dist,xt=xt,yt=yt))
}

#impute the x trace into the timeframe defined by targ
TraceImpute <- function(x,targ=NULL)
{
	if(is.null(targ))
	{
		x[,1] <- x[,1]-min(x[,1])
		mint <- max(x[,1])
		targ <- data.frame(Time=seq(0,mint,by=1/120))
	}
	#impute into the targ time frame
	xt <- targ
	for(i in names(x[,-1]))
	{
		xspan <- 5/nrow(x)
		xloe <- loess(substitute(i ~ Time,list(i=as.name(i))),span=xspan,data=x)
		xt[i] <- predict(xloe,newdata=targ)
		xt[i] <- xt[,i]-min(xt[,i])
	}
	return(xt)
}

keep2func <- function(x)
{
	ox <- sort(x,decreasing=T)[2]
	x[x<ox | is.na(x)] <- 0
	return(x)
}

TraceImpute.2 <- function(x,ts=seq(0,(length(x)-1))*.03,xspan=5/length(x))
{
	targ <- data.frame(ts=seq(min(ts),max(ts),by=1/120))
	
	xloe <- loess(x ~ ts,span=xspan)
	xp <- predict(xloe,newdata=targ)
	return(xp)
}

#baseline disturbances
#compute the maximum correlation shifting one against the other
#return the start indices and the cor value and pvalue

ShiftCor <- function(x,y)
{
	x.n <- length(x)
	y.n <- length(y)
	mid.x <- as.integer(x.n/2)	
	mid.y <- as.integer(y.n/2)
	lp <- min(mid.x,mid.y)
	rr <- 1:min(x.n,y.n)
	cor.t <- cor.test(x[rr],y[rr])
	cor.tab <- data.frame(x.i = 1,y.i=1,r=cor.t$estimate,p=cor.t$p.value)
	cnt.i <- 1
	for(i in 1:lp)
	{
		r1 <- i:min((x.n-1),y.n)
		r2 <- 1:length(r1)
		cor.t <- cor.test(x[r1],y[r2])
		cnt.i <- cnt.i+1
		cor.tab[cnt.i,] <- c(x.i = i,y.i=1,r=cor.t$estimate,p=cor.t$p.value)		
		r2 <- i:min(x.n,(y.n-1))
		r1 <- 1:length(r2)
		cor.t <- cor.test(x[r1],y[r2])
		cnt.i <- cnt.i+1
		cor.tab[cnt.i,] <- c(x.i = 1,y.i=i,r=cor.t$estimate,p=cor.t$p.value)		
		
	}
	return(cor.tab)
		
}



#must have a criteria for predicting 
#Should use raw trace data in specifically defined windows.
#polynomial fits within these windows should be generated before hand
#coefficients of the polynomials 
#the window regions must be sync'd  SHOULD DO THAT FIRST
#review typing from multiple experiments.
ScoreReviewMulti <- function(lookat,predvar,p.names)
{	
	PredFunc <- function()
	{
		glt <- glm(lookat[,predvar] ~ as.matrix(lookat[,p.names]),weights=(1-lookat[,"drop"]),family="binomial")
		preds <- predict(glt)
#		preds <- log((preds-min(preds))+1)
		z <- 1
		pmean <- mean(preds)
		psd <- sd(preds)
		ptail <- (preds-pmean)/psd
		x1 <- (rank(preds[ptail < -z])/sum(ptail < -z))-(z+1)
		x1 <- pmean + x1*psd
		preds[ptail < -z] <- x1
		preds[ptail > z] <- ((rank(preds[ptail > z])/sum(ptail > z))+z)*psd + pmean
#		preds <- log((preds-min(preds))+1)
#		preds <- rank(preds) + preds
#		preds <- log((preds-min(preds))+1)
		preds <- (preds-min(preds))/(max(preds)-min(preds))
		return(preds)
	}
	lookat["x"] <- PredFunc()
	lookat["y"] <- jitter(as.integer(as.factor(lookat[,"rd.name"])),amount=.4)
	lookat["col"] <- "black"
	lookat[lookat[,predvar]==1,"col"] <- "red"
	lookat["pch"] <- 16
	lookat[lookat[,"drop"]==1,"pch"] <- 1
	lookat["selected"] <- 0
	buttons.x <- rep(-.02,5)
	buttons.y <- seq(from=min(lookat[,"y"]),to=max(lookat[,"y"])/2,length.out=5)
	buttons.col <- rep("green",5)
	buttons.pch <- rep(16,5)
	buttons.text <- c("Finish","Drop","Toggle","Recalc","Clear")
	
	dev.new(height=4,width=14)
	rr.dev <- dev.cur()
	dev.new(height=4,width=14)
	xs <- c(lookat[,"x"],buttons.x)
	ys <- c(lookat[,"y"],buttons.y)
	cols <- c(lookat[,"col"],buttons.col)
	pchs <- c(lookat[,"pch"],buttons.pch)
	plot(xs,ys,col=cols,pch=pchs,main=predvar,xlab="",ylab="")
	text(buttons.x,buttons.y,buttons.text,pos=2,cex=.5)
	i <- identify(xs,ys,n=1,plot=F)
	my.dev <- dev.cur()
	while(length(i) > 0)
	{
		if(i > nrow(lookat))#command button press
		{
			i <- i-nrow(lookat)
			if(i == 1)
			{
				return(lookat)
			}
			if(i ==2)
			{
				xi <- row.names(lookat)[lookat[,"selected"]==1]
				if(length(xi) > 0)
				{
					lookat[xi,"drop"] <- 1
					lookat[xi,"pch"] <- 1
					points(lookat[xi,"x"],lookat[xi,"y"],pch=8,cex=1.3,col="white")
					points(lookat[xi,"x"],lookat[xi,"y"],pch=lookat[xi,"pch"],col=lookat[xi,"col"])
					lookat[xi,"selected"] <- 0
				}					
			}
			if(i == 3)
			{
				xi <- row.names(lookat)[lookat[,"selected"]==1]
				if(length(xi) > 0)
				{
					x0 <- xi[lookat[xi,predvar]==0]
					x1 <- xi[lookat[xi,predvar]==1]
					if(length(x0) > 0)
					{
					lookat[x0,predvar] <- 1
					lookat[x0,"col"] <- "red"
					}
					if(length(x1) > 0)
					{
					lookat[x1,predvar] <- 0
					lookat[x1,"col"] <- "black"
					}
					points(lookat[xi,"x"],lookat[xi,"y"],pch=8,cex=1.3,col="white")
					points(lookat[xi,"x"],lookat[xi,"y"],pch=lookat[xi,"pch"],col=lookat[xi,"col"])
					lookat[xi,"selected"] <- 0					
				}									
			}			
			if(i == 4)
			{
				xi <- row.names(lookat)[lookat[,"selected"]==1]
				if(length(xi) > 0)
				{
					points(lookat[xi,"x"],lookat[xi,"y"],pch=8,cex=1.4,col="white")
					points(lookat[xi,"x"],lookat[xi,"y"],pch=lookat[xi,"pch"],col=lookat[xi,"col"])
					lookat[xi,"selected"] <- 0
				}
				lookat["x"] <- PredFunc()
				xs <- c(lookat[,"x"],buttons.x)
				ys <- c(lookat[,"y"],buttons.y)
				cols <- c(lookat[,"col"],buttons.col)
				pchs <- c(lookat[,"pch"],buttons.pch)
				plot(xs,ys,col=cols,pch=pchs,main=predvar,xlab="",ylab="")
				text(buttons.x,buttons.y,buttons.text,pos=2,cex=.5)
			}
			if(i == 5)
			{
				xi <- row.names(lookat)[lookat[,"selected"]==1]
				if(length(xi) > 0)
				{
					points(lookat[xi,"x"],lookat[xi,"y"],pch=8,cex=1.3,col="white")
					points(lookat[xi,"x"],lookat[xi,"y"],pch=lookat[xi,"pch"],col=lookat[xi,"col"])
					lookat[xi,"selected"] <- 0
				}					
				
			}
			i <- identify(xs,ys,n=1,plot=F)							
		}
		else
		{
			#tag selected and show trace TOGGLE
			if(lookat[i,"selected"]==1)
			{
				points(lookat[i,"x"],lookat[i,"y"],pch=8,cex=1.1,col="white")				
				lookat[i,"selected"]<-0
				points(lookat[i,"x"],lookat[i,"y"],pch=lookat[i,"pch"],col=lookat[i,"col"])
			}
			else
			{
				lookat[i,"selected"]<-1
				points(lookat[i,"x"],lookat[i,"y"],pch=8,col="blue")				
			}	
			
		x.names <- lookat[i,"trace.id"]
		tmp <- get(lookat[i,"rd.name"])
		levs <- unique(tmp$w.dat[,"wr1"])
		lmain <- paste(i,lookat[i,"rd.name"])
		dev.set(which=rr.dev)
		PeakFunc2(tmp$t.dat,x.names,3,30,TRUE,tmp$w.dat[,"wr1"],lmain=lookat[i,"rd.name"])
		dev.set(which=my.dev)
		i <- identify(xs,ys,n=1,plot=F)
		}
	}
}

#jitter scatter with marked medians
#x is a factor y quantitative

#resample to get CI for median
MedianBS <- function(x,CI=.95,nrs=1000)
{
	x <- x[!is.na(x)]
	x.md <- median(x)
	bx.md <- sort(replicate(nrs,median(sample(x,length(x),replace=T))))
	low.i <- 1
	high.i <- 1
	ex <- (1-CI)/2
	low.i <- floor(nrs*ex)
	high.i <- ceiling(nrs*(1-ex))
	if(low.i < 1){low.i <- 1;warning("increase nrs in MedianBS")}
	if(high.i > nrs){high.i <- nrs;warning("increase nrs in MedianBS")}
	return(list(lowCI=bx.md[low.i],med=median(bx.md),highCI=bx.md[high.i]))
}

#calculate means and sems for all c.names of dat
#divided by the levels of fac.name
#make a bargraph of these
#augment this to indicate magnitudes of effect.
MeanSemGraph <- function(dat,c.names,fac.name,t.cols=NULL,ylab=NULL,main.lab=NULL,x.labs=NULL,bt=.1,lgc="topleft",ylim=NULL,den=NULL,cex.lab=1.25,yaxt="s")
{
	semfunc <- function(x)
	{
		n <- sum(!is.na(x))
		if(n < 3){return(NA)}
		return(sd(x,na.rm=T)/sqrt(n))
	}
	if(is.null(ylab)){ylab <- "Response Frequency"}
	x <- as.factor(dat[,fac.name])
	x.levs <- levels(x)
	if(1/length(x.levs) < bt){bt <- 1/length(x.levs)}
	sem.levs <- paste(x.levs,"sem",sep=".")
	x.res <- data.frame(apply(dat[x==x.levs[1],c.names,drop=F],2,mean,na.rm=T))
	for(i in x.levs)
	{
		x.res[i] <- apply(dat[x==i,c.names,drop=F],2,mean,na.rm=T)
		x.res[paste(i,"sem",sep=".")] <- apply(dat[x==i,c.names,drop=F],2,semfunc)
	}
	xlim <- c(1,length(c.names)+length(x.levs)*bt)
	if(is.null(ylim)){ylim <- c(0,max(x.res[,x.levs]+x.res[,sem.levs]*2)*1.2)}
	
	if(is.null(t.cols)){t.cols <- rainbow(length(x.levs));names(t.cols) <- x.levs}
	plot(x.res[,x.levs[1]],xlim=xlim,ylim=ylim,type="n",xaxt="n",xlab="",ylab=ylab,main=main.lab,yaxt=yaxt)
	for(i in 1:length(x.levs))
	{
		x1 <- seq(1,length(c.names))+(i-1)*bt
		y1 <- x.res[,x.levs[i]]
		rect(x1,rep(0,length(x1)),x1+bt,y1,col=t.cols[x.levs[i]],density=den[x.levs[i]])
	}
	for(i in 1:length(x.levs))
	{
		x1 <- seq(1,length(c.names))+(i-1)*bt+(bt)/2
		y1 <- x.res[,x.levs[i]] + x.res[,sem.levs[i]]*2
		y2 <- x.res[,x.levs[i]] - x.res[,sem.levs[i]]*2		
		arrows(x1,y2,x1,y1,angle=90,col="black",length=bt*.25,code=3)
	}
	
	if(is.null(x.labs)){x.labs <- row.names(x.res)}
	center.off <- (length(x.levs)*bt)/2 
#	text(seq(1,(length(x.labs)))+center.off,rep(-.02,length(x.labs)),x.labs,pos=4,cex=1,offset= -nchar(x.labs)/10)
	mtext(x.labs,at=seq(1,(length(x.labs)))+center.off,cex=cex.lab,side=1,line=1)
	if(!is.na(lgc)){legend(lgc,fill=t.cols,names(t.cols),cex=cex.lab,density=den[x.levs])}
	
	return(x.res[,-1])
}

MeanSemGraph2 <- function(dat,c.names,fac.name,t.cols=NULL,ylab=NULL,main.lab=NULL,x.labs=NULL,bt=.1,lgc="topleft",ylim=NULL,den=NULL,cex.lab=1.25,yaxt="s",xgap=.01)
{
	semfunc <- function(x)
	{
		n <- sum(!is.na(x))
		if(n < 3){return(NA)}
		return(sd(x,na.rm=T)/sqrt(n))
	}
	if(is.null(ylab)){ylab <- "Response Frequency"}
	x <- as.factor(dat[,fac.name])
	x.levs <- levels(x)
	if(1/length(x.levs) < bt){bt <- 1/length(x.levs)}
	sem.levs <- paste(x.levs,"sem",sep=".")
	x.res <- data.frame(apply(dat[x==x.levs[1],c.names,drop=F],2,mean,na.rm=T))
	for(i in x.levs)
	{
		x.res[i] <- apply(dat[x==i,c.names,drop=F],2,mean,na.rm=T)
		x.res[paste(i,"sem",sep=".")] <- apply(dat[x==i,c.names,drop=F],2,semfunc)
	}
	xlim <- c(1,length(c.names)+length(x.levs)*bt)
	if(is.null(ylim)){ylim <- c(0,max(x.res[,x.levs]+x.res[,sem.levs]*2)*1.2)}
	
	if(is.null(t.cols)){t.cols <- rainbow(length(x.levs));names(t.cols) <- x.levs}
	plot(x.res[,x.levs[1]],xlim=xlim,ylim=ylim,type="n",xaxt="n",xlab="",ylab="",main=main.lab,yaxt="n",cex.main=cex.lab)
	axis(2,at=seq(ylim[1],ylim[2],by=.1),labels=paste(seq(ylim[1],ylim[2],by=.1)*100,"%",sep=""),las=T,font=2,cex.axis=cex.lab*.75)	
	for(hx in seq(ylim[1],ylim[2],by=.1)){abline(h=hx,col="lightgrey",lty=2)}
	for(i in 1:length(x.levs))
	{
		x1 <- seq(1,length(c.names))+(i-1)*bt
		y1 <- x.res[,x.levs[i]]
		if(is.na(den[x.levs[i]]))
		{
			rect(x1,rep(0,length(x1)),x1+(bt-xgap),y1,col=t.cols[x.levs[i]])
		}
		else
		{
			rect(x1,rep(0,length(x1)),x1+(bt-xgap),y1,col="white")			
			rect(x1,rep(0,length(x1)),x1+(bt-xgap),y1,col=t.cols[x.levs[i]],density=den[x.levs[i]])			
		}	
	}
	for(i in 1:length(x.levs))
	{
		x1 <- seq(1,length(c.names))+(i-1)*bt+(bt)/2
		y1 <- x.res[,x.levs[i]] + x.res[,sem.levs[i]]*2
		y2 <- x.res[,x.levs[i]] - x.res[,sem.levs[i]]*2		
		arrows(x1-xgap/2,y2,x1-xgap/2,y1,angle=90,col="black",length=bt*.25,code=3)
	}
	
	if(is.null(x.labs)){x.labs <- row.names(x.res)}
	center.off <- (length(x.levs)*bt)/2 
#	text(seq(1,(length(x.labs)))+center.off,rep(-.02,length(x.labs)),x.labs,pos=4,cex=1,offset= -nchar(x.labs)/10)
	mtext(x.labs,at=seq(1,(length(x.labs)))+center.off-xgap/2,cex=cex.lab,side=1,line=1)
	if(!is.na(lgc))
	{
		lg.tab <- table(dat[,fac.name])
		leg.text <- paste(names(lg.tab),"\n","(n=",lg.tab, ")",sep="")
		legend(lgc,fill=t.cols[x.levs],legend=leg.text,cex=cex.lab,density=den[x.levs]*2,ncol=length(lg.tab),bg="white")
	}

	return(x.res[,-1])
}



#two columns for each 
#col 1 is freq/mean
#col 2 is sem
BarChartSE <- function(dat,c.names=NULL,bt=1)
{
	if(is.null(c.names)){c.names <- row.names(dat)}
	se1 <- dat[,1]-2*dat[,2]
	se2 <- dat[,1]+2*dat[,2]
	y1 <- dat[,1]
#	x1 <- seq(1:nrow(dat))
	x1 <- seq(1,length(c.names))+(seq(1,length(c.names))-1)*bt+(bt)/2	
	barplot(y1,width=bt)
	for(i in 1:length(y1))
	{
		arrows(x1,se1,x1,se2,angle=90,code=3)
	}

}

#borrowed from R-blogger
scatterhist = function(x, y, xlab="", ylab=""){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  plot(x,y)
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
    at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
    at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}


#Cheenu's suggested binning no NAs
Size.Bin2 <- function(x)
{
	x.bin <- rep("S",length(x))
	x.bin[x > 350] <- "M"
	x.bin[x > 500] <- "L"
	x.bin[x > 650] <- "H"	
	return(x.bin) 
}
#bin the staining data 
#based on Kp and log transformed mean intensity.
Stain.Bin <- function(tmp,stain.name=NULL,k.name=NULL,plotit=F,k.n=3)
{
	require(cluster)
	if(is.null(stain.name)){stain.name <- grep("tritc",names(tmp$c.dat),ignore.case=T,value=T)[1]}

	stain.val <- tmp$c.dat[,stain.name]
	stain.val <- log((stain.val-min(stain.val)+1))
	pm1 <- pam(stain.val,k=k.n)
	pclust <- pm1$clustering
	pmed <- tapply(stain.val,pclust,median)
	pmed[] <- rank(pmed)
	ret.val <- pmed[pclust]	
	if(plotit==T)
	{
		pcol <- pmed
		pcol[] <- rainbow(length(pcol))
		plot(stain.val,rank(stain.val)/length(stain.val),type="n",xlab=stain.name,ylab="rank",main=paste("Stain.Bin for",stain.name))
		for(i in 1:k.n)
		{
			x <- stain.val[ret.val==i]
			points(x,(rank(x)-1)/(length(x)-1),pch=16,col=pcol[i])
		}
		legend("bottomright",names(pcol),pch=16,col=pcol)
	}
	return(ret.val)	
}


#write a function that plots all the deltas and let's you click to adjust the points.
WindowClick <- function(tmp,wr.i="wr1",plevs=NULL)
{
	if(is.null(plevs)){plevs <- unique(tmp$w.dat[,wr.i]);plevs <- plevs[plevs != ""]}
	wt <- tmp$t.dat
	wtd <- wt[-1,]-wt[-nrow(wt),]
	wtd <- sweep(wtd[,-1],1,wtd[,1],'/')
	pam5 <- row.names(pam(t(wtd),k=7)$medoids)
	wtd[wtd < 0] <- 0
	xseq <- wt[-nrow(wt),1]
	dev.new(width=14,height=4)
	par(mar=c(0,0,0,0),xaxt="n",yaxt="n")
	plot(xseq,wtd[,1],ylim=c(-1,5),type="n")
	apply(wtd,2,lines,x=xseq,lwd=.1)
	for(i in 1:7){lines(xseq,wt[-nrow(wt),pam5[i]]-min(wt[-nrow(wt),pam5[i]])+i/2,col="white",lwd=3)}
	for(i in 1:7){lines(xseq,wt[-nrow(wt),pam5[i]]-min(wt[-nrow(wt),pam5[i]])+i/2,col="orange",lwd=2)}
	
	text(0,5,"Done")
	tmp.wr <- tmp$w.dat[-nrow(tmp$w.dat),wr.i]
	yi <- rep(-.1,2*length(plevs))
	xi <- NULL
	for(i in plevs){xi <- c(xi,min(xseq[tmp.wr==i]),max(xseq[tmp.wr==i]))}
	s1 <- seq(1,length(xi),by=2)
	s2 <- seq(2,length(xi),by=2)
#	segments(xi[s1],yi[s1],xi[s2],yi[s2],col="red",lwd=3)
	text(xi[s1],yi[s1],plevs,pos=1,cex=.5)
	xc <- 0
	yc <- 0
	while(yc < 4.5)
	{
		arrows(xi[s1],yi[s1],xi[s2],yi[s2],col="red",lwd=3,angle=90,code=3,length=.1)
		#text(xi[s1],yi[s1],plevs,pos=1,cex=.5)
		aclick <- locator(n=1)
		xc <- aclick$x
		yc <- aclick$y
		if(yc < 4.5)
		{
			click.i <- which.min(abs(xi-xc))
			arrows(xi[s1],yi[s1],xi[s2],yi[s2],col="white",lwd=3.2,angle=90,code=3,length=.1)
			#text(xi[s1],yi[s1],plevs,pos=1,cex=.65,col="white")
			
			xi[click.i] <- xc
		}
	}
	wr.out <- tmp$w.dat[,wr.i]
	wr.out[] <- ""
	for(i in 1:length(xi)){xi[i] <- which.min(abs(xi[i]-xseq))}
	for(i in 1:length(plevs)){wr.out[xi[s1[i]]:xi[s2[i]]] <- plevs[i]}
	return(wr.out)					
			
	#just click and redraw
}

ScaleRange <- function(x,min.x,max.x)
{
	xv <- as.vector(as.matrix(x))
	xv <- xv-min(x)
	xf <- (max.x-min.x)/max(xv)
	xv <- xv*xf+min.x
	xv <- matrix(xv,nrow=nrow(x))
	dimnames(xv)[[2]] <- names(x)
	dimnames(xv)[[1]] <- row.names(x)
	return(xv)
}

PlotMeanGroup <- function(tdat,grp,min.n=4,xval=NULL,col="black",pch=16)
{
	rmean <- function(x){tapply(x,as.factor(grp),mean)}
	retval <- t(apply(tdat,1,rmean))
	n.tab <- table(grp)
	usem <- names(n.tab)[n.tab > min.n]
	retval <- retval[,usem]
	col <- rainbow(length(n.tab))
	names(col) <- names(n.tab)
	if(is.null(xval)){xval <- seq(1,nrow(tdat))}

	for(i in usem)
	{
		yval <- retval[,i]
		points(xval,yval,col=col[i],pch=pch)
		lines(xval,yval,col=col[i],lwd=.5)
	}	
	text(rep(max(xval),length(usem)),retval[nrow(retval),],n.tab,pos=2)
	return(retval)
}



CalcMeanGroup <- function(tdat,grp,max.cnt=10)
{
	#means most common groups up to max.cnt
	#0 is the special collect all group
	rmean <- function(x){tapply(x,as.factor(grp),mean)}
	n.tab <- table(grp)
	if(length(n.tab) > max.cnt)
	{
		usem <- names(sort(n.tab,decreasing=T))[1:max.cnt]
		xem <- setdiff(names(n.tab),usem)
		grp[is.element(grp,xem)] <- "0"
		usem <- c(usem,"0")
	}
	retval <- t(apply(tdat,1,rmean))	
	#retval <- ScaleRange(retval,min(tdat),max(tdat))
	return(list(retval=retval,grp=grp))
}

#only keep the lowest X and highest X 
#mask out the others.
KeepX <- function(x,xlow=5,xhigh=5,mask=0)
{
	if(sum(is.na(x))>0){stop("no NAs in KeepX, please")}
	o.i <- order(x)
	x[o.i][(xlow+1):max(1,(length(x)-xhigh))] <- mask
	return(x)
}

#wt has "Time" and traces
ClusterDendroView <- function(wt,xsec=NULL,wdist=NULL,d50=NULL)
{
	require(cluster)
	linesfunc <- function(x,wt=NULL)
	{
		da.dev <- dev.cur()
		dev.new(width=14,height=8)
		
		if(length(x) > 1)
		{
			plot(xsec,wt[,names(x)[1]],ylim=range(wt[,names(x)]),type="n")
			apply(wt[,names(x),drop=F],2,lines,x=xsec,lwd=.3)
		}
		else
		{
			plot(xsec,wt[,names(x)[1]],ylim=range(wt[,names(x)]),type="l")
		}	
		dev.set(da.dev)		
	}
	if(is.null(xsec)){xsec=seq(1,nrow(wt))}
	if(is.null(w.dist)){w.dist <- dist(t(wt))}
	if(is.null(d50)){d50 <- median(w.dist)}
	da.dev <- dev.cur()
	w.hclust <- hclust(w.dist,method="complete")
	dev.new()
	plot(w.hclust)
	identify(w.hclust,FUN=linesfunc,wt=wt)
	dev.set(da.dev)	

}

#given a trace file with one challenge
#lasting ptime find the best window
#time must be the first column 
#and a decimal representation of minutes
getPrime <- function(wt,ptime=1,start.lag=1)
{
	wtd <- wt[-1,]-wt[-nrow(wt),]
	wtd <- sweep(wtd[,-1],1,wtd[,1],'/')
	wts <- wt[,-1]
	x.cnt <- ceiling(ptime/((max(wt[,1])-min(wt[,1]))/nrow(wt)))
	hab.vote <- apply(wts,2,KeepX,xlow=0,xhigh=x.cnt,mask=NA)	
	dt.vote <-  abs(apply(wtd,2,KeepX,xlow=ceiling(x.cnt/2),xhigh=ceiling(x.cnt/2),mask=NA))
	hab.sum <- apply(hab.vote,1,sum,na.rm=T)
	hab.sum[is.na(hab.sum)] <- 0
	dt.sum <- apply(dt.vote,1,sum,na.rm=T)
	dt.sum[is.na(dt.sum)] <- 0
	hab.sum <- KeepX(hab.sum,xlow=0,xhigh=x.cnt)
	dt.sum <- KeepX(dt.sum,xlow=0,xhigh=x.cnt)
	ignore <- rep(0,nrow(wt))	
	tzed <- wt[,1]-wt[1,1]
	i <- which.min(abs(tzed - start.lag))
	
	#ignoring above let's go rogue

	wt.sd <- apply(wt[,-1],1,sd)
	wt.sd.delta <- wt.sd[-1]-wt.sd[-length(wt.sd)]
	i <- which.max(wt.sd.delta)
	i2 <- which.min(abs((tzed[i]+ptime)-tzed)) 
	ignore[i:i2] <- 1	
	#select window of widtch ptime with highest vote count
	# vote.max <- 0
	# last.i <- which.min(abs((max(wt[,1])-ptime)-wt[,1]))
	# for(i in 1:last.i)
	# {
		# i2 <- which.min(abs((wt[i,1]+ptime)-wt[-nrow(wt),1]))
		# vote.sum <- sum(c(hab.sum[i:i2],dt.sum[i:i2]))
		# if(vote.sum > vote.max){vote.max <- vote.sum;vote.i <- i}
	# }
	# i <- vote.i
	# i2 <- which.min(abs((wt[i,1]+ptime)-wt[,1]))
	# ignore[i:i2] <- 1	
	# print(vote.max)
	#select highest votes and remove singles	
	#ignore <- rep(1,nrow(wt))
	# ignore[hab.sum==0 & (c(dt.sum,0)==0) ] <- 0
	# ig.sing <- ignore + c(ignore[-1],1)+(c(1,ignore[-length(ignore)]))
	# ignore[ig.sing == 1 & ignore==1] <- 0
	if(sum(ignore)==0){warning("No prime region for clustering");ignore[] <- 1}
	return(ignore)
}



BackgroundRaster <- function(wt,ht,wd,col50,xlim,ylim)
{
	require(png)
	png("tmp.png",width=wd,height=ht,units="in",res=72,bg="transparent",type="cairo")
	par(bty="n",ann=F,fig=c(0,1,0,1),mar=c(0,0,0,0),mgp=c(0,0,0),xpd=NA,xaxs="i",yaxs="i")
	plot(wt[,1],wt[,2],xaxt="n",yaxt="n",ylim=ylim,xlim=xlim,type="n")
	for(i in 2:ncol(wt)){lines(wt[,1],wt[,i],lwd=1,col=col50[i-1])}
	dev.off()
	tmp.png <- readPNG("tmp.png")
	return(tmp.png)			
}

	
ScoreLev <- function(wt,bin=NULL,ignore=NULL,col1=NULL,main.lab=NULL,bin.cat=3,rd.name="",wh=14,hh=8)
{
	#score 1 region
	#wt is trace data with Time in the first column
	#ignore is an indicator of rows to ignore

	require(RColorBrewer)
		
	graphics.off()
	orig.time<-round(wt[,1],digits=2)
	if(wt[1,1] > 0){wt[,1] <- wt[,1]-wt[1,1]} #start time at 0
	if(is.null(bin)){bin <- rep(0,ncol(wt)-1)}
	names(bin) <- names(wt)[-1]
	names(col1) <- names(wt)[-1]
	if(is.null(col1)){col1=rep("gray70",length(bin))}
	xsec <- round(wt[,1]*60,1)
	tcol50 <- col1
	colfunc <- function(x){return(rgb(red=x[1],green=x[2],blue=x[3],alpha=64,maxColorValue=255))}
	tcol50 <- apply(col2rgb(col1),2,colfunc)	
	if(is.null(ignore)|(sum(ignore)==0))
	{
		if(sum(bin)< 2)
			{ignore <- getPrime(wt,bin)}
		else
			{ignore <- getPrime(wt[,c("Time",names(bin[bin==1]))])}
	}
	#select 20 seconds of max activity
	#cluster on those for both HAB and D
	#assuming baseline corrected data.

	#despike and blc data for clustering this should only be done
	#once.

	wts <- wt#SpikeTrim(wt)#[,-1]#ScaleRange(wt[,-1],0,1)
	wtd <- wts[-1,]-wts[-nrow(wts),]
	wtd <- sweep(wtd[,-1],1,wtd[,1],'/')
	wts <- wts[,-1]##test
	
	sub.names <- names(wt[,-1])
	if(bin.cat != 3)
	{
		sub.names <- sub.names[bin[sub.names]==bin.cat]
	}
	if(length(sub.names)==0){sub.names <- names(wt[,-1])[1]}
	my.list <- list()
	my.list[[1]] <- sub.names
	list.i <- 1
	first.time=TRUE
	not.done <- TRUE
	click.i <- 1
	redraw = T
	reclust = T
	while(not.done)
	{
		if(redraw==T)
		{
			redraw <- FALSE
			sub.names <- my.list[[list.i]]
			if(bin.cat!=3)
			{
				sub.names <- sub.names[bin[sub.names]==bin.cat]
			}
			if(length(sub.names)==0){sub.names <- names(wt[,-1])[1]}			
			if(reclust == TRUE)
			{
				reclust <- FALSE
				if(length(sub.names) > 10)
				{
					wts.min <- apply(wts[ignore==1,sub.names],2,min)
					wts.blc <- sweep(wts[ignore==1,sub.names],2,wts.min,'-')
					wtc <- rbind(wtd[ignore==1,sub.names],wts.blc)
					#wtc <- wtd[ignore==1,sub.names]
					w.dist <- dist(t(wtc),method="manhattan")
					d50 <- quantile(w.dist,probs=seq(.1,.9,by=.1))			
					w.hclust <- hclust(w.dist,method="complete")
				
					for(dval in d50)
					{
						cut1 <- cutree(w.hclust,h=dval)
						if(length(table(cut1)) < 20){break}
					}
					cmg <- CalcMeanGroup(wts[,sub.names,drop=F],cut1[sub.names],max.cnt=10)
					ret.group <- cmg[["retval"]]
					cut1 <- cmg[["grp"]]				
				}
				else
				{
					cut1 <- seq(1,length(sub.names))
					names(cut1) <- sub.names
					ret.group <- wts[,sub.names,drop=F]
					names(ret.group) <- seq(1,length(sub.names))

				}
				ylim <- range(wts)
				ret.max <- max(ret.group)
				ret.spread <- max(ylim)-ret.max
				ret.sum <- apply(ret.group,2,sum)
				ret.ord <- order(ret.sum)
				ret.seq <- seq(0,ret.spread,length.out=ncol(ret.group))
				ret.group <- sweep(ret.group,2,ret.seq[ret.ord],'+')
				
#				cmg <- CalcMeanGroup(wts[,sub.names,drop=F],cut1[sub.names],max.cnt=10)
#				ret.group <- cmg[["retval"]]
#				cut1 <- cmg[["grp"]]
				
				grp.tab <- data.frame(tot=apply(ret.group,2,sum))
				grp.tab["sd"] <- apply(ret.group,2,sum)
				grp.tab["n"] <- tapply(bin[sub.names],cut1[sub.names],length)[row.names(grp.tab)]	
				grp.tab["bin"] <- tapply(bin[sub.names],cut1[sub.names],mean)[row.names(grp.tab)]

			}
					
			ng <- ncol(ret.group)
			i1s <- grep("1",ignore)
			if(length(sub.names) > 10 & length(i1s) > 10)
			{
				mi <- i1s[which.max(apply(wt[,sub.names,drop=F],1,var)[grep("1",ignore)[1:10]])]
			}
			else
			{
				mi <- i1s[min(2,length(i1s))]
			}
	
			ylim <- range(wts)
			ygap <- ylim[2]-ylim[1]
			colMean <- function(x){xcol <- apply(col2rgb(x),1,mean);return(rgb(xcol[1],xcol[2],xcol[3],alpha=255,maxColor=255))}
			grp.tab["color"] <- tapply(tcol50[sub.names],cut1[sub.names],colMean)[row.names(grp.tab)]#t.cols[rank(grp.tab[,"tot"],ties="first")]
			grp.tab["m2x"] <- rep(xsec[mi],ng)
			grp.tab["m2y"] <- as.vector(unlist(ret.group[mi,]))
			grp.tab["m1x"] <- rep(xsec[1]-xsec[2],ng)
			yseq <- seq(ylim[1],ylim[2],length.out=ng+2)[2:(ng+1)]
			grp.tab["m1y"] <- yseq[rank(grp.tab[,"m2y"],ties="first")]			
			#grp.tab["m1y"] <- seq(ylim[1]*(1+1/ng),ylim[2]*((ng-1)/ng),length.out=ng)[rank(grp.tab[,"m2y"],ties="first")]
		
			col3 <- format(grp.tab[,"n"],width=max(nchar(grp.tab[,"n"])),justify="right")
			col2 <- format(round(grp.tab[,"bin"],2),width=4,justify="right")
			col1 <- format(row.names(grp.tab),width=max(nchar(row.names(grp.tab))),justify="right")
			grp.tab["col.text"] <- paste(col1,col2,col3)
			if(first.time==TRUE)
			{
				first.time <- FALSE
				dev.new(width=wh,height=hh,family="mono",canvas="black") 
				par(fg="white",col.axis="white",col.lab="grey",col.main="grey",mar=c(3,2,2,0))
				plot(xsec,wt[,2],xlim=c((xsec[1]-xsec[2]),max(xsec)),ylim=ylim,xlab="Time (seconds)",ylab="Ratio",type="n",xaxt="n")	
				smax <- max(strwidth(paste(grp.tab[,"col.text"]," "),family="mono"))
				xlim = range(xsec)
				wt[,1] <- wt[,1]*60
				tmp.png <- BackgroundRaster(wt,8,14,tcol50,xlim,ylim)
				wt[,1] <- wt[,1]/60						
			}
			xlim <- c((xsec[1]-xsec[2])-smax,max(xsec))
			plot(xsec,wt[,2],ylim=ylim,xlim=xlim,xlab="Time (seconds)",ylab="Ratio",type="n",xaxt="n")		
			xlim = range(xsec)
			rasterImage(tmp.png,xlim[1],ylim[1],xlim[2],ylim[2])
			
			require(RColorBrewer)
			(color.number<-length(dimnames(ret.group)[[2]]))
			cols <-brewer.pal(8,"Dark2")
			cols <- rep(cols,ceiling(color.number/length(cols)))
			cols <- cols[1:color.number]
			names(cols)<-dimnames(ret.group)[[2]]

			
			for(i in dimnames(ret.group)[[2]]){lines(xsec,t(ret.group[,i]),col="white",lwd=5);lines(xsec,t(ret.group[,i]),col=cols[i],lwd=4)}
				axis(1,at=xsec, labels=FALSE)
				text(x = xsec,par("usr")[3]-.02,labels=orig.time, srt = 50, pos = 1, xpd = TRUE, col="white", cex=.7)
			
			if(length(sub.names)==1)
				{title(main=paste(main.lab,rd.name,sub.names))}
			else
				{title(main=paste(main.lab,rd.name,sep=" "))}
			for(i in row.names(grp.tab))
			{
				lines(grp.tab[i,c("m1x","m2x")],grp.tab[i,c("m1y","m2y")],col="white",lwd=1.8,type="b",cex=2)
				lines(grp.tab[i,c("m1x","m2x")],grp.tab[i,c("m1y","m2y")],col=cols[i],lwd=2,lty=2,type="b",cex=2)
				text(grp.tab[i,c("m1x","m1y")],grp.tab[i,"col.text"],col="grey",pos=2,family="mono")#grp.tab[i,"color"]
			}
			points(xsec,rep(ylim[1],length(ignore)),pch=c(2,17)[ignore+1],col=c("grey","red")[ignore+1],cex=ignore+.1)
			points(grp.tab[,"m1x"]-smax*1.3,grp.tab[,"m1y"],pch="+",cex=2)
			points(grp.tab[,"m1x"]-smax*1.1,grp.tab[,"m1y"],pch="-",cex=2)
			o.x <- rep(max(xsec),7)
			o.y <- seq(ylim[2]*.80,ylim[2],length.out=9)[2:8]
			points(o.x,o.y,pch=16,cex=2)
			text(o.x,o.y,c("Top","Done","Back","All","1s","0s","png.out"),pos=2)
	
			loc.x <- c(grp.tab[,"m1x"]-smax*1.3,grp.tab[,"m1x"]-smax*1.1,grp.tab[,"m1x"],xsec,o.x)
			loc.y <- c(grp.tab[,"m1y"],grp.tab[,"m1y"],grp.tab[,"m1y"],rep(ylim[1],length(ignore)),o.y)
		
			
		}
		
		ptx <- locator(n=1)
		click.i <- which.min(sqrt((ptx$x-loc.x)^2 + (ptx$y-loc.y)^2))
		if(click.i <= ng) #negative selected for a group
		{
			x.names <- names(cut1)[cut1==row.names(grp.tab)[click.i]]
			bin[x.names] <- 1
				grp.tab["bin"] <- tapply(bin[sub.names],cut1[sub.names],mean)[row.names(grp.tab)]
			text(grp.tab[click.i,c("m1x","m1y")],grp.tab[click.i,"col.text"],col="black",pos=2,family="mono")
			col3 <- format(grp.tab[,"n"],width=max(nchar(grp.tab[,"n"])),justify="right")
			col2 <- format(round(grp.tab[,"bin"],2),width=4,justify="right")
			col1 <- format(row.names(grp.tab),width=max(nchar(row.names(grp.tab))),justify="right")
			grp.tab["col.text"] <- paste(col1,col2,col3)
			text(grp.tab[click.i,c("m1x","m1y")],grp.tab[click.i,"col.text"],col="grey",pos=2,family="mono")	#grp.tab[click.i,"color"]						

		}
		if(click.i > ng & click.i <= ng*2) #positive selected for a group
		{
			x.names <- names(cut1)[cut1==row.names(grp.tab)[click.i-ng]]
			bin[x.names] <- 0
			grp.tab["bin"] <- tapply(bin[sub.names],cut1[sub.names],mean)[row.names(grp.tab)]
			text(grp.tab[click.i-ng,c("m1x","m1y")],grp.tab[click.i-ng,"col.text"],col="black",pos=2,family="mono")
			col3 <- format(grp.tab[,"n"],width=max(nchar(grp.tab[,"n"])),justify="right")
			col2 <- format(round(grp.tab[,"bin"],2),width=4,justify="right")
			col1 <- format(row.names(grp.tab),width=max(nchar(row.names(grp.tab))),justify="right")
			grp.tab["col.text"] <- paste(col1,col2,col3)
			text(grp.tab[click.i-ng,c("m1x","m1y")],grp.tab[click.i-ng,"col.text"],col="grey",pos=2,family="mono") #grp.tab[click.i-ng,"color"]							

		}
		if(click.i > ng*2 & click.i <= ng*3) #group selected
		{
			i.names <- names(cut1)[cut1==row.names(grp.tab)[click.i-(ng*2)]]
			list.i <- list.i+1
			my.list[[list.i]] <- i.names
			reclust <- TRUE
			redraw <- TRUE
		}
		if(click.i > ng*3 & click.i <= (ng*3+length(ignore))) #toggle ignore points
		{
			ig.i <- click.i-ng*3 
			if(ignore[ig.i] == 0)
			{points(loc.x[click.i],loc.y[click.i],pch=17,col="red",cex=1.1)}
			else
			{points(loc.x[click.i],loc.y[click.i],pch=17,col="black",cex=1.1);points(loc.x[click.i],loc.y[click.i],pch=2,col="grey",cex=.1)}
			ignore[ig.i] <- abs(ignore[ig.i]-1)			
		}
		else
		{
			o.i <- click.i-(ng*3+length(ignore))
			if(o.i==1) #reset
			{
				#print("reset")
				list.i <- 1
				reclust <- TRUE
				redraw <- TRUE
			}
			if(o.i==2) #done
			{
				#print("done")
				not.done=FALSE
				sub.names <- names(wt[,-1])
				
				wtc <- rbind(wtd[ignore==1,sub.names],wts[ignore==1,sub.names])
				w.dist <- dist(t(wtc),method="manhattan")
				d50 <- quantile(w.dist,probs=seq(.1,.9,by=.1))			
				w.hclust <- hclust(w.dist,method="complete")
				for(dval in d50)
				{
					cut1 <- cutree(w.hclust,h=dval)
					if(length(table(cut1)) < 10){break}
				}
				cmg <- CalcMeanGroup(wts[,sub.names,drop=F],cut1[sub.names],max.cnt=6)
				ret.group <- cmg[["retval"]]
				cut1 <- cmg[["grp"]]	
				ret.max <- max(ret.group)
				ret.spread <- max(ylim)-ret.max
				ret.sum <- apply(ret.group,2,sum)
				ret.ord <- order(ret.sum)
				ret.seq <- seq(0,ret.spread,length.out=ncol(ret.group))
				ret.group <- sweep(ret.group,2,ret.seq[ret.ord],'+')
							
			}
			if(o.i==3) #back
			{
				#print("back")
				list.i <- max(list.i-1,1)
				reclust <- TRUE
				redraw <- TRUE

			}	
			if(o.i==4) #all
			{
				#print("all")
				bin.cat <- 3
				reclust <- TRUE
				redraw <- TRUE

			}
			if(o.i==5) #1s
			{
				#print("1s")
				bin.cat <- 1
				reclust <- TRUE
				redraw <- TRUE
				
			}
			if(o.i==6) #0s
			{			
				#print("0s")
				bin.cat <- 0	
				reclust <- TRUE
				redraw <- TRUE
				
			}
			if(o.i==7) #png.out things are getting messy
			{

				png(paste(rd.name,main.lab,bin.cat,"out.png",sep="."),width=1200,height=800,family="mono") 
				par(fg="black",col.axis="white",col.lab="grey",col.main="grey",mar=c(3,3,3,0),cex.main=2)
				plot(xsec,wt[,2],xlim=c((xsec[1]-xsec[2]),max(xsec)),ylim=ylim,xlab="Time (seconds)",ylab="Ratio",type="n",xaxt="n")	
			
				xlim <- c((xsec[1]-xsec[2])-smax,max(xsec))
				plot(xsec,wt[,2],ylim=ylim,xlim=xlim,xlab="Time (seconds)",ylab="Ratio",type="n",xaxt="n")		
				xlim = range(xsec)
				
				for(i in dimnames(ret.group)[[2]]){lines(xsec,t(ret.group[,i]),col="black",lwd=6);lines(xsec,t(ret.group[,i]),col=grp.tab[i,"color"],lwd=5)}
				axis(1,at=xsec, labels=FALSE)
				text(x = xsec,par("usr")[3]-.02,labels=orig.time, srt = 45, pos = 1, xpd = TRUE, col="white", cex=.75)

				if(length(sub.names)==1)
					{title(main=paste(main.lab,,rd.name,sub.names))}
				else
					{title(main=paste(main.lab,rd.name,sep=" "))}
				for(i in row.names(grp.tab))
				{
					lines(grp.tab[i,c("m1x","m2x")],grp.tab[i,c("m1y","m2y")],col="white",lwd=1.8,type="b",cex=2)
					lines(grp.tab[i,c("m1x","m2x")],grp.tab[i,c("m1y","m2y")],col=grp.tab[i,"color"],lwd=2,lty=2,type="b",cex=2)
					text(grp.tab[i,c("m1x","m1y")],as.character(grp.tab[i,"n"]),col="black",pos=2,family="mono",cex=3)#grp.tab[i,"color"]
				}
				points(xsec,rep(ylim[1],length(ignore)),pch=c(2,17)[ignore+1],col=c("grey","red")[ignore+1],cex=ignore+.1)
				dev.off()	
				
			}
		}
	}
	graphics.off()

	return(list(ignore=ignore,bin=bin,groups=cut1))	
	#return(list(means=ret.group,groups=cut1,ignore=ignore,bin=bin,resp.dets=resp.dets))
}

SetShades <- function(tmp,trans=.1)
{
	gfp.frac <- NULL
	ib4.frac <- NULL
	gfp.name <- grep("gfp",names(tmp$c.dat),ignore.case=T,value=T)
	if(length(gfp.name) > 0){gfp.name <- gfp.name[[1]]}
	ib4.name <- grep("ib4",names(tmp$c.dat),ignore.case=T,value=T)
	if(length(ib4.name) > 0){ib4.name <- ib4.name[[1]]}	
	roi.name <- grep("area",names(tmp$c.dat),ignore.case=T,value=T)
	if(length(roi.name) > 0){roi.name <- roi.name[[1]]}
	cy <- grep("center.*y$",names(tmp$c.dat),ignore.case=T,value=T)
	if(length(cy) > 0){cy <- cy[[1]]}
	cx <- grep("center.*x$",names(tmp$c.dat),ignore.case=T,value=T)
	if(length(cy) > 0){cy <- cy[[1]]}
	if(length(c(cx,cy,gfp.name))==3)
	{
		y <- log(tmp$c.dat[,gfp.name]+1)
		if(sum(is.na(y)) > 0){y <- runif(nrow(tmp$c.dat))}
		x.pos <- tmp$c.dat[,cx]
		y.pos <- tmp$c.dat[,cy]
		cent.x <- max(x.pos)-c(480,640,1024)
		cent.x[cent.x < 0] <- 480
		cent.x <- c(240,320,512)[which.min(cent.x)]
		cent.y <- max(y.pos)-c(480,640,1024)
		cent.y[cent.y < 0] <- 480
		cent.y <- c(240,320,512)[which.min(cent.y)]
		
		glt <- glm(y ~ sqrt((x.pos-cent.x)^2+(y.pos-cent.y)^2) + x.pos + y.pos)
		gfp.frac <- residuals(glt)
		gfp.frac <- (rank(gfp.frac)-1)/(length(gfp.frac)-1)
	}
	if(length(c(cx,cy,ib4.name))==3)
	{
		y <- log(tmp$c.dat[,ib4.name]+1)
		if(sum(is.na(y)) > 0){y <- runif(nrow(tmp$c.dat))}
		x.pos <- tmp$c.dat[,cx]
		y.pos <- tmp$c.dat[,cy]
		cent.x <- max(x.pos)-c(480,640,1024)
		cent.x[cent.x < 0] <- 480
		cent.x <- c(240,320,512)[which.min(cent.x)]
		cent.y <- max(y.pos)-c(480,640,1024)
		cent.y[cent.y < 0] <- 480
		cent.y <- c(240,320,512)[which.min(cent.y)]
		
		glt <- glm(y ~ sqrt((x.pos-cent.x)^2+(y.pos-cent.y)^2) + x.pos + y.pos)
		ib4.frac <- residuals(glt)
		ib4.frac <- (rank(ib4.frac)-1)/(length(ib4.frac)-1)
	}
			
# # 	glt <- glm(log(mean.tritc+1) ~ sqrt((center.x-512)^2+(center.y-512)^2) + center.x + center.y,data=tmp$c.dat)
	# ib4.frac <- residuals(glt)
	# ib4.frac <- (rank(ib4.frac)-1)/(length(ib4.frac)-1)

	# glt <- glm(log(mean.gfp.1+1) ~ sqrt((center.x-512)^2+(center.y-512)^2) + center.x + center.y,data=tmp$c.dat)
	# gfp.frac <- residuals(glt)
	# gfp.frac <- (rank(gfp.frac)-1)/(length(gfp.frac)-1)

	size.frac <- tmp$c.dat[,roi.name]
	size.frac <- (rank(size.frac)-1)/(length(size.frac)-1)
	if(is.null(size.frac)){size.frac <- seq(1,nrow(tmp$c.dat))/nrow(tmp$c.dat)}	
	names(size.frac) <- row.names(tmp$bin)
#	names(gfp.frac) <- row.names(tmp$bin)
#	names(ib4.frac) <- row.names(tmp$bin)

	#size.frac <- size.frac/2+.5

	if(is.null(gfp.frac)){gfp.frac <- size.frac*.5}
	if(is.null(ib4.frac)){ib4.frac <- size.frac*.5}
	col50 <- rgb(ib4.frac,gfp.frac,size.frac,alpha=trans) 
	names(col50) <- row.names(tmp$bin)
	return(col50)
	
}

#tmp is an RD object
#wr.i is the window region definition.
#rd.name is the name of the RD object (used for png.out)
#c.i is an alternate list of cells to review (instead of all)
#rscale is a boolean for rescaling the data
#wh is the window height
#hh is the window width (why the hell did I name it hh?)
RDView <- function(tmp,wr.i="wr1",rd.name=NULL,c.i=NULL,rscale=F,wh=14,hh=8)
{
	if(!is.element("bin",names(tmp))){stop("No bin ")}
	if(!is.element("drop",names(tmp$bin))){tmp$bin[,"drop"] <- 0}
	col50 <- SetShades(tmp,trans=200/(nrow(tmp$bin)^1.5))
	tlevs <- unique(tmp$w.dat[,wr.i])
 	tlevs <- tlevs[tlevs!=""]
 	tlevs <- intersect(tlevs,names(tmp$bin))
	if(!is.element("ignore",names(tmp$w.dat)))
	{
		tmp$w.dat[,"ignore"] <- 0
	}
	sel.i <- 1
#	if(!is.element("grps",names(tmp))){tmp$grps <- tmp$bin[,tlevs];tmp$grps[] <- NA}

	while(sel.i != 0)
	{
		sel.i <- menu(tlevs,,title="Select Window To Review")
		if(sel.i > 0)
		{
			r.i <- tmp$w.dat[,wr.i]==tlevs[sel.i]
			if(is.null(c.i)){c.i <- row.names(tmp$bin)[tmp$bin[,"drop"]==0]}
			if(is.element("mp",names(tmp)))
				{wt <- tmp$mp[r.i,c("Time",c.i)]}
			else
				{wt <- tmp$t.dat[r.i,c("Time",c.i)]}
			if(rscale){wt <- sweep(wt,2,apply(wt,2,min),"-")}
			wt.sl <- ScoreLev(wt,tmp$bin[c.i,tlevs[sel.i]],tmp$w.dat[r.i,"ignore"],col50,main.lab=tlevs[sel.i],rd.name=rd.name,wh=wh,hh=hh)
			tmp$bin[c.i,tlevs[sel.i]] <- wt.sl$bin
			tmp$w.dat[r.i,"ignore"] <- wt.sl$ignore
			tmp$grps <- wt.sl[["groups"]]
		}
	}		
	return(tmp)
}

#off loaded from above this should calculate the relevant features of a response
#presumably wts is spike trimed data wtd is the delta data. ignore is a vector of points to use
#xsec is the time in seconds
DescribeResponse <- function(wts,wtd,ignore,xsec)
{
	#calculate max delta up max delta down max total rise max total fall
	resp.dets <- data.frame(dmean=apply(wtd,2,mean))
	resp.dets["max.up"] <- apply(wtd[ignore[-length(ignore)]==1,],2,max)
	resp.dets["max.dn"] <- apply(wtd[ignore[-length(ignore)]==1,],2,min)
	resp.dets["n.up"] <- apply(wtd[ignore[-length(ignore)]==1,]>0,2,sum)
	resp.dets["n.dn"] <- apply(wtd[ignore[-length(ignore)]==1,]<0,2,sum)
	
#
	upfunc <- function(x)
	{
		xdiffunc <- function(i){x[i[2]]-x[i[1]]}	
		dn.x <- abs(min(apply(combn(1:length(x),2),2,xdiffunc)))
		up.x <- abs(max(apply(combn(1:length(x),2),2,xdiffunc)))
		return(c(dn.x,up.x))
	}


	#change values in the wts to indicate the baseline as the first point
	#this removes the "sipping trough"
	#ignore should not have anything left of the sipping trough
	xsec <- xsec-min(xsec)
	if(sum(xsec < 20 & ignore==0) > 5)
	{bl.med <- apply(wts[xsec < 20 & ignore==0,row.names(resp.dets)],2,median)}
	else
	{bl.med <- apply(wts[,row.names(resp.dets)],2,min)}
	wi <- max(match(1,ignore)-1,1)
	wts[wi,row.names(resp.dets)] <- as.vector(bl.med)
	ignore[wi] <- 1

	t20.1 <- min(xsec[ignore==1]) + 20
	dn.up20 <- apply(wts[ignore==1 & xsec < t20.1,],2,upfunc)
	resp.dets["max.rise20"] <- dn.up20[2,]
	dn.up <- apply(wts[ignore==1,],2,upfunc)
	resp.dets["max.fall"] <- dn.up[1,]
	resp.dets["max.rise"] <- dn.up[2,]
	resp.dets["poly0"] <- NA
	resp.dets["poly1"] <- NA	
	resp.dets["poly2"] <- NA	
	resp.dets["poly3"] <- NA	

	for(i in row.names(resp.dets))
	{
		lt <- lm(wts[ignore==1,i] ~ poly(xsec[ignore==1],3))
		resp.dets[i,c("poly0","poly1","poly2","poly3")] <- coefficients(lt)
		#plot(wts[ignore==1,i] ~ xsec[ignore==1])
		#lines(xsec[ignore==1],predict(lt))
	}
	#above baseline area
	wti <- wts[ignore==1,]
	wti <- sweep(wti,2,apply(wti,2,min))
	resp.dets["tsum"] <- apply(wti,2,sum) * (max(xsec[ignore==1])-min(xsec[ignore==1])) #should use m%*%n FIGURE IT OUT!
	return(resp.dets)	
	
}

#
DescribeAll <- function(tmp,t.levs=NULL,wr="wr1",mpyes=FALSE)
{
	if(!is.element("mp",names(tmp))){mpyes <- FALSE}
	if(is.null(t.levs)){t.levs <- setdiff(unique(tmp$w.dat[,wr]),c(""))}
	for(t in t.levs)
	{
		if(mpyes)
		{wt <- tmp$mp[tmp$w.dat[,wr]==t,]}
		else
		{wt <- tmp$t.dat[tmp$w.dat[,wr]==t,]}
		ignore <- tmp$w.dat[tmp$w.dat[,wr]==t,"ignore"]
		
		wts <- wt
		wtd <- wts[-1,]-wts[-nrow(wts),]
		wtd <- sweep(wtd[,-1],1,wtd[,1],'/')
		wts <- wts[,-1]
		xsec <- round(wt[,1]*60,0)		
		if(sum(ignore) > 5)
		{
			resp.dets <- DescribeResponse(wts,wtd,ignore,xsec)
			names(resp.dets) <- paste(t,names(resp.dets),sep=".")
			for(rp.i in names(resp.dets)){tmp$scp[,rp.i] <- resp.dets[row.names(tmp$scp),rp.i]}		
		}
	}
	return(tmp)
}

#
DescribeAll_lee <- function(tmp,t.levs=NULL,wr="wr1",blc_yes=FALSE)
{
	if(!is.element("blc",names(tmp))){blc_yes <- FALSE}
	if(is.null(t.levs)){t.levs <- setdiff(unique(tmp$w.dat[,wr]),c(""))}
	for(t in t.levs)
	{
		if(blc_yes)
		{wt <- tmp$blc[tmp$w.dat[,wr]==t,]}
		else
		{wt <- tmp$t.dat[tmp$w.dat[,wr]==t,]}
		
		ignore <- tmp$w.dat[tmp$w.dat[,wr]==t,"ignore"]
		wts <- wt
		wtd <- wts[-1,]-wts[-nrow(wts),]
		wtd <- sweep(wtd[,-1],1,wtd[,1],'/')
		wts <- wts[,-1]
		xsec <- round(wt[,1]*60,0)		
		if(sum(ignore) > 5)
		{
			resp.dets <- DescribeResponse(wts,wtd,ignore,xsec)
			names(resp.dets) <- paste(t,names(resp.dets),sep=".")
			for(rp.i in names(resp.dets)){tmp$scp[,rp.i] <- resp.dets[row.names(tmp$scp),rp.i]}		
		}
	}
	return(tmp)
}


#Display the analysis of a single trace 
#dat is the trace dataframe with "Time" in the first column and cell trace intensities in subsequent columns
#i is the index column to be analyzed and displayed.
#shws is the smoothing half window size
#Plotit is a flag indicating that the results should be ploted or not.
#wr is the response window factor 
#SNR.lim is the signal to noise ratio limit for peak detection
#bl.meth is the method for baseline correction.
PlotFunc <- function(dat,i,shws=2,phws=20,Plotit=F,wr=NULL,SNR.lim=2,bl.meth="TopHat",lmain=NULL)
{
    library("MALDIquant")
    xt <- dat[,"Time"]
    xd.med <- median(xt[-1]-xt[-length(xt)])
    x.pad <- (seq(1:5)*xd.med)+max(xt)
    xt <- c(xt,x.pad)
    xi <- c(seq(1,nrow(dat)),1:5)
    y <- dat[xi,i]
    s1 <- createMassSpectrum(xt,y)
    if(shws > 1)
        s3 <- smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=shws)
    else
        s3 <- s1
    if(Plotit)
    {
        bSnip <- estimateBaseline(s3, method="SNIP")
        bTopHat <- estimateBaseline(s3, method="TopHat")
        bMed <- estimateBaseline(s3, method="median",halfWindowSize=50)
        bConvec <- estimateBaseline(s3, method="ConvexHull")
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














