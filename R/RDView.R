#' RDView 
#' This is a function that allows you to correct the binary
#' scoring for the application of interest. This uses advanced
#' heirarchal clustering to guide you through score correction.
#' 
#' @param tmp is an RD object
#' @param wr.i is the window region definition.
#' @param rd.name is the name of the RD object (used for png.out)
#' @param c.i is an alternate list of cells to review (instead of all)
#' @param rscale is a boolean for rescaling the data
#' @param wh is the window height
#' @param hh is the window width (why the hell did I name it hh?)
#' @export 
RDView <- function(tmp,c.i=NULL,wr.i="wr1",rd.name=NULL,rscale=F,wh=14,hh=8){
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
    #if(!is.element("grps",names(tmp))){tmp$grps <- tmp$bin[,tlevs];tmp$grps[] <- NA}

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

SetShades <- function(tmp,trans=.1){
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
			
    #glt <- glm(log(mean.tritc+1) ~ sqrt((center.x-512)^2+(center.y-512)^2) + center.x + center.y,data=tmp$c.dat)
	# ib4.frac <- residuals(glt)
	# ib4.frac <- (rank(ib4.frac)-1)/(length(ib4.frac)-1)

	# glt <- glm(log(mean.gfp.1+1) ~ sqrt((center.x-512)^2+(center.y-512)^2) + center.x + center.y,data=tmp$c.dat)
	# gfp.frac <- residuals(glt)
	# gfp.frac <- (rank(gfp.frac)-1)/(length(gfp.frac)-1)

	size.frac <- tmp$c.dat[,roi.name]
	size.frac <- (rank(size.frac)-1)/(length(size.frac)-1)
	if(is.null(size.frac)){size.frac <- seq(1,nrow(tmp$c.dat))/nrow(tmp$c.dat)}	
	names(size.frac) <- row.names(tmp$bin)
    #names(gfp.frac) <- row.names(tmp$bin)
    #names(ib4.frac) <- row.names(tmp$bin)

	#size.frac <- size.frac/2+.5

	if(is.null(gfp.frac)){gfp.frac <- size.frac*.5}
	if(is.null(ib4.frac)){ib4.frac <- size.frac*.5}
	col50 <- rgb(ib4.frac,gfp.frac,size.frac,alpha=trans) 
	names(col50) <- row.names(tmp$bin)
	return(col50)
	
}
	
ScoreLev <- function(wt,bin=NULL,ignore=NULL,col1=NULL,main.lab=NULL,bin.cat=3,rd.name="",wh=14,hh=8){
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
				
			#cmg <- CalcMeanGroup(wts[,sub.names,drop=F],cut1[sub.names],max.cnt=10)
			#ret.group <- cmg[["retval"]]
			#cut1 <- cmg[["grp"]]
				
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

#given a trace file with one challenge
#lasting ptime find the best window
#time must be the first column 
#and a decimal representation of minutes
getPrime <- function(wt,ptime=1,start.lag=1){
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

#only keep the lowest X and highest X 
#mask out the others.
KeepX <- function(x,xlow=5,xhigh=5,mask=0){
	if(sum(is.na(x))>0){stop("no NAs in KeepX, please")}
	o.i <- order(x)
	x[o.i][(xlow+1):max(1,(length(x)-xhigh))] <- mask
	return(x)
}

CalcMeanGroup <- function(tdat,grp,max.cnt=10){
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

BackgroundRaster <- function(wt,ht,wd,col50,xlim,ylim){
	require(png)
	png("tmp.png",width=wd,height=ht,units="in",res=72,bg="transparent",type="cairo")
	par(bty="n",ann=F,fig=c(0,1,0,1),mar=c(0,0,0,0),mgp=c(0,0,0),xpd=NA,xaxs="i",yaxs="i")
	plot(wt[,1],wt[,2],xaxt="n",yaxt="n",ylim=ylim,xlim=xlim,type="n")
	for(i in 2:ncol(wt)){lines(wt[,1],wt[,i],lwd=1,col=col50[i-1])}
	dev.off()
	tmp.png <- readPNG("tmp.png")
	return(tmp.png)			
}



