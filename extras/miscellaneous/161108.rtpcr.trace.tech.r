#####
#rtpcr analysis


#Plot Traces with rtpcr data next to trace.
#melt.plot, adds melt curve analysis
rtpcr.plotter<-function(dat,pdf=F,bcex=1,melt.plot=T){

	primer.selection<-select.list(names(dat$rtpcr$pcr.quant.sum)) # select which primer to assess
	pcr.quant.select<-dat$rtpcr$pcr.quant.sum[[primer.selection]] # make that the data frame to reference
	select.cells<-select.list(names(pcr.quant.select), multiple=T) #
	if(length(grep("blank",select.cells))>=1){x.names<-select.cells[-grep("blank",select.cells)]
	}else{x.names<-select.cells}
	for(i in x.names){
	if(melt.plot){width=6;height=5}else{width=10; height=3}
	
	if(pdf){pdf(paste(i,".pdf"),width=width, height=height)
		}else{dev.new(width=width, height=height)}
		
	if(melt.plot){
		layout(rbind(matrix(1, 1, 3),matrix(c(2,3,4),1,3)), widths=c(2,2,2), heights=c(2.5,2.5))
		}else{layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(14,6), heights=6)}	
	PeakFunc7(dat, i,bcex=bcex)
	par(mar=c(4,4,4,4) ,cex=.8)
	plot(pcr.quant.select$Cycle, pcr.quant.select[,2], ylim=c(-10, max(pcr.quant.select)), type="l", col=NULL, xlab="Cycle", ylab="Signal")
	matlines(pcr.quant.select$Cycle, pcr.quant.select[,select.cells], col="black", lwd=2, lty=1)
	par(xpd=T, bty="l")
	for(j in select.cells){
		text(max(pcr.quant.select$Cycle)*1.05,max(pcr.quant.select[,j]), labels=j, cex=.9)
	}
	
	lines(pcr.quant.select$Cycle, pcr.quant.select[,i], col="red", lwd=2, lty=1)
	text(max(pcr.quant.select$Cycle)*1.05,max(pcr.quant.select[,i]), labels=i, cex=.9, col="red")
	par(xpd=F)
	
	
	if(melt.plot){
		melt.curve<-dat$rtpcr$melt.curve[[primer.selection]]
		plot(melt.curve$Temperature, melt.curve[,i], ylim=c(-10, max(melt.curve)),lwd=2, type="l", col="black", xlab="Cycle", ylab="Temperature")
	
		melt.deriv<-dat$rtpcr$melt.deriv[[primer.selection]]
		plot(melt.deriv$Temperature, melt.deriv[,i], ylim=c(-10, max(melt.deriv)),lwd=2, type="l", col="black", xlab="Cycle", ylab="DF/DT")
	}

	if(pdf){dev.off()}
	
	}
}


# This is how i plot rtpcr data when there are more than one primer
#updated for melt plotter
rtpcr.multi.plotter<-function(dat,x.names=NULL, primer=NULL,pdf=F,bcex=1, melt.plot=T, plot.new=T,lmain=NULL){

	if(is.null(primer)){
	primer.selection<-select.list(names(dat$rtpcr$pcr.quant.sum),multiple=T)
	}else{primer.selection<-primer} # select which primer to assess
	pcr.quant.select<-dat$rtpcr$pcr.quant.sum[[primer.selection[1]]] # make that the data frame to reference
	if(is.null(x.names)){select.cells<-select.list(names(pcr.quant.select), multiple=T)}else{select.cells<-x.names}
	#if(length(grep("blank",select.cells))>=1){x.names<-select.cells[-grep("blank",select.cells)]
	#}else{x.names<-select.cells}
	x.names<-select.cells

		
	for(i in x.names){
	if(melt.plot){width=6;height=5}else{width=10; height=3}
	
	if(pdf){pdf(paste(i,".pdf"),width=width, height=height)
	}else{if(plot.new){dev.new(width=width, height=height)}}
	
	if(melt.plot){
	layout(rbind(matrix(1, 1, 3),matrix(c(2,3,4),1,3)), widths=c(2,2,2), heights=c(2.5,2.5))}else{layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(14,6), heights=6)}	
	
	if(length(grep("blank|extracellular",i))>=1){
		plot(0,0, xlim=c(0,10),ylim=c(0,10))
		text(5,5, labels=i)
		}else{
		#if(is.null(lmain)){lmain<-i}else{lmain<-lmain}
		PeakFunc7(dat, i,bcex=bcex)}
	
	par(mar=c(4,4,4,4) ,cex=.8)
	plot(pcr.quant.select$Cycle, pcr.quant.select[,i], ylim=c(-10, max(pcr.quant.select)), type="l", col=NULL, xlab="Cycle", ylab="Signal")
	cols<-c("black","red")
	for(j in 1:length(primer.selection)){
		lines(dat$rtpcr$pcr.quant.sum[[j]][,"Cycle"], dat$rtpcr$pcr.quant.sum[[j]][,i], col=cols[j], lwd=2, lty=1)
		#text(max(dat$rtpcr$pcr.quant.sum[[j]][,"Cycle"])*1.05,max(dat$rtpcr$pcr.quant.sum[[j]][i]), labels=names(dat$rtpcr$pcr.quant.sum)[j], cex=.9, cols=cols[j])
		}
	legend("topright", primer.selection, lwd=2, col=cols, bty="n", cex=.6)
	par(xpd=T, bty="l")

	#lines(pcr.quant.select$Cycle, pcr.quant.select[,i], col="red", lwd=2, lty=1)
	#text(max(pcr.quant.select$Cycle)*1.05,max(pcr.quant.select[,i]), labels=i, cex=.9, col="red")
	par(xpd=F)
	
	if(melt.plot){
	
		melt.curve<-dat$rtpcr$melt.curve[[primer.selection[1]]]
		plot(melt.curve$Temperature, melt.curve[,i], ylim=c(-10, max(melt.curve)), type="l", col=NULL, xlab="Temp", ylab="Signal", lwd=2)
		par(xpd=T)
		for(j in 1:length(primer.selection))
		{
			lines(dat$rtpcr$melt.curve[[j]][,"Temperature"], dat$rtpcr$melt.curve[[j]][,i], col=cols[j], lwd=2, lty=1)
			}
		legend("topright", primer.selection, lwd=2, col=cols, bty="n", cex=.6)
			
		melt.deriv<-dat$rtpcr$melt.deriv[[primer.selection[1]]]
		plot(melt.deriv$Temperature, melt.deriv[,i], ylim=c(-10, max(melt.deriv)), type="l", col=NULL, xlab="Temp", ylab="DF/DT",lwd=2)
		par(xpd=T)
		for(j in 1:length(primer.selection)){
			lines(dat$rtpcr$melt.deriv[[j]][,"Temperature"], dat$rtpcr$melt.deriv[[j]][,i], col=cols[j], lwd=2, lty=1)
			}
		}
			legend("topright", primer.selection, lwd=2, col=cols, bty="n", cex=.6)


	if(pdf){dev.off()}
	}
}


# How to create a cutoff for a barplot comparison
pcr.x.cutter<-function(dat){

	primer.selection<-select.list(names(dat$rtpcr$pcr.quant.sum)) # select which primer to assess
	pcr.quant.select<-dat$rtpcr$pcr.quant.sum[[primer.selection]] # make that the data frame to reference
	select.cells<-select.list(names(pcr.quant.select), multiple=T) #select the cells you would like to observe
	x.names<-select.cells[-grep("blank",select.cells)]#create a clean set to plot traces
	
		#cols <- rainbow(length(m.names),start=.55)
		require(RColorBrewer)
		cols <-brewer.pal(8,"Dark2")
        cols <- rep(cols,ceiling(length(select.cells)/length(cols)))
        cols <- cols[1:length(select.cells)]
	LinesEvery.5(dat, select.cells, cols=cols,plot.new=T, img=dat$img1)
	
	dev.new(width=8, height=4)
	par(mfrow=c(1,2))
	par(mar=c(4,4,4,4) ,cex=.8,bty="l")
	plot(pcr.quant.select$Cycle, pcr.quant.select[,2], ylim=c(-10, max(pcr.quant.select)), type="l", col=NULL, xlab="Cycle", ylab="Signal")
	matlines(pcr.quant.select$Cycle, pcr.quant.select[,select.cells], col=cols, lwd=2, lty=1)
	par(xpd=T, bty="l")
	for(j in 1:length(select.cells)){
		text(max(pcr.quant.select$Cycle)*1.05,max(pcr.quant.select[,select.cells[j]]), labels=select.cells[j], cex=.9, col=cols[j])
	}

	cutoff<-locator(n=1)
	cutoff<-ceiling(cutoff$x)
	par(xpd=F)
	abline(v=cutoff, col="red")
	
	bob<-t(pcr.quant.select[row.names(pcr.quant.select)[pcr.quant.select["Cycle"]==cutoff],select.cells])
	
	barplot(bob, beside=T, col=cols, names.arg=select.cells, las=2)
}
	
pcr.y.cutter<-function(dat, cutoff=NULL){

	primer.selection<-select.list(names(dat$rtpcr$pcr.quant.sum)) # select which primer to assess
	pcr.quant.select<-dat$rtpcr$pcr.quant.sum[[primer.selection]] # make that the data frame to reference
	select.cells<-select.list(names(pcr.quant.select), multiple=T) #select the cells you would like to observe
	x.names<-select.cells[-grep("blank",select.cells)]#create a clean set to plot traces
	
		#cols <- rainbow(length(m.names),start=.55)
		require(RColorBrewer)
		cols <-brewer.pal(8,"Dark2")
        cols <- rep(cols,ceiling(length(select.cells)/length(cols)))
        cols <- cols[1:length(select.cells)]
	LinesEvery.5(dat, select.cells, cols=cols,plot.new=T, img=dat$img1)
	
	dev.new(width=8, height=4)
	par(mfrow=c(1,2))
	par(mar=c(4,4,4,4) ,cex=.8,bty="l")
	plot(pcr.quant.select$Cycle, pcr.quant.select[,2], ylim=c(-10, max(pcr.quant.select)), type="l", col=NULL, xlab="Cycle", ylab="Signal")
	matlines(pcr.quant.select$Cycle, pcr.quant.select[,select.cells], col=cols, lwd=2, lty=1)
	par(xpd=T, bty="l")
	for(j in 1:length(select.cells)){
		text(max(pcr.quant.select$Cycle)*1.05,max(pcr.quant.select[,select.cells[j]]), labels=select.cells[j], cex=.9, col=cols[j])
	}

	
	if(is.null(cutoff)){
		cutoff<-locator(n=1)
		cutoff<-cutoff$y
		par(xpd=F)
		abline(h=cutoff, col="red")
	}
	
		#for(i in select.cells){
	#	thresh[i]<-min(which(pcr.quant.select[,i]>=cutoff, arr.ind=T))
	#	if(thresh[i]=="Inf"){thresh[i]<-NA}
	#	}

	thresh<-vector()
	for(i in select.cells){
	cell<-pcr.quant.select[,i]
	Cycle<-pcr.quant.select[,1]

	cell<-cell[25:length(cell)]
	Cycle<-Cycle[25:length(Cycle)]
	#plot(cell~Cycle, main=i,ylim=c(0,max(Reduce(cbind,dat$rtpcr$pcr.quant.sum))))

	model<-tryCatch(nls(cell~SSfpl(Cycle,a,b,c,d)),error=function(e) NULL)
	if(!is.null(model)){

		xv<-seq(min(Cycle),max(Cycle), .05)
		yv<-predict(model, list(Cycle=xv))
		xyv<-cbind(xv,yv)
		thresh[i]<-tryCatch(xyv[min(which(yv>=cutoff, arr.ind=T)),1],error=function(e) NULL)

		#lines(xv,yv)
		summary(model)
	}else{thresh[i]<-NA}
}

	#for(i in select.cells){
	#	thresh[i]<-min(which(pcr.quant.select[,i]>=cutoff, arr.ind=T))
	#	if(thresh[i]=="Inf"){thresh[i]<-NA}
	#	}
	
	barplot(thresh, beside=T, col=cols, names.arg=select.cells, las=2, log="y",ylab="cq")
	return(thresh)
	return(cutoff)
}
	
PeakFunc.rtpcr <- function(dat,n.names,t.type=NULL,select.trace=F,Plotit.trace=T,Plotit.both=F, info=T,lmain=NULL, bcex=1, ylim.max=NULL)
{
	par(xpd=FALSE)
	if(select.trace==TRUE){
		dat.select<-menu(names(dat))
		dat.t<-dat[[dat.select]]
	}else{
	
		if(is.null(t.type)){dat.t<-dat$mp}else{dat.t<-dat[[t.type]]}
		}
	if(is.null(ylim.max)){ylim.max<-1.4}else{ylim.max<-ylim.max}
	xlim <- range(dat.t[,1]) # use same xlim on all plots for better comparison
	if(Plotit.trace){ylim <- c(.2,ylim.max)}
	if(Plotit.both){ylim <- c(-.5,ylim.max)}
	
	levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
	par(mar=c(9,4.5,3.5,9),xpd=T, bty="l")
	plot(dat.t[,n.names]~dat.t[,1], main=paste(lmain,n.names),xlim=xlim,ylim=ylim,xlab="", ylab="",type="n", cex=bcex)
	#axis(1, at=seq(0, length(dat.t[,1]), 5),tick=TRUE )  
	par(xpd=F)
	
	# Tool for labeling the binary score
	wr<-dat$w.dat[,"wr1"]
	levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
	x1s <- tapply(dat.t[,"Time"],as.factor(wr),min)[levs]
	#abline(v=x1s,col="black")
	levs <- setdiff(unique(wr),"")
	text(dat.t[match(levs,wr),"Time"],c(min(ylim), min(ylim)+.1),levs,pos=4,offset=0,cex=bcex)
	
	# Tool for labeling cellular aspects, gfp.1, gfp.2, tritc, area
	legend(x=max(xlim)*.75, y=ylim.max+yinch(.3), xpd=TRUE, inset=c(0,-.14), 
	legend=c(
		#if(!is.null(dat$c.dat[n.names, "CGRP"])){paste("CGRP","",round(dat$c.dat[n.names,"CGRP"],digits=0))},
		#if(!is.null(dat$c.dat[n.names, "mean.gfp"])){paste("GFP","",round(dat$c.dat[n.names,"mean.gfp"],digits=0))},
		#if(!is.null(dat$c.dat[n.names, "mean.gfp.1"])){paste("GFP.1","",round(dat$c.dat[n.names,"mean.gfp.1"],digits=0))},
		#if(!is.null(dat$c.dat[n.names, "mean.gfp.2"])){paste("GFP.2","",round(dat$c.dat[n.names,"mean.gfp.2"],digits=0))},
		#if(!is.null(dat$c.dat[n.names, "mean.dapi"])){paste("DAPI","",round(dat$c.dat[n.names,"mean.dapi"],digits=0))},
		if(!is.null(dat$c.dat[n.names, "IB4"])){paste("IB4","",round(dat$c.dat[n.names,"IB4"],digits=0))},
		if(!is.null(dat$c.dat[n.names, "mean.tritc"])){paste("IB4","",round(dat$c.dat[n.names, "mean.tritc"], digits=0))}, 
		if(!is.null(dat$c.dat[n.names, "area"])){paste("Area","", round(dat$c.dat[n.names, "area"], digits=0))},
		if(!is.null(dat$c.dat[n.names, "ROI.Area"])){paste("Area","", round(dat$c.dat[n.names, "ROI.Area"], digits=0))},
		if(!is.null(dat$c.dat[n.names, "perimeter"])){paste("Perimeter","", round(dat$c.dat[n.names, "perimeter"], digits=0))},
		if(!is.null(dat$c.dat[n.names, "circularity"])){paste("Circularity","", round(dat$c.dat[n.names, "circularity"], digits=3))}
		)
	,bty="n", cex=bcex)

	 
	par(xpd=FALSE)
	if(Plotit.both){
	par(xpd=T)
		if(!is.null(dat$der)){lines(dat$der[,n.names]~dat.t[-1,1], lwd=.01, col="paleturquoise4")}
		abline(h=0)
		lines(dat.t[,n.names]~dat.t[,1])
		points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
		par(xpd=F)
	}
	
	if(Plotit.trace){
		par(xpd=T)
		lines(dat.t[,n.names]~dat.t[,1])
		points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
		par(xpd=F)
	}
	
	if(info){
		x.name<-n.names
		levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])), "")
		levs.loc<-tapply(dat.t[,"Time"],as.factor(wr),mean)[levs]
		mtext(c("Peak Height","Area"), side=1, at=-max(dat.t[,1])*.05, line=c(1, 2), cex=bcex*.8)
		for(i in levs){
			max.name<-paste(i,".max", sep="")
			max.val<-round(dat$scp[x.name, max.name], digits=3)
			mtext(max.val, side=1, at=levs.loc[i], line=1, cex=bcex*.8)
			
			tot.name<-paste(i,".tot", sep="")
			tot.val<-round(dat$scp[x.name, tot.name], digits=3)
			mtext(tot.val, side=1, at=levs.loc[i], line=2, cex=bcex*.8)
		}?mtext
	}
	
	## Tool for adding rasterImages to plot
	img.dim<-dim(dat$img1)[1]
	zf<-20
	x<-dat$c.dat[n.names,"center.x"]
	left<-x-zf
	if(left<=0){left=0; right=2*zf}
	right<-x+zf
	if(right>=img.dim){left=img.dim-(2*zf);right=img.dim}
	
	y<-dat$c.dat[n.names,"center.y"]
	top<-y-zf
	if(top<=0){top=0; bottom=2*zf}
	bottom<-y+zf
	if(bottom>=img.dim){top=img.dim-(2*zf);bottom=img.dim}
	
	par(xpd=TRUE)
	
	#ymax<-max(dat.t[,n.names])*1.05
	#ymin<-min(dat.t[,n.names])*.95
	ymax<-max(ylim)
	ymin<-min(ylim)
	yrange<-ymax-ymin
	
	if(!is.null(dat$img1)){
		xleft<-min(dat.t[,1])-xinch(.2)
		xright<-min(dat.t[,1])-xinch(.2)
		ytop<-ymax-(yrange*.05)
		ybottom<-ymax-(yrange*.35)
		rasterImage(dat$img1[top:bottom,left:right,],xleft,ybottom,xright,ytop)
	}
	if(!is.null(dat$img2)){
		xleft<-max(dat.t[,1])*1.15
		xright<-max(dat.t[,1])*1.23
		ytop<-ymax-(yrange*.05)
		ybottom<-ymax-(yrange*.35)
		
		rasterImage(dat$img2[top:bottom,left:right,],xleft,ybottom,xright,ytop)
	}
	
	if(!is.null(dat$img3)){
		xleft<-max(dat.t[,1])*1.05
		xright<-max(dat.t[,1])*1.13
		ytop<-ymax-(yrange*.45)
		ybottom<-ymax-(yrange*.75)
		rasterImage(dat$img3[top:bottom,left:right,],xleft,ybottom,xright,ytop)
	}
	
	if(!is.null(dat$img4)){
		xleft<-max(dat.t[,1])*1.15
		xright<-max(dat.t[,1])*1.23
		ytop<-ymax-(yrange*.45)
		ybottom<-ymax-(yrange*.75)
		rasterImage(dat$img4[top:bottom,left:right,],xleft,ybottom,xright,ytop)
	}
		
	
}	
	

	
	
	
	