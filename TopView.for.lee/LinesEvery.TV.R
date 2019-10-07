LinesEvery.TV <- function(dat,m.names, img=dat$img1,pic.plot=TRUE, multi.pic=T,zf=NULL, t.type="mp", snr=NULL,lmain="",cols=NULL, levs=NULL, levs.cols="grey90",  m.order=NULL,rtag=NULL,rtag2=NULL,rtag3=NULL,plot.new=T,sf=1,lw=2,bcex=.6,p.ht=7,p.wd=10, lns=T, pts=F)
{
	require(png)
	if(class(t.type)=="character"){t.dat<-dat[[t.type]]}# if trace type is empty select the data, you would like your trace to be
	else{t.type<-menu(names(dat));t.dat<-dat[[t.type]]}
	wr<-dat$w.dat[,"wr1"]
	if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")}
	else{levs<-levs}
    m.names <- intersect(m.names,names(t.dat))
    xseq <- t.dat[,1]
	hbc <- length(m.names)*sf+max(t.dat[,m.names])
	hb <- ceiling(hbc)
	library(RColorBrewer)
    
	if(length(m.names) > 0)
    {
		if(!is.null(m.order))
		{	
			tmp<-dat$c.dat[m.names,]
			n.order<-tmp[order(tmp[,m.order]),]
			m.names <- row.names(n.order)
		}
### Picture Plotting!		
	#if(XY.plot==T){cell.zoom.2048(dat, cell=m.names,img=img, cols="white",zoom=F, plot.new=T)}
	## Tool for color labeleing
		if(is.null(cols)){
			cols <-brewer.pal(8,"Dark2")
			cols <- rep(cols,ceiling(length(m.names)/length(cols)))
			cols <- cols[1:length(m.names)]
		} 
	## Tool for single color labeling
		else {cols<-cols
			cols <- rep(cols,ceiling(length(m.names)/length(cols)))
			cols <- cols[1:length(m.names)]
		}
		
		if(multi.pic){
			if(plot.new){
				if(length(m.names)>10){dev.new(width=16,height=10);layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(10,6), heights=c(6,6))}
				else(dev.new(width=12,height=8))
			}
			else{
				if(length(m.names)>10){layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(10,6), heights=c(6,6))}
				}
		}else{dev.new(width=12,height=8)}
		par(xpd=TRUE,mar=c(4,3,2,2), bty="l")
        plot(xseq,t.dat[,m.names[1]],ylim=c(0,hbc),xlab="Time (min)",main=lmain,type="n", xaxt="n",yaxt="n",xlim=c(min(xseq)-1.5,max(xseq)))#-sf
		bob<-dev.cur()
        axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
	    #axis(2, 1.4, )
		text(rep(0,length(m.names)),seq(1,length(m.names))*sf+t.dat[1,m.names],m.names, cex=.5,col=cols,pos=2)

	## Tool for adding window region labeling
		if(length(wr) > 0){
            #levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
            x1s <- tapply(dat$w.dat[,1],as.factor(wr),min)[levs]
            x2s <- tapply(dat$w.dat[,1],as.factor(wr),max)[levs]
            y1s <- rep(-.3,length(x1s))
            y2s <- rep(hbc+.2,length(x1s))
            rect(x1s,y1s,x2s,y2s,col=levs.cols,border="black")
            cpx <- xseq[match(levs,wr)+round(table(wr)[levs]/2,0)]
            offs <- nchar(levs)*.5
#			par(xpd=TRUE)
            text(dat$t.dat[match(levs,wr),"Time"],rep(c((sf*.7)*.5,(sf*.7),(sf*.7)/5),length=length(levs)),levs,pos=4,offset=0,cex=bcex*.8)#,offset=-offs}
			#par(xpd=FALSE)
		}
	
	## Tool for adding line, point and picture to the plot
		for(i in 1:length(m.names)){
			ypos<-t.dat[,m.names[i]]+i*sf
			if(lns){lines(xseq,ypos, lty=1,col=cols[i],lwd=lw)}
			if(pts){points(xseq,ypos,pch=16,col=cols[i],cex=.3)}
			if(!is.null(snr)){
				pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
				pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
				points(xseq[pp1],t.dat[pp1,m.names[i]]+i/10,pch=1,col=cols[i])
				points(xseq[pp2],t.dat[pp2,m.names[i]]+i/10,pch=0,col=cols[i])
			}
		}

	if(is.null(img)){
	img.p<-dat[[select.list(grep("img",names(dat), value=T))]]
	if(is.null(img.p)){img.p<-dat$img1}
	}else{img.p<-img}
	if(is.null(zf)){zf<-20}else{zf<-zf}
	
	#if(pic.plot==TRUE & length(m.names)<=10){
	if(pic.plot==TRUE){
	if(length(m.names)<=100){

		pic.pos<-list()
		for(i in 1:length(m.names)){
			ypos<-t.dat[1,m.names[i]]+i*sf
			pic.pos[[i]]<-ypos}

			for(i in 1:length(m.names)){		
				#if(dat$bin[m.names[1],"mean.gfp.bin"]!=1 & dat$bin[m.names[1],"mean.tritc.bin"]!=1){img.p<-dat$img.gtd #if the cell is neither red or green, then make the img to plot img.gtd
				#}else{img.p<-img} 
				#img.p<-img
				
				img.dim<-dim(dat$img1)[1]

				x<-dat$c.dat[m.names[i],"center.x"]
				left <- x-zf
				if(left<=0){left=0; right=2*zf}
				right<- x+zf
				if(right>=img.dim){left=img.dim-(2*zf);right=img.dim}
				
				y<-dat$c.dat[m.names[i],"center.y"]
				top<-y-zf
				if(top<=0){top=0; bottom=2*zf}
				bottom<-y+zf
				if(bottom>=img.dim){top=img.dim-(2*zf);bottom=img.dim}
				
				#par(xpd=TRUE)
				xleft<-min(dat$t.dat[,1])-xinch(1)
				xright<-min(dat$t.dat[,1])-xinch(.5)
				ytop<-pic.pos[[i]]+yinch(.25)
				ybottom<-pic.pos[[i]]-yinch(.25)
				
				tryCatch(rasterImage(img.p[top:bottom,left:right,],xleft,ybottom,xright,ytop),error=function(e) rasterImage(img.p[top:bottom,left:right],xleft,ybottom,xright,ytop))
			}
		}
	else{
		par(mar=c(0,0,0,0))
		plot(0,0,xlim=c(0,6), ylim=c(0,6), xaxs="i",yaxs="i", xaxt='n', yaxt='n')
		tmp.img<-multi.pic.zoom.2(dat, m.names,img=img.p, labs=T, zf=zf, cols=cols)
		dev.set(bob) # FUCK THIS!
		rasterImage(tmp.img, 0,0,6,6)
	}
	}
	}
		#if(!is.null(pdf.name))
        #{dev.off()}

#return(pic.pos)
}
