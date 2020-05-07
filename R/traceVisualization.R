#This peak func allows for multiple t.types to be plotted
#170515: added pts and lns: (logical)
#added dat.n for insertation of the name for the rd file
PeakFunc7 <- function(dat,n.names,t.type="t.dat",Plotit.trace=T,Plotit.both=F, info=T,lmain=NULL, bcex=.7, yvar=T, ylim.max=NULL, zf=40, pts=T, lns=T, levs=NULL, underline=T, dat.n=""){
    dat.name<-deparse(substitute(dat))
    if(dat.name=="dat"){dat.name<-dat.n
    }else{dat.name<-dat.name}	
    
    if(is.null(lmain)){
        lmain=n.names
    }else{lmain=lmain}
    if(class(t.type)=="character")
    {
        dat.select<-t.type
        dat.t<-dat[[dat.select]]
    }else{
        dat.select<-select.list(names(dat), multiple=T)
        dat.t<-dat[[dat.select]]
    }

    if(yvar){
        ymax<-max(dat.t[,n.names])*1.05
        ymin<-min(dat.t[,n.names])*.95
        yrange<-ymax-ymin
    }else{		
        if(is.null(ylim.max)){ylim.max<-1.4}else{ylim.max<-ylim.max}
        if(Plotit.trace){ylim <- c(-.1,ylim.max)}
        if(Plotit.both){ylim <- c(-.5,ylim.max)}
        ymin<-min(ylim)
        ymax<-max(ylim)
        yrange<-ymax-ymin
    }

    if(Plotit.trace){ylim <- c(ymin,ymax)}
    if(Plotit.both){ymin<- -.5;ylim <- c(ymin,ymax)}
    par(xpd=FALSE)
    xlim <- range(dat.t[,1]) # use same xlim on all plots for better comparison
    
    #   ylim <- range(intensity(s1))
    if(is.null(levs)){levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    }else{levs<-levs}
    par(mar=c(9,6.2,3.5,13), bty="n")
    plot(dat.t[,n.names]~dat.t[,1], main=lmain,xlim=xlim,ylim=ylim,xlab="", ylab="",pch="", cex=.5)

    #axis(3,tick=TRUE, outer=F )
    axis(1, at= seq(0, max(dat.t[,1]),10), tick=TRUE)
    
    # Tool for labeling window regions
    wr<-dat$w.dat[,"wr1"]
    #levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs]
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs]
    y1s <- rep(par("usr")[4],length(x1s))
    y2s <- rep(par("usr")[3],length(x1s))
    rect(x1s,y1s,x2s,y2s,col="grey95")
    
    
    # Tool for labeling cellular aspects, gfp.1, gfp.2, tritc, area 
    legend(x=par("usr")[1]-xinch(1.45), y=par("usr")[3]-yinch(.25), xpd=TRUE, inset=c(0,-.14),bty="n", cex=.7, legend=c(
        if(!is.null(dat$c.dat[n.names, "CGRP"])){paste("CGRP","",round(dat$c.dat[n.names,"CGRP"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp"])){paste("GFP","",round(dat$c.dat[n.names,"mean.gfp"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.1"])){paste("GFP.1","",round(dat$c.dat[n.names,"mean.gfp.1"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.2"])){paste("GFP.2","",round(dat$c.dat[n.names,"mean.gfp.2"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.start"])){paste("mean.gfp.start","",round(dat$c.dat[n.names,"mean.gfp.start"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.end"])){paste("mean.gfp.end","",round(dat$c.dat[n.names,"mean.gfp.end"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.immuno"])){paste("CGRP immunostain","",round(dat$c.dat[n.names,"mean.gfp.immuno"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.dapi"])){paste("DAPI","",round(dat$c.dat[n.names,"mean.dapi"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "IB4"])){paste("IB4","",round(dat$c.dat[n.names,"IB4"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.tritc"])){paste("IB4","",round(dat$c.dat[n.names, "mean.tritc"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.tritc.start"])){paste("IB4.start","",round(dat$c.dat[n.names, "mean.tritc.start"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.tritc.end"])){paste("IB4.end","",round(dat$c.dat[n.names, "mean.tritc.end"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.tritc.immuno"])){paste("NF200 immunostain","",round(dat$c.dat[n.names, "mean.tritc.immuno"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.cy5.start"])){paste("IB4.start","",round(dat$c.dat[n.names, "mean.cy5.start"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.cy5.end"])){paste("IB4.end","",round(dat$c.dat[n.names, "mean.cy5.end"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "area"])){paste("area","", round(dat$c.dat[n.names, "area"], digits=4))},
        if(!is.null(dat$c.dat[n.names, "ROI.Area"])){paste("area","", round(dat$c.dat[n.names, "ROI.Area"], digits=4))},
        #if(!is.null(dat$c.dat[n.names, "perimeter"])){paste("perimeter","", round(dat$c.dat[n.names, "perimeter"], digits=0))},
        if(!is.null(dat$c.dat[n.names, "circularity"])){paste("circularity","", round(dat$c.dat[n.names, "circularity"], digits=4))}
        )
    )
    
    legend(x=par("usr")[2]+xinch(.8), y=par("usr")[3]-yinch(.9), xpd=TRUE, inset=c(0,-.14), bty="n", cex=.7, legend=dat.name)

    
    #Adding binary scoring for labeling to plot
    par(xpd=TRUE)
    if(!is.null(dat$bin[n.names, "gfp.bin"])){text(y=par("usr")[4]+yinch(.5), x=par("usr")[2]+xinch(1.8), paste("GFP:",dat$bin[n.names,"gfp.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "tritc.bin"])){text(y=par("usr")[4]+yinch(.25), x=par("usr")[2]+xinch(1.8), paste("IB4 :",dat$bin[n.names,"tritc.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "cy5.bin"])){text(y=par("usr")[4]+yinch(.25), x=par("usr")[2]+xinch(1.8), paste("IB4 :",dat$bin[n.names,"cy5.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "drop"])){text(y=par("usr")[4]+yinch(0), x=par("usr")[2]+xinch(1.8), paste("Drop :",dat$bin[n.names,"drop"]), cex=.7)}


    # Tool for lableing window region information
    levs.loc<-tapply(dat$w.dat[,"Time"],as.factor(wr),mean)[levs]
    if(info){
        x.name<-n.names
        #levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])), "")
        mtext(c("max","snr"), side=3, at=-max(dat.t[,1])*.05, line=c(0, .7), cex=.6)
        for(i in 1:length(levs)){
            max.name<-paste(levs[i],".max", sep="")
            max.val<-round(dat$scp[x.name, max.name], digits=3)
            mtext(max.val, side=3, at=levs.loc[ levs[i] ], line=0, cex=.6)
            
            tot.name<-paste(levs[i],".snr", sep="")
            tot.val<-round(dat$scp[x.name, tot.name], digits=3)
            mtext(tot.val, side=3, at=levs.loc[ levs[i] ], line=.7, cex=.6)
        }
        
    # Tool for labeling the binary score
        #levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
        z<-t(dat$bin[n.names,levs])
        zz<-z==1
        zi<-attributes(zz)
        zzz<-which(zz, arr.ind=T)
        #levs<-zi$dimnames[[2]][zzz[,2]]
        levs1<-unique(as.character(row.names(zzz)))
        x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs1]
        x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs1]
        y1s <- rep(par("usr")[4],length(x1s))
        y2s <- rep(par("usr")[3],length(x1s))
        rect(x1s,y1s,x2s,y2s,col="grey80")
        #levs <- setdiff(unique(wr),"")
    }
    
    #text(dat.t[match(levs,wr),"Time"],c(ymin, ymin+(yrange*.2)),levs,pos=4,offset=0,cex=bcex)	
    #text(dat.t[match(levs,wr),"Time"],par("usr")[3],levs,pos=3,offset=-4.2,cex=bcex, srt=90)    
    levs_cex <- nchar(levs)
	levs_cex[ levs_cex <= 12*1.3  ] <- 1
	levs_cex[ levs_cex > 12*1.3  ] <- 12/levs_cex[ levs_cex>12*1.3  ]*1.3

    text(levs.loc,par("usr")[3],levs,pos=3,offset=-4.3,cex=levs_cex, srt=90)	

    if(Plotit.both){
        if(!is.null(dat$der)){lines(dat$der[,n.names]~dat.t[-1,1], lwd=.01, col="paleturquoise4")}
        par(xpd=T)
        abline(h=0)
        if(lns){lines(dat.t[,n.names]~dat.t[,1])
        }else{}
        if(pts){points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
        }else{}
        par(xpd=F)
    }
    
    if(Plotit.trace){
        par(xpd=T)
        if(lns){lines(dat.t[,n.names]~dat.t[,1])
        }else{}
        
        if(pts){points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
        }else{}
        
        par(xpd=F)
    }
    
    ##Tool for adding underline to plot
    if(underline){
        par(xpd=F)
        abline(h=min(dat.t[,n.names]), col="black")
        par(xpd=T)
    }else{}

    
    ## Tool for adding rasterImages to plot
    
    ###Finding the picture loaction of the cells
    if(!is.null(dat$img1)){
        if(is.null(zf)){zf<-20
        }else{zf<-zf}

        img.dim<-dim(dat$img1)[1]
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
    }
    ### Where to plot pictures
    #ymax<-max(dat.t[,n.names])*1.05
    #ymin<-min(dat.t[,n.names])*.95
    #yrange<-ymax-ymin
    

    
    ymax<-par("usr")[4]
    xmax<-par("usr")[2]
    if(!is.null(dat$img1)){
        img1<-dat$img1
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax+yinch(.8)
        ybottom<-ymax
        tryCatch(
            rasterImage(img1[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img1[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img2)){
        img2<-dat$img2
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax+yinch(.8)
        ybottom<-ymax
        tryCatch(
            rasterImage(img2[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img2[top:bottom,left:right],xleft,ybottom,xright,ytop))

    }

    if(!is.null(dat$img3)){
        img3<-dat$img3
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax
        ybottom<-ymax-yinch(.8)
        tryCatch(
            rasterImage(img3[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img3[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img4)){
        img4<-dat$img4
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax
        ybottom<-ymax-yinch(.8)
        tryCatch(
            rasterImage(img4[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img4[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }

    if(!is.null(dat$img5)){
        img5<-dat$img5
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax-yinch(.8)
        ybottom<-ymax-yinch(1.6)
        tryCatch(
            rasterImage(img5[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img5[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img6)){
        img6<-dat$img6
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax-yinch(.8)
        ybottom<-ymax-yinch(1.6)
        tryCatch(
            rasterImage(img6[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img6[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img7)){
        img7<-dat$img7
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax-yinch(1.6)
        ybottom<-ymax-yinch(2.4)
        tryCatch(
            rasterImage(img7[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img7[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img8)){
        img8<-dat$img8
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax-yinch(1.6)
        ybottom<-ymax-yinch(2.4)
        tryCatch(
            rasterImage(img8[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img8[top:bottom,left:right],xleft,ybottom,xright,ytop))

    }
    

    
    
}	

#' LinesEvery same as .4 but has image at begining of trace and moves to pic plot at >10 
#' multipi does not work on this.  Instead, if greater than 10, the traces are plotted as a portrait orientation
#' Also window labels are rotated on axis and place on the bottom of the plot
#' I am also adding two more images to the left side of the plot
#' 171009 added underline.  Helps to show irreversibility
#' 171031 added dat.n for the name fotthe experiment
#' @export
LinesEvery.5 <- function(dat,m.names, img="img1",channel=NULL,pic.plot=TRUE,zf=NULL, t.type="mp.1", snr=NULL,lmain="",cols="black", levs=NULL, levs.cols="grey90", values=NULL,plot.new=T,sf=1,lw=1,bcex=1,p.ht=7,p.wd=10, lns=T, pts=F, underline=T,dat.n=NULL){
    #require(RColorBrewer)
    dat.name<-deparse(substitute(dat))
    if(dat.name=="dat" | dat.name == "tmp.rd" | dat.name == "tmp_rd"){
        dat.name<-dat.n
    }else{
        dat.name<-dat.name
    }
    
    #Trace Selector if t.type is empty.  t.type must be character input
    if(class(t.type)=="character"){
        t.dat<-dat[[t.type]]# if trace type is empty select the data, you would like your trace to be
    }else{
        t.type<-menu(names(dat));t.dat<-dat[[t.type]]
    }
    
    m.names <- intersect(m.names,names(t.dat))
    xseq <- t.dat[,1]
    #upper ylimit 
    hbc <- length(m.names)*sf+max(t.dat[,m.names])

    #Selecting multiple images
    if(is.null(img)){
        img.l<-select.list(grep("img",names(dat), value=T), multiple=T)
    }else{
        img.l<-img
    }

    if(length(m.names) > 0){
        #For pdf output
            #if(is.null(pdf.name))
            #	{dev.new(width=14,height=8)}
            #else
                #{if(length(grep("\\.pdf",pdf.name))>0){pdf(pdf.name,width=p.wd,height=p.ht)}else{png(pdf.name,width=1200,height=600)}}
        ## Tool for addind value tags displayed on the right side of trace
        #See line 3016 for where values come into play
        #values<-c("area", "mean.gfp.start", "mean.gfp.end" "mean.tritc.start", "mean.tritc.end")
            if(is.null(values)){
                values<-c("area")
            }else{values<-values}

        ## Tool for color labeleing
        ## Tool for single color labeling
            if(cols=="brew.pal"){
                #cols <- rainbow(length(m.names),start=.55)
                require(RColorBrewer)
                cols <-brewer.pal(8,"Dark2")
                cols <- rep(cols,ceiling(length(m.names)/length(cols)))
                cols <- cols[1:length(m.names)]
            }
            if(cols=="rainbow"){
                cols<-rainbow(length(m.names),start=.7,end=.1)
            }
            if(cols=="topo"){
                cols<-topo.colors(length(m.names))
            }else{		
                cols<-cols
                cols <- rep(cols,ceiling(length(m.names)/length(cols)))
                cols <- cols[1:length(m.names)]
            }
            if(plot.new){
                if(length(m.names)>10){dev.new(width=10+length(img)+(length(values)*.6),height=12)}
                else(dev.new(width=10+length(img)+(length(values)*.6),height=8))
            }
 
            xinch(length(img))
            par(xpd=FALSE,mai=c(2,.5+(.5*length(img.l)), 1, 0.6*length(values)), bty="l")
            plot(xseq,t.dat[,m.names[1]],ylim=c(0,hbc),xlab="",main=lmain,type="n",yaxt="n",xlim=c(min(xseq)-1.5,max(xseq)), ylab="")#-sf
            #axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
            #axis(2, 1.4, )
            #Label cell names
            text(rep(0,length(m.names))-xinch(.1),seq(1,length(m.names))*sf+t.dat[1,m.names],m.names, cex=.5,col=cols,pos=3)

        ## Tool for adding window region labeling
            if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
            }else{levs<-levs}
            wr<-dat$w.dat$wr1
            if(length(wr) > 0){
                x1s <- tapply(dat$w.dat[,1],as.factor(wr),min)[levs]
                x2s <- tapply(dat$w.dat[,1],as.factor(wr),max)[levs]
                y1s <- rep(par("usr")[3],length(x1s))
                y2s <- rep(par("usr")[4],length(x1s))
                rect(x1s,y1s,x2s,y2s,col=levs.cols,border="black")
                par(xpd=TRUE)
                #text(x1s-xinch(.1),par("usr")[3]-yinch(1),levs,cex=.8*bcex, srt=90)	
                #dat$t.dat[match(levs,wr),"Time"]
                levs.loc<-tapply(dat$w.dat[,"Time"],as.factor(wr),mean)[levs]
				levs_cex <- nchar(levs)
				levs_cex[ levs_cex <= 12*1.3 ] <- 1
				levs_cex[ levs_cex > (12*1.3) ] <- 12/levs_cex[ levs_cex>(12*1.3) ]*1.3
                text(levs.loc,par("usr")[3],levs,pos=3,offset=-4.3,cex=levs_cex, srt=90)	
                par(xpd=FALSE)
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
                if(underline){abline(h=min(ypos), col="black")}else{}
            }
            par(xpd=TRUE)
            
        ## Tool for adding Value info on right side of trace
            placement<-seq(0,length(values),.5)
            digits<-c(0,rep(4,length(values)))
            text(max(xseq)+xinch(placement[1:length(values)]), par("usr")[4]+yinch(.2), pos=4, values,cex=bcex*.75, srt=30)
            
            for(i in 1:length(values)){
                if(!is.null(dat$c.dat[m.names, values[i]])){
                rtag<-values[i]
                rtag <- round(dat$c.dat[m.names,rtag], digits=digits[i])
                text(
                    rep(max(xseq)+xinch(placement[i]),length(m.names)),
                    seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],
                    paste(rtag),
                    cex=.65*bcex,
                    col=cols,
                    pos=4)
                }
            }
            
        ##Tool for adding images to the left side of the plot
        if(is.null(zf)){
            zf<-20
        }else{zf<-zf}
        
        pic.pos<-list()
        for(i in 1:length(m.names)){
            ypos<-t.dat[1,m.names[i]]+i*sf
            pic.pos[[i]]<-ypos
        }

        xinchseq1<-seq(1,5,.5)
        xinchseq2<-seq(.5,5,.5)
        
        if(is.null(channel)){channel<-rep(list(c(1:3)),length(img.l))
        }else{channel<-channel}

        for(j in 1:length(img.l)){
            for(i in 1:length(m.names)){		
                img.dim<-dim(dat$img1)[1]
                x<-dat$c.dat[m.names[i],"center.x"]
                left<-x-zf
                if(left<=0){
                    left=0
                    right=2*zf
                }
                
                right<-x+zf
                if(right>=img.dim){
                    left=img.dim-(2*zf)
                    right=img.dim
                }
                    
                y<-dat$c.dat[m.names[i],"center.y"]
                top<-y-zf
                if(top<=0){
                    top=0
                    bottom=2*zf
                }
                bottom<-y+zf
                
                if(bottom>=img.dim){
                    top=img.dim-(2*zf)
                    bottom=img.dim
                }
                
                par(xpd=TRUE)
                xleft<-min(dat$t.dat[,1])-xinch(xinchseq1[j])
                xright<-min(dat$t.dat[,1])-xinch(xinchseq2[j])
                ytop<-pic.pos[[i]]+yinch(.25)
                ybottom<-pic.pos[[i]]-yinch(.25)
                
                tryCatch(
                    rasterImage(dat[[img.l[j]]][top:bottom,left:right,channel[[j]]],xleft,ybottom,xright,ytop),
                    error=function(e) rasterImage(dat[[img.l[j]]][top:bottom,left:right],xleft,ybottom,xright,ytop)
                )
            }
        }
    }

        tryCatch(
			legend(x=par("usr")[2]-xinch(1.2), y=par("usr")[3]-yinch(1.6), xpd=TRUE, inset=c(0,-.14), bty="n", cex=.7, legend=dat.name),
		error=function(e) NULL)


        #if(!is.null(pdf.name))
        #{dev.off()}

    #return(pic.pos)
}


#' intereact: LOGICAL; 
#' TRUE select cell groups to work though and return list of groups of cells
#' FALSE only plot out the groups, and dont return group of cells

#' #region.group: Select a region to group the cells around.  Brings up option to select region to group around
#' 170403 bp logical: lets you choose whether to boxplot

#' 170508 Allows to select the trace you would like to use for grouping with option:
#' t.type:input character

#' 170605:  Adding a drop function to this.  It will automatically update the RD.file.  I need something to drop cells much faster
#' @export
LinesStack.2<- function(dat,m.names=NULL,t.type=NULL,lmain="", interact=T, region.group=T,levs=NULL, plot.new=TRUE,bcex=.7, sf=1.1, subset.n=NULL, img=NULL,bp.param=NULL, bp=F, bp.pts=F){
    #graphics.off()
    #
    if(is.null(img)){img<-dat$img1}
    
    # If a list of cells is not input, then look at all cells
    if(is.null(m.names)){
        dropped.cells<-cellzand(dat$bin, "drop",1)
        m.names<-setdiff(dat$c.dat$id, dropped.cells)
    }else{
        dropped.cells<-cellzand(dat$bin, "drop",1)
        m.names<-setdiff(m.names, dropped.cells)
    }
        
    if(is.null(subset.n)){subset.n<-as.numeric(select.list(as.character(c(5,10,15,20,25,30))))}
    if(plot.new){
        if(subset.n>=10){
            dev.new(width=14,height=10)
            }
        else{dev.new(width=14,height=10)}
        
        linesstack.win<-dev.cur()
    }
    if(length(m.names)>subset.n){
        
        if(is.null(t.type)){t.dat<-dat$t.dat}
        else{t.dat<-dat[[t.type]]}
        
        wr<-dat$w.dat[,2]
        if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")}
        else{levs<-levs}
        m.names <- intersect(m.names,names(t.dat))
        hbc <-(max(t.dat[,m.names])+subset.n)*sf
        #hbc <- (subset.n*(.8*sf))+max(t.dat[,m.names])
        
        xseq <- t.dat[,1]
        library(RColorBrewer)
        par(mai=c(2,1,1,1))
        
        ylim<-c(-.1,hbc)
        #ylim <- c(-.1,2.5)
        
        plot(xseq,
            t.dat[,m.names[1]],
            ylim=ylim,
            xlab="",
            ylab='',
            main=lmain,
            type="n",
            xlim=c(min(xseq)-1.5,max(xseq)+10),
            bty='n'
        )#-sf
        #axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
        
        ## Tool for adding window region labeling
        if(length(wr) > 0){
            #levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
            x1s <- tapply(xseq,as.factor(wr),min)[levs]
            x2s <- tapply(xseq,as.factor(wr),max)[levs]
            y1s <- par('usr')[3]
            y2s <- par('usr')[4]
            rect(x1s,y1s,x2s,y2s,col="grey90",border="black")
            levs.loc<-tapply(dat$w.dat[,"Time"],as.factor(wr),mean)[levs]
            par(xpd=T)
            text(levs.loc, par("usr")[3] - yinch(.5),levs,cex=bcex, srt=90)	
        }
        ## Tool for adding line and point plot for all lines
            #matlines(xseq, blc[,m.names], col=rgb(0,0,0,3, maxColorValue=100), lwd=.01)
            #matpoints(xseq, blc[,m.names], col=rgb(0,0,0,3, maxColorValue=100), pch=16, cex=.03)
        
        #cols <- rainbow(length(m.names),start=.55)
     
        library(cluster)	
        ## To select data within the experiment to group around
        if(region.group){
            dev.new(width=15, height=10)
            LinesEvery.5(dat, sample(dat$c.dat$id)[1:5], plot.new=F, lmain="Click to Select region to Groups Cells", t.type="t.dat", img="img1")
            #LinesEvery.4(dat, sample(row.names(dat$c.dat)[1:5]), plot.new=F, lmain="Click to Select region to Groups Cells", img=dat$img1)
            b.xseq<-locator(n=2, type="o", pch=15, col="red")$x
            dev.off()
            x.min<-which(abs(t.dat$Time-b.xseq[1])==min(abs(t.dat$Time-b.xseq[1])))
            x.max<-which(abs(t.dat$Time-b.xseq[2])==min(abs(t.dat$Time-b.xseq[2])))
            
            pam5 <- pam(t(t.dat[x.min:x.max,m.names]),k=subset.n)
            s.names <- row.names(pam5$medoids)
            pam5.tab <- table(pam5$clustering)
            #tags <- paste(paste("#",names(pam5.tab),sep=""),as.vector(pam5.tab),sep=":")
            group.means<-list()
            group.names<-list()
            for(i in 1:subset.n){
                x.names<-names(which(pam5$clustering==i, arr.ind=T))
                group.names[[i]]<-x.names
                group.means[i]<-paste(
                tryCatch(round(mean(dat$c.dat[x.names, "area"]), digits=0),error=function(e) NULL),
                "\u00b1",
                tryCatch(round(sd(dat$c.dat[x.names, "area"]), digits=1),error=function(e) NULL))#," : ",
                #tryCatch(round(mean(dat$c.dat[x.names, "mean.gfp"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.gfp"]), digits=0),error=function(e) NULL)," : ",	
                #tryCatch(round(mean(dat$c.dat[x.names, "mean.tritc"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.tritc"]), digits=0),error=function(e) NULL), sep="")
                #tryCatch(round(sd(dat$c.dat[x.names, "area"]), digits=0),"\u00b1",error=function(e) NULL)
            }
    
        }else{
            library(cluster)
            pam5 <- pam(t(t.dat[,m.names]),k=subset.n)
            s.names <- row.names(pam5$medoids)
            pam5.tab <- table(pam5$clustering)
            #tags <- paste(paste("#",names(pam5.tab),sep=""),as.vector(pam5.tab),sep=":")
            group.means<-list()
            group.names<-list()
            for(i in 1:subset.n){
                x.names<-names(which(pam5$clustering==i, arr.ind=T))
                group.names[[i]]<-x.names
                group.means[i]<-paste(
                tryCatch(round(mean(dat$c.dat[x.names, "area"]), digits=0),error=function(e) NULL),
                "\u00b1",
                tryCatch(round(sd(dat$c.dat[x.names, "area"]), digits=1),error=function(e) NULL))				
                #round(mean(dat$c.dat[x.names, "mean.gfp"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.gfp"]), digits=0)," : ",	
                #round(mean(dat$c.dat[x.names, "mean.tritc"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.tritc"]), digits=0)
                #adding standard deviation
                #"\u00b1",round(sd(dat$c.dat[x.names, "area"]), digits=0), sep="")
            }
        }			
        tags <- paste(as.vector(pam5.tab),":",group.means)
        info<-pam5$clustering
        
        ## Tool For adding color to selected Traces
        cols <-brewer.pal(8,"Dark2")
        cols <- rep(cols,ceiling(length(s.names)/length(cols)))
        cols <- cols[1:length(s.names)]

        ## Tool for adding labeling for single line within stacked traces
        par(xpd=T)
        dev.set(which=linesstack.win)
        for(i in 1:length(s.names)){
            if(length(group.names[[i]])>=2){
                matlines(xseq, (t.dat[,group.names[[i]]]+i)*sf, col=rgb(0,0,0,20, maxColorValue=100), lwd=.01)
                lines(xseq, apply(t.dat[,group.names[[i]]],1,mean)+i*sf, col=cols[i], lwd=.2)
                points(xseq, apply(t.dat[,group.names[[i]]],1,mean)+i*sf, col=cols[i], pch=16, cex=.02)
                text(x=min(t.dat[,1]), y=t.dat[1,s.names[i]]+i*sf, labels=i, col=cols[i], pos=2, cex=bcex)
                text(x=max(t.dat[,1]), y=t.dat[nrow(dat$t.dat),s.names[i]]+i*sf, labels=tags[i], col=cols[i], pos=4, cex=bcex)
            }else{
                lines(xseq, t.dat[,group.names[[i]]]+i*sf, col=cols[i], lwd=.2)
                points(xseq, t.dat[,group.names[[i]]]+i*sf, col=cols[i], pch=16, cex=.02)
                text(x=min(t.dat[,1]), y=t.dat[1,s.names[i]]+i*sf, labels=i, col=cols[i], pos=2, cex=bcex)
                text(x=max(t.dat[,1]), y=t.dat[nrow(dat$t.dat),s.names[i]]+i*sf, labels=tags[i], col=cols[i], pos=4, cex=bcex)
                }
        }
        if(region.group){
            text(mean(b.xseq),
                par('usr')[4]+yinch(3),
                "Cluster Region",
                col='blue',
                font=2,
                cex=2
            )            
            par(xpd=F)
            abline(v=b.xseq, col="blue", lwd=2)
        }else{}

        par(xpd=F)

    #### Tool for adding boxplot
        if(bp){
                par(xpd=T)
                dev.current<-dev.cur()
                if(is.null(bp.param)){
                    #dat.select<-"c.dat"
                    #bp.param<-c(
                    #grep("area",names(dat$c.dat),value=T),
                    ##tryCatch(grep("mean.gfp",names(dat$c.dat)),error=function(e) NULL),
                    #grep("mean.gfp",names(dat$c.dat),value=T),
                    #grep("mean.tritc",names(dat$c.dat),value=T))
                    
                    #cols<-c("blue", "darkgreen","red")
                #}else{
                    dat.select<-select.list(names(dat))
                    bp.param<-as.character(select.list(names(dat[[dat.select]]), multiple=T))
                    cols<-NULL
                }else{
                dat.select<-"c.dat"
                bp.param<-bp.param
                }
                
                #for(i in 1:length(group.names)){
                    #if(length(group.names[[i]])>5){
                    #	xleft<-max(t.dat[,1])+xinch(.3)
                    #	xright<-xleft+xinch(1)*length(bp.param)
                    #	y<-(apply(t.dat[nrow(t.dat),group.names[[i]]],1,mean)+i)*sf
                    #	ybottom<- y-yinch(.5)
                    #	ytop<-y+yinch(.5)

                    #	bp.img<-bpfunc.3(dat,group.names[[i]],dat.select, bp.param, print.out=T, cols=cols, bcex=bcex)
                    #	dev.set(dev.current)
                    #	rasterImage(bp.img,xleft, ybottom, xright, ytop)
                    #}else{}
                    
                    #170509 How to create a new window with these boxplots
                dev.new(width=length(bp.param), height=subset.n)
                bp.win<-dev.cur()
                par(mfrow=c(subset.n,1))
                group.names.rev<-rev(group.names)
                    
                for(i in 1:length(group.names.rev)){
                        par(mar=c(0,0,0,0))
                        plot(0,0)
                        dim<-par("usr")
                        xleft<-par("usr")[1]
                        xright<-par("usr")[2]
                        ybottom<- par("usr")[3]
                        ytop<-par("usr")[4]
                        bp.img<-bpfunc.3(dat,group.names.rev[[i]],dat.select, bp.param, print.out=T, cols=cols, bcex=bcex, bp.pts=bp.pts)
                        dev.set(bp.win)
                        rasterImage(bp.img,xleft, ybottom, xright, ytop)
                        text(xleft+xinch(.1), 0, subset.n-i+1, cex=2)
                    }
            
            }
        }
        
            if(interact){
                continue<-select.list(c("yes", "no"))
                if(continue=="yes"){
                    while(i!=length(s.names)+1){
                        i<-scan(n=1)
                        if(i>length(s.names)| i==0){i<-length(s.names)+1}
                        else{
                            assesment.selection<-select.list(c("Trace.Click","LinesEvery","LinesStack", "drop"))
                            
                            if(assesment.selection=="Trace.Click"){
                                Trace.Click.dev(dat,names(which(info==i, arr.ind=T)))
                            }
                            
                            if(assesment.selection=="LinesEvery"){
                                number.to.display<-as.numeric(select.list(as.character(c(5,10,20))))
                                LinesEvery.5(dat,sample(names(which(info==i, arr.ind=T)))[1:number.to.display], img, pic.plot=T, lmain=i, plot.new=T, col="black")
                            }
                            
                            if(assesment.selection=="LinesStack"){
                            
                                LinesStack.2(dat,names(which(info==i, arr.ind=T)),bp=F,lmain=i, interact=T, region.group=T,levs=NULL, plot.new=TRUE,bcex=.7, img=dat$img1, t.type="mp.1")
                            }
                            if(assesment.selection=="drop"){
                                rd.namels2 <- as.character(substitute(dat))
                                dat$bin[names(which(info==i, arr.ind=T)), "drop"]<-1
                                assign(rd.namels2, dat, envir=.GlobalEnv)
                                print(paste("You Dropped Group",i))
                            }
                        }
                    }
                }
                #return(pam5$clustering)	

            }
        
    #dev.off(which=linesstack.win)
    return(group.names)
}

