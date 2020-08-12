#' This peak func allows for multiple t.types to be plotted
#' 170515: added pts and lns: (logical)
#' added dat.n for insertation of the name for the rd file
#' @export
PeakFunc7 <- function(dat,cell,t.type="t.dat",Plotit.trace=T, info=T,lmain=NULL, bcex=.7, yvar=T, ylim.max=NULL, zf=40, pts=T, lns=T, levs=NULL, underline=T, dat.n=""){
    # Plot ontop always unless you dont need it
    par(xpd = T)

    dat.name<-deparse(substitute(dat))
    if(dat.name=="dat"){
        dat.name<-dat.n
    }else{
        dat.name<-dat.name
    }	
    
    if(is.null(lmain)){
        lmain=cell
    }else{
        lmain=lmain
    }
    
    if(class(t.type)=="character"){
        dat.select<-t.type
        dat.t<-dat[[dat.select]]
    }else{
        dat.select<-select.list(names(dat), multiple=T)
        dat.t<-dat[[dat.select]]
    }

    # Set the y limits for the plotting
    if(yvar){
        ymax<-max(dat.t[,cell])*1.05
        ymin<-min(dat.t[,cell])*.95
    }else{		
        if(is.null(ylim.max)){
            ylim.max <- 1.4
        }
        ylim <- c(-0.1, ylim.max)
        ymin <- min(ylim)
        ymax <- max(ylim)
    }
    
    ylim <- c(ymin,ymax)
    xlim <- range(dat.t[,1]) # use same xlim on all plots for better comparison
    
    if(is.null(levs)){
        levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    }

    # Plot it up
    par(mar=c(9,6.2,3.5,13), bty="n")
    plot(
        1, 
        type='n',
        main=lmain,
        xlim=xlim,
        ylim=ylim,
        xlab="", 
        ylab="",
        xaxt = 'n'
    )

    axis(1, at= seq(0, max(dat.t[,1]),10), tick=TRUE)
    
    # Tool for labeling window regions
    wr <- dat$w.dat[,"wr1"]
    #levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs]
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs]
    y1s <- rep(par("usr")[4],length(x1s))
    y2s <- rep(par("usr")[3],length(x1s))

    # Colors the windows regions based on if the window is scored as a 1 or a 0
    if(info){
        cols <- ifelse(dat$bin[cell,levs] == 1, "grey80", "grey95")
    }else{
        cols <- "grey95"
    }
    rect(
        x1s,
        y1s,
        x2s,
        y2s,
        col=cols
    )
    
    # Tool for labeling cellular aspects, gfp.1, gfp.2, tritc, area 
    cellLabs <- c(
        "CGRP",                 "CGRP",
        "mean.gfp",             "GFP",
        "mean.gfp.1",           "GFP.1",
        "mean.gfp.2",           "GFP.2",
        "mean.gfp.start",       "mean.gfp.start",
        "mean.gfp.end",         "mean.gfp.end",
        "mean.gfp.immuno",      "CGRP immunostain",
        "mean.dapi",            "DAPI",
        "IB4",                  "IB4",
        "mean.tritc",           "IB4",
        "mean.tritc.start",     "IB4.start",
        "mean.tritc.end",       "IB4.end",
        "mean.tritc.immuno",    "NF200 immunostain",
        "mean.cy5.start",       "IB4.start",
        "mean.cy5.end",         "IB4.end",
        "area",                 "area",
        "ROI.Area",             "area",
        "perimeter",            "perimeter",
        "circularity",          "circularity")

    dim(cellLabs) <- c(2,length(cellLabs)/2)
    cellLabs <- t(cellLabs)
    
    legendAdds <- c()
    for(i in 1:dim(cellLabs)[1]){
        if( !is.null(dat$c.dat[cell, cellLabs[i,1] ]) ){
            toAdd <- paste0(cellLabs[i,2],": ",round(dat$c.dat[ cell, cellLabs[i,1] ],digits=4))
            legendAdds <- c(legendAdds, toAdd)
        }
    }

    legend(
        x=par("usr")[1]-xinch(1.45), 
        y=par("usr")[3]-yinch(.25), 
        xpd=TRUE, 
        inset=c(0,-.14),
        bty="n", 
        cex=.7, 
        legend=legendAdds
    )
    
    # Adds the rd name
    legend(
        x=par("usr")[2]+xinch(.8), 
        y=par("usr")[3]-yinch(.9), 
        xpd=TRUE, 
        inset=c(0,-.14), 
        bty="n", 
        cex=.7, 
        legend=dat.name
    )

    #Adding binary scoring for labeling to plot
    par(xpd=TRUE)
    binLabs <- c(
        "gfp.bin",      "GFP",
        "tritc.bin",    "IB4",
        "cy5.bin",      "IB4",
        "drop",         "Drop"
    )
    dim(binLabs) <- c(2, length(binLabs)/2)
    binLabs <- t(binLabs)

    binLabsToAdd <- c()
    for(i in 1:dim(binLabs)[1]){
        if( !is.null(dat$bin[cell, binLabs[i,1] ]) ){
            toAdd <- paste0(binLabs[i,2],": ", dat$bin[cell, binLabs[i,1] ])
            binLabsToAdd <- c(binLabsToAdd, toAdd)
        }
    }

    legend(
        x=par("usr")[2]+xinch(1.8), 
        y=par("usr")[4]+yinch(.5), 
        xpd=TRUE, 
        inset=c(0,-.14), 
        bty="n", 
        cex=.7, 
        legend=binLabsToAdd
    )

    # Tool for lableing window region information
    levs.loc <- tapply(dat$w.dat[,"Time"], as.factor(wr), mean)[levs]
    if(info){
        xLoc <- par('usr')[1] - xinch(.01)
        yLoc <- par('usr')[4]
        toAdd <- c('max', 'tot', 'snr')
        bcex <- 0.55
        for(j in 1:length(toAdd)){
            yLocJ <- yLoc + (yinch(.09) * j)
            text(
                xLoc, 
                yLocJ,
                toAdd[j],
                cex = bcex
            )

            valNames <- paste(levs, toAdd[j], sep=".")
            val <- apply(dat$scp[cell, valNames], 1, round, digits=2)

            text(
                x = levs.loc[ levs ],
                y = rep(yLocJ, length(val)),
                val,
                cex = bcex
            )
        }    
    }
    
    # Adding the window labels
    shrinkFactor <- 1
    maxLength <- 12
    levs_cex <- nchar(levs)
	levs_cex[ levs_cex <= maxLength * shrinkFactor ] <- 1
	levs_cex[ levs_cex > maxLength * shrinkFactor  ] <- maxLength/levs_cex[ levs_cex> maxLength * shrinkFactor  ] * shrinkFactor
    
    text(
        levs.loc,
        par("usr")[3] - yinch(.38),
        levs,
        adj = 1,
        cex=levs_cex, 
        srt=90
    )

    # Adding the uncertainty  to the windows
    tryCatch({
        val <- apply(dat$uncMat[cell, levs], 1, round, digits=2)
        
        text(
            x = levs.loc[ levs ],
            y = rep(yLoc, length(val)),
            val,
            cex = bcex
        )

        yLoc <- par('usr')[3] + (yinch(.09))
        text(
            xLoc, 
            yLoc,
            'unc',
            cex = bcex
        )
    },error=function(e) NULL)

    # Add lines or points
    if(lns){
        lines(dat.t[,cell]~dat.t[,1])
    }
    
    if(pts){
        points(dat.t[,cell]~dat.t[,1], pch=16, cex=.3)
    }

    ##Tool for adding underline to plot
    if(underline){
        par(xpd=F)
        abline(h=min(dat.t[,cell]), col="black")
        par(xpd=T)
    }

    ## Tool for adding rasterImages to plot
    for(i in 1:9){
        if( !is.null(dat[[paste0('img',i)]]) ){
            break
        }
    }
    # If the loops never breaks, there are no images to raster
    if(i != 9){
        # This gathers the cells location
        imgToUse <- paste0('img', i)
        if(!is.null(dat[[imgToUse]])){
            img.dim <- dim(dat$img1)[1]
            x<-dat$c.dat[cell,"center.x"]
            y <- dat$c.dat[cell, "center.y"]
            if(is.null(zf)){
                zf<-20
            }else{
                zf<-zf
            }

            left <- x-zf
            right <- x+zf
            top <- y-zf
            bottom<-y+zf

            if(left <= 0){
                left = 0
                right = 2*zf
            }
            
            if(right >= img.dim){
                left = img.dim-(2*zf)
                right=img.dim
            }
            
            if(top<=0){
                top=0
                bottom=2*zf
            }

            if(bottom >= img.dim){
                top=img.dim-(2*zf)
                bottom=img.dim
            }       
        }
        
        desDims <- c(
            1,1,
            1,2,
            2,1,
            2,2,
            3,1,
            3,2,
            4,1,
            4,2
        )
    
        dim(desDims) <- c(2, length(desDims)/2)
        desDims <- t(desDims)
        desDims <- desDims - 1
        
        yMax <- par("usr")[4]
        xMax <- par("usr")[2]
        
        for(i in 1:8){
            tryCatch({
                imgName <- paste0('img', i)
                xLeft <- xmax + (xinch(0.8) * desDims[i,2])
                xRight <- xLeft + xinch(0.8)
                
                yBottom <- yMax + (yinch(0.8) * -desDims[i,1])
                yTop <- yBottom + yinch(0.8)

                tryCatch(
                    rasterImage(
                        dat[[imgName]][top:bottom,left:right,],
                        xLeft,
                        yBottom,
                        xRight,
                        yTop
                    )
                    ,error=function(e){
                        rasterImage(
                            dat[[imgName]][top:bottom,left:right],
                            xLeft,
                            yBottom,
                            xRight,
                            yTop
                        )
                    }
                )
            },error = function(e) 'uh')
        }
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

