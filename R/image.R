#' Add images to the RD.experiment file
#' 
#' @param dat RD.experiment input
#' @export 
ImageFiller<-function(dat){
    require(png)
    potential.images<-list.files(pattern='png')
    print(potential.images)
    tryCatch(bringToTop(-1), error=function(e)NULL)
    print("######################")
    print("These are the images you have the option of selecting")
    print("Now select the images to fill in for image 1 to 8")
    
    for(i in 1:8){
        image.to.add<-select.list(list.files(pattern='png'),title=paste('img',i,sep=''))
        
        if(image.to.add==""){dat[[paste('img',i,sep='')]]<-NULL
        }else{
            dat[[paste('img',i,sep='')]]<-png::readPNG(image.to.add)
        }
    }
    return(dat)
}

#' Add images to the RD.experiment file
#' 
#' @param dat RD.experiment input
#' @export 
ImageFillerv2 <- function(dat, img_name_vec = NULL){

    if( is.null(img_name_vec) ){
        img_name_vec<-c(
            "bf.gfp.tritc.start.png",
            "gfp.tritc.start.ci.ltl.rs.png",
            "tritc.start.ci.ltl.rs.png",
            "gfp.start.ci.ltl.rs.png",
            "bf.start.lab.png",
            "fura2.png",
            "fura2.divide.start.png",
            "roi.img.png")
        image_question <- T
    }else{
        image_question <- F
    }
    # Add images
    img_list<-list()
    for( j in 1:length(img_name_vec) ){
        dat[[ paste0("img",j) ]] <- tryCatch(png::readPNG(img_name_vec[j]), error=function(e)NULL)
    }
    if(image_question == T){
        cat('\nThese are the images I have attempted to load for you\nIf any are NULL, and want to add different images say yes to the \nnext question. You will be asked to select a png image for each loaction.\n\n')
        cat(img_name_vec, sep='\n')
        cat(str(img_list))

        cat('\nDO YOU WANT DIFFERENT IMAGES[y,n]?\n')
        img_reselect <- scan(n=1,what='character')
        if( img_reselect=='y' ){
            cat("\nAlright buddy I am going to give you options if you don't\nwant any image there just go ahead and press 0\n\n")
            png_imgs <- list.files(pattern='png')
            for( j in 1:8 ){
                cat("\nWhat do you want for image ", j, '\n')
                selection <- menu(png_imgs)
                if(selection==0){
                    dat[[paste0("img",j)]] <- NULL
                }else{
                    dat[[paste0("img",j)]] <- png::readPNG(png_imgs[selection])
                }
                cat('\nI have added ', png_imgs[selection],' to position ',j,'\n')
            }
            }
    }
    return(dat)    
}

#' Function to View all cells in grid, zoomed in.
#' @export 
cell.view <- function(dat, cell=NULL,img=NULL,  zoom=TRUE, cols=NULL,lmain="", bcex=.8, labs=T, plot.new=T, cell.name=T){
    if(plot.new){dev.new()}
    require(png)
    par(mar=c(0,0,1,0))
    x<-dat$c.dat[,"center.x"]
    y<-dat$c.dat[,"center.y"]
    cell.x<-dat$c.dat[cell,"center.x"]
    cell.y<-dat$c.dat[cell,"center.y"]

    if(is.null(img)){img<-dat$img1}
    else{img<-img}
    if(is.null(cols)){cols="white"}
    else{cols=cols}
    img.dimy<-dim(img)[1]
    img.dimx<-dim(img)[2]

    plot(0, 0, xlim=c(0,img.dimx),ylim=c(img.dimy,0), main=lmain,xaxs="i", yaxs="i", xlab="Pixels", ylab="Pixels")
    rasterImage(img, 0, img.dimy, img.dimx, 0)

    if(labs){
        if(!is.null(cell)){
            points(cell.x, cell.y, col=cols, pch=0, cex=2)
            text(cell.x, cell.y, labels=cell, col=cols, pos=2, cex=bcex)
        }
        else{
            points(x, y, col=cols, pch=4, cex=2)
            text(x, y, labels=dat$c.dat[,1], col=cols, pch=0, pos=2, cex=bcex)
        }
    }

    if(zoom==TRUE & length(cell)>1){
        cell.1<-row.names(dat$c.dat[order(dat$c.dat$center.x),])
        cell<-intersect(cell,cell.1)
        multi.pic.zoom(dat,cell,img) 
    }
    }

    cell.zoom.640.480<-function(dat, img=NULL, cell=NULL, zoom=NULL, cols=NULL, labs=T, plot.new=T, cell.name=T)
    {
    if(plot.new){dev.new()}
    require(png)
    require(zoom)
    par(mar=c(0,0,0,0))
    x<-dat$c.dat[,"center.x"]
    y<-dat$c.dat[,"center.y"]
    cell.x<-dat$c.dat[cell,"center.x"]
    cell.y<-dat$c.dat[cell,"center.y"]

    if(is.null(img)){img<-dat$img1}
    else{img<-img}
    if(is.null(cols)){cols="white"}
    else{cols=cols}

    plot(0, 0, xlim=c(0,640),ylim=c(480,0), xaxs="i", yaxs="i", xlab="Pixels", ylab="Pixels")
    rasterImage(img, 0, 480, 640, 0)

    if(labs){
    if(!is.null(cell)){
        points(cell.x, cell.y, col=cols )
        text(cell.x, cell.y, labels=cell, col=cols, pos=2, cex=.8)
        }
    else{
        points(x, y, col=cols)
        text(x, y, labels=dat$c.dat[,1], col=cols, pos=2, cex=.5)
    }}

    if(!is.null(zoom)){
    zoomplot.zoom(x=cell.x, y=cell.y, fact=zoom)
    }
    else{zm()}
}

#' View Individual cell picture
#' @export 
multi.pic.zoom<-function(dat, m.names, img, labs=T,plot.new=T, zf=20){
    col.row<-ceiling(sqrt(length(m.names)))
    
    if(plot.new)
    {
        dev.new()
        par(mfrow=c(col.row, col.row))
        par(mar=c(0,0,0,0))
    }
    else{
        par(mfrow=c(col.row, col.row))
        par(mar=c(0,0,0,0))
    }	
    
    m.names<-rev(m.names)
    for(i in 1:length(m.names)){		

        img.dim<-dim(img)
        x<-dat$c.dat[m.names[i],"center.x"]
        left<-x-zf
        if(left<=20){left=0; right=zf}
        right<-x+zf
        if(right>=img.dim-zf){left=img.dim-zf;right=img.dim}
                    
        y<-dat$c.dat[m.names[i],"center.y"]
        top<-y-zf
        if(top<=20){top=0; bottom=zf}
        bottom<-y+zf
        if(bottom>=img.dim-zf){top=img.dim-zf;bottom=img.dim}

        par(xpd=TRUE)
        xleft<-0
        xright<-20
        ytop<-0
        ybottom<-20
        plot(c(xright, xleft), c(ytop, ybottom), ylim=c(20,0) ,xaxs="i", yaxs="i", axes=F)
        rasterImage(img[top:bottom,left:right,],xleft,ybottom,xright,ytop)
        text(4,1.5, m.names[i], col="white", cex=.8)
        box(lty = 1, col = "white",lwd=2)
        text(16.5, 2, labels=dat$c.dat[m.names[i], "area"], col="white")

        if(labs){
            points(x=10,y=10, type="p", pch=3, cex=2,col="white")
            text(16.5, 2, labels=dat$c.dat[m.names[i], "ROI.Area"], col="white")
            text(16.5, 3.5, labels=dat$c.dat[m.names[i], "mean.gfp.1"], col="green")
            text(16.5, 3.5, labels=dat$c.dat[m.names[i], "mean.gfp"], col="green")
            text(16.5, 3.5, labels=dat$c.dat[m.names[i], "CGRP"], col="green")
            text(16.5, 5, labels=dat$c.dat[m.names[i], "mean.tritc"], col="red")
            text(16.5, 5, labels=dat$c.dat[m.names[i], "IB4"], col="red")
            text(16.5, 6.5, labels=dat$c.dat[m.names[i], "mean.dapi"], col="blue")
        }
    }
        
}
