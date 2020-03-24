#tmp is an RD object, x.names are the cell ids to investiage
#pad is the extra amount of image to select around the cell e.g. 1 = at cell bondaries 1.05 = 5% extra
#stain.name is the stain to display ("tritc","gfp","dapi") anything else defaults to yellow ROI boundaries
#title1 will be the title of the grid selection window.

#160819
# Select grid anything that isn't ib4 or gfp is ordered according to what the user wants
# also, soom window has c.dat values addeds
# problem: cell data does not match cell selected.
# correct name of collumn saved to in bin table
# previously you had t saving to c.dat, since it is a binary score i added to the bin table,
# I also changed the name of the collumns, to what i have a few functions dependent upon; mean.gfp.bin, mean.tritc.bin

#160829
# roireview now has mutliple options for dropping, and the function allows you to drop cells 
# as many rounds as needed


#161114
#added kevins subset.ne to ROIreview to allow subset of cells to add
#updated selectgrid to recognize the dimensions of the images
SelectGrid <- function(tmp,x.names,pad=1.05,stain.name="tritc",title1="SelectRed",window.h=7,window.w=7,l.col="red")
{
	#img1 is all colors
	#img2 is blue and green
	#img3 is blue and red
	#img4 has yellow roi lines

	imgs <- grep("img",names(tmp),value=T)	
	imgs.yes <- rep(F,length(imgs))
	for(i in 1:length(imgs)){imgs.yes[i] <- length(dim(tmp[[imgs[i]]])) == 3}
	imgs <- imgs[imgs.yes]
	if(length(imgs) < 1){stop("no image data")}	
	imgs.yes <- rep(F,length(imgs))
	for(i in 1:length(imgs)){imgs.yes[i] <- dim(tmp[[imgs[i]]])[3] == 3}
	imgs <- imgs[imgs.yes]
	
	if(length(imgs) < 1){stop("no image data")}
	img.rgb <- data.frame(name=imgs)
	img.rgb["r"] <- 0
	img.rgb["g"] <- 0
	img.rgb["b"] <- 0
	
	for(j in 1:nrow(img.rgb))
	{
		img.rgb[j,"r"] <- mean(tmp[[imgs[j]]][,,1])
		img.rgb[j,"g"] <- mean(tmp[[imgs[j]]][,,2])
		img.rgb[j,"b"] <- mean(tmp[[imgs[j]]][,,3])			
	}
	#set the channel to use and subtract the others. red=1, green=2, blue=3
	#also select the best image.
	img.red <- imgs[which.max(img.rgb[,"r"]-img.rgb[,"g"]-img.rgb[,"b"])]
	img.green <- imgs[which.max(img.rgb[,"g"]-img.rgb[,"r"]-img.rgb[,"b"])]
	img.blue <- imgs[which.max(img.rgb[,"b"]-img.rgb[,"r"]-img.rgb[,"g"])]
	img.yellow <- imgs[which.max(img.rgb[,"b"]+img.rgb[,"r"]-img.rgb[,"g"])]
	if(is.element(stain.name,c("tritc","gfp","dapi")))
	{
		sn <- grep(stain.name,names(tmp$c.dat),ignore.case=T,value=T)[1]
		if(is.null(sn)){stop("no stain value data")}
		x.names <- x.names[order(tmp$c.dat[x.names,sn])]
		if(stain.name=="tritc")
		{
			img.name <- imgs[which.max(img.rgb[,"r"]-img.rgb[,"g"]-img.rgb[,"b"])]
			chn <- 1
		}
		if(stain.name=="gfp")
		{
			img.name <- imgs[which.max(img.rgb[,"g"]-img.rgb[,"r"]-img.rgb[,"b"])]
		  	chn <- 2		
		}
		if(stain.name=="dapi")
		{
			img.name <- imgs[which.max(img.rgb[,"b"]-img.rgb[,"r"]-img.rgb[,"g"])]
			chn <- 3
		}
		
		
		img <- tmp[[img.name]]
		img.dat <- img[,,chn]
		for(i in setdiff(c(1,2,3),chn)){gt.mat <- img.dat < img[,,i];img.dat[gt.mat] <- 0} 
		#single.img <- tmp$img1
	}else
	{
		img.name <- intersect("img.bl",imgs)
		if(is.null(img.name)){img.name <- imgs[which.max(img.rgb[,"b"]+img.rgb[,"r"]-img.rgb[,"g"])]}
		sn <- select.list(names(tmp$c.dat), title="Select variable to sort") # choose variable to sort
		sort.dir<-select.list(c("TRUE", "FALSE"), title="Decreasing?")
		x.names <- x.names[order(tmp$c.dat[x.names,sn], decreasing=sort.dir)]
		img <- tmp[[img.name]]
		img.dat <- (img[,,1]+img[,,2])/2
		med.r <- .99
		med.b <- .99
		if(sum(as.vector(img[,,1]) > med.r)==0){med.r <- quantile(as.vector(img[,,1]),probs=c(.95))[1]}
		if(sum(as.vector(img[,,2]) > med.b)==0){med.b <- quantile(as.vector(img[,,2]),probs=c(.95))[1]}
		img.dat[img[,,1] < med.r] <- 0
		img.dat[img[,,2] < med.b] <- 0		

		#single.img <- tmp$img4
	}
	
	#set up two devices
	graphics.off()
	dev.new(height=window.h,width=window.w,canvas="black",title="SingleCell")
	dev.single <- dev.cur()
	op <- par(mar=c(0,0,0,0))	
	plot(c(0,1),c(0,1),xaxt="n",yaxt="n",type="n",ylab="",xlab="")	
	
	dev.new(height=window.w,width=window.h,canvas="black",title=title1)
	dev.grid <- dev.cur()
	op <- par(mar=c(0,0,0,0))	
	plot(c(0,1),c(0,1),xaxt="n",yaxt="n",type="n",ylab="",xlab="")	
	xn <- length(x.names)
	num.grid <- xn+3
	nr <- floor(sqrt(num.grid))
	nc <- ceiling((num.grid)/nr)
	mtx <- max(nr,nc)
	dx <- seq(0,1,length.out=(mtx+1))[-1]
	sl <- (dx[2]-dx[1])/2
	dx <- dx-sl
	all.x <- as.vector(matrix(rep(dx,mtx),byrow=F,ncol=mtx))
	all.y <- as.vector(matrix(rep(dx,mtx),nrow=mtx,byrow=T))
	
	#zf<-(sqrt(tmp$c.dat[x.names,"area"])/pi)*pad
	zf<-rep(5,length(x.names)) #maintain consitent window sizing
	x <- tmp$c.dat[x.names,"center.x"]
	y <- tmp$c.dat[x.names,"center.y"]
	
	img.dim<-dim(img.dat)[1]
	img.dim.max<-img.dim-1
	
	zf[zf > x] <- x[zf > x]
	zf[zf > y] <- y[zf > y]
	zf[x+zf > img.dim] <- img.dim-x[x+zf > img.dim]
	zf[y+zf > img.dim] <- img.dim-y[y+zf > img.dim]
	
	img.left<-x-zf
	img.left[img.left < 1] <- 1
	img.right<-x+zf
	img.right[img.right > img.dim] <- img.dim
	img.top<-y-zf
	img.top[img.top < 1] <- 1
	img.bottom<-y+zf
	img.bottom[img.bottom > img.dim] <- img.dim

	img.bottom[img.top >= img.bottom & img.top < img.dim] <- img.top[img.top >= img.bottom]+1
	img.right[img.left >= img.right & img.left < img.dim] <- img.left[img.left >= img.right]+1

	img.top[img.top == img.dim] <- img.dim.max
	img.left[img.left == img.dim] <- img.dim.max
		
	for(i in 1:xn)
	{
		xl <- all.x[i]-sl*.9
		xr <- all.x[i]+sl*.9
		xt <- all.y[i]-sl*.9
		xb <- all.y[i]+sl*.9
		#rasterImage(tmp$img1[img.bottom[i]:img.top[i],img.left[i]:img.right[i],],xl,xb,xr,xt)
		rasterImage(img.dat[img.bottom[i]:img.top[i],img.left[i]:img.right[i]],xl,xb,xr,xt)
	}
	fg <- rep("black",length(all.x))
	fg[1:xn] <- "grey"
	cexr <- sl/.04
	symbols(all.x,all.y,squares=rep(sl*1.9,length(all.x)),add=T,inches=F,fg=fg,lwd=cexr)
	text(all.x[xn+1],all.y[xn+1],"Done",col="white",cex= cexr)
	text(all.x[xn+2],all.y[xn+2],"All",col="white",cex= cexr)
	text(all.x[xn+3],all.y[xn+3],"None",col="white",cex= cexr)

	#first click defines the split
	all.sel <- rep(0,xn)
	names(all.sel) <- x.names	
	not.done=TRUE
	click1 <- locator(n=1)
	dist <- sqrt((click1$x[[1]]-all.x)^2 + (click1$y[[1]]-all.y)^2)
	sel.i <- which.min(dist)
	if(sel.i == xn+1){not.done=FALSE;return(all.sel)}
	if(sel.i == xn+2){all.sel[1:xn] <- 1;fg[1:xn] <- l.col}
	if(sel.i == xn+3){all.sel[1:xn] <- 0;fg[1:xn] <- "grey"}
	if(sel.i <= xn)
	{
	dev.set(which=dev.single)
#	rasterImage(single.img[img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],0,0,1,1,interpolate=F)
	rasterImage(tmp[[img.red]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],0,0,.5,.5,interpolate=F)
	rasterImage(tmp[[img.green]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],0,.5,.5,1,interpolate=F)
	rasterImage(tmp[[img.blue]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],.5,0,1,.5,interpolate=F)
	rasterImage(tmp[[img.yellow]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],.5,.5,1,1,interpolate=F)
	abline(h=.5,col="grey")
	abline(v=.5,col="grey")
	
	dev.set(which=dev.grid)	
	neg.i <- 1:max((sel.i-1),1) 
	all.sel[neg.i] <- 0
	pos.i <- sel.i:xn	
	all.sel[pos.i] <- 1
	fg[neg.i] <- "grey"
	fg[pos.i] <- l.col
	}
	while(not.done)
	{
		symbols(all.x,all.y,squares=rep(sl*1.9,length(all.x)),add=T,inches=F,fg=fg,lwd=cexr)
		click1 <- locator(n=1)
		dist <- sqrt((click1$x[[1]]-all.x)^2 + (click1$y[[1]]-all.y)^2)
		sel.i <- which.min(dist)
		if(sel.i == xn+1){not.done=FALSE;return(all.sel)}
		if(sel.i == xn+2){all.sel[1:xn] <- 1;fg[1:xn] <- l.col}
		if(sel.i == xn+3){all.sel[1:xn] <- 0;fg[1:xn] <- "grey"}
		if(sel.i <= xn)
		{
		dev.set(which=dev.single)
#		rasterImage(single.img[img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],0,0,1,1,interpolate=F)
		rasterImage(tmp[[img.red]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],0,0,.5,.5,interpolate=F)
		rasterImage(tmp[[img.green]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],0,.5,.5,1,interpolate=F)
		rasterImage(tmp[[img.blue]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],.5,0,1,.5,interpolate=F)
		rasterImage(tmp[[img.yellow]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],.5,.5,1,1,interpolate=F)
		abline(h=.5,col="grey")
		abline(v=.5,col="grey")

		dev.set(which=dev.grid)	
		if(all.sel[sel.i] ==0)
		{
			all.sel[sel.i] <- 1
			fg[sel.i] <- l.col
		}
		else
		{
			all.sel[sel.i] <- 0
			fg[sel.i] <- "grey"
		}
		}
	}		
	
}

#three tests Drop (confirm), Red (confirm) and Green (confirm)
#return and RD object with the changes made to c.dat and bin
ROIreview <- function(tmp,x.names=NULL,pad=2,wh=7,hh=7, subset.n=NA)
{
	dice <- function(x, n,min.n=10)
	{
		x.lst <- split(x, as.integer((seq_along(x) - 1) / n))
		x.i <- length(x.lst)
		if(length(x.lst[x.i]) < min.n & x.i > 1)
		{
			x.lst[[x.i-1]] <- c(x.lst[[x.i-1]],x.lst[[x.i]])
			x.lst <- x.lst[1:(x.i-1)]
		}
		return(x.lst)
	}
	
if(is.null(x.names)){x.names <- row.names(tmp$c.dat)}
x.names <- x.names[tmp$bin[x.names,"drop"]==0]
if(is.na(subset.n) | subset.n > length(x.names)){subset.n=length(x.names)}
subset.list <- dice(x.names,subset.n,subset.n/4)
for(x.names in subset.list){

	continue<-"Yes"	
	while(continue!="No"){
		d.names <- SelectGrid(tmp,x.names,pad,"circularity","SelectDrops",window.h=hh,window.w=wh)
		d1.names <- names(d.names[d.names==1])
		if(length(d1.names) > 5)
		{
			ad.names<-row.names(tmp$bin)[tmp$bin[,"drop"]==1]# find already dropped cells
			d1.names<-setdiff(d1.names, ad.names)# only observe cells that haven't been dropped
			d1.names <- SelectGrid(tmp,d1.names,pad,"area","ConfirmDrops",window.h=hh,window.w=wh) 
			d1.names <- names(d1.names)[d1.names==1]
			if(length(d1.names) > 0){tmp$bin[d1.names,"drop"] <- 1}
		}
		continue<-select.list(c("Yes", "No"), title="Continue Dropping Cells?")	
		}
		
		
		r.names <- SelectGrid(tmp,x.names,pad,"tritc","SelectRed",window.h=hh,window.w=wh)
		r1.names <- names(r.names[r.names==1])
		q1 <- 1:floor(length(r1.names)*.25)
		r2.names <- r1.names[q1]
		if(length(r2.names) > 5)
		{
			r2.names <- SelectGrid(tmp,r2.names,pad*2,"tritc","ConfirmRed",window.h=hh,window.w=wh)
			r.names[names(r2.names)] <- r2.names
		}
		tmp$bin[names(r.names),"mean.tritc.bin"] <- r.names

		r.names <- SelectGrid(tmp,x.names,pad,"gfp","SelectGreen",window.h=hh,window.w=wh)
		r1.names <- names(r.names[r.names==1])
		q1 <- 1:floor(length(r1.names)*.25)
		r2.names <- r1.names[q1]
		if(length(r2.names) > 5)
		{
			r2.names <- SelectGrid(tmp,r2.names,pad*2,"gfp","ConfirmGreen",window.h=hh,window.w=wh)
			r.names[names(r2.names)] <- r2.names
		}
				tmp$bin[names(r.names),"mean.gfp.bin"] <- r.names
		}

		levs<-setdiff(unique(as.character(tmp$w.dat$wr1)),"")
		cells<-list()
		neuron.response<-select.list(levs, title="What defines Neurons?", multiple=T)
		neurons<-cellz(tmp$bin,neuron.response, 1)
		drop<-cellz(tmp$bin, "drop", 1)
		neurons<-setdiff(neurons,drop)
		pf<-apply(tmp$bin[,c("mean.gfp.bin", "mean.tritc.bin")],1,paste, collapse="")
		tmp$bin["lab.pf"]<-as.factor(pf)
		lab.groups<-unique(tmp$bin$lab.pf)
			
		cells<-list()
		for(i in lab.groups){
			x.names<-cellz(tmp$bin[neurons,], "lab.pf", i)
			cells[[i]]<-x.names
		}
			
		glia.response<-select.list(c(levs, "none"), title="What defines glia?", multiple=T)
		if(glia.response!="none"){
			drop<-cellz(tmp$bin, "drop", 1)
			glia<-cellz(tmp$bin,glia.response, 1)
			glia<-setdiff(glia,drop)
			cells[["000"]]<-setdiff(glia, neurons)
			} 
			else {cells[["000"]]<-setdiff(row.names(tmp$c.dat), neurons)}
			tmp$cells<-cells
		
		return(tmp)			
			
	}

