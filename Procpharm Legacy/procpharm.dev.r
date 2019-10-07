#XYtrace improvements
# this function needs regional selection
# there should be two options for this
# 1. Self selection
# 2. regional selection

# Self selection
#open and name all windows
dat<-RD.170104.22.m.p1
img<-dat$img1
cols<-"red"


dev.new(height=4,width=12)
trace.window<-dev.cur()
	
dev.new(width=8, height=8)
pic.window<-dev.cur()
	
dev.new(height=8,width=12)
lines.window<-dev.cur()
	
#plot image in the window
if(is.null(img)){img<-dat$img1}
img.dim.x<-dim(img)[1]	
img.dim.y<-dim(img)[2]
dev.set(which=pic.window)
par(mar=c(0,0,0,0))
plot(0, 0, xlim=c(0,img.dim.x),ylim=c(img.dim.y,0),xaxs="i", yaxs="i",col=cols,pch=".")
if(!is.null(img)){
	rasterImage(img, 0, img.dim.y, img.dim.x, 0)
}
#points(cell.coor[,1],cell.coor[,2],col=cols,pch=0)}


x.sel<-locator(n=2, type="l", col="Red")$x
y.sel<-locator(n=2, type="l", col="Red")$y

rect(x.sel[1],y.sel[2],x.sel[2],y.sel[1])

# now i need to clsoe the window and open a new one with the same type of selection
x.size<-abs(x.sel[1]-x.sel[2])
y.size<-abs(y.sel[1]-y.sel[2])

#width vs height ratio
w.v.h<-y.size/x.size
h<-8*w.v.h
dev.off()
dev.new(width=8, height=h)
pic.window<-dev.cur()

plot(0, 0, xlim=c(x.sel[1],x.sel[2]),ylim=c(y.sel[1],y.sel[2]),xaxs="i", yaxs="i",col=cols,pch=".")
rasterImage(img[x.sel[1]:x.sel[2],y.sel[1]:y.sel[2], ], x.sel[1], y.sel[2], x.sel[2], y.sel[1])






if(is.null(img)){img<-dat$img1}
img.dim.x<-dim(img)[1]	
img.dim.y<-dim(img)[2]

	
if(is.null(img)){img<-dat$img1}
img.dim.x<-dim(img)[1]	
img.dim.y<-dim(img)[2]
dev.set(which=pic.)
par(mar=c(0,0,0,0))


plot(0, 0, xlim=c(0,img.dim.x),ylim=c(img.dim.y,0),xaxs="i", yaxs="i",col=cols,pch=".")
if(!is.null(img)){
	rasterImage(img, 0, img.dim.y, img.dim.x, 0)
	points(cell.coor[,1],cell.coor[,2],col=cols,pch=0)
	}
else{
points(cell.coor[,1],cell.coor[,2], col=cols, cex=dat$c.dat[,area]/200)
points(cell.coor[,1],cell.coor[,2],col=cols, pch=4)}

# working with
dat<-RD.170104.22.m.p1





# Function allows for selection and deselection of cells to build stacked traces
XYtrace <- function(dat, cell=NULL, img=NULL, cols=NULL, labs=F, y.var=T)
{
	graphics.off()
	x.coor<-grep("\\.x",names(dat$c.dat), value=T, ignore.case=T)
	if(length(x.coor)>1){x.coor<-"center.x"}
	y.coor<-grep("\\.y",names(dat$c.dat), value=T, ignore.case=T)
	if(length(x.coor)>1){y.coor<-"center.y"}
	area<-grep("area",names(dat$c.dat), value=T, ignore.case=T)
	
	lab1<-grep("cgrp",names(dat$c.dat), value=T, ignore.case=T)
	if(length(lab1)==0){lab1<-grep("gfp.1",names(dat$c.dat), value=T, ignore.case=T)}
	
	lab1.1<-grep("cgrp",names(dat$c.dat), value=T, ignore.case=T)
	if(length(lab1.1)==0){lab1.1<-grep("gfp.2",names(dat$c.dat), value=T, ignore.case=T)}
	
	lab2<-grep("ib4",names(dat$c.dat), value=T, ignore.case=T)
	if(length(lab2)==0){lab2<-grep("tritc",names(dat$c.dat), value=T, ignore.case=T)}
	
	if(is.null(cell)){cell<-row.names(dat$c.dat)}
	else{cell<-cell}
	cell.coor<-dat$c.dat[cell,c(x.coor, y.coor)]

	
	
	# select the names of the collumns containing coordinates
	levs <- unique(dat$w.dat[,"wr1"])
	levs<-setdiff(levs, "")
	if(labs==TRUE){
	if(is.null(cols)){cols="orangered1"} else{cols=cols}}
	pch=16
	
	dev.new(height=4,width=12)
	trace.window<-dev.cur()
	
	dev.new(width=8, height=8)
	pic.window<-dev.cur()
	
	dev.new(height=8,width=12)
	lines.window<-dev.cur()
	
	lmain<-"XY ROI"
	
	
	if(is.null(img)){img<-dat$img1}
	img.dim.x<-dim(img)[1]	
	img.dim.y<-dim(img)[2]
	dev.set(dev.list()[2])
	par(mar=c(0,0,0,0))
	plot(0, 0, xlim=c(0,img.dim.x),ylim=c(img.dim.y,0),xaxs="i", yaxs="i",col=cols,pch=".")
	if(!is.null(img)){rasterImage(img, 0, img.dim.y, img.dim.x, 0);points(cell.coor[,1],cell.coor[,2],col=cols,pch=0)}
	else{
	points(cell.coor[,1],cell.coor[,2], col=cols, cex=dat$c.dat[,area]/200)
	points(cell.coor[,1],cell.coor[,2],col=cols, pch=4)}

	i <- identify(cell.coor[,1],cell.coor[,2],n=1,plot=F, col=NA, tolerance=0.05)
	i.names<-row.names(dat$c.dat[cell,])[i]
	while(length(i) > 0)
	{	#selected name of cell
		s.names <- row.names(dat$c.dat[cell,])[i]
		dev.set(dev.list()[1])
		if(y.var){PeakFunc6(dat,s.names, Plotit.both=F)}
		else{PeakFunc5(dat,s.names, Plotit.both=T)}

		dev.set(dev.list()[2])
		# If a cell is selected, that has already been selected, 
		# then remove that cell from the list
		if(length(intersect(i.names,s.names))==1){
			i.names<-setdiff(i.names,s.names)
			points(cell.coor[s.names,1],cell.coor[s.names,2],col="gray70",pch=0,cex=2.4)
			points(cell.coor[i.names,1],cell.coor[i.names,2],col="red",pch=0,cex=2.4)	
		}
		# If it han't been selected, then add it to the list
		else{i.names<-union(i.names,s.names)
		points(cell.coor[i.names,1],cell.coor[i.names,2],col="red",pch=0,cex=2.4)}
		
		if(length(i.names)>=1){
			dev.set(dev.list()[3])
			LinesEvery.5(dat,m.names=i.names, plot.new=F, img=img, cols="black",sf=.2)}				
			dev.set(dev.list()[2])
			i <- identify(cell.coor[,1],cell.coor[,2],labels=dat$c.dat[cell,1],n=1,plot=T, pch=0,col="white", tolerance=0.05)
		}
	dev.off()
	graphics.off()
	return(row.names(dat$c.dat[i.names,]))
	   
}








