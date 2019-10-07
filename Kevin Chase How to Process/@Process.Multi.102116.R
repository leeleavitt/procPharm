# Standard multiexperiment analysis.

require(cluster)

#1 get a list of all the RD files in the current directory
c.fn <- list.files(pattern="^RD.*\\.Rdata$")
#load all the rd files
for(i in c.fn){load(i)}
#2 convert file names to rd names
c.rd <- sub("^.*RD","RD",sub("\\.Rdata","",c.fn));setdiff(c.rd,ls())

#3 despike each rd file and save the trace data as $mp in the RD 
for(k in c.rd)
{	
	tmp <- get(k)
	wts <- tmp$t.dat
	for(i in 1:5) #run the despike 5 times.
	{
		wt.mn3 <- Mean3(wts)
		wts <- SpikeTrim2(wts,1,-1)
		print(sum(is.na(wts))) #this prints out the number of points removed should be close to 0 after 5 loops.
		wts[is.na(wts)] <- wt.mn3[is.na(wts)]
	}
	tmp$mp <- wts
	assign(k,tmp)
}


#3a. #brief look at windows.
#not quite so extensive and slow as RoughReview
for(k in c.rd)
{
	tmp <- get(k)
	if(is.element("wr2",names(tmp$w.dat))){tmp$w.dat[,"wr1"] <- tmp$w.dat[,"wr2"]}
	k.name <- grep("^kc",names(tmp$bin),value=T,ignore.case=T)[1]
	x.names <- sample(row.names(tmp$c.dat[tmp$bin[,k.name]==1,]),min(50,nrow(tmp$c.dat)))
	LinesSome(tmp$mp,,x.names,as.character(tmp$w.dat[,"wr1"]),unique(as.character(tmp$w.dat[,"wr1"])),lmain=k,subset.n=20)
}

#3d check for data spacing errors
for(k in c.rd)
{
	tmp <- get(k)
	t <- tmp$t.dat[,1]*60
	dt <- t[-1]-t[-length(t)]
	dev.new();
	plot(t[-1],dt,main=k)
}
	

#4a. expand the window regions a bit. 5 ticks left and 5 ticks right
for(k in c.rd)
{
	tmp <- get(k)
	tmp.wr <- tmp$w.dat[,"wr1"]
	tmp.wr.1 <- c(tmp.wr[-(1:5)],rep("",5))
	tmp.wr[tmp.wr==""] <- tmp.wr.1[tmp.wr==""]
	tmp.wr.1 <- c(rep("",5),tmp.wr)[1:length(tmp.wr)]
	tmp.wr[tmp.wr==""] <- tmp.wr.1[tmp.wr==""]
	tmp$w.dat[,"wr1"] <- tmp.wr
	assign(k,tmp)
}



#5 set the thresholds for scoring and run the automatic scoring
sm <- 2 #smooth window size set
ws <- 30 #window peak size
snr.lim <- 4 #signal to noise threshold
hab.lim <- .05 #height above baseline threshold
for(i in c.rd)
{
	tmp <- get(i)
	t.dat <- tmp$t.dat
	tmp$t.dat <- tmp$mp #use the despiked data
	tlevs <- setdiff(unique(tmp$w.dat[,"wr1"]),"")
	tmp.pcp <- ProcConstPharm(tmp,sm,ws,"SNIP") #note I use "SNIP" baseline now instead of TopHat
	tmp.scp <- ScoreConstPharm(tmp$t.dat,tmp.pcp$blc,tmp.pcp$snr,snr.lim,hab.lim,tmp$w.dat[,"wr1"],sm)
	tmp.bin <- bScore(tmp.pcp$blc,tmp.pcp$snr,snr.lim,hab.lim,tlevs,tmp$w.dat[,"wr1"])
	tmp$bin <- tmp.bin[,tlevs]
	tmp$bin["drop"] <- 0
	tmp$scp <- tmp.scp
	tmp$t.dat <- t.dat
	assign(i,tmp)
}

#6 select traces to drop from analysis if necessary
#if you have a control pulse or a wash pulse
#name it here if not leave it NULL
#for each RD you will be shown a set of "bad" traces.
#select the ones you want to drop and click "finish"
#once you click "finish" without any traces selected it will
#go on to the next RD.
wash.name <- NULL
#wash.name <- "X3.Washes..DRG.Obs."
for(i in c.rd)
{
	tmp <- get(i)
	eq.mat <- tmp$mp!=tmp$t.dat
	eq.cnt <- apply(eq.mat,2,sum)
	x.names <- names(sort(eq.cnt,decreasing=T))[1:20] #check the 20 most spike-fixed cells.
	if(length(x.names) > 0)
	{x.names <- TraceSelectLarge(tmp$t.dat,,x.names,tmp$w.dat[,"wr1"],unique(tmp$w.dat[,"wr1"]),lmain="SpikeReplace")}
	if(length(x.names) > 0){tmp$bin[x.names,"drop"] <- 1}
	if(!is.null(wash.name))
	{
		x.names <- row.names(tmp$bin)[tmp$bin[,wash.name]==1]	
		if(length(x.names)>0)
		{x.names <- TraceSelectLarge(tmp$t.dat,,x.names,tmp$w.dat[,"wr1"],unique(tmp$w.dat[,"wr1"]),lmain="Wash Response")}
		if(length(x.names) > 0){tmp$bin[x.names,"drop"] <- 1}
	}
		
	tmp <- DropTestList(tmp)
	assign(i,tmp)	
}

#6b display all images for each rd
#alternate form for Lee's new cellprofiler data
for(i in c.rd)
{
	tmp <- get(i)
	img.names <- grep("^img",names(tmp),value=T)
	j.cnt <- length(img.names)
	p.name <- paste("img",i,"png",sep=".")
	png(p.name,512*j.cnt,512,bg="black")
	op <- par(mar=c(0,0,0,0))	
	plot(c(0,1),c(0,1),xaxt="n",yaxt="n",type="n",ylab="",xlab="")	
	for(j in 1:j.cnt)
	{
		x.i <- (j-1)*(1/j.cnt)
		rasterImage(tmp[[img.names[j]]],x.i,0,x.i+1/j.cnt,1)
		text(x.i+.1,.1,img.names[j],cex=4,col="white")
	}	
	dev.off()
}



#may have to scale cell size measures
#scale to median of 300
for(i in c.rd)
{
	tmp <- get(i)
	tmp$c.dat[,"area"] <- tmp$c.dat[,"area"] * (300/median(tmp$c.dat[,"area"]))
	assign(i,tmp)
}

#6c
#source('SelectGrid.R', chdir = TRUE)

for(i in c.rd)
{
#note that 'pad' shows extra pixels around the ROI
#increase this to show more pixels.
#the second input, if left blank, shows all the rois
#you can put a list of roi names in there if you want to 
#look at a subset of your rois.	
	tmp<- get(i)
	tmp <- ROIreview(tmp,,pad=2,7,7) 
	assign(i,tmp)	
}

#6d add bins for stain data
for(i in c.rd)
{
	tmp <- get(i)
	#tmp$c.dat["gfp5.bin"] <- Stain.Bin(tmp,"mean.gfp",plotit=T,k.n=5)
	tmp$c.dat["IB45.bin"] <- Stain.Bin(tmp,"mean.tritc",plotit=T,k.n=5)
	assign(i,tmp)	
}		



#7 save all these changes to disk
for(i in 1:length(c.rd)){save(list=c.rd[i],file=c.fn[i])}

#8 I use the new function RDView
#this function allows you to review the scoring 
#for all cells 1 pulse at a time.  You should pay
#close attention to setting the red diamonds to define 
#the data points of interest.  The function will 
#attempt to find the best 1 minute of interest but is often wrong.
#the first point should be the first point that is clearly showing an effect

	
for(i in c.rd)
{
#    i <- c.rd[1] # replace [1] with [5],[6] etc.
	tmp <- get(i)
	tmp <- RDView(tmp) 
	assign(i,tmp)
}


#9 DescribeAll calculates curve and response parameters for each window region
#named in b.names I usually omit the long incubations
#you MUST have reviewed with RDView and set the regions to use (red triangles).
for(i in c.rd)
{
	tmp <- get(i)
	b.names <- setdiff(unique(tmp$w.dat[,"wr1"]),c("","hoe.10nM","r715.1uM")) # it is to remove those window regions( hoe and r715)
	t.dat <- tmp$t.dat
	tmp$t.dat <- tmp$mp #again using despiked data.
	tmp <- DescribeAll(tmp,b.names,"wr1")
	tmp$t.dat <- t.dat
	assign(i,tmp)		
}


#10 save all these changes to disk again
for(i in 1:length(c.rd)){save(list=c.rd[i],file=c.fn[i])}

#now you can collect all the data in summary dataframe


#11
c.tot <- CollectMulti(c.fn)

#and run all sorts of analysis on the combined data. e.g.
#12 Show the response frequencies by rd.name
#b.names <- c("kcl.30mM","ach.1mM","atp.20uM","brady.10uM","menthol.400uM","aitc.100uM")
b.names <- sub("\\.snr","",grep("snr$",names(c.tot),value=T))
k.names <- grep("^kcl",b.names,ignore.case=T,value=T)
k.names
kp.tab <- apply(c.tot[,k.names],1,sum)
table(c.tot[,"rd.name"],kp.tab)
Kp <- kp.tab > 0 & c.tot[,"drop"]==0
MeanSemGraph(c.tot,b.names,"rd.name")
MeanSemGraph(c.tot[Kp,],c(b.names,"gfp.bin","IB4.bin"),"rd.name")

MeanSemGraph(c.tot,b.names,"ach.1mM") # find percent within a given pulse (all cells), here ach example.
pdf("kfreq.rd.pdf",width=16,height=8);MeanSemGraph(c.tot[Kp,],c.names,"rd.name",main.lab=paste("KCl Dose Response",sum(Kp),"K+ cells"));dev.off()

#12 Broad and Narrow cell classes
#note that names must match
c.tot[,"broad"] <- apply(c.tot[,c("gfp.bin","IB4.bin")],1,paste,collapse="")
c.tot[,"broad"] <- factor(c.tot[,"broad"],labels=c("GFP-, IB4-","GFP-, IB4+","GFP+, IB4-","GFP+, IB4+","unknown")[1:length(unique(c.tot[,"broad"]))])
c.tot[,"broad"] <- as.character(c.tot[,"broad"])
table(c.tot[,"broad"])
broad.col <- c("darkgrey","red","forestgreen","purple","white")
names(broad.col) <- c("GFP-, IB4-","GFP-, IB4+","GFP+, IB4-","GFP+, IB4+","unkown")

#
#note that names must match.
c.tot[,"narrow"] <- apply(c.tot[,c("gfp.bin","IB4.bin","AITC.100uM","Capsacin")],1,paste,collapse="")
levels(factor(c.tot[,"narrow"]))
jp <- is.element(c.tot[,"narrow"],c("1001","0110"))
c.tot[!jp,"narrow"] <- NA
jp <- !is.na(c.tot[,"narrow"])
c.tot[,"narrow"] <- factor(c.tot[,"narrow"],labels=c("GFP-, IB4+, AITC+, Cap-","GFP+, IB4-, AITC-, Cap+"))
c.tot[,"narrow"] <- as.character(c.tot[,"narrow"])
narrow.col <- table(c.tot[,"narrow"])
narrow.col[] <- c("darkred","green")

#12b cell size classes
c.tot[,"size.b"] <- as.factor(Size.Bin2(c.tot[,"area"]))
c.tot[,"size.b"] <- relevel(c.tot[,"size.b"],"H")
c.tot[,"size.b"] <- relevel(c.tot[,"size.b"],"L")
c.tot[,"size.b"] <- relevel(c.tot[,"size.b"],"M")
c.tot[,"size.b"] <- relevel(c.tot[,"size.b"],"S")
levels(c.tot[,"size.b"]) <- c("area < 350","350 - 500","500 - 650","> 650")

#
#
#
#DOSE RESPONSE
#
tmp <- get(unique(c.tot[,"rd.name"])[1])
unique(tmp$w.dat[,"wr1"])
k.names <- unique(tmp$w.dat[,"wr1"])[c(2,3,4,5,6,8)]
k.names
x <- c(10,20,30,40,50,0)
cbind(k.names,x)
#k.names <- grep("^kc",names(tmp$bin),value=T,ignore.case=T)
m20.names <- paste(k.names,"max.rise20",sep=".")
bm.names <- m20.names
setdiff(m20.names,names(c.tot))
phi.names <- c("phi1","phi2","phi3")
for(i in phi.names){c.tot[,i] <- NA}
for(i in c.rd)
{	
	tmp <- get(i)
	for(j in phi.names){tmp$scp[,j] <- NA}
	x.names <- c.tot[c.tot[,"rd.name"]==i & Kp,"trace.id"]
	for(j in x.names)
	{
		y <- as.vector(unlist(tmp$scp[j,bm.names]))
		dr1 <- DoseResponse(x,y,plotit=F)
		c.tot[c.tot[,"rd.name"]==i & c.tot[,"trace.id"]==j,phi.names] <- dr1
		tmp$scp[j,phi.names] <- dr1		
	}
	assign(i,tmp)
}		

#look at KCl Dose response
d.names <- PointTrace(c.tot[Kp,],x.trt="ec50",y.trt="rd.jit")
library(cluster)
figs <- matrix(c(rep(c(0,.8,.1,1),10),rep(c(.8,1,.1,1),10)),ncol=4,byrow=T)
figs[,3] <- rep(seq(1,0,by= -.1)[-1],2)
figs[,4] <- rep(seq(1,0,by= -.1)[-11],2)


y.names <- d.names #row.names(pm1$medoids) 
x.names <- c.tot[y.names,"trace.id"]

pdf(paste("KCl.ec50.h.fig.pdf",sep=""),height=40,width=20)
par(mar=c(0,0,0,0),family="sans",font=2)

#dev.new()
split.screen(figs)
for(i in 1:10)
{
	j <- c.tot[y.names[i],"rd.name"]
	tmp <- get(c.tot[y.names[i],"rd.name"])
	tmp$t.dat <- tmp$mp
	screen(i)
	par(mar=c(3,3,1,1),mgp=c(2,.25,0),tcl=NA,lwd=2)
	PeakFunc6(tmp,x.names[i],lmain=paste(x.names[i],"  ",j))	
}

for(i in 1:10)
{
	tmp <- get(c.tot[y.names[i],"rd.name"])
	j <- x.names[i]
	screen(i+10)
	par(mar=c(3,3,1,1),mgp=c(2,.25,0),tcl=NA,lwd=2)	
	y <- as.vector(unlist(tmp$scp[j,bm.names]))
	plot(x,y,xlab="KCl mM",ylab="Max",main="",type="p",pch=16,ylim=c(0,1))
	dr1 <- DoseResponse(x,y,plotit=T)
}
dev.off()
