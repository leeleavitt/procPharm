
#1 get a list of all the RD files in the current directory
c.fn <- list.files(pattern="^RD.*\\.Rdata$")
#load all the rd files
for(i in c.fn){load(i)}
#2 convert file names to rd names
c.rd <- sub("^.*RD","RD",sub("\\.Rdata","",c.fn));setdiff(c.rd,ls())

#3 despike each rd file and save the trace data as $mp in the RD 
for(k in exp.names)
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
	



#ROIreview <- function(tmp,x.names=NULL,pad=2,wh=7,hh=7,subset.n=NA)
#three tests Drop (confirm), Red (confirm) and Green (confirm)
#return and RD object with the changes made to c.dat and bin
#tmp is an RD object with images, "tritc.mean" and "gfp.mean" in c.dat
#x.names is a list of specific cells to review
#pad is the expansion factor about the center of the cell.
#subset.n is number of cells to review at once instead of all at once.

for(i in c.rd)
{
	tmp <- get(i)
	tmp <- ROIreview(tmp,subset.n=500)
	assign(i,tmp)
}


#I use the new function RDView
#this function allows you to review the scoring 
#for all cells 1 pulse at a time.  You should pay
#close attention to setting the red diamonds to define 
#the data points of interest.  The function will 
#attempt to find the best 1 minute of interest but is often wrong.
#the first point should be the first point that is clearly showing an effect

#RDView <- function(tmp,wr.i="wr1",rd.name=NULL,c.i=NULL,rscale=F,wh=14,hh=8)
#tmp is an RD object
#wr.i is the window region definition.
#rd.name is the name of the RD object (used for png.out)
#c.i is an alternate list of cells to review (instead of all)
#rscale is a boolean for rescaling the data (set this to T for wild baseline shifts)
#wh is the window height
#hh is the window width (why the hell did I name it hh?)

	
for(i in c.rd)
{
	tmp <- get(i)
	tmp <- RDView(tmp) 
	assign(i,tmp)
}

# DescribeAll calculates curve and response parameters for each window region
#named in b.names I usually omit the long incubations
#and replace the t.dat with the spike trimmed mp
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

