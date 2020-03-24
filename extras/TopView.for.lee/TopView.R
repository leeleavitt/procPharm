#tmp is an RD object
#x.names defines the cells to display
#wt = device window width
#ht = device window height
#scale.var is the variable used to scale each cell.  If not in the RD object defaults to a log scale transform then each row scale 0-1 min to max.
#aux.var is a list of auxillary variables to be displayed to the right of the traces.  If there are missing values, no variation or variables not in the RD object they are not shown.
#img is the img name sent to LinesEvery.TV
#t.type is the trace type data sent to LinesEvery.TV
#title is the device window title.
TopView <- function(tmp,x.names=NULL,wt=7,ht=4,scale.var="mean.sm.sd",aux.var=c("diameter","IB45.bin","gfp5.bin"),img="img1",t.type="t.dat",title="TOPVIEW")
{
	#vet the vars
	m.tot <- CollectMulti(rd.names=c(deparse(substitute(tmp))))
	scale.var <- intersect(scale.var,names(m.tot))
	aux.var <- intersect(aux.var,names(m.tot))
	if(is.null(x.names)){x.names <- row.names(tmp$bin)}

	blc.s <- as.matrix(t(tmp$blc[,x.names]))
	name2 <- make.names(tmp$w.dat[,"wr1"],unique=F)
	name2[is.element(name2,c("X","epad"))] <- "gap"	
	dimnames(blc.s)[[2]] <- name2

	if(length(scale.var)==1)
	{
		blc.s <- sweep(blc.s,1,m.tot[x.names,scale.var],'/')
	}
	else
	{
		blc.s <- log10(blc.s+1)
		med <- apply(blc.s,1,min)
		blc.s <- sweep(blc.s,1,med,'-')
		blc.s[blc.s<0] <- 0
	}
	blc.s <- blc.s-min(blc.s)
	blc.s <- blc.s/max(blc.s)	

	if(length(aux.var) > 0)
	{
		aux.mat <- NULL
		n <- ceiling(nrow(tmp$w.dat)/50)
		for(i in aux.var)
		{
			mat1 <- matrix(rep(m.tot[x.names,i],n),ncol=n)
			dimnames(mat1)[[1]] <- x.names
			dimnames(mat1)[[2]] <- rep(i,n)
			aux.mat <- cbind(aux.mat,mat1)
		}
		aux.min <- apply(aux.mat,2,min)
		aux.mat <- sweep(aux.mat,2,aux.min,'-')
		aux.max <- apply(aux.mat,2,max)
		aux.mat <- sweep(aux.mat,2,aux.max,'/')
		aux.mean <- apply(aux.mat,2,mean)
		aux.mat <- aux.mat[,!is.na(aux.mean)]
		blc.s <- cbind(blc.s,aux.mat)
	}

	name2 <- dimnames(blc.s)[[2]]	
	name2[!is.element(name2,names(m.tot))] <- "gap"
	seqi <- seq(0,(ncol(blc.s)-1))/(ncol(blc.s)-1)
	click.id <- data.frame(x=tapply(seqi,as.factor(dimnames(blc.s)[[2]]),median))
	click.id[,"y"] <- 1
	click.id <- click.id[intersect(row.names(click.id),names(m.tot)),]
	click.vals <- m.tot[x.names,row.names(click.id)]
#	click.vals[is.na(click.vals)] <- 0
	PlotHeatMat(blc.s,wt=wt,ht=ht,title=title)
	xy.click <- list(x=1,y=1)
	while(xy.click$x > 0 | xy.click$y > 0)
	{
		xy.click <- locator(n=1,type="n")
		#print(xy.click)
		if(xy.click$y > 1)
		{
			sort.trt <- row.names(click.id)[which.min(abs(xy.click$x-click.id[,"x"]))]
			sval <- click.vals[x.names,sort.trt]
			x.names <- x.names[order(sval)]
			blc.s <- blc.s[x.names,]
			PlotHeatMat(blc.s,new.dev=F,wt=wt,ht=ht)
		}
		if(xy.click$y < 1)
		{
			len1 <- length(x.names)-1
			y.i <- abs(seq(0,len1)/len1 - (1-xy.click$y))
			names(y.i) <- x.names			
			sort.i <- names(sort(y.i)[1:10])
			di <- dev.cur()
			LinesEvery.TV(tmp,sort.i,lw=3,levs.cols=grey(.95),img=tmp[[img]],t.type=t.type,m.order <- seq(1,length(sort.i)),rtag="diameter",rtag2="gfp5.bin",rtag3="IB45.bin",zf=15,cols="black")
			dev.set(di)
		}
		
	}
	
}

