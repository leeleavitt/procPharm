#matrix image
PlotHeatMat <- function(mat,wt=14,ht=10,new.dev=T,title="TOPVIEW")
{
#	mat <- log10(mat)
	mat <- mat-min(mat)
	mat <- mat/max(mat)
#	gy <- !grepl("blank",dimnames(mat)[[1]])	
#	pdf("tmp.pdf",width=40,height=8,family="mono",colormodel="grey")
	if(new.dev)
	{
		dev.new(width=wt,height=ht,family="mono",canvas="black",title=title)
		par(fg="darkgrey",col.axis="white",col.lab="grey",col.main="grey",mar=c(1,3,6,1))
		plot(c(0,1),c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",type="n")
	}
	
	rasterImage(mat,0,0,1,1,interpolate=F)
	gx <- !grepl("gap",dimnames(mat)[[2]])
	ux <- unique(dimnames(mat[,gx])[[2]])
	gi <- match(ux,dimnames(mat)[[2]])/ncol(mat)
	xi <- seq(0,(ncol(mat)-1))/(ncol(mat)-1)
	xi.mat <- data.frame(min=tapply(xi,dimnames(mat)[[2]],min))
	xi.mat[,"max"] <- tapply(xi,dimnames(mat)[[2]],max)
	xi.mat[,"med"] <- (xi.mat[,"max"]+xi.mat[,"min"])/2
	xi.lab <- row.names(xi.mat)
	xi.lab[grep("gap",xi.lab)] <- ""
	#xi.mat <- xi.mat[!grepl("gap",row.names(xi.mat)),]
	axis(side=3,at=xi.mat[xi.lab != "","med"],labels=xi.lab[xi.lab != ""],las=2,col.axis="darkgrey")
	axis(side=3,at=xi.mat[xi.lab=="","min"],labels=NA,tck=1,lwd=.2)
	axis(side=3,at=xi.mat[xi.lab=="","max"],labels=NA,tck=1,lwd=.2)
	yi <- seq(0,(nrow(mat)-1))/(nrow(mat)-1)
	y.names <- dimnames(mat)[[1]]
	y.names <- sub("^w[0987654321]*\\.","",y.names)
	yi.mat <- data.frame(min=tapply(yi,y.names,min))
	yi.mat[,"max"] <- tapply(yi,y.names,max)
	yi.mat[,"med"] <- (yi.mat[,"min"]+yi.mat[,"max"])/2
	yi.mat <- yi.mat[!grepl("blank",row.names(yi.mat)),]
	axis(side=2,at=1-yi.mat[,"med"],labels=row.names(yi.mat),las=3,cex.axis=.5,col.axis="darkgrey")
	#dev.off()
}
