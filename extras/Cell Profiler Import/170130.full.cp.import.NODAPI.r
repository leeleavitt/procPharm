# Do this for multiple folders.  So i am going to select.list ll the folders
#in my current directory that hass all my experiments listed.

#first select the folder that contains ALL experiments
main.dir<-"E:/Data/Lee Leavitt"
setwd(main.dir)
exp.dir<-select.list(list.dirs(), multiple=T)

#name the experiments you will assess today.
rd.names<-c(
"RD.170128.f.47.p1",
"RD.170128.f.47.p2",
"RD.170128.f.47.p3",
"RD.170128.f.47.p4",
"RD.170128.f.47.p5",
"RD.170128.f.47.p6")





for(i in 1:length(rd.names)){
setwd(exp.dir[i])


# read in cell data.  This is a big data frame, but we only need a 
# few collumns for now.
c.dat<-read.delim("celldatacells_filtered.txt", header=T, sep="\t")

#These are the collumns needed for analysis
green<-grep("CGRP", names(c.dat))
if(length(green)>1){
	c.dat.names<-c(
	"ObjectNumber", 
	"AreaShape_Area", 
	"AreaShape_Center_X", 
	"AreaShape_Center_Y", 
	"AreaShape_FormFactor", 
	"AreaShape_Perimeter",
	"Intensity_MeanIntensity_CGRP_1",
	"Intensity_MeanIntensity_IB4_1",
	"Location_Center_X",
	"Location_Center_Y")

	c.dat.1<-c.dat[,c.dat.names] 

	# i like to rename the collumns for my procpharm
names(c.dat.1)<-c(
	"id",
	"area",
	"center.x.simplified", 
	"center.y.simplified",
	"circularity",
	"perimeter",
	"mean.gfp",
	"mean.tritc",
	"center.x",
	"center.y")
c.dat<-cbind(c.dat.1, c.dat)	
}else{
c.dat.names<-c(
	"ObjectNumber", 
	"AreaShape_Area", 
	"AreaShape_Center_X", 
	"AreaShape_Center_Y", 
	"AreaShape_FormFactor", 
	"AreaShape_Perimeter",
	"Intensity_MeanIntensity_IB4_1",
	"Intensity_MeanIntensity_DNA_1",
	"Location_Center_X",
	"Location_Center_Y")

c.dat.1<-c.dat[,c.dat.names]

	# i like to rename the collumns for my procpharm
names(c.dat.1)<-c(
	"id",
	"area",
	"center.x.simplified", 
	"center.y.simplified",
	"circularity",
	"perimeter",
	"mean.tritc",
	"mean.dapi",
	"center.y")
c.dat<-cbind(c.dat.1, c.dat)	

}



# Now read in video Data
require(data.table)
f2.img<-fread("video_dataROI.txt")
#cell names

t.340<-xtabs(Intensity_MeanIntensity_f2_340~ImageNumber+ObjectNumber, f2.img)
t.380<-xtabs(Intensity_MeanIntensity_f2_380~ImageNumber+ObjectNumber, f2.img)
t.e.340<-xtabs(Intensity_MeanIntensityEdge_f2_340~ImageNumber+ObjectNumber, f2.img)
t.e.380<-xtabs(Intensity_MeanIntensityEdge_f2_380~ImageNumber+ObjectNumber, f2.img)
t.dat<-t.340/t.380
t.e.dat<-t.e.340/t.e.380

# add time info from file exported on nis veiwer
time.info<-read.delim("time.info.txt", sep="\t", fileEncoding="UCS-2LE")
if(length(which(is.na(time.info[3]), arr.ind=T)[,1])>=1){
	time.info<-time.info[-which(is.na(time.info[3]), arr.ind=T)[,1],] 
}else{time.info<-time.info}#remove any nas
time.min<-round(time.info["Time..s."]/60, digits=3)

# Create row.names
cell.names<-paste("X.",c.dat$id,sep="")
c.dat[,"id"]<-cell.names
row.names(c.dat)<-cell.names

t.340<-cbind(time.min[,1], t.340)
t.340<- t.340[unique(row.names(t.340)),]#incase of duplicated rownames
colnames(t.340)<- c("Time",cell.names)
row.names(t.340)<-time.min[,1]
t.340<-as.data.frame(t.340)

t.e.340<-cbind(time.min[,1], t.e.340)
t.e.340<- t.e.340[unique(row.names(t.e.340)),]#incase of duplicated rownames
colnames(t.e.340)<- c("Time",cell.names)
row.names(t.e.340)<-time.min[,1]
t.e.340<-as.data.frame(t.e.340)

t.380<-cbind(time.min[,1], t.380)
t.380<- t.380[unique(row.names(t.380)),]#incase of duplicated rownames
colnames(t.380)<- c("Time",cell.names)
row.names(t.380)<-time.min[,1]
t.380<-as.data.frame(t.380)

t.e.380<-cbind(time.min[,1], t.e.380)
t.e.380<- t.e.380[unique(row.names(t.e.380)),]#incase of duplicated rownames
colnames(t.e.380)<- c("Time",cell.names)
row.names(t.e.380)<-time.min[,1]
t.e.380<-as.data.frame(t.e.380)


t.dat<-cbind(time.min[,1], t.dat)
t.dat<- t.dat[unique(row.names(t.dat)),]#incase of duplicated rownames
colnames(t.dat) <- c("Time",cell.names)
row.names(t.dat)<-time.min[,1]
t.dat<-as.data.frame(t.dat)

t.e.dat<-cbind(time.min[,1], t.e.dat)
t.e.dat<- t.e.dat[unique(row.names(t.e.dat)),]#incase of duplicated rownames
colnames(t.e.dat) <- c("Time",cell.names)
row.names(t.e.dat)<-time.min[,1]
t.e.dat<-as.data.frame(t.e.dat)



## Proof my dataframes are already sorting in the correct way.
#first create a dataframe from the video file x and y location data
#t.loc.x<-matrix(f2.img[,"Location_Center_X"],byrow=FALSE, nrow=max(cell.names))[,1]
#t.loc.y<-matrix(f2.img[,"Location_Center_Y"],byrow=FALSE, nrow=max(cell.names))[,1]
#t.loc<-cbind(t.loc.x, t.loc.y)
#t.loc<-round(t.loc, digits=3) #Then round dataframe to three digits

## Create datatrame from c.dat(c created above)
#c.xy<-c.dat[,c("center.x", "center.y")]
#c.xy<-round(c.xy, digits=3) # also round to three digits
#setdiff(t.loc[,2], c.xy[,2]) #There should be no difference between the two data frames

#setequal(t.loc[,1], c.xy[,1])#There should be no difference between the two data frame
#setequal(t.loc[,2], c.xy[,2])#There should be no difference between the two data frame

#example
#a<-seq(from=1,to=10,by=10)
#b<-seq(from=1,to=10,by=10)

############ wr1 import
wrdef<-"wr1.csv"
if(!is.null(wrdef))
	{
		wr <- ReadResponseWindowFile(wrdef)
		Wr<-length(wr[,1])#complete and revise this section
		if(length(colnames(wr))<2){w.dat<-WrMultiplex(t.dat,wr,n=Wr)}
		else{w.dat <- MakeWr(t.dat,wr)}
		}

# Initial and simple Data processing
tmp.rd <- list(t.dat=t.dat,w.dat=w.dat,c.dat=c.dat)
levs<-setdiff(unique(as.character(w.dat[,2])),"")
snr.lim=4;hab.lim=.05;sm=3;ws=30;blc="SNIP"
pcp <- ProcConstPharm(tmp.rd,sm,ws,blc)
scp <- ScoreConstPharm(tmp.rd,pcp$blc,pcp$snr,pcp$der,snr.lim,hab.lim,sm)
bin <- bScore(pcp$blc,pcp$snr,snr.lim,hab.lim,levs,tmp.rd$w.dat[,"wr1"])
bin <- bin[,levs]
bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
bin<-pf.function(bin,levs)

# Add images
require(png)
img1<-png::readPNG("bf.png")
img2<-png::readPNG("gfp.png")
img3<-png::readPNG("tritc.png")
img4<-png::readPNG("bf.lab.png")
img5<-png::readPNG("roi.img.png")
img6<-png::readPNG("fura2.png")
img7<-NULL
img8<-NULL


#img.t<-png::readPNG("tritc.png")
#img.b<-png::readPNG("bf.png")
#img.bl<-png::readPNG("bf.lab.png")
#img.f<-png::readPNG("fura2.png")

tmp.rd <- list(t.dat=t.dat,t.340=t.340,t.380=t.380,te.340=t.e.340,te.380=t.e.380, te.dat=t.e.dat,
w.dat=w.dat,c.dat=c.dat, bin=bin, scp=scp, snr=pcp$snr, blc=pcp$blc, der=pcp$der, 
img1=img1,img2=img2,img3=img3,img4=img4,img5=img5,img6=img6,img7=img7)
#img.gtd=img.gtd, img.g=img.g, img.t=img.t, img.b=img.b,img.bl=img.bl, img.f=img.f)


## Check for duplicated rows	
if(length(which(duplicated(row.names(t.dat))))>=1){
dup<-which(duplicated(row.names(t.dat)))
paste(dup)
t.dat<-t.dat[-dup,]
w.dat<-w.dat[-dup,]
}

##DESPIKE
wts <- tmp.rd$t.dat
for(j in 1:5) #run the despike 5 times.
{
	wt.mn3 <- Mean3(wts)
	wts <- SpikeTrim2(wts,1,-1)
	print(sum(is.na(wts))) #this prints out the number of points removed should be close to 0 after 5 loops.
	wts[is.na(wts)] <- wt.mn3[is.na(wts)]
}
tmp.rd$mp <- wts

#170127
# Take the despiked data, subtract the minimum value from the trace, then divide by the maximun value
# to create traces that are all on the same 0 to 1 scale
tmp.dat<-tmp.rd$mp

for(j in 1:length(colnames(tmp.rd$mp))){

	tmp.dat[,j]<-tmp.rd$mp[,j]-min(tmp.rd$mp[,j])
	tmp.dat[,j]<-tmp.dat[,j]/max(tmp.dat[,j])

}
tmp.dat[,1]<-tmp.rd$t.dat[,1]

tmp.rd$mp.1<-tmp.dat


rd.name<-rd.names[i]
f.name <- paste(rd.name,".Rdata",sep="")
assign(rd.name,tmp.rd)
save(list=rd.name,file=f.name)

rm(f2.img)
rm(rd.name)
rm(tmp.rd)
setwd(main.dir)
}



















