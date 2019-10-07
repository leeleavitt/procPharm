# Do this for multiple folders.  So i am going to select.list ll the folders
#in my current directory that hass all my experiments listed.
main.dir<-"E:/DRG Profiling/Short Protocols"
setwd(main.dir)
exp.dir<-select.list(list.dirs(), multiple=T)

rd.names<-c(
"RD.161011.50.f.p1",
"RD.161011.50.f.p2"
)



for(i in 1:length(rd.names)){
setwd(exp.dir[i])


# read in cell data.  This is a big data frame, but we only need a 
# few collumns for now.
c.dat<-read.delim("celldatacells_filtered.txt", header=T, sep="\t")

#These are the collumns needed for analysis
c.dat.names<-c(
"ObjectNumber", 
"AreaShape_Area", 
"AreaShape_Center_X", 
"AreaShape_Center_Y", 
"AreaShape_FormFactor", 
"AreaShape_Perimeter",
"Intensity_MeanIntensity_CGRP",
"Intensity_MeanIntensity_IB4_1",
"Intensity_MeanIntensity_DNA_2",
"Location_Center_X",
"Location_Center_Y")

c.dat<-c.dat[,c.dat.names]

# i like to rename the collumns for my procpharm
names(c.dat)<-c(
"id",
"area",
"center.x.simplified", 
"center.y.simplified",
"circularity",
"perimeter",
"mean.gfp",
"mean.tritc",
"mean.dapi",
"center.x",
"center.y")



# Now read in video Data
# this takes a long time and needs to be optimzed outside of R
require(data.table)
f2.img<-fread("video_dataROI.txt")
f2.img<-read.delim("video_dataROI.txt", se)
f2.img2<-f2.img[-drop.rows]
#cell names
cell.names<-unique(f2.img[,"ObjectNumber"])
#table of frames per cell

row.names(f2.img)<-NULL

cell.names<-unique(f2.img[,"ObjectNumber"])

x.tab<-table(f2.img[,"ObjectNumber"])

x.row<-max(x.tab)

t.340<-xtabs(Intensity_MeanIntensity_f2_340~ImageNumber+ObjectNumber, f2.img)
t.380<-xtabs(Intensity_MeanIntensity_f2_380~ImageNumber+ObjectNumber, f2.img)

t.dat<-t.340/t.380

time.info<-read.delim("time.info.txt", sep="\t", fileEncoding="UCS-2LE")
if(length(which(is.na(time.info[3]), arr.ind=T)[,1])>=1){
	time.info<-time.info[-which(is.na(time.info[3]), arr.ind=T)[,1],] 
}else{time.info<-time.info}#remove any nas
time.min<-round(time.info["Time..s."]/60, digits=3)

# Create row.names
cell.names<-paste("X.",c.dat$id,sep="")
c.dat[,"id"]<-cell.names
row.names(c.dat)<-cell.names
t.dat <- as.data.frame(t.dat)

t.dat["Time"]<-time.min
t.dat<- t.dat[unique(row.names(t.dat)),]#incase of duplicated rownames

names(t.dat) <- c("Time",cell.names)
t.340["Time"]<-time.min
t.340<-as.data.frame(t.340)
row.names(t.340)<-time.min

t.380["Time"]<-time.min
t.380<-as.data.frame(t.380)
row.names(t.380)<-time.min



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

# Initial Data processing
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
img.gtd<-readPNG("gfp.tritc.dapi.png")
img.g<-readPNG("gfp.png")
img.t<-readPNG("tritc.png")
img.b<-readPNG("bf.png")
img.bl<-readPNG("bf.lab.png")
img.r<-readPNG("roi.img.png")
img.f<-readPNG("fura2.png")

tmp.rd <- list(t.dat=t.dat,t.340=t.340, t.380=t.380,w.dat=w.dat,c.dat=c.dat, bin=bin, scp=scp, snr=pcp$snr, blc=pcp$blc, der=pcp$der, 
img.gtd=img.gtd, img.g=img.g, img.t=img.t, img.b=img.b,img.bl=img.bl,img.r=img.r, img.f=img.f)

## Check for duplicated rows	
if(length(which(duplicated(row.names(t.dat))))>=1){
dup<-which(duplicated(row.names(t.dat)))
paste(dup)
t.dat<-t.dat[-dup,]
w.dat<-w.dat[-dup,]
}

rd.name<-rd.names[i]
f.name <- paste(rd.name,".Rdata",sep="")
assign(rd.name,tmp.rd)
save(list=rd.name,file=f.name)

setwd(main.dir)
rm(f2.img)
}



















