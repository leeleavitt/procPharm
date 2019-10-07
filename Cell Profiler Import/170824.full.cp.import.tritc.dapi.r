# Do this for multiple folders.  So i am going to select.list ll the folders
#in my current directory that hass all my experiments listed.

#first select the folder that contains ALL experiments
main.dir<-"E:/Data/Abby/"
setwd(main.dir)
# select the folder for each experimetn from the list
exp.dir<-select.list(list.dirs(), multiple=T)

#Create names for the experiments you will assess today.
rd.names<-paste("RD.",sub("rotocol ","",sub("/",".",sub("./","",exp.dir))), sep="")
#rd.names<-paste("RD.",sub("/",".",sub("./","",exp.dir)), sep="")


#select your um for pixel conversion
area.conversion<-1.625

#Input the names of the files that CP created
cell.data.name<-"celldatacells_filtered.txt"
video.data.name<-"video_dataROI.txt"

#Names of the files you want loaded into the RD file
img1.1<-"bf.tritc.dapi.start.png"
img2.1<-"fura2.png"
img3.1<-"tritc.start.ci.ltl.rs.png"
img4.1<-"tritc.end.ci.ltl.rs.png"
img5.1<-"tritc.dapi.start.ci.ltl.rs.png"
img6.1<-"tritc.dapi.end.ci.ltl.rs.png"
img7.1<-"bf.end.lab.png"
img8.1<-"dapi.end.lab.png"

setwd(main.dir)

for(i in 1:length(rd.names)){
setwd(exp.dir[i])


# read in cell data.  This is a big data frame, but we only need a 
# few collumns for now.
c.dat<-read.delim(cell.data.name, header=T, sep="\t")

#These are the collumns needed for analysis
c.dat.names<-c(
	"ObjectNumber", 
	"AreaShape_Area", 
	"AreaShape_Center_X", 
	"AreaShape_Center_Y", 
	"AreaShape_FormFactor", 
	"AreaShape_Perimeter",
	"Intensity_MeanIntensity_IB4_start_ci",
	"Intensity_MeanIntensity_IB4_end_ci",
	"Intensity_MeanIntensity_DAPI_ci",
	"Intensity_MeanIntensity_BF_start",
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
	"mean.tritc.start",
	"mean.tritc.end",
	"mean.dapi",
	"mean.bf",
	"center.x",
	"center.y")
c.dat<-cbind(c.dat.1, c.dat)	


# Now read in video Data
require(data.table)
f2.img<-fread(video.data.name)
#cell names

##Mean Trace
t.340<-xtabs(Intensity_MeanIntensity_f2_340~ImageNumber+ObjectNumber, f2.img)
t.380<-xtabs(Intensity_MeanIntensity_f2_380~ImageNumber+ObjectNumber, f2.img)
t.dat<-t.340/t.380

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

t.380<-cbind(time.min[,1], t.380)
t.380<- t.380[unique(row.names(t.380)),]#incase of duplicated rownames
colnames(t.380)<- c("Time",cell.names)
row.names(t.380)<-time.min[,1]
t.380<-as.data.frame(t.380)

t.dat<-cbind(time.min[,1], t.dat)
t.dat<- t.dat[unique(row.names(t.dat)),]#incase of duplicated rownames
colnames(t.dat) <- c("Time",cell.names)
row.names(t.dat)<-time.min[,1]
t.dat<-as.data.frame(t.dat)

##LowerQuartile Trace
#tlq.340<-xtabs(Intensity_LowerQuartileIntensity_f2_340~ImageNumber+ObjectNumber, f2.img)
#tlq.380<-xtabs(Intensity_LowerQuartileIntensity_f2_380~ImageNumber+ObjectNumber, f2.img)
#tlq.dat<-tlq.340/tlq.380

#tlq.340<-cbind(time.min[,1], tlq.340)
#tlq.340<- tlq.340[unique(row.names(tlq.340)),]#incase of duplicated rownames
#colnames(tlq.340)<- c("Time",cell.names)
#row.names(tlq.340)<-time.min[,1]
#tlq.340<-as.data.frame(tlq.340)

#tlq.380<-cbind(time.min[,1], tlq.380)
#tlq.380<- tlq.380[unique(row.names(tlq.380)),]#incase of duplicated rownames
#colnames(tlq.380)<- c("Time",cell.names)
#row.names(tlq.380)<-time.min[,1]
#tlq.380<-as.data.frame(tlq.380)

#tlq.dat<-cbind(time.min[,1], tlq.dat)
#tlq.dat<- tlq.dat[unique(row.names(tlq.dat)),]#incase of duplicated rownames
#colnames(tlq.dat) <- c("Time",cell.names)
#row.names(tlq.dat)<-time.min[,1]
#tlq.dat<-as.data.frame(tlq.dat)

##upperQuartile Trace
#tuq.340<-xtabs(Intensity_UpperQuartileIntensity_f2_340~ImageNumber+ObjectNumber, f2.img)
#tuq.380<-xtabs(Intensity_UpperQuartileIntensity_f2_380~ImageNumber+ObjectNumber, f2.img)
#tuq.dat<-tuq.340/tuq.380

#tuq.340<-cbind(time.min[,1], tuq.340)
#tuq.340<- tuq.340[unique(row.names(tuq.340)),]#incase of duplicated rownames
#colnames(tuq.340)<- c("Time",cell.names)
#row.names(tuq.340)<-time.min[,1]
#tuq.340<-as.data.frame(tuq.340)

#tuq.380<-cbind(time.min[,1], tuq.380)
#tuq.380<- tuq.380[unique(row.names(tuq.380)),]#incase of duplicated rownames
#colnames(tuq.380)<- c("Time",cell.names)
#row.names(tuq.380)<-time.min[,1]
#tuq.380<-as.data.frame(tuq.380)

#tuq.dat<-cbind(time.min[,1], tuq.dat)
#tuq.dat<- tuq.dat[unique(row.names(tuq.dat)),]#incase of duplicated rownames
#colnames(tuq.dat) <- c("Time",cell.names)
#row.names(tuq.dat)<-time.min[,1]
#tuq.dat<-as.data.frame(tuq.dat)

##edge Trace
#te.340<-xtabs(Intensity_MeanIntensityEdge_f2_340~ImageNumber+ObjectNumber, f2.img)
#te.380<-xtabs(Intensity_MeanIntensityEdge_f2_380~ImageNumber+ObjectNumber, f2.img)
#te.dat<-te.340/te.380

#te.340<-cbind(time.min[,1], te.340)
#te.340<- te.340[unique(row.names(te.340)),]#incase of duplicated rownames
#colnames(te.340)<- c("Time",cell.names)
#row.names(te.340)<-time.min[,1]
#te.340<-as.data.frame(te.340)

#te.380<-cbind(time.min[,1], te.380)
#te.380<- te.380[unique(row.names(te.380)),]#incase of duplicated rownames
#colnames(te.380)<- c("Time",cell.names)
#row.names(te.380)<-time.min[,1]
#te.380<-as.data.frame(te.380)

#te.dat<-cbind(time.min[,1], te.dat)
#te.dat<- te.dat[unique(row.names(te.dat)),]#incase of duplicated rownames
#colnames(te.dat) <- c("Time",cell.names)
#row.names(te.dat)<-time.min[,1]
#te.dat<-as.data.frame(te.dat)

##MaxIntenistyEdge Trace
#tme.340<-xtabs(Intensity_MaxIntensityEdge_f2_340~ImageNumber+ObjectNumber, f2.img)
#tme.380<-xtabs(Intensity_MaxIntensityEdge_f2_380~ImageNumber+ObjectNumber, f2.img)
#tme.dat<-tme.340/tme.380

#tme.340<-cbind(time.min[,1], tme.340)
#tme.340<- tme.340[unique(row.names(tme.340)),]#incase of duplicated rownames
#colnames(tme.340)<- c("Time",cell.names)
#row.names(tme.340)<-time.min[,1]
#tme.340<-as.data.frame(tme.340)

#tme.380<-cbind(time.min[,1], tme.380)
#tme.380<- tme.380[unique(row.names(tme.380)),]#incase of duplicated rownames
#colnames(tme.380)<- c("Time",cell.names)
#row.names(tme.380)<-time.min[,1]
#tme.380<-as.data.frame(tme.380)

#tme.dat<-cbind(time.min[,1], tme.dat)
#tme.dat<- tme.dat[unique(row.names(tme.dat)),]#incase of duplicated rownames
#colnames(tme.dat) <- c("Time",cell.names)
#row.names(tme.dat)<-time.min[,1]
#tme.dat<-as.data.frame(tme.dat)


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
	snr.lim=5;hab.lim=.05;sm=2;ws=3;blc="SNIP"
pcp <- ProcConstPharm(tmp.rd,sm,ws,blc)
scp <- ScoreConstPharm(tmp.rd,pcp$blc,pcp$snr,pcp$der,snr.lim,hab.lim,sm)
bin <- bScore(pcp$blc,pcp$snr,snr.lim,hab.lim,levs,tmp.rd$w.dat[,"wr1"])
bin <- bin[,levs]
bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
bin<-pf.function(bin,levs)

# Add images
require(png)
img1<-tryCatch(readPNG(img1.1), error=function(e)NULL)
img2<-tryCatch(readPNG(img2.1), error=function(e)NULL)
img3<-tryCatch(readPNG(img3.1), error=function(e)NULL)
img4<-tryCatch(readPNG(img4.1), error=function(e)NULL)
img5<-tryCatch(readPNG(img5.1), error=function(e)NULL)
img6<-tryCatch(readPNG(img6.1), error=function(e)NULL)
img7<-tryCatch(readPNG(img7.1), error=function(e)NULL)
img8<-tryCatch(readPNG(img8.1), error=function(e)NULL)
#readPNG("bf.gfp.tritc.dapi.png")

#img.t<-readPNG("tritc.png")
#img.b<-readPNG("bf.png")
#img.bl<-readPNG("bf.lab.png")
#img.f<-readPNG("fura2.png")

tmp.rd <- list(t.dat=t.dat,t.340=t.340,t.380=t.380, 
#tlq.dat=tlq.dat, tlq.340=tlq.340,tlq.380=tlq.380,
#tuq.dat=tuq.dat,tuq.340=tuq.340,tuq.380=tuq.380,
#te.dat=te.dat,te.340=te.340,te.380=te.380,
#tme.dat=tme.dat,tme.340=tme.340,tme.380=tme.380,
w.dat=w.dat,c.dat=c.dat, bin=bin, scp=scp, snr=pcp$snr, blc=pcp$blc, der=pcp$der, 
img1=img1,img2=img2,img3=img3,img4=img4,img5=img5,img6=img6,img7=img7, img8=img8)
#img.gtd=img.gtd, img.g=img.g, img.t=img.t, img.b=img.b,img.bl=img.bl, img.f=img.f)

tmp.rd$c.dat[,"area"]<-tmp.rd$c.dat$area*area.conversion

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

for(k in 1:length(colnames(tmp.rd$mp))){

	tmp.dat[,k]<-tmp.rd$mp[,k]-min(tmp.rd$mp[,k])
	tmp.dat[,k]<-tmp.dat[,k]/max(tmp.dat[,k])

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
rm(c.dat)
rm(cell.names)
setwd(main.dir)
}



















