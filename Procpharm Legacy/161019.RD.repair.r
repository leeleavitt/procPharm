### The RD files already have improper t data, there for t.dat 
# needs to be repaired and w.dat, scp, bin,blc,snr all need to be redone
main.dir<-"E:/DRG Profiling/Short Protocols"
setwd(main.dir)
exp.dir<-select.list(list.dirs(), multiple=T)

for(i in 1:length(exp.dir)){


setwd(exp.dir[i])
load(list.files(pattern="RD"))
dat<-get(ls(pattern="RD"))

#dapi.info<-dat$c.dat$mean.tritc
#tritc.info<-dat$c.dat$mean.dapi
#dat$c.dat["mean.tritc"]<-tritc.info
#dat$c.dat["mean.dapi"]<-dapi.info

time.info<-read.delim("time.info.txt", sep="\t", fileEncoding="UCS-2LE")
if(length(which(is.na(time.info[3]), arr.ind=T)[,1])>=1){
	time.info<-time.info[-which(is.na(time.info[3]), arr.ind=T)[,1],] 
}else{time.info<-time.info}#remove any nas
time.min<-round(time.info["Time..s."]/60, digits=3)
dat$t.dat[1]<-time.min
dat$t.340["Time"]<-time.min
dat$t.380["Time"]<-time.min




t.dat<-dat$t.dat
c.dat<-dat$c.dat
t.340<-dat$t.340
t.380<-dat$t.380


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



rd.name<-ls(pattern="RD")
f.name <- paste(rd.name,".Rdata",sep="")
assign(rd.name,tmp.rd)
save(list=rd.name,file=f.name)
rm(list=ls(pattern="RD"))
setwd(main.dir)
}

