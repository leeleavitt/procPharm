ce#WAZZZUPPPPP

#this is how to use ROIreview
# assign whatever you do to a temporary object
#then assign that temprary object back to RD
# then save
#dat is the RD.blah.blah.blah
tmp<-cell.creator(dat, score=T, subset.n=300)

dat<-tmp

save(dat, file="dat.rdata")

#To use cellz i have changed it up a bit
# there are now 2 cellz functions
# you can enter either the data frame you are interested
#or you can enter the entire RDfile and then select


#how do it one way
# I'm working through cellzand so you can input a list of cells, to select from only that list of cells
selected.cells<-cellzand(dat,,,1)
#-or-
selected.cells<-cellzand(dat$bin,,,1)

#The difference between 'and and 'or'
#is, 

#AND
#IF you select more than one variable, the cell MUST be positive for BOTH

#OR
#If you select more than one variable, the cell can be positive for either

#170109
#create csv for each group

dat<-RD.161221.60.m.p1

for(i in 1:length(dat$cells)){

write.csv(
	dat$cells[[i]], 
	file=paste(names(dat$cells[i]),".csv")
)

}

# How to Use use LinesStack.2.1 and LinesStack.2.2
# New LinesStack.2
# interact: LOGICAL, True allows you to search through groups
# region.group: LOGICAL, TRUE allows you to select a region for PAM clustering
# levs is an option to select which window regions to display on your plot
# subset.n: NUMERICAL allows you to select how many groups to PAM cluster into
LinesStack.2(dat,m.names,lmain="", interact=T, region.group=T,levs=NULL, plot.new=TRUE,bcex=.7, sf=.7, subset.n=5, img=NULL)




#170403
source("E:\\Data\\Lee Leavitt\\Box Sync\\Box Sync\\procpharm\\procPharm 170210.r")
Trace.Click(RD.170125.45.f.p20
Trace.Click(RD.170125.45.f.p2, c.sort.2(RD.170125.45.f.p2))
Trace.Click(RD.170125.45.f.p2, c.sort.2(RD.170125.45.f.p2))
tmp.rd<-ROIreview(RD.170125.45.f.p2,subset.n=500)
graphics.off()
source("E:\\Data\\Lee Leavitt\\Box Sync\\Box Sync\\procpharm\\170103.doodles.lee.R")
tmp.rd<-RDView(RD.170125.45.f.p2)
RD.170125.45.f.p2<-tmp.rd
savehistory("mario.r")

LinesStack.2(dat,m.names,lmain="", interact=T, region.group=T,levs=NULL, plot.new=TRUE,bcex=.7, sf=.7, subset.n=5, img=NULL)

af5.1<-.Last.value
write.csv(plxiva, file="plxiva.csv")

riiijdirect.cells<-read.csv(file="riiijdirect.csv")
(riiijdirect.cells<-paste("X.",riiijdirect.cells[,1],sep=""))


name<-last.value

Trace.Click(RD.170210.60.f.af5.1.w2, cells=c("X.1", "X.35"))

colSums(RD.170331.m.50.af5.1.w1.updated$bin[1:10])


#170612 census

#How to find cell.types
dat$cell.types$FR

#how to use census brewer
tmp.rd<-RD.
tmp.rd<-census.brewer(RD.)






> RD.170526.53.m.p1$census$LU$r3j.100nm.de


> FR<-cellzand(tmp.rd,,1,tmp.rd$cell.types$FR)
> FR.ide<-FR
> Trace.Click.dev(RD.170519.46.m.p1, FR.ide)


#########################################################
#Quick census
neurons<-cellzand(tmp.rd$bin, , 1)
length(neurons)
tmp.rd$cell.types$neurons<-neurons

tmp.rd<-census.brewer.2(tmp.rd)
#SAVE

#now remove the census that you have just created like below
tmp.rd$census<-NULL
#and then do census.brewer.2 again
tmp.rd<-census.brewer.2(tmp.rd)


#How to fix
  





boxplot(RD.171103.66.f.p2$c.dat[aitc.only, "mean.cy5.start"], 
RD.171103.66.f.p2$c.dat[cap.cells,"mean.cy5.start"], 
RD.171103.66.f.p2$c.dat[his.brady,"mean.cy5.start"], 
names=c("AITC only (NP)", "Capsaicin cells", "His responsive"), ylab="CY5 Intensity", lwd=2, col=c("darkred","purple","red"))

stripchart(list(RD.171103.66.f.p2$c.dat[aitc.only, "mean.cy5.start"], 
RD.171103.66.f.p2$c.dat[cap.cells,"mean.cy5.start"],
RD.171103.66.f.p2$c.dat[his.brady,"mean.cy5.start"]
),add=T, vertical=T, method="jitter", jitter=.2, lwd=2, cex=.3)

# How to Clean Up window regions
# The new function belows will fix the window regiuons of desire, or redo 
# The entire window region process.  It is important to look at you data
# before moving into RDView or ROIreview since this process will completely 
# erase all work provied by the two mentioned functions

RD.180104.mg.snt.c.l3l4l5.p3$img8<-readPNG("nf200.immuno.ci.rs.png")

write.csv(tmp.rd$c.dat[LU,"id"] ,file="LUid.csv")

Trace.Click.dev(tmp.rd, tmp.rd$cell.types$LU)
LU<-.Last.value
write.csv(tmp.rd$c.dat[LU[[1]],"id"] ,file="LUpropid.csv")
write.csv(tmp.rd$c.dat[LU[[1]],"area"] ,file="LUproparea.csv")
write.csv(tmp.rd$c.dat[LU[[2]],"area"] ,file="LUjagarea.csv")
write.csv(tmp.rd$c.dat[LU[[2]],"id"] ,file="LUjagid.csv")
write.csv(tmp.rd$c.dat[LU[[3]],"id"] ,file="LUideid.csv")
write.csv(tmp.rd$c.dat[LU[[3]],"area"] ,file="LUidearea.csv")
write.csv(tmp.rd$c.dat[LU[[4]],"area"] ,file="LUnearea.csv")
write.csv(tmp.rd$c.dat[LU[[4]],"id"] ,file="LUneid.csv")

Trace.Click.dev(tmp.rd, tmp.rd$cell.types$FR)
notFR<-.Last.value[[1]]
FR<-setdiff(tmp.rd$cell.types$G.0, notFR)
write.csv(tmp.rd$c.dat[FR,"area"] ,file="FRarea.csv")
write.csv(tmp.rd$c.dat[FR,"id"] ,file="FR.csv")

large<-union(LU$g.names, FR)
write.csv(tmp.rd$c.dat[large,"area"] ,file="largearea.csv")
write.csv(tmp.rd$c.dat[large,"id"] ,file="largeid.csv")

#######################################################################################33
write.csv(tmp.rd$c.dat[tmp.rd$cell.types$UL.1, "id"], file="LUpropid.csv") 
write.csv(tmp.rd$c.dat[tmp.rd$cell.types$UL.1, "area"], file="LUproparea.csv")
write.csv(tmp.rd$c.dat[tmp.rd$cell.types$UL.2, "id"], file="LUjagid.csv")
write.csv(tmp.rd$c.dat[tmp.rd$cell.types$UL.2, "area"], file="LUjagarea.csv")
write.csv(tmp.rd$c.dat[tmp.rd$cell.types$UL.3, "id"], file="LUideid.csv")
write.csv(tmp.rd$c.dat[tmp.rd$cell.types$UL.3, "area"], file="LUidearea.csv")
write.csv(tmp.rd$c.dat[tmp.rd$cell.types$UL.4, "id"], file="LUneid.csv")
write.csv(tmp.rd$c.dat[tmp.rd$cell.types$UL.4, "area"], file="LUnearea.csv")

> large<-union(FR, tmp.rd$cell.types$UL.1)
> large<-union(large, tmp.rd$cell.types$UL.2)
> large<-union(large, tmp.rd$cell.types$UL.3)
> large<-union(large, tmp.rd$cell.types$UL.4)

write.csv(tmp.rd$c.dat[peptidergic, "id"], file="peptidergic.id.csv") 
write.csv(tmp.rd$c.dat[peptidergic, "area"], file="peptidergicarea.csv")
write.csv(tmp.rd$c.dat[nonpep, "id"], file="nonpepid.csv") 
write.csv(tmp.rd$c.dat[nonpep, "area"], file="nonpeparea.csv")
write.csv(tmp.rd$c.dat[SU, "id"], file="SUid.csv") 
write.csv(tmp.rd$c.dat[SU, "area"], file="SUarea.csv")

tmp.rd<-RD.170503.41.f.p2
peptidergic<-tmp.rd$cell.types$PN
nonpep<-tmp.rd$cell.types$NP
SU<-tmp.rd$cell.types$SU
SU<-union(US, tmp.rd$cell.types$SU)
 
SU<-tmp.rd$cell.types$US.A
SU<-union(SU, tmp.rd$cell.types$US.C)
SU<-union(SU, tmp.rd$cell.types$US.0)
SU<-union(SU, tmp.rd$cell.types$thermos)
peptidergic<-tmp.rd$cell.types$P.C
peptidergic<-union(peptidergic, tmp.rd$cell.types$P.A)
peptidergic<-union(peptidergic, tmp.rd$cell.types$P.C.A)
peptidergic<-union(peptidergic, tmp.rd$cell.types$P.0)
peptidergic<-union(peptidergic, tmp.rd$cell.types$P.M)
nonpep<-tmp.rd$cell.types$N.A
nonpep<-union(nonpep, tmp.rd$cell.types$N.C)
nonpep<-union(nonpep, tmp.rd$cell.types$N.other)
 
 
 