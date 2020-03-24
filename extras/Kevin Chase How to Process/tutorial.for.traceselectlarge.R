#source the file.  You'll have to change the path.
source('~/Box Sync/Rfunctions/procPharm.0.2.9.R', chdir = TRUE)
#load an RD file (again with the path)
load('~/Box Sync/KCreview/Time course experiments/1 day post surgery snl/170224/RD.170224.d1.snl3.ipsi.l5l6.w2.m.Rdata')
#work with a temporary variable.
tmp <- get("RD.170224.d1.snl3.ipsi.l5l6.w2.m")
#criteria for cells
gi <- tmp$bin[,"kcl"]==1 & tmp$bin[,"drop"]==0 & tmp$bin[,"ach"]==1
#get the trace names
x.names <- row.names(tmp$bin)[gi]
#select only the traces that match
y.names <- TraceSelectLarge(tmp$mp,,x.names,tmp$w.dat[,"wr1"],unique(tmp$w.dat[,"wr1"]))
#now y.names has a list of the selected traces.
#you can store this list as a separate list.
#if you want to also keep track of the rd.file
ach.group1 <- NULL
ach.group1[["RD.170224.d1.snl3.ipsi.l5l6.w2.m"]] <- y.names

#now you can genralize this to run through a set of rd files.
#set the working directory
setwd("~/Box Sync/KCreview/Time course experiments/1 day post surgery snl/170224")
#get a list of all the RD files in that directory
c.fn <- list.files(pattern="RD.*\\.Rdata$",recursive=T)
#load all the rd files
for(i in c.fn){load(i)}
#2 convert file names to rd names
c.rd <- sub("^.*RD","RD",sub("\\.Rdata","",c.fn));setdiff(c.rd,ls())

#set up the list variable to store the cell ids do this only once
#if you do it again it will erase the previous variable.
ach.group1 <- NULL
#now c.rd has a list of loaded RD files
i <- c.rd[1] #change the number to access the separate rd files
tmp <- get(i)
#criteria for cells
gi <- tmp$bin[,"kcl"]==1 & tmp$bin[,"drop"]==0 & tmp$bin[,"ach"]==1
#get the trace names
x.names <- row.names(tmp$bin)[gi]
#select only the traces that match
y.names <- TraceSelectLarge(tmp$mp,,x.names,tmp$w.dat[,"wr1"],unique(tmp$w.dat[,"wr1"]))
#now y.names has a list of the selected traces.
#you can store this list as a separate list.
#if you want to also keep track of the rd.file

ach.group1[[i]] <- y.names

#after running through the rd files save the list variable.
save(ach.group1,file="ach.group1.Rdata")

