
#read in this data, for next export, clean up the name please
#f2.v<-fura2 video
f2.frames<-read.delim("video_dataImage.txt", header=T, sep="/t")
#these are all the names from the data set.  as you can see we
#have some added types of data.
all.names<-names(f2.v)

#will need to convert image number to a time value
# based on what the user specifies
##Or
#based on a 1.5s/frame conversion
time.name<-"ImageNumber"

#will need to align objects by XY location
# need to think about this
id.name<-"ObjectNumber"

# Unique measurements.  If memory allows, we should have 6 initial
# trace dataframes within our list of dataframes, RD

#fura2 meanintensity from each object
#f2.mean.340
f2.mean.340<-"Intensity_MeanIntensity_f2_340"
#f2.mean.380
m2.mean.380<-"Intensity_MeanIntensity_f2_380"

#fura2 meanintensityedge from the single pixel on the edge of the ROI
#f2.meanedge.380
f2.meanedge.340<-"Intensity_MeanIntensityEdge_f2_380"
#f2.meanedge.380
f2.meanedge.380<-"Intensity_MeanIntensityEdge_f2_380"

# to maintain a constant sorting order i will order the data frame first
# by y cooridinates, then by x cooridinates
# first thing

cell.names<-unique(f2_v[,id.name])
cell.names
x.tab<-table(f2_v[,id.name])
x.tab[1:5,1:5]
head(x.tab)
max(x.tab)
min(x.tab)
(x.row<-max(x.tab))
t.dat<-matrix(f2_v[,f2_340],byrow=FALSE, nrow=x.row)
t.dat[1:5,1:5]
savehistory()
