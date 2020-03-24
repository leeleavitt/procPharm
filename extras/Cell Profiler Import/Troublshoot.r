

#I need to ensure the ROI from my video is the same as the ROI from my ce data file
roi.dat<-read.delim("ROI Data.txt", header=T, sep="\t", fileEncoding="UCS-2LE")
vid.dat<-read.delim("ROI video Data.txt", header=T, sep="\t", fileEncoding="UCS-2LE")

#create coordinate dataframe for picture file
centerxy.pic<-roi.dat[,c("RoiID","CentreXpx", "CentreYpx")]
centerxy.pic<-centerxy.pic[order(centerxy.pic[,"RoiID"]),]

#create coordinate dataframe for video file
centerxy.vid<-vid.dat[,c("RoiID","CentreXpx", "CentreYpx")]
#put it into coridenates that should match pic
centerxy.vid[c("CentreXpx","CentreYpx")]<-centerxy.vid[c("CentreXpx","CentreYpx")]*2
centerxy.vid<-centerxy.vid[order(centerxy.vid[,"RoiID"]),]

#Round out both
#centerxy.vid<-round(centerxy.vid, digits=0)
#centerxy.pic<-round(centerxy.pic, digits=0)

# Kevin did this assesment
quantile(centerxy.pic[,2] - centerxy.vid[,2])

quantile(centerxy.pic[,3] - centerxy.vid[,3])
	
#see if they are identical sorted by RoiID  They are not
identical(centerxy.vid,centerxy.vid) #this is true
identical(centerxy.vid,centerxy.pic) #this is false

#Alright, sort both dataframes by CentreXpx, and see if the rois are identical
centerxy.pic<-centerxy.pic[order(centerxy.pic[,"CentreXpx"]),]
centerxy.vid<-centerxy.vid[order(centerxy.vid[,"CentreXpx"]),]

identical(centerxy.pic[,"RoiID"],centerxy.vid[,"RoiID"])
#shows false See below
centerxy.pic[,"RoiID"]==centerxy.vid[,"RoiID"]


#Alright, sort both dataframes by CentreYpx, and see if the rois are identical
centerxy.pic<-centerxy.pic[order(centerxy.pic[,c("CentreYpx","CentreXpx")]),]
centerxy.vid<-centerxy.vid[order(centerxy.vid[,"CentreYpx"]),]

identical(centerxy.pic[,"RoiID"],centerxy.vid[,"RoiID"])
#shows false
centerxy.pic[,"RoiID"]==centerxy.vid[,"RoiID"]

# Over lunch though of this.  Add CentreXpx with CentreYpx sory by that value and see if RoiID matches up
centerxy.pic["XplusY"]<-centerxy.pic[,"CentreYpx"]+centerxy.pic[,"CentreXpx"]
centerxy.vid["XplusY"]<-centerxy.vid[,"CentreYpx"]+centerxy.vid[,"CentreXpx"]

#Still dont match up
centerxy.pic[order(centerxy.pic[,"XplusY"]),]
centerxy.vid[order(centerxy.vid[,"XplusY"]),]

centerxy.pic[order(centerxy.pic[,"XplusY"]),"XplusY"]
centerxy.vid[order(centerxy.vid[,"XplusY"]),"XplusY"]

centerxy.pic[order(centerxy.pic[,"XplusY"]),"XplusY"]=centerxy.vid[order(centerxy.vid[,"XplusY"]),"XplusY"]
centerxy.pic[order(centerxy.pic[,"XplusY"]),"XplusY"]-centerxy.vid[order(centerxy.vid[,"XplusY"]),"XplusY"]

cbind(centerxy.pic[order(centerxy.pic[,"XplusY"]),"RoiID"],centerxy.vid[order(centerxy.vid[,"XplusY"]),"RoiID"])

