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