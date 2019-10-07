#Select the t.dat, blc and, w.dat
t_to_view<-select.list(names(RD.181130.44.f.m1.p1), multiple=T)
#create both experiments in a tmp formate
dat1<-RD.181130.44.f.m1.p1
dat2<-RD.181130.44.f.m1.p2

#How long were your experiments seperated?
#I think at least 7 min is mandatory
time_buffer<-0
#look at the max time value in the t.dat on the first experiemnt
#add it to the time_buffer
time_to_change<-max(dat1$t.dat[1])+time_buffer
#now increase every time value in the second experiment by the 
#value above
changed_time<-dat2$t.dat[,1]+time_to_change

#Now we need make unique names of rhte windew regions
names_to_change<-setdiff(unique(dat2$w.dat[,"wr1"]),"")
new_names<-paste(names_to_change,"_2", sep="")
dat2$w.dat["ignore"]<-0
for(i in 1:length(names_to_change)){
    dat2$w.dat[ dat2$w.dat["wr1"]==names_to_change[i] ,2]<-new_names[i]
}

##now merg ethe datasets
for(i in 1:length(t_to_view)){
    #change the row names first
    row.names( dat2[[ t_to_view[i] ]] )<-as.character(changed_time)
    #change the first collumn value
    dat2[[ t_to_view[i] ]][1]<-changed_time
    #combine the experiments trace dataframes together
    dat1[[ t_to_view[i] ]]<-rbind(dat1[[ t_to_view[i] ]], dat2[[ t_to_view[i] ]])
}

#And reprocess the data
levs<-setdiff(unique(as.character(dat1$w.dat[,2])),"")
	snr.lim=5;hab.lim=.05;sm=2;ws=3;blc="SNIP"
pcp <- ProcConstPharm(dat1,sm,ws,blc)
scp <- ScoreConstPharm(dat1,pcp$blc,pcp$snr,pcp$der,snr.lim,hab.lim,sm)
bin <- bScore(pcp$blc,pcp$snr,snr.lim,hab.lim,levs,dat1$w.dat[,"wr1"])

dat1$scp<-scp
dat1$snr<-pcp$snr
dat1$blc<-pcp$blc
dat1$bin<-bin

alarm()
alarm()
alarm()

RD.181130.44.f.m1.p1<-dat1

#The w.dat is giving us a headach 
# as is the bin and scp