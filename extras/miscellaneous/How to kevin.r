#tcd is the new Trace.Click.dev
#steps
#1.Trace Click
#2. ROI review
#3. RDView
#4. Cell Type



#How to re-name
#example > (escape)
RD.170331.m.50.af5.1.w1.updated> RD.170331.50.m.p1<-RD.170331.m.50.af5.1.w1.updated
save(RD.170331.50.m.p1,file="RD.170331.50.m.p1.Rdata")
#check to see if the celltype is in the new one($cell.type) before deleting the old one

#How to categorize and save
Response<-.Last.value$g.names1
save(Response, file="Response.Rdata")

de<-.Last.value$g.names2                                    #put all the de into g.names2
save(de, file="de.Rdata")                                   #save de

#How to list
list(Response, de)                                          #creating a list of Response and de, unamed
list(Response =Response, de = de)                           #creating a list of Response and de, named

#How to add file into RD. (drag them into the main RD.)
load("C:\\Users\\Administrator\\Box Sync\\Calcium Imaging\\Voltage Gated Channel experiments\\Af5.1\\170110.57.f.w1.direct response to Af5.1,Ach, atp, AITC, capsaicin\\Response.Rdata")
load("C:\\Users\\Administrator\\Box Sync\\Calcium Imaging\\Voltage Gated Channel experiments\\Af5.1\\170110.57.f.w1.direct response to Af5.1,Ach, atp, AITC, capsaicin\\de.Rdata")
de
Response

#How to add mp.1

Trace.Click.dev(RD.170110.57.f.w1, Response)                #to see the response traces
tmp.rd<-mp.brewer(RD.170110.57.f.w1)                        #assign mp.brewer to tmp.rd
Trace.Click.dev(tmp.rd.170110.57.f.w1, Response)
Trace.Click.dev(tmp.rd, Response)
RD.170110.57.f.w1<-tmp.rd                                   #assign tmp.rd back to RD.
save(RD.170110.57.f.w1, file="RD.170110.57.f.w1.Rdata")     #save RD. with mp.1
rm(tmp.rd)                                                  #remove tmp.rd
tmp.rd<-RD.170110.57.f.w1
class(tmp.rd)                                               #what is tmp.rd
RD.170110.57.f.w1<-mp.brewer(RD.170110.57.f.w1)             #a faster way to assign mp.brewer to RD.

#How to names, class, dim, head, tail
names(RD.170110.57.f.w1)                                    #what is in the experiment
Trace.Click.dev(tmp.rd, Response)
names(RD.170110.57.f.w1$bin)                                #what is in the list bin
class(RD.170110.57.f.w1$bin)                                #what is bin, bin is a list
class(RD.170110.57.f.w1$scp)                                #what is scp, scp is a data frame
names(RD.170110.57.f.w1$scp)                                #what is in scp
dim(RD.170110.57.f.w1$bin)                                  #what is the dimension of bin
dim(RD.170110.57.f.w1$scp)                                  #what is the dimension of scp
dim(RD.170110.57.f.w1$t.dat)                                #what is the dimension of t.dat
dim(RD.170110.57.f.w1$blc)                                  #what is the dimension of blc
head(RD.170110.57.f.w1$bin)                                 #show the beginning of the list bin
tail(RD.170110.57.f.w1$bin)                                 #show the end of the list bin

#How to delete
RD.170110.57.f.w1[[14]]<-NULL                               #delete 14 from the experiment, which was mp.1

#How to add mp.1 part2
mp.brewer(RD.170110.57.f.w1)                                #a faster way to add mp.1 into the experiment, not recommended
names(RD.170110.57.f.w1)
RD.170110.57.f.w1<-mp.brewer(RD.170110.57.f.w1)             #the traditional way of adding mp.1 to the experiment, recommended

#Example of a list that's everywhere (not good)
ScoreConstPharm(RD.170110.57.f.w1)
require(MALDIquant)
ScoreConstPharm(RD.170110.57.f.w1)
class(.Last.value)
bob<-ScoreConstPharm(RD.170110.57.f.w1)
head(bob)
RD.170110.57.f.w1$bob<-ScoreConstPharm(RD.170110.57.f.w1)
identical(RD.170110.57.f.w1$scp, RD.170110.57.f.w1$bob)     #are bob and scp identical?

#How to cellzand
neurons<-cellzand(
neurons<-cellzand(RD.180302.59.m.p1)   #using cellzand(RD.) to find neurons, select bin, then capsasin and K.40mM
length(neurons)                        #there are 1252 neurons
neurons<-cellzand(RD.180302.59.m.p1$bin,"K.40mM",parameter=1)        #using cellzand to find neurons
length(neurons)                                                      #the amount of neurons
bp.selector(RD.180302.59.m.p1, neurons)                              #use bp.selector with neurons to have cleaner data

#How to window repair
#in the wr1
#1. treatment column: only need volume for the first K
#control between every K
#duration: 1 or 1.5 minute for K (prefer 1.5) 
#"at" and "duration" add up need to be less than "at" in the next row
#K+ 40mM, ATP and Capsaicin usually are duration of 4

tmp.rd<-WindowRepair(RD.170210.60.f.k.w1,complete=T)                 #save to tmp.rd, T or TRUE only when treatment was edit, use F or FALSE if not and select what were edit
Trace.Click.dev(tmp.rd)                                              #check to see if the window matches
RD.170210.60.f.k.w1<-tmp.rd                                          #assign tmp.rd to RD.
save(RD.170210.60.f.k.w1, file="RD.170210.60.f.k.w1.Rdata")          #save RD.

#####################################################
tmp.rd<-RD.170210.60.f.k.w1

#Do RDView for capsaicin, 40mM K+, and maybe atp
tmp.rd<-RDView(tmp.rd)

#Do ROIreview(tmp.rd) for drops only, weird shaped cells only
tmp.rd<-ROIreview(tmp.rd)

neurons<-cellzand(tmp.rd$bin)                        #capsaicin and K+40mM
drops<-cellzand(tmp.rd$bin,'drop',1)

?setdiff

neurons<-setdiff(neurons, drops)

# Assign back to and save new

#How to TraceBrewer
tmp.rd<-TraceBrewer(RD.170503.41.f.p1)
RD.170503.41.f.p1<-tmp.rd

#How to RDView
tmp.rd<-RDView(RD.170503.41.f.p1)        #Do K+40mM, ATP, and Capsaicin

######################################################
#How to locate
locator()
#######################################################
#How to binary column
load("C:\\Users\\Administrator\\Box Sync\\Calcium Imaging\\Voltage Gated Channel experiments\\Af5.1\\top 5\\170421.29.f.c.p2.direct response to Af5.1, cce9a, plexiva, atp, capsaicin\\Response.Rdata")
#drag file, in this case, Response, to the main R.data
Response
tmp.rd<-RD.170421.29.f.c.p2          #assign RD to tmp.rd
names(tmp.rd$bin)
tmp.rd$bin['af5.1.amp']=0            #assign the non response to 0
tmp.rd$bin['af5.1.amp']              
sum(tmp.rd$bin['af5.1.amp'])
tmp.rd$bin[Response,'af5.1.amp']=1   #assign response to 1
sum(tmp.rd$bin['af5.1.amp'])         #this should match with length(Response)
length(Response)                     #match with sum(tmp.rd$bin['af5.1.amp])
tmp.rd$bin['af5.1.amp']              #make sure the column looks right
RD.<-tmp.rd                          #assign back to RD and save
save
####################################################
#How to colSums
cols.to.sum<-select.list(names(tmp.rd$bin),multiple=T)     #select.list is select items from a list, select K+40mM, ATP, Capsaicin, af5.1amp, block, and de
cols.to.sum                                                #show what you select
colSums(tmp.rd$bin[cols.to.sum])                           #form row and column sums of the selected items

######################################################
#How to table in excel for ALL
tmp.rd<-RD.170421.29.f.c.p2
colSums(tmp.rd$bin[,collumns])
bin.summary<-colSums(tmp.rd$bin[,collumns])
write.csv(bin.summary,'all.csv')
#####################################################
#How to table in excel for neurons
tmp.rd<-RD.
neurons<-cellzand(tmp.rd$bin)                           #pick K.40mM and Capsaicin
drops<-cellzand(tmp.rd$bin, 'drop')
neurons<-setdiff(neurons,drops)
collumns<-select.list(names(tmp.rd$bin),multiple=T)     #pick atp, capsaicin, K.40, amp, blocks, de
colSums(tmp.rd$bin[neurons,collumns])
bin.summary<-colSums(tmp.rd$bin[neurons,collumns])
write.csv(bin.summary,'neurons.csv')

#################################################
#better version of table in excel
tmp.rd<-RD.170421.29.f.c.p2
collumns<-select.list(names(tmp.rd$bin),multiple=T)
colSums(tmp.rd$bin[neurons,collumns])
bin.summary<-colSums(tmp.rd$bin[,collumns])
###############################################
#ROIreview for IB4, gfp
tmp.rd<-RD.
tmp.rd<-ROIreview(tmp.rd, neurons)
tmp.rd$bin[,c('tritc.bin','gfp.bin','drop')]
tmp.rd$bin[is.na(tmp.rd$bin)]<-0
tmp.rd$bin[,c('tritc.bin','gfp.bin','drop')]
tmp.rd$bin[neurons,]
###############################################
#How to add images for ROIreview
tmp.rd<-ImageFiller(tmp.rd)
 [1] "bf.gfp.IB4"   
 [2] "gfp.IB4"       
 [3] "gfp"
 [4] "IB4"            
 [5] "fura"
 [6] "ROI"
 [7] "bf start"
 [8] "bf end"

###############################################
#How to add IB4 and GFP to .ALL
tmp.rd<-RD.170421.53.m.p1
neurons<-cellzand(tmp.rd$bin)
drops<-cellzand(tmp.rd$bin, 'drop')
neurons<-setdiff(neurons,drops)
tmp.rd<-ROIreview(tmp.rd, neurons)
tmp.rd$bin[,c('tritc.bin','gfp.bin','drop')]
tmp.rd$bin[is.na(tmp.rd$bin)]<-0
tmp.rd$bin[,c('tritc.bin','gfp.bin','drop')]
collumns<-select.list(names(tmp.rd$bin),multiple=T)
colSums(tmp.rd$bin[,collumns])
bin.summary<-colSums(tmp.rd$bin[,collumns])
write.csv(bin.summary,'all.csv')

#How to table and bar graph
collumns.to.sum<-select.list(names(RD.170421.29.f.c.p2$bin),multiple=T)   #to select multiple
bob<-colSums(RD.170421.29.f.c.p2$bin[neurons,collumns.to.sum])            #assign it (this is for neurons), leave the bracket blank for all

#neurons stuff
write.csv(bob, 'all.csv')                                                 #turns it to a excel
names(RD.170421.29.f.c.p2$bin)
RD.170421.29.f.c.p2$cell.types<-list()
neurons<-cellzand(RD.170421.29.f.c.p2$bin)
neurons
class(RD.170421.29.f.c.p2$cell.types)
RD.170421.29.f.c.p2$cell.types$neurons<-neurons
drops<-cellzand(RD.170421.29.f.c.p2$bin)
RD.170421.29.f.c.p2$cell.types$drops<-drops
RD.170421.29.f.c.p2$cell.types[['drops']]<-drops
RD.170421.29.f.c.p2$cell.types$neurons
RD.170421.29.f.c.p2$bin[RD.170421.29.f.c.p2$cell.types$neurons,]
######################################################
read.csv('age.comparison.csv')

age.comp<-read.csv('age.comparison.csv')
barplot(age.comp)
class(age.comp)
age.comp<-as.matrix(age.comp)
barplot(age.comp)
row.names(age.comp)<-age.comp[,1]
age.comp
age.comp<-age.comp[-1]
age.comp
age.comp<-read.csv('age.comparison.csv')
age.comp<-as.matrix(age.comp)
row.names(age.comp)<-age.comp[,1]
age.comp
age.comp<-read.csv('age.comparison.csv')
row.names(age.comp)<-age.comp[,1]
age.comp[,-1]
age.comp<-age.comp[,-1]
age.comp<-as.matrix(age.comp)

barplot(age.comp)
barplot(t(age.comp), beside=T)
legend('topright', names(t(age.comp)), fill=c('blue',red'), bty='n')
legend('topright', names(t(age.comp)), fill=c('blue',red'))
legend('topright', names(t(age.comp)), fill=c('blue','red'))
legend('topright', legend=names(t(age.comp)), fill=c('blue','red'))
names(t(age.comp))
names(age.comp)
age.comp

#######################################################################################
#How to graph
age.comp<-as.matrix(read.table('age.comparison.csv', header=T, row.names=1,sep=','))
graph.title<-"All Cells"                       					#can be change to other such as 'neurons'
age.comp<-t(age.comp)
age.comp.names<-select.list(colnames(age.comp),multiple=T)
require(RColorBrewer)
#my colors is based on the number of row(which is my experiments) to 
#color my graph
number.of.color<-length(row.names(age.comp))
cols<-brewer.pal(number.of.color,"Dark2")
dev.new(width=8,height=6)
par(mai=c(2,1,1,1),bg='gray70')
barplot(age.comp[,age.comp.names], beside=T, col=cols, font=2,las=2,main=graph.title)
legend('topright', legend=row.names(age.comp), fill=cols, cex=.7,bty='n')

#######################################################################################
#How to for loop

main.dir<-"C:/Users/Administrator/Box Sync/Calcium Imaging/Voltage Gated Channel experiments/Af5.1/top 5"

Kevins.barplotter<-function(main.dir,graph.title="yo"){
	setwd(main.dir)
	# select the folder for each experiment from the list
	exp.dir<-select.list(list.dirs(), multiple=T)

	setwd(main.dir)
	neurons.df.list<-list()
	for(i in 1:length(exp.dir)){
		#change the folder
		setwd(exp.dir[i])
		#find the experiment name
		exp.name<-list.files(pattern='RD.')
		exp.name<-sub('.Rdata','',exp.name)
		#find the new neuron dataframe
		csv.import<-list.files(pattern='neurons.csv')
		#import that dataframe
		neurons.csv<-read.csv(csv.import,col.names=c('responses',exp.name))
		#store that dataframe
		neurons.df.list[[exp.name]]<-neurons.csv
		#reset the working directory
		setwd(main.dir)
		}

	pulse.type<-select.list(colnames(neuron.df),multiple=T)
	#merge function; all.x=T and all.y=T to show a supertable and assign it to neuron.df
	neuron.df<-Reduce(function(x,y) merge(x,y,by='responses',all.x=T,all.y=T),neurons.df.list)
	#give row names and assign
	row.names(neuron.df)<-neuron.df[,1]
	#assign it to neuron.df with [-1] to get rid of the extra collumn
	neuron.df<-neuron.df[-1]
	#transpose
	neuron.df<-t(neuron.df)
	#do simmple math to find the percentage of neurons responses, etc.
	neuron.df.percent<-neuron.df/neuron.df[,'K..40mM']*100
	#my colors is based on the number of row(which is my experiments) to 
	#color my graph
	number.of.color<-length(row.names(neuron.df.percent[,pulse.type]))
	cols<-brewer.pal(number.of.color,"Dark2")
	dev.new(width=8,height=6)
	par(mai=c(2,1,1,1),bg='gray70')
	barplot(neuron.df.percent[,pulse.type],beside=T, col=cols, las=2,main=graph.title)
	legend('topright', legend=row.names(neuron.df.percent[,pulse.type]), fill=cols, cex=.7,bty='n')
}


##############################################################################################

legend('topright', legend=colnames(age.comp), fill=cols)
barplot(t(age.comp), beside=T, col=cols, las=2)
par(mai=c(0,0,0,0)
)
barplot(t(age.comp), beside=T, col=c('red','blue'), las=2)
par(mai=c(3,0,0,0))
barplot(t(age.comp), beside=T, col=c('red','blue'), las=2)

barplot(t(age.comp), beside=T, col=c('red','blue','yellow','green','purple'), las=2)
par(mai=c(3,1,1,1),bg='gray70', font=2)
barplot(t(age.comp), beside=T, col=c('red','blue'), las=2)
barplot(t(age.comp), beside=T, col=c('red','blue','yellow','green','purple'), las=2, border=1)
barplot(t(age.comp), beside=T, col=c('red','blue','yellow','green','purple'), las=2, main='Age.Comparison.ALL')
savehistory('lees.table.prep.r')

###################################################################################################
#cell type script
#Start with this if you haven't already done it. If you have, move past the hashtag line to cell type part 2

tmp.rd<-RD.170421.29.f.c.p2
Trace.Click.dev(tmp.rd)

#If window regions are off, repair the wr1.csv file, and use the function below
tmp.rd<-WindowRepair(tmp.rd, complete=T)

#if a region of trace data needs to be remove use region deleter
tmp.rd<-RegionDeleter(tmp.rd)

#RDview allows for an accurate scoring of a pulse
tmp.rd<-RDView(tmp.rd)
#ROIreview score visual aspects of the cells

tmp.rd<-ROIreview(tmp.rd, subset.n=500)

RD.170421.29.f.c.p2<-tmp.rd
save(RD.170421.29.f.c.p2, file="RD.170421.29.f.c.p2.Rdata")
##########################################################################################################
#cell type part 2

dropped<-cellzand(tmp.rd$bin,"drop",1)#selected bin and dropped
neurons<-cellzand(tmp.rd$bin, , 1)#select capsaicin and K.40mM
neurons<-setdiff(neurons, dropped)
greeR.Cells<-cellzand(tmp.rd$bin,"gfp.bin" ,1) #selected bin then gfp.bin
red.cells<-cellzand(tmp.rd$bin,"cy5.bin" ,1)
caG.Cells<-cellzand(tmp.rd$bin,  grep("caps",names(tmp.rd$bin),ignore.case=T, value=T), 1)

glia<-setdiff(tmp.rd$c.dat$id, neurons)
glia<-setdiff(glia, dropped)

peptidergic<-greeR.Cells
not.cgrp<-setdiff(neurons, greeR.Cells)

G.C<-intersect(greeR.Cells, caG.Cells)
G.0<-setdiff(greeR.Cells, G.C)

R.A<-intersect(not.cgrp, red.cells)

R.C<-intersect(R.A, caG.Cells)
R.A<-setdiff(R.A, R.C)

unlabeled<-setdiff(not.cgrp, R.A)
unlabeled<-setdiff(unlabeled, R.C)

UL<-intersect(large.cells.330,unlabeled)
UL<-intersect(UL,neurons)
UL<-setdiff(UL, caG.Cells)

US<-setdiff(unlabeled, UL)
US<-intersect(US,neurons)

#go through the unlabeled small cells in trace.click and mark as "1" any that have a dipping baseline.
thermos.low<-Trace.Click.dev(tmp.rd, US)
thermos.low<-thermos.low[[1]]

US<-setdiff(US, thermos.low)

US.C<-intersect(US, caG.Cells)
US<-setdiff(US, US.C)

cell.types<-named.list(
	neurons,
	glia,
	UL,
	G.C,
	G.0,
	R.A,
	R.C,
	US,
	US.C,
	thermos.low)
	
tmp.rd$cell.types<-cell.types
RD.170421.29.f.c.p2<-tmp.rd

save(RD.170421.29.f.c.p2, file="RD.170421.29.f.c.p2.Rdata")

##############################################################################################**data analysis**
#analysis into mario's cell class

#Sorting into Mario's cell broad cell classes
#Start out the same as a typical cell profiling-RDview and ROIreview
#Do these three things to each data folder before continuing
#Before ventruiong into RDview and ROireivew, use Trace.Click
#to amke usre you data looks good
## rename your data to tmp.rd for ease of use

tmp.rd<-RD.180614.58.m.p2

tcd(tmp.rd)
tmp.rd<-TraceBrewer(tmp.rd)
#If window regions are off, repair the wr1.csv file, and use the function below
#tmp.rd<-WindowRepair(tmp.rd, complete=T)

#if a region of trace data needs to be remove use region deleter
#tmp.rd<-RegionDeleter(tmp.rd)

#RDview allows for an accurate scoring of a pulse
tmp.rd<-RDView(tmp.rd)
#ROIreview score visual aspects of the cells
neurons<-cellzand(tmp.rd$bin, , 1)
tmp.rd<-ROIreview(tmp.rd,neurons,subset.n=500)
#This is how to save
RD.180614.58.m.p2<-tmp.rd
save(RD.180614.58.m.p2, file="RD.180614.58.m.p2.Rdata")

##############################################################

dropped<-cellzand(tmp.rd$bin,"drop",1)#selected bin and dropped
neurons<-cellzand(tmp.rd$bin, , 1)
neurons<-setdiff(neurons, dropped)
greeR.Cells<-cellzand(tmp.rd$bin,"gfp.bin" ,1) #selected bin then gfp.bin
red.cells<-cellzand(tmp.rd$bin,"cy5.bin" ,1)
caG.Cells<-cellzand(tmp.rd$bin,  grep("caps",names(tmp.rd$bin),ignore.case=T, value=T), 1)
aitc.cells<-cellzand(tmp.rd$bin,  grep("aitc",names(tmp.rd$bin),ignore.case=T, value=T), 1)
menth.cells<-cellzand(tmp.rd$bin, grep("menth",names(tmp.rd$bin),ignore.case=T, value=T), 1)
menth.only<-setdiff(menth.cells, aitc.cells)
large.cells.330<-cellzand(tmp.rd$c.dat,"area" ,330)#selected c.dat then area

glia<-setdiff(tmp.rd$c.dat$id, neurons)
glia<-setdiff(glia, dropped)

peptidergic<-greeR.Cells
not.cgrp<-setdiff(neurons, greeR.Cells)

#Sort green cells first by capsaicin then aitc
G.C<-intersect(greeR.Cells, caG.Cells)
G.0<-setdiff(greeR.Cells, G.C)

G.C.A<-intersect(G.C, aitc.cells)
G.C<-setdiff(G.C, G.C.A)

G.A<-intersect(G.0, aitc.cells)
G.A<-setdiff(G.A, G.C.A)
G.0<-setdiff(G.0, G.A)
G.M<-intersect(G.0, menth.cells)
G.0<-setdiff(G.0, G.M)

#This gives us G.C, G.C.A, G.A, G.M, and G.0 (zero) under the green cells
#next we seperate red from unlabeled

nonpep<-intersect(not.cgrp, red.cells)
nonpep<-setdiff(nonpep, menth.only)
unlabeled<-setdiff(not.cgrp, nonpep)

#Chase down the red classes, the two that are pretty unambiguous are R.A, R.C and R.other

R.A<-intersect(nonpep, aitc.cells)
R.A<-setdiff(R.A, caG.Cells)

R.other<-setdiff(nonpep, R.A)
R.C<-intersect(R.other, caG.Cells)
R.C<-setdiff(R.C, aitc.cells)

R.other<-setdiff(R.other, R.C)

#This gives us our red groups: R.A, R.C and R.other
#Finally we chase down our unlabeled groups (unlabeled)

thermos<-menth.only
thermos<-intersect(thermos, neurons)
unlabeled<-setdiff(unlabeled, thermos)
UL<-intersect(large.cells.330,unlabeled)
UL<-intersect(UL,neurons)
UL<-setdiff(UL, caG.Cells)
UL<-setdiff(UL, aitc.cells)

US<-setdiff(unlabeled, UL)
US<-intersect(US,neurons)

US.A<-intersect(US, aitc.cells)
US.A<-setdiff(US.A, caG.Cells)
US.C<-intersect(US, caG.Cells)
US.C<-setdiff(US.C, US.A)
US.0<-setdiff(US, US.A)
US.0<-setdiff(US.0, US.C)

#review the autosorted cell classes. Remove the cells that are not part of each class and put into discard pile, press "1" to move cells to the discard pile

#Remove anything that responds to capsaicin, menthol or aitc. Also remove cells labeled with IB4 or CGRP-GFP
#No UL for this experiment
UL.edit<-tcd(tmp.rd, c(UL))
	discard<-UL.edit[[1]]
UL<-setdiff(UL, discard)

#Remove cells that are not green or do not respond to menthol. Also remove red cells and cells responding to AITC	
#No UL for this experiment
G.M.edit<-tcd(tmp.rd, c(G.M))
	discard1<-G.M.edit[[1]]
G.M<-setdiff(G.M, discard1)
discard<-union(discard, discard1)

#Remove cells that respond to menthol, capsaicin or aitc, cells that are not green and cells that are red.
G.0.edit<-tcd(tmp.rd, c(G.0))
	discard1<-G.0.edit[[1]]
G.0<-setdiff(G.0, discard1)
discard<-union(discard, discard1)

#Remove cells that respond to Capsaicin or have "noisy" variable baseline, or cells that AITC response does not return to baseline
G.A.edit<-tcd(tmp.rd, c(G.A))
	discard1<-G.A.edit[[1]]
G.A<-setdiff(G.A, discard1)
discard<-union(discard, discard1)

G.C.edit<-tcd(tmp.rd, c(G.C))
	discard1<-G.C.edit[[1]]
G.C<-setdiff(G.C, discard1)
discard<-union(discard, discard1)

G.C.A.edit<-tcd(tmp.rd, c(G.C.A))
	discard1<-G.C.A.edit[[1]]
G.C.A<-setdiff(G.C.A, discard1)
discard<-union(discard, discard1)

R.A.edit<-tcd(tmp.rd, c(R.A))
	discard1<-R.A.edit[[1]]
R.A<-setdiff(R.A, discard1)
discard<-union(discard, discard1)

#Remove cells that respond to aitc or are green
R.C.edit<-tcd(tmp.rd, c(R.C))
	discard1<-R.C.edit[[1]]
R.C<-setdiff(R.C, discard1)
discard<-union(discard, discard1)

R.other.edit<-tcd(tmp.rd, c(R.other))
	discard1<-R.other.edit[[1]]
R.other<-setdiff(R.other, discard1)
discard<-union(discard, discard1)

#Remove cells that either do not respond to menthol or have an aitc response larger than the menthol response
thermos.edit<-tcd(tmp.rd, c(thermos))
	discard1<-thermos.edit[[1]]
thermos<-setdiff(thermos, discard1)
discard<-union(discard, discard1)

US.A.edit<-tcd(tmp.rd, c(US.A))
	discard1<-US.A.edit[[1]]
US.A<-setdiff(US.A, discard1)
discard<-union(discard, discard1)

US.C.edit<-tcd(tmp.rd, c(US.C))
	discard1<-US.C.edit[[1]]
US.C<-setdiff(US.C, discard1)
discard<-union(discard, discard1)

US.0.edit<-tcd(tmp.rd, c(US.0))
	discard1<-US.0.edit[[1]]
US.0<-setdiff(US.0, discard1)
discard<-union(discard, discard1)

#Sort the discard pile
hand.sort<-tcd(tmp.rd, c(discard))
#1 UL
#2 G.M
#3 G.0
#4 G.A
#5 G.C
#6 G.A.C
#7 R.A
#8 R.C
#9 R.other
#10 thermos
#11 US.A
#12 US.C

#Union between hand sorted and edited autosort. Sort UL into 4 groups based on R3J response (1=prop, 2=jagged, 3=ide, 4=NE), Create a few thermos groups
UL<-union(UL, hand.sort[[1]])
#UL.groups<-tcd(tmp.rd, c(UL))
	#UL.1<-UL.groups[[1]] #proprioceptor	
	#UL.2<-UL.groups[[2]] #jagged
	#UL.3<-UL.groups[[3]] #IDE only
	#UL.4<-UL.groups[[4]] #no effect
G.M<-union(G.M, hand.sort[[2]])
G.0<-union(G.0, hand.sort[[3]])
G.A<-union(G.A, hand.sort[[4]])
G.C<-union(G.C, hand.sort[[5]])
G.C.A<-union(G.C.A, hand.sort[[6]])
R.A<-union(R.A, hand.sort[[7]])
R.C<-union(R.C, hand.sort[[8]])
R.other<-union(R.other, hand.sort[[9]])
thermos<-union(thermos, hand.sort[[10]])
#thermos.groups<-tcd(tmp.rd, c(thermos))
	#thermos.high<-thermos.groups[[1]]
	#thermos.low<-thermos.groups[[2]]
	thermos.C<-intersect(thermos, caG.Cells)
US.A<-union(US.A, hand.sort[[11]])
US.C<-union(US.C, hand.sort[[12]])

cell.types<-named.list(
	neurons, 
	glia,
	#UL.1,
	#UL.2,
	#UL.3,
	#UL.4,
	UL,
	G.M,
	G.0,
	G.A,
	G.C,
	G.C.A,
	R.A,
	R.C,
	R.other,
	thermos,
	#thermos.high,
	#thermos.low,
	thermos.C,
	US.A, 
	US.C,
	US.0
	)

tmp.rd$cell.types<-cell.types

RD.180614.58.m.p2<-tmp.rd

save(RD.180614.58.m.p2, file="RD.180614.58.m.p2.Rdata")

tmp.rd<-census.brewer(tmp.rd)
RD.180614.58.m.p2<-tmp.rd

save(RD.180614.58.m.p2, file="RD.180614.58.m.p2.Rdata")
















