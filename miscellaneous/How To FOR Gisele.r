# First look at the data with Trace.Click.dev()
Trace.Click.dev(RD.180126.60.m.p1)

# If window region are off use, WindowRepair
tmp.rd<-WindowRepair(RD.180126.60.m.p1)
# After making correct changes to the wr.csv file, we doubled checked
#the work of the Windowrepair function with Trace.Click.dev
#We then double check and preeat the process until it is correct.

#If a region needs to be removed from the data use,
tmp.rd<-RegionDeleter(RD.180126.60.m.p1)

# Once we know the changes we made were correct, we reassign the name back to the
#original name.
#### This is how to save your work.
RD.180126.60.m.p1<-tmp.rd
save(RD.180126.60.m.p1, file="RD.180126.60.m.p1.Rdata")

# The first you need to do is use ROIreview to drop weird shaped cell, and score
#green and red.
tmp.rd<-ROIreview(RD.180126.60.m.p1)
RD.180126.60.m.p1<-tmp.rd
save(RD.180126.60.m.p1, file="RD.180126.60.m.p1.Rdata")

#Now score the cells with RDView.  This function will help you work
# DO NOT exit out of the window with the x button, only click DONE1
# through each response type to eventually be used in TableBrewer.
# Enter 0 when finshed!
tmp.rd<-RDView(RD.180126.60.m.p1)
RD.180126.60.m.p1<-tmp.rd
save(RD.180126.60.m.p1, file="RD.180126.60.m.p1.Rdata")

##############################################################################
#first find dropped cells
dropped.cells<-cellzand(RD.180126.60.m.p1$bin, "drop", 1)
#now select all cells
all.cells<-RD.180126.60.m.p1$c.dat$id

clean.cells<-setdiff(all.cells, dropped.cells)		

cells.of.interest<-Trace.Click.dev(RD.180126.60.m.p1, clean.cells)

# If you have quantified a group of cells (for example g.names4), 
# we need to collect this information in the bin data.frame

#rename this group so it is easier to deal with
b.d<-cells.of.interest$g.names4

#create a blank collumn
RD.180126.60.m.p1$bin["ro51.2mM.d.b"]<-0

#Now make all the rows from above, score in the collumn created above as 1
RD.180126.60.m.p1$bin[b.d,"ro51.2mM.d.b"]<-1

#Create 2 groups for the TableBrewer
dropped<-cellzand(RD.180126.60.m.p1$bin,"drop",1)#selected bin and dropped
neurons<-cellzand(RD.180126.60.m.p1$bin, , 1)
neurons<-setdiff(neurons, dropped)
glia<-setdiff(RD.180126.60.m.p1$c.dat$id, neurons)
glia<-setdiff(glia, dropped)
cell.types<-named.list(neurons, glia)

RD.180126.60.m.p1$cell.types<-cell.types

#######################################
TableBrewer(RD.180126.60.m.p1)

###############################################



















