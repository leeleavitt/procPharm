#The question that can and has been raised is, 
#Are the effects you see independent or overlaping

#Using tcd one can create a collection of cells manually to add to
# The bin dataframe

#STEP 1: using tcd create your cell groups. Be sure to name with with 'r', and
#when you press q to exit be sure to save your cells
tmp.rd<-tcd(tmp.rd)

#STEP 2: Now add these groups to the bin dataframe using ColBinner
tmp.rd<-col_binner(tmp.rd, gt.names)

#STEP 3: Now using this function we create in a sense a barcode for us to see
#combinatorial effects.
#This function adds this information to the scp dataframe.
tmp.rd <- combinner(tmp.rd)

#Here we us pf.function to create a barcode out of select binary collumns,
#this function automatically updated the tmp, thus not requiring a "<-"
pf.function_2(tmp)
#this is how we look at the new collumn created with pf.function, we will call this response.classes
(col_to_summarize<-menu(names(tmp$bin)))

(response.classes<-sort(summary(tmp$bin[,col_to_summarize]),decreasing=T))

#pf.summary, allows you to take these classes to create new collumns in your binary data.frame
tmp<-pf.summary(tmp, response.classes)

#Now that the binary dataframe has been created, we can creat a table  to look at
TableBrewer(tmp, RD.181120.41.f.m1.p2$cell_types)

#This function will go into the new bin dataframes and create groups of
#cells that are scored as a 1 for a specific collumn/collumns
#select the "000, 101,etc..."
( newest_groups<-bin_to_group(tmp) )

(jasmine.cc_and_newest_groups<-c(jasmine.cc,newest_groups))
#Now that we a new updated and complicated census, this will allow you to dissect and look at it
census_viewer(tmp, jasmine.cc_and_newest_groups)
