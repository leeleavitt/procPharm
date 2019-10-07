#This program will take in each cell type then allow you to 
#1 discard these cells into a general discard bin, or overall drop
#2 Circulate the dropped cells, and then allow for these cells to be re classified
source("Y:/Box Sync/procpharm/procPharm 170210.r")

#we will for start in one of halens experiments
setwd("Y:/Halen/CalciumImaging/190912.44.f.m3.p1 GDNF hour incubation M2 with nicotine/")
#load the rdata
expName <- list.files(pattern='^RD.*.Rdata')
load(expName)

#Now lets dive into the cell types
rdTmp <- get(ls(pattern='^RD.'))
#Need to check for any cell types names
cellTypeId <- grep('^cell',names(rdTmp), value=T)
rdTmp[cellTypeId]
#we need to reduce the scope to just cell types, so remove neurons and glia from this
cellTypeNames <- names(rdTmp[[cellTypeId]])
cellTypesToClean <- setdiff(cellTypeNames,c('neurons','glia'))

#FOR EACH CELL TYPE
i = 3
cat('\n Please scan cell type ',names(rdTmp[[cellTypeId]][i]),'\n' )
cat('Place the cells to drop into key 1, You can always redo this by pressing shift+ that number\n')
drops <- c()
forTcd <- list(drops, rdTmp[[cellTypeId]][[i]])
names(forTcd) <- c('drops', names(rdTmp[[cellTypeId]][i]))


