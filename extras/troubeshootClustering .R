dat <- get(load("Y:/Jackson/THP experiments/200803.M.31.lavage.LPS.m3.p1/RD.200803.M.31.lavage.LPS.m3.p1.Rdata"))

dat$blc<- dat$blc[1:dim(dat$t.dat)[1],]

names(pam5)

pam5$clusinfo

objs <- c()
asw <- c()
for(i in 1:20){
    print(i)
    pam5 <- pam(t(t.dat[x.min:x.max,m.names]),k=i)
    objs[i] <- pam5$objective[1]

    asw[i] <- pam5$silinfo$avg.width
}

plot(objs)

k.best <- which.max(asw)
dev.new()
plot(asw, type= "b", main = "pam() clustering assessment",
     xlab= "k  (# clusters)", ylab = "average silhouette width")
axis(1, k.best, paste("best",k.best,sep="\n"), col = "red", col.axis = "red")


main_dir <- "Y:/Dermaxon/200901.23.f.p1/"
setwd(main_dir)
tmpRD <- get(load(list.files(pattern ="^RD.*[.]Rdata$")))

source("C:\\Users\\leele\\Documents\\procPharm\\R\\import.R")
source("C:\\Users\\leele\\Documents\\procPharm\\R\\tracePrep.R")
source("C:\\Users\\leele\\Documents\\procPharm\\R\\traceStats.R")


dat <- tmpRD

.libPaths()
.libPaths(.libPaths()[3] )

#install
devtools::document()
devtools::install()

## quickSource
tmpRD <- get(load("Y:/NP-Nehls/nAChR Research/Other data/223L/200812.30.m.m3.p1 K10_K15_K20_K40_mc_dose/RD.200812.30.m.m3.p1.Rdata"))

files.sources = list.files("C:/Users/leeLeavittLapTop/Documents/procPharm/procPharm/R/",full.names = T)
sapply(files.sources, source)

tmpRD <- drStatMaker(tmpRD)
names(tmpRD$scp)

tcd(tmpRD)


