require(procPharm)

tmpRD <- get(load("./extras/RD.200309.30.m.m3.p1.Rdata"))

# Import the python packages
require(reticulate)
#pyPharm <-import_from_path("python_pharmer","C:/Users/leele/Documents/procPharm/python_packages/python_pharmer/python_pharmer/")
pyPharm <- import('python_pharmer')

# Find where the AITC is located
minWin <- min( grep("^AITC.*", tmpRD$w.dat$wr1) )
maxWin <- minWin + 119

testTraces <- as.data.frame(t(tmpRD$blc[minWin:maxWin,-1]))

pyPharm$featureMaker(testTraces, 20)

