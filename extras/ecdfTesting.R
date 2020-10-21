getwd()
formals(select.list)$graphics <- T

source("./procPharm/R/traceStats.R")
source("./procPharm/R/visualization.R")
tmpRD <- get(load("Z:/Changshan Niu/gabaResearch/20201009.32m.m3.p5 gaba340_APN20_6_40uM/RD.20201009.32m.m3.p5.Rdata"))

# THis experiment needs to go through the de stat maker

# We need to determine automatically which repeat pulse is the testing pulse
#' This is a function that returns a regular expression to find the testPulse
testPulseFinder <- function(dat){
    # Find windows that start with k
    uniqueNames <-  grep("^K",names(dat$bin), value = T)
    # discover the min length of these repeat pulses
    toId <- min(sapply(strsplit(uniqueNames, "[.]"), function(x){length(x)}))

    #Only show all pulse by this length.
    pulseTypes <- sapply(uniqueNames, function(x){
        paste0(strsplit(x, "[.]")[[1]][1:toId], collapse = '.')
    })
    pulseTypesSumm <- summary(as.factor(pulseTypes))
    majorPulseTest <- names(which(pulseTypesSumm == max(pulseTypesSumm)))
    pusleTestSplit <- strsplit(majorPulseTest, "[.]")[[1]]

    # Now create the regulr expression that will get me to select the 
    pulseTestRegex <- paste0("^", pusleTestSplit[1], "[.]", pusleTestSplit[2])
    return(pulseTestRegex)
}

# Do the stats exists?
totalStat <- length(grep("ide[.]mmnorm", names(tmpRD$scp), value = T))

if(totalStat < 1){
    pulseTestRegex <- testPulseFinder(tmpRD)
    tmpRD <- ideStatMaker(tmpRD, testPulseNames = pulseTestRegex)
}

stat <- 'ide'

# Now lets plot the ecdf
allNames <- grep(paste0('[.]', stat, ".", "mmnorm"), names(tmpRD$scp), value= T)
controlNames <- grep("control", allNames, value = T)
testNames <- setdiff(allNames, controlNames)

controlChoices <- select.list(controlNames, multiple = T, title = "Select the controls to view")
testChoices <- select.list(testNames, multiple = T, title = "Select the test to view")

# Here we select the names of the collumns to observer
dev.new(
    width = 5,
    height = 15
)

ecdfPlotter(tmpRD, controlChoices, testChoices, legendSep = .3)



