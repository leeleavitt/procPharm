# #' This function creates the binary scores for AITC, menthol, capsaicin and K40
# #' @export
# traceProbMaker <- function(dat){
#     pyPharm <- reticulate::import('python_pharmer')
#     pulsesWithNN <- c("^AITC.*", "^[cC]aps.*", "^[mM]enth.*", "[kK][.]40.*")
#     nnNames <- c('aitc','capsaicin', 'menthol', 'k40')

#     for( i in 1:length(pulsesWithNN)){        
#         # Make sure the pulse exists
#         pulse <- grep(pulsesWithNN[i], dat$w.dat$wr1)
        
#         if(length(pulse) > 0){
#             # Grab the model
#             model <- pyPharm$modelLoader(nnNames[i])
            
#             # Snag the pulse for all cells
#             minWin <- min( pulse )
#             maxWin <- minWin + 119

#             pulseToScore <- as.data.frame(t(dat$blc[minWin:maxWin,-1]))


#             # Now use the python score all the responses of interest
#             tryCatch({
#                 featureFrame <- pyPharm$featureMaker(pulseToScore, 10)
#                 probs <- model$predict(featureFrame)
#                 colnames(probs) <- c(0,1)
#                 # Transfer these scoring to the binary dataframe
#                 pulseName <- grep(pulsesWithNN[i], names(dat$bin), value=T)[1]
#                 dat[['probs']][[pulseName]] <- probs
#                 dat$bin[,1] <- model$predict_classes(featureFrame)
#                 }
#                 , error=function(e) print(paste("Could not score", pulsesWithNN[i]))
#             )
#         }
#     }
#     return(dat)
# }

#' This function creates the binary scores for AITC, menthol, capsaicin and K40
#' @export
#' @param dat this is the RD.experiment
#' @param minute This is a question for collecting the window region if it is minute then the window will be prepared such that from the input window 4 minutes following. This is useful when the data points are not normallly collected or if there was a pause.
#' @param pulsesToScore This is a vector of regular expressions to select the nn for the pulses of AITC, Menthol, Capsaicin, and K.40mM
traceProbMaker <- function(dat, minute = TRUE, pulsesToScore = NA){
    pyPharm <- reticulate::import('python_pharmer')

    nnNames <- c('aitc', 'menthol', 'capsaicin', 'k40')
    if(is.na(pulsesToScore)){
        pulsesWithNN <- c("^[aA][iI][Tt][Cc]","^[mM][eE][nN][tT][hH]", "^[cC][aA][pP][sS]","^[kK].*40")    
    }else{
        pulsesWithNN <- pulsesToScore
        nnNames <- unlist(sapply(pulsesWithNN, function(x) grep(x, nnNames, value = T)))
    }

    # Make an uncertainty matrix, or collect it.
    if(is.null(dat$uncMat)){
        uncertainMat <- data.frame(matrix(nrow = dim(dat$c.dat)[1], ncol = 0))
        row.names(uncertainMat) <- row.names(dat$c.dat)
    }else{
        uncertainMat <- dat$uncMat
    }
    
    for( i in 1:length(pulsesWithNN)){        
        # Make sure the pulse exists
        window <- grep(pulsesWithNN[i], unique(dat$w.dat$wr1), value=T)[1]
        pulseWindow <- dat$w.dat[dat$w.dat$wr1 == window,]
        
        if(dim(pulseWindow)[1] > 0){
            # Grab the model
            model <- pyPharm$modelLoader(nnNames[i])
            
            # Snag the pulse for all cells
            if(!minute){
                minWin <- min( pulseWindow$Time )
                maxWin <- minWin + 119
                pulseToScore <- as.data.frame(t(dat$blc[minWin:maxWin,-1]))
            }else{
                windowStart <- min(pulseWindow$Time)
                windowEnd <- windowStart + 4
                windowLogic <- dat$w.dat$Time > windowStart & dat$w.dat$Time < windowEnd
                pulseToScore <- t(dat$blc[windowLogic,-1])
            }

            # Now use the python score all the responses of interest
            tryCatch({
                if(!minute){
                    featureFrame <- pyPharm$featureMaker(pulseToScore, 10)
                }else{
                    featureFrame <- pyPharm$featureMaker2(pulseToScore, 12)
                }
                # Transfer these scoring to the binary dataframe
                pulseName <- grep(pulsesWithNN[i], names(dat$bin), value=T)[1]
                dat$bin[,pulseName] <- model$predict_classes(featureFrame)

                # add to the uncMat
                probs <- model$predict(featureFrame)
                uncertainty <- sqrt(probs[,1]^2 + probs[,2]^2)
                uncertainMat[pulseName] <- uncertainty
                }
                , error=function(e) print(paste("Could not score", pulsesWithNN[i]))
            )
        }else{
            print(paste("Could not score", pulsesWithNN[i])) 
        }
    }

    dat$uncMat <- uncertainMat

    return(dat)
}

#' Function to convert the full image into individual cell images.
#' @param dat is the RD.experiment to load in. You should have images loaded into this list
#' @param image is the character name of the image within your list
#' @param channel is the rgb channel of your image to extrace
#' @param range is the side in all directions from center to define your view
#' @export
imageExtractor <- function(dat, image = 'img3', channel = 1, range = 20){
    slice <- (range * 2) + 1

    # Gather the center x and y data
    cellXValue <- dat$c.dat$center.x.simplified
    cellYValue <- dat$c.dat$center.y.simplified

    # This checks to see if we will be able to observe the cell
    canView <- 
        0 < (cellXValue - range) & 
        (cellXValue + range) < dim(dat[[ image ]] )[2] &
        0 < (cellYValue - range) &
        (cellYValue + range) < dim(dat[[ image ]])[1] 

    imageArray <- array(
        dim = c(
            dim(dat$c.dat[canView, ])[1], 
            slice, 
            slice, 
            length(channel)
        ) 
    )

    for( j in 1:dim(dat$c.dat[canView, ])[1] ){
        xLeft <- cellXValue[canView][j] - range
        xRight <- cellXValue[canView][j] + range

        yTop <- cellYValue[canView][j] - range
        yBottom <- cellYValue[canView][j] + range

        imageArray[j, , , ] <- dat[[ image ]][yTop:yBottom,xLeft:xRight, channel]
    }

    image <- list()
    image[['imageArray']] <- imageArray
    image[['cellNames']] <- dat$c.dat$id[canView]

    return(image)
}

#' Function to gather images used inside cell_type_modeler. Unlike imageExtractor, this function
#' can collect all cells. Where ImageExtractor struggles with cells near the border by simply ignoring them
#' this function pads these images with zeros in these regions.
#' @param dat RD.experiment
#' @param range is the side in all directions from center to define your view
#' @param image is the character name of the image within your list
#' @param channel is the rgb channel of your image to extrace
#' @export
imageExtractorAll <- function(dat, image = 'img3', channel = 1, range = 20){
    slice <- (range * 2) + 1
    mainImageArray <- array(dim = c(0, slice, slice, length(channel)))
    # Gather the center x and y data
    cellXValue <- dat$c.dat$center.x.simplified
    cellYValue <- dat$c.dat$center.y.simplified

    subImageArray <- array(
        dim = c(
            dim(dat$c.dat)[1], 
            slice, 
            slice, 
            length(channel)
        ) 
    )


    for( j in 1:dim(dat$c.dat)[1] ){
        xLeft <- cellXValue[j] - range
        xRight <- cellXValue[j] + range
        xDims <- xLeft:xRight

        yTop <- cellYValue[j] - range
        yBottom <- cellYValue[j] + range
        yDims <- yTop:yBottom

        yLogic <-   yDims > 0 &
                    yDims < dim(dat[[ image ]] )[1]

        xLogic <-   xDims > 0 &
                    xDims < dim(dat[[ image ]] )[2]
        
        # Make a temporary array full of zeros
        tmpArray <- array(dim = c(slice, slice, 3))
        tmpArray[is.na(tmpArray)] <- 0
        
        # This is how i obtain the region of the image that is within the boundaries.
        tmpArray[yLogic, xLogic, channel] <-  dat[[ image ]][yDims[yLogic], xDims[xLogic], channel]
                
        subImageArray[j, , , ] <- tmpArray[,,channel]
    }

    mainImageArray <- abind::abind(mainImageArray, subImageArray, along = 1)

    image <- list()
    image[['imageArray']] <- mainImageArray
    image[['cellNames']] <- dat$c.dat$id

    return(image)
}

#' Function to create the cy5.bin and the gfp.bin collumns in the binary
#' data frame. This also adds probability scores for the images
#' @export 
imageProbMaker <- function(dat, verbose = T){
    if(verbose){
        cat("\nThis function scores the cy5,gfp, and drops\nMake sure that:\nimg3: cy5 only\nimg4: gfp only\nimg8: dapi.lab.png\n")
    }
    # Load in the data
    tryCatch({
        pyPharm <- reticulate::import('python_pharmer')
        image <- imageExtractorAll(dat, 'img3', 1)
        # grab correct model
        model <- pyPharm$modelLoader('cy5')
        # predict classes
        predictedClasses <- model$predict_classes(image$imageArray)
        # assign classes
        dat$bin[image$cellNames, 'cy5.bin'] <- predictedClasses
        
        # Grab the probabilities
        predictedClassProbs <- model$predict(image$imageArray)
        predictedClassProbsDF <- dat$bin[, c(1,2)]
        predictedClassProbsDF[,c(1,2)] <- NA
        colnames(predictedClassProbsDF) <- c(0,1)
        predictedClassProbsDF[image$cellNames,] <- predictedClassProbs

        dat[['probs']][['cy5']] <- predictedClassProbsDF
    }, error=function(e) print("Could not score IB4, you are most likely missing image 3"))

    tryCatch({
        image <- imageExtractorAll(dat, 'img4', 2)
        # grab correct model
        model <- pyPharm$modelLoader('gfp')
        # predict classes
        predictedClasses <- model$predict_classes(image$imageArray)
        # assign classes
        dat$bin[image$cellNames, 'gfp.bin'] <- predictedClasses
        
        # Grab the probabilities
        predictedClassProbs <- model$predict(image$imageArray)
        predictedClassProbsDF <- dat$bin[, c(1,2)]
        predictedClassProbsDF[,c(1,2)] <- NA
        colnames(predictedClassProbsDF) <- c(0,1)
        predictedClassProbsDF[image$cellNames,] <- predictedClassProbs

        dat[['probs']][['gfp']] <- predictedClassProbsDF
    }, error=function(e)print("Could not score GFP, you are most likely missing image 4"))

    tryCatch({
        image <- imageExtractorAll(dat, 'img8', c(1,2,3))
        # grab correct model
        model <- pyPharm$modelLoader('drop')
        # predict classes
        predictedClasses <- model$predict_classes(image$imageArray)

        if(sum(predictedClasses) < (.6 * length(dat$c.dat$id)) ){
            # assign classes
            dat$bin[image$cellNames, 'drop'] <- predictedClasses
            
            # Grab the probabilities
            predictedClassProbs <- model$predict(image$imageArray)
            predictedClassProbsDF <- dat$bin[, c(1,2)]
            predictedClassProbsDF[,c(1,2)] <- NA
            colnames(predictedClassProbsDF) <- c(0,1)
            predictedClassProbsDF[image$cellNames,] <- predictedClassProbs

            dat[['probs']][['drop']] <- predictedClassProbsDF
        }else{
            dat$bin[image$cellNames, 'drop'] <- 0
            cat("\nThe drop model dropped too many cells.\nImage 8 is most likely the wrong image.")
        }
    }, error=function(e)print("Could not score drop, you are most likely missing image 8"))

    return(dat)
}

#' This function calculates uncertainty for each probability data frame in the 
#' RD.experiment. Use probMaker and imageProbMaker prior to this or you will
#' get an error.
#' @param dat this is the experiment list to add.
#' @export
uncertaintyMaker <- function(dat){
    tryCatch({
        for(i in 1:length(dat$probs)){
            dat$probs[[i]][is.na(dat$probs[[i]])] <- .5
            uncertainty <- sqrt(dat$probs[[i]][1]^2 + dat$probs[[i]][2]^2)
            dat$scp[paste0(names(dat$probs),'.','unc')] <- uncertainty
        }
        return(dat)
    }, error = function(e)
        cat("\n Probabilities have not bee created. Plese do both, \n imageProbMaker(RD.exerpiment) \n and \n traceProbMaker(RD.experiment)\n")
    )
}

labelBinder <- function(dat, windowSizeMin = 3, tType = "t.dat", winMin = 1.2, winMax = 5, model = "./VRC.h5", winSel = T, featureWindows = 24, window = NA){
    keras <- reticulate::import('keras')
    
    if(all(class(model) == 'character')){
        model <- keras::load_model_hdf5(model)
    }
    pyPharm <- reticulate::import('python_pharmer')

    # Lets create the window region
    wr <- dat$w.dat$wr1
    levs <- setdiff(unique(dat$w.dat$wr1), "")

    # Lets find the window sizes to include,
    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)

    # Perform a test to include windows of the right sizes
    windowLen <- x2s - x1s
    winCol <- names(windowLen[windowLen < winMax & windowLen > winMin])

    # Provides the user with control over the windows includesd
    if(winSel){
        cat('\nSelect windows to add to the labeled data\n')
        winCol <- select.list(levs, preselect = winCol,  multiple = T, title = "Select windows", graphics = T)
    }else{
        if(is.na(window)){
            winCol <- grep(window, levs, T, value = T)
        }else{
            winCol <- levs
        }
    }
    winStart <- sort(x1s[winCol])
    winEnd <- winStart + windowSizeMin

    # Test the window size since we are losing some data. due to differing window sizes
    windowSize <- c()
    for(i in 1:length(winStart)){
        windowLogic <-  dat$w.dat$Time > winStart[i] & 
                        dat$w.dat$Time < winEnd[i]

        windowSize[i] <- dim(dat$t.dat[-1][windowLogic, ])[1]
    }

    windowSizeMax <- max(windowSize)
    windowLogic <- (windowSizeMax - windowSize) < 10
    if(length(winStart[!windowLogic])> 0 ){
        cat("\nThese windows were not scored\n")
        print(names(winStart[!windowLogic]))
        cat("\nAnd these windows are defined to not be scored\n")
        print(setdiff(levs, names(winStart[windowLogic])))
    }

    # This is the way we set the maximun windo size to the minimum window size
    # for all window regions.
    winStart <- winStart[windowLogic]
    windowSizeMini <- min(windowSize)

    startPoints <- Reduce(c,lapply(winStart, function(x) which(dat$t.dat$Time == x, arr.ind=T) ))
    endPoints <- startPoints + windowSizeMini

    # Create the uncertainty matrix and add the binary scoring to the data frame
    uncertainMat <- data.frame(matrix(nrow = dim(dat$c.dat)[1], ncol = length(winStart)))
    for(i in 1:length(winStart)){
        pulseToScore <- as.data.frame(t(dat[[tType]][-1][startPoints[i]:endPoints[i],]))
        featureFrame <- pyPharm$featureMaker2(pulseToScore, featureWindows)
        
        probs <- model$predict(featureFrame)
        classes <- model$predict_classes(featureFrame)

        uncertainty <- sqrt(probs[,1]^2 + probs[,2]^2)
        uncertainMat[,i] <- uncertainty

        dat$bin[,names(winStart)[i]] <- classes
    }

    names(uncertainMat) <- names(winStart)
    row.names(uncertainMat) <- row.names(dat$c.dat)
    dat$uncMat <- uncertainMat

    # A plot to show the uncertainty per response type.
    dev.new(width = 8, height = 4)
    par(bty = 'l', mar = c(7,6,4,2))
    bpDim <- boxplot(
        uncertainMat,
        xaxt = "n",
        border = NA,
        boxfill = 'gray90',
        medcol = 'black',
        medlty = 1,
        medlwd =2,
        whisklty = 1,
        whisklwd = 2,
        whiskcol = 'black',
        staplelty = 1, 
        staplelwd = 2,
        staplecol = 'black',
        main = 'Uncertainty per response',
        ylab = 'Class probabilities\nsqrt(Sum of squares)'
    )

    # stripchart(
    #     uncertainMat, 
    #     vertical = T,
    #     jitter = .2,
    #     method = 'jitter',
    #     pch = 18,
    #     cex = .2,
    #     col = rgb(0,0,0,.2),
    #     add = T
    #     )

    par(xpd = T)
    text(
        seq(1,length(winStart), by=1),
        par('usr')[3] - yinch(.1),
        names(winStart),
        adj = 1,
        srt = 90,
        cex = .7
    )

    return(dat)
}

#' Using marios cell type models I will create the cell_type model list. This is BETA version,
#' Simply inputing a RD.experiment which contains AMCK, R3J, and and image of IB4 only at
#' img3 and an image of GFP only at img4, this function returns a list of model output matrices.
#' 
#' @export
cell_type_modeler <- function(dat, modelDir = NA){
    if(is.na(modelDir)){
        modelDir <- "Y:/Computer Setup/R/models/"
    }

    pyPharm <- reticulate::import("python_pharmer")
    driveLoc <- strsplit(getwd(), "/")[[1]][1]

    # R12 is very dirty lots of N15 reassignment
    # R13 vs N14 is very dirty, so I am going to try and clean it up using neuralNets
    aitcModel <- paste0(modelDir,"aitc.h5")
    menthModel <- paste0(modelDir,"menth.h5")
    capsModel <- paste0(modelDir,"caps.h5")
    k40Model <- paste0(modelDir,"k40.h5")
    r3jModel <- paste0(modelDir,"r3j.h5")

    # ImageModels
    ib4Model <- paste0(modelDir,"ib4.h5")
    gfpModel <- paste0(modelDir,"gfp.h5")

    models <- c(aitcModel, menthModel, capsModel, k40Model, r3jModel, ib4Model, gfpModel)

    modelNames <- c('aitc', 'menth', 'caps', 'k40', 'r3J', 'gfp', 'ib4')
    predictedClasses <- list()

    # Traces
    pulsesWithNN <- c("^[aA][iI][Tt][Cc]","^[mM][eE][nN][tT][hH]", "^[cC][aA][pP][sS]","^[kK].*40", '[rR](3|[I]{3})[jJ]')

    predictedClasses <- list()
    rowNames <- dat$c.dat$id
    for(i in 1:4){
        model <- invisible(keras::load_model_hdf5(models[i]))
        window <- grep(pulsesWithNN[i], unique(dat$w.dat$wr1), value=T)
        pulseWindow <- dat$w.dat[dat$w.dat$wr1 == window,]
        
        windowStart <- min(pulseWindow$Time)
        windowEnd <- windowStart + 4
        windowLogic <- dat$w.dat$Time > windowStart & dat$w.dat$Time < windowEnd
        pulseToScore <- t(dat$blc[windowLogic,-1])

        featureFrame <- pyPharm$featureMaker2(pulseToScore, 12)
        predictedClasses[[ modelNames[i] ]] <- model$predict(featureFrame)
        row.names(predictedClasses[[ modelNames[i] ]]) <- rowNames
    }

    # R3J
    tryCatch({
        i = 5
        model <- invisible(keras::load_model_hdf5(models[i]))
        winRegs <- unique(dat$w.dat$wr1)
        r3JLoc <- grep('[rR](3|[I]{3})[jJ]', winRegs)[1]
        kBeforeR3JLoc <- r3JLoc - 1
        kAfterR3JLoc <- r3JLoc + 1
        r3JStart <- min(which(dat$w.dat$wr1 == winRegs[kBeforeR3JLoc] , arr.ind=T))
        r3JEnd <- max(which(dat$w.dat$wr1 == winRegs[kAfterR3JLoc] , arr.ind=T))

        pulseToScore <- t(dat$blc[r3JStart:r3JEnd, -1])

        featureFrame <- pyPharm$featureMaker2(pulseToScore, 22)
        predictedClasses[[ modelNames[i] ]] <- model$predict(featureFrame)
        row.names(predictedClasses[[ modelNames[i] ]]) <- rowNames
    }, error = function(e) NULL)

    # Ib4
    i = 6
    model <- invisible(keras::load_model_hdf5(models[i]))
    image <- imageExtractorAll(dat, 'img3', c(1))[[1]]
    predictedClasses[[ modelNames[i] ]] <- model$predict(image)
    row.names(predictedClasses[[ modelNames[i] ]]) <- rowNames

    # GFP
    i = 7
    model <- invisible(keras::load_model_hdf5(models[i]))
    image <- imageExtractorAll(dat, 'img4', c(2))[[1]]
    predictedClasses[[ modelNames[i] ]] <- model$predict(image)
    row.names(predictedClasses[[ modelNames[i] ]]) <- rowNames

    # Unpack the predicted classes into a list of model frames
    modelFrames <- list()
    classNames <- c('L1', "L2", "L3", "L4", "L5", "L6", "G7", "G8", "G9", 'G10', 'R11', 'R12', 'R13', 'N14', 'N15', 'N16')
    for( i in 1:length(dat$c.dat$id)){
        # Create the data frame of all collected models
        modelFrame <- data.frame(
            matrix(
                nrow = length(predictedClasses), 
                ncol = dim(predictedClasses[[1]])[2]
        ))
        row.names(modelFrame) <- names(predictedClasses)
        colnames(modelFrame) <- classNames

        # Loop through each model and collect the cells scores
        for(j in 1:length(predictedClasses)){
            modelFrame[j, ] <- predictedClasses[[j]][dat$c.dat$id[i],]
        }
        
        modelFrames[[ dat$c.dat$id[i] ]] <- modelFrame
    }

    dat$cellTypeModel <- modelFrames
    
    # Add the kurtosis to the experiment
    kurtosisinfo <- lapply(dat$cellTypeModel, function(x){kurtosis(apply(x,2,sum))})
    dat$scp['cell_type_kurtosis'] <- Reduce(c,kurtosisinfo)

    return(dat)
}

#' Add the original cell types to the RD.experiment
#' This fucntion will look for double classified cells and return 
#' what is double classified
#' @export
cellTypeAdder <- function(dat){
    cell_type_id <- grep("^cell([_]|[.])types$",names(dat),value = T )[1]

    #selectedCT <- select.list(names(dat$cell_types), multiple = T)
    selectedCT <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'G7', 'G8', 'G9', 'G10', 'R11', 'R12', 'R13', 'N14', 'N15', 'N16', 'UC')

    cell_types <- data.frame(matrix(nrow = dim(dat$c.dat)[1], ncol = length(dat[[cell_type_id]][selectedCT])))
    colnames(cell_types) <- names(dat[[cell_type_id]][selectedCT])
    row.names(cell_types) <- dat$c.dat$id
    cell_types[is.na(cell_types)] <- 0

    dat$c.dat['cell_types'] <- NA
    for(i in 1:length(dat[[cell_type_id]][selectedCT])){
        tryCatch({
            cell_types[ dat[[cell_type_id]][[ selectedCT[i] ]], names(dat[[cell_type_id]][selectedCT])[i] ] <- 1
            dat$c.dat[dat[[cell_type_id]][selectedCT][[i]], 'cell_types'] <- names(dat[[cell_type_id]][selectedCT])[i]
        }, error = function(e) NULL)
    }

    # How many cells are double classified?
    ctSum <- apply(cell_types, 1, sum)
    multiClassCells <- names(ctSum[ctSum > 1])
    if(length(multiClassCells) > 0){
        cat("\nThese cells were double classified.\nThey will be referred to as the later class\n")
        for(i in 1:length(multiClassCells)){
            cat(multiClassCells[i],": ")
            ctLogic <- as.logical(cell_types[multiClassCells[i],])
            cat(names(cell_types[multiClassCells[i],][ctLogic]), sep=" ")
            cat("\n")
        }
    }

    return(dat)
}

#' Function taken from e1071 package
#' This observes the tails of a distribution
#' big kurtosis means small tails: certain
#' small kurtosis means long tails: uncertain
kurtosis <- function(x, na.rm = FALSE, type = 3){
    x <- x[!is.na(x)]
    
    if(!(type %in% (1 : 3)))
       stop("Invalid 'type' argument.")
    
    n <- length(x)
    x <- x - mean(x)
    r <- n * sum(x ^ 4) / (sum(x ^ 2) ^ 2)
    y <- if(type == 1)
        r - 3
    else if(type == 2) {
        if(n < 4){
            stop("Need at least 4 complete observations.")
        }
        ((n + 1) * (r - 3) + 6) * (n - 1) / ((n - 2) * (n - 3))
    }else{
        r * (1 - 1 / n) ^ 2 - 3
    }

    return(y)
}

#' Function to view the cell type models. Pretty good vis
modelViewer <- function(dat, cell, plot.new = T){
    cellTypes <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'G7', 'G8', 'G9', 'G10', 'R11', 'R12', 'R13', 'N14', 'N15', 'N16', 'UC')
    cols <- RColorBrewer::brewer.pal(n = dim(dat$cellTypeModel[[1]])[1], 'Dark2')
    if(plot.new){
        dev.new(width = 8, height = 3)
    }
    par(mar=c(5,4,4,5))
    kurtStat <- round(dat$scp[cell,"cell_type_kurtosis"], digits = 3)
    
    ctName <- ifelse(is.na(dat$c.dat[cell,'cell_types']), "Not classified", dat$c.dat[cell,'cell_types'])

    bpDims <- barplot(
        as.matrix(dat$cellTypeModel[[ cell ]]), 
        beside = F, 
        col = cols,
        ylim = c(0,10),
        xaxt = 'n',
        border = NA,
    )
    
    mainName <- paste(ctName, " : ", kurtStat)
    text(
        par('usr')[1] - xinch(.5),
        par('usr')[4] + yinch(.5),
        "Cell Name : Certainty",
        cex = 1.5, 
        font =2, adj = 0
    )

    tryCatch({
        text(
            par('usr')[1] - xinch(.5),
            par('usr')[4] + yinch(.25),
            mainName,
            cex = 1, 
            font =2, 
            adj = 0
        )
    })

    legend(par('usr')[2] + xinch(0),
        par('usr')[4] + yinch(0.1),
        legend = row.names(dat$cellTypeModel[[1]]),
        fill = cols,
        horiz = F,
        bty='n',
        border = NA)

    par(xpd=T)
    text(
        apply(as.matrix(bpDims), 1, mean),
        par('usr')[3] -yinch(.2),
        colnames(dat$cellTypeModel[[1]]),
        col = ifelse(dat$c.dat[cell,'cell_types'] == cellTypes, 'red', 'black'),
        font = ifelse(dat$c.dat[cell,'cell_types'] == cellTypes, 2, 1),
        cex = ifelse(dat$c.dat[cell,'cell_types'] == cellTypes, 1.5, 1)
        )
    
    text(
        par('usr')[1] - xinch(.5),
        par('usr')[3] - yinch(.65),
        "Mario: Neuralnetworks",
        cex = .5,
        adj = 0
    )


}

#' Function to apply a general purpose potassium model which was trained on a potassium 
#' dose response experiment. This will apply the model to any potassium window
#' except the first k40 window, which has its own model to use
#' @export
potassiumModeler <- function(dat, modelLoc = NA){
    if(is.na(modelLoc)){
        model <- keras::load_model_hdf5("Y:/Computer Setup/R/models/kdr.h5")
    }else{
        model <- keras::load_model_hdf5(modelLoc)
    }
    pyPharm <- reticulate::import('python_pharmer')

    windowSize <- 3
    featureWindows <- 12
    tType <- c("blc")

    # Here we are collecting the windows that were scored
    # These are all the smaller windows that were ensured to be correct
    wr <- dat$w.dat$wr1
    levs <- setdiff(unique(dat$w.dat$wr1), "")
    kWins <- grep("^[kK]", levs, value = T)
    cellTypeK40 <- grep("^[kK][.]40", levs, value = T)[1]
    kWins <- setdiff(kWins, cellTypeK40)

    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)
    windowLen <- x2s - x1s
    windowLen <- windowLen[kWins]

    winCol <- names(windowLen[windowLen < 2 & windowLen > 1.2])

    winStart <- sort(x1s[winCol])
    winEnd <- winStart + windowSize

    # Make an uncertainty matrix, or collect it.
    if(is.null(dat$uncMat)){
        uncertainMat <- data.frame(matrix(nrow = dim(dat$c.dat)[1], ncol = 0))
        row.names(uncertainMat) <- row.names(dat$c.dat)
    }else{
        uncertainMat <- dat$uncMat
    }

    #Apply the model to everything to see how well it performs
    for(i in 1:length(winStart)){
        windowLogic <-  dat$w.dat$Time > winStart[i] & 
                        dat$w.dat$Time < winEnd[i]

        pulseToScore <- as.data.frame(t(dat[[tType]][windowLogic,-1]))
        featureFrame <- pyPharm$featureMaker2(pulseToScore, featureWindows)
        
        dat$bin[,names(winStart)[i]] <- model$predict_classes(featureFrame)

        # add to the uncMat
        probs <- model$predict(featureFrame)
        uncertainty <- sqrt(probs[,1]^2 + probs[,2]^2)
        uncertainMat[names(winStart)[i]] <- uncertainty
        
    }

    dat$uncMat <- uncertainMat

    return(dat)
}


#' Function to apply a general purpose potassium model which was trained on a potassium 
#' dose response experiment. This will apply the model to any potassium window
#' except the first k40 window, which has its own model to use
#' @export
lungModeler <- function(dat){
    model <- keras::load_model_hdf5( "Y:/Computer Setup/R/models/lungMulti.h5")
    pyPharm <- reticulate::import('python_pharmer')

    windowSize <- 4
    featureWindows <- 12
    tType <- c("blc")

    # Here we are collecting the windows that were scored
    # These are all the smaller windows that were ensured to be correct
    wr <- dat$w.dat$wr1
    levs <- setdiff(unique(dat$w.dat$wr1), c("", 'epad'))

    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)
    windowLen <- x2s - x1s
    windowLen <- windowLen[levs]

    winStart <- sort(x1s[levs])
    winEnd <- winStart + windowSize

    # Make an uncertainty matrix, or collect it.
    if(is.null(dat$uncMat)){
        uncertainMat <- data.frame(matrix(nrow = dim(dat$c.dat)[1], ncol = 0))
        row.names(uncertainMat) <- row.names(dat$c.dat)
    }else{
        uncertainMat <- dat$uncMat
    }

    #Apply the model to everything to see how well it performs
    for(i in 1:length(winStart)){
        windowLogic <-  dat$w.dat$Time > winStart[i] & 
                        dat$w.dat$Time < winEnd[i]

        pulseToScore <- as.data.frame(t(dat[[tType]][windowLogic,-1]))
        featureFrame <- pyPharm$featureMaker2(pulseToScore, featureWindows)
        
        dat$bin[,names(winStart)[i]] <- model$predict_classes(featureFrame)

        # add to the uncMat
        probs <- model$predict(featureFrame)
        uncertainty <- sqrt(probs[,1]^2 + probs[,2]^2)
        uncertainMat[names(winStart)[i]] <- uncertainty
        
    }

    dat$uncMat <- uncertainMat

    return(dat)
}

#' This function allows you to choose models to apply to your applications
#' This will allow you to choose variouse models and apply it to the applications you desire
# load("Y:/Lee Leavitt/nicotinic receptor/201103.36.m.m1.p2 mc_mc.am2_PNU/RD.201103.36.m.m1.p2.Rdata")
#' @export
modelChooser <- function(dat, levPat = NA, modelName = NA, windowSize = 3){
    pyPharm <- reticulate::import('python_pharmer')
    cat("####################################################\n\n####################################################\n\n####################################################\nTHIS WILL WIPE THE SCORES OF YOUR PUSLES\n\nello, welcome to modelChooser, select a model to apply to your applications\nIt might break, so try again.")
    
    # Find and load up models based on user input
    modelLocations <- "Y:/Computer Setup/R/models/"
    if(is.na(modelName)){
        modelNames <- list.files(modelLocations)
        cat("As a suggestion try\nkdr.h5, trained on potassium dose response\nlungMulti.h5, trained on lung profiling experiments\n")
        selectedModelName <- select.list(modelNames)
    }else{
        selectedModelName <- modelName
    }
    model <- keras::load_model_hdf5(paste0(modelLocations, selectedModelName))

    # Define the window size to send into the neural networks
    if(is.na(windowSize)){
        cat("How many minutes per pulse should we feed in minutes?\n 3 min or 4 min is standard\n\nEnter an integer ")
        windowSize <- scan(n =1, what = integer())
    }

    # Define other model parameters
    featureWindows <- 12
    tType <- c("blc")

    # Select the windows to send through the models
    wr <- dat$w.dat$wr1
    levs <- setdiff(unique(wr), "")

    if(is.na(levPat)){
        cat("####################################################\n\n####################################################\n\n####################################################\n\nSelect the windows to place into the neural networks.\n")
        selectedLevs <- select.list(levs, multiple = T, title = "Select Windows")
    }else{
        selectedLevs <- grep(levPat, levs, TRUE, value = T)
    }

    # Obtain the start and end times of the selected windows
    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)
    windowLen <- x2s - x1s
    windowLen <- windowLen[selectedLevs]
    #make sure the window regions are larger than two minutes
    winCol <- names(windowLen[windowLen < 2])
    winStart <- sort(x1s[winCol])
    winEnd <- winStart + windowSize

     # Make an uncertainty matrix, or collect it.
    if(is.null(dat$uncMat)){
        uncertainMat <- data.frame(matrix(nrow = dim(dat$c.dat)[1], ncol = 0))
        row.names(uncertainMat) <- row.names(dat$c.dat)
    }else{
        uncertainMat <- dat$uncMat
    }

    # Apply the model to everything to see how well it performs
    for(i in 1:length(winStart)){
        windowLogic <-  dat$w.dat$Time > winStart[i] & 
                        dat$w.dat$Time < winEnd[i]

        pulseToScore <- as.data.frame(t(dat[[tType]][windowLogic,-1]))
        featureFrame <- pyPharm$featureMaker2(pulseToScore, featureWindows)
        
        dat$bin[,names(winStart)[i]] <- model$predict_classes(featureFrame)

        # add to the uncMat
        probs <- model$predict(featureFrame)
        uncertainty <- sqrt(probs[,1]^2 + probs[,2]^2)
        uncertainMat[names(winStart)[i]] <- uncertainty
        
    }

    dat$uncMat <- uncertainMat

    return(dat)



}
