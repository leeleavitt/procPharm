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
        window <- grep(pulsesWithNN[i], unique(dat$w.dat$wr1), value=T)
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

#' Function to create the cy5.bin and the gfp.bin collumns in the binary
#' data frame. This also adds probability scores for the images
#' @export 
imageProbMaker <- function(dat, verbose = T){
    if(verbose){
        cat("\nThis function scores the cy5,gfp, and drops\nMake sure that:\nimg3: cy5 only\nimg4: gfp only\nimg8: dapi.lab.png")
    }
    # Load in the data
    tryCatch({
        pyPharm <- reticulate::import('python_pharmer')
        image <- imageExtractor(dat, 'img3', 1)
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
    },
    error=function(e) print("Could not score IB4, you are most likely missing image 3"))

    tryCatch({
        image <- imageExtractor(dat, 'img4', 2)
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
        image <- imageExtractor(dat, 'img8', c(1,2,3))
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


