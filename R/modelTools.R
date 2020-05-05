#' This function creates the binary scores for AITC, menthol, capsaicin and K40
#' @export
traceProbMaker <- function(dat){
    pyPharm <- reticulate::import('python_pharmer')
    pulsesWithNN <- c("^AITC.*", "^[cC]aps.*", "^[mM]enth.*", "[kK][.]40.*")
    nnNames <- c('aitc','capsaicin', 'menthol', 'k40')

    for( i in 1:length(pulsesWithNN)){        
        # Make sure the pulse exists
        pulse <- grep(pulsesWithNN[i], dat$w.dat$wr1)
        
        if(length(pulse) > 0){
            # Grab the model
            model <- pyPharm$modelLoader(nnNames[i])
            
            # Snag the pulse for all cells
            minWin <- min( pulse )
            maxWin <- minWin + 119

            pulseToScore <- as.data.frame(t(dat$blc[minWin:maxWin,-1]))


            # Now use the python score all the responses of interest
            tryCatch({
                featureFrame <- pyPharm$featureMaker(pulseToScore, 10)
                probs <- model$predict(featureFrame)
                colnames(probs) <- c(0,1)
                # Transfer these scoring to the binary dataframe
                pulseName <- grep(pulsesWithNN[i], names(dat$bin), value=T)
                dat[['probs']][[pulseName]] <- probs
                dat$bin[,pulseName] <- model$predict_classes(featureFrame)
                }
                , error=function(e) print(paste("Could not score", pulsesWithNN[i]))
            )
        }
    }
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
imageProbMaker <- function(dat){
    # Load in the data
    tryCatch({
        pyPharm <- reticulate::import('python_pharmer')
        image <- imageExtractor(dat, 'img3', 1)
        # grab correct model
        model <- pyPharm$modelLoader('cy5')
        # predict classes
        predictedClasses <- model$predict_classes(image$imageArray)
        # assign classes
        dat$bin[image$cellNames, 'drop'] <- predictedClasses
        
        # Grab the probabilities
        predictedClassProbs <- model$predict(image$imageArray)
        predictedClassProbsDF <- dat$bin[, c(1,2)]
        predictedClassProbsDF[,c(1,2)] <- NA
        colnames(predictedClassProbsDF) <- c(0,1)
        predictedClassProbsDF[image$cellNames,] <- predictedClassProbs

        dat[['probs']][['drop']] <- predictedClassProbsDF
    },
    error=function(e) print("Could not score IB4, you are most likely missing image 3"))

    tryCatch({
        image <- imageExtractor(dat, 'img4', 2)
        # grab correct model
        model <- pyPharm$modelLoader('gfp')
        # predict classes
        predictedClasses <- model$predict_classes(image$imageArray)
        # assign classes
        dat$bin[image$cellNames, 'drop'] <- predictedClasses
        
        # Grab the probabilities
        predictedClassProbs <- model$predict(image$imageArray)
        predictedClassProbsDF <- dat$bin[, c(1,2)]
        predictedClassProbsDF[,c(1,2)] <- NA
        colnames(predictedClassProbsDF) <- c(0,1)
        predictedClassProbsDF[image$cellNames,] <- predictedClassProbs

        dat[['probs']][['drop']] <- predictedClassProbsDF
    }, error=function(e)print("Could not score GFP, you are most likely missing image 4"))

    tryCatch({
        image <- imageExtractor(dat, 'img8', c(1,2,3))
        # grab correct model
        model <- pyPharm$modelLoader('drop')
        # predict classes
        predictedClasses <- model$predict_classes(image$imageArray)
        # assign classes
        dat$bin[image$cellNames, 'drop'] <- predictedClasses
        
        # Grab the probabilities
        predictedClassProbs <- model$predict(image$imageArray)
        predictedClassProbsDF <- dat$bin[, c(1,2)]
        predictedClassProbsDF[,c(1,2)] <- NA
        colnames(predictedClassProbsDF) <- c(0,1)
        predictedClassProbsDF[image$cellNames,] <- predictedClassProbs

        dat[['probs']][['drop']] <- predictedClassProbsDF
    }, error=function(e)print("Could not score drop, you are most likely missing image 8"))

    return(dat)
}

#' This function calculates uncertainty for each probability data frame in the 
#' RD.experiment. Use probMaker and imageProbMaker prior to this or you will
#' get an error.
#' @param RD.experiment
#' @export
uncertaintyMaker <- function(dat){
    for(i in 1:length(dat$probs)){
        dat$probs[[i]][is.na(dat$probs[[i]])] <- .5
        uncertainty <- sqrt(dat$probs[[i]][1]^2 + dat$probs[[i]][2]^2)
        dat$scp[paste0(names(dat$probs),'.','unc')] <- uncertainty
    }
    return(dat)
}
