#' This function plots the stat(data.frame form) of all cells, and individual
#' cell type densities.
#' @param dat RD data
#' @param cells Total group
#' @param cell_types How to seperate the groups
#' @param stat premade statistic in data.frame formate where row names are cell.names
#' @param xlim_top this is the maximun xlim value to display
#' @param xlim_bottom this is the minimun xlim value to display
#' @param overlay will plot the density plot ontop of the cells density plot
#' @param dens_sep This will plot out the densitys on seperate plots
#' @param plot_new Will create a new window for this plot
#' @param abline_loc where to display the added line to help display data better
#' @export
density_ct_plotter<-function(dat, cells, cell_types = NA, stat=dat$c.dat["area"],xlim_top=NULL, xlim_bottom=NULL,overlay=T,dense_sep=T,plot_new=T,env=NULL,dat.name=NULL, abline_loc=0){
    stat[is.na(stat)]<-0

    par(xpd=F, bty='n')
    if(is.null(dat.name)){
        dat.name<-deparse(substitute(dat))
    }else{dat.name<-dat.name}

    if(plot_new & dense_sep){
        dev.new(width=10,height=10)
        density_window<-dev.cur()
        #density_plot<-dev.cur()
    }

    if(plot_new & dense_sep==F){
        dev.new(width=5,height=5)
        density_window<-dev.cur()
        #density_plot<-dev.cur()
    }
    #Now add a density plot per cell type to show the distribution of cell type effects
    color <- c(
        "dodgerblue2", "#E31A1C", # red
        "green4",
        "#6A3D9A", # purple
        "#FF7F00", # orange
        "black", "gold1",
        "skyblue2", "#FB9A99", # lt pink
        "palegreen2",
        "#CAB2D6", # lt purple
        "#FDBF6F", # lt orange
        "gray70", "khaki2",
        "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
        "darkturquoise", "green1", "yellow4", "yellow3",
        "darkorange4", "brown"
    )

    all.cells.density<-density(stat[,1])
    #Overlay plot used in bp.selector
    if(is.na(cell_types)){
        if(!is.null(dat$cell_types)){
            cell_types <- dat$cell_types
        }else{
            overlay=F
            dense_sep=F
        }
    }else{
        # tryCatch(bringToTop(-1), error=function(e)NULL)
        # print("which Cell Types would you like to view on the plotter")
        # selected_cell_types<-select.list(names(cell_types), multiple=T)
        # cell_types<-cell_types[selected_cell_types]
    }

    if(dense_sep==T){
        plot_sep<-ceiling(sqrt(length(cell_types)+2))
        par(mfrow=c(plot_sep,plot_sep),mai=c(.25,.25,.25,.25))
    }

    if( is.null(xlim_top) ){
        xlim_top<-max(stat[,1])
    }else{
        xlim_top<-xlim_top
    }

    if( is.null(xlim_bottom) ){
        xlim_bottom<-min(stat[,1])
    }else{
        xlim_bottom<-xlim_bottom
    }

    density_window<-dev.cur()
    xlim <- c(xlim_bottom, xlim_top)
    dev.set(density_window)
    plot(
        all.cells.density,
        xlab="",
        xlim=xlim,
        ylim=c(0,max(all.cells.density$y)*1.5),
        pch="",lwd=3, col="black",
        main=names(stat),
        xaxt = 'n'
        #main = ""
    )

    if(max(xlim) > 10){
        axisSeq <- ceiling(seq(xlim[1], xlim[2], length.out = 15))
    }else{
        axisSeq <- round(seq(xlim[1], xlim[2], length.out = 17), digits = 2)
    }
    axis(1, axisSeq)
    polygon(all.cells.density,col=rgb(1, 0, 0 ,0.2),lwd=1)

    #Provide density plots with lines overla
    par(xpd=T)
    if(overlay==T){
        for(i in 1:length(cell_types)){
            if(length(cell_types[[i]])>2){
                tryCatch({
                    cell_type_density <- density(stat[cell_types[[i]],])
                    lines(cell_type_density, col="black", lwd=5)
                    lines(cell_type_density, col=color[i], lwd=2)
                }, error = function(e) NULL)
            }
        }
    legend("topleft",legend=names(cell_types), fill=color, cex=.6,box.col="Black")
    }

    if(dense_sep){
        for(i in 1:length(cell_types)){
            if(length(cell_types[[i]])>2){
                cell_type_density<-density(stat[cell_types[[i]],])

                plot(
                    cell_type_density, col="black", lwd=5,
                    xlim=xlim,
                    ylim=c(0,max(all.cells.density$y)*1.5),
                    main=paste(names(cell_types[i])," n=",length(cell_types[[i]])),
                    bty="l",
                    xlab = ''
                )
                abline(v=abline_loc,col="red")
                lines(cell_type_density, col=color[i], lwd=2)
            }else{
                plot(0,0,pch="",main=paste(names(cell_types[i])," n=",length(cell_types[[i]])),bty="l", xlab='')
            }
        }
        plot(0,0,main=NULL,xlab=NULL,ylab=NULL,xaxt=NULL,yaxt=NULL,bty="n",pch="")
        #legend("topleft",legend=names(cell_types), fill=color, cex=.8,bg="gray70")
        text(0+xinch(.2),0,dat.name, cex=1.1)

    }
}

#' Display census in a barplot
#' @param dat this can either be and experiment RD.experiment or a vector of file paths see example
#' @param cols this is the colo.brewer pallette to use. \href{https://www.datanovia.com/en/wp-content/uploads/dn-tutorials/ggplot2/figures/0101-rcolorbrewer-palette-rcolorbrewer-palettes-colorblind-friendly-1.png}{click here}
#' @param selectCT logical T or F for cleaning up cell types
#' @param horiz logical T or F for how to plot the barplot
#' @examples
#' \dontrun{
#' # First step make single csv files from each experiment. Follow this general pattern
#' RD.experiment <- census_to_table(RD.experiment)
#' # Only select a single thing. For example only select amp
#' TableBrewer(RD.experiment)
#' #
#' #
#' # Once it is saved, open the cvs and manually rename the amp/first row something like "CNF EP1 1uM"
#' # additionally you have full control over the data in the cvs files. I warn against changing any numbers
#' # unless you are summing together other tables.
#' #
#' # Once it is saved, open the cvs and manually rename the amp something like "CNF EP1 1uM"
#' #
#' # Set the working directory to have the correct location. This will be where the files are located
#' # In the example you can see that i've set the working directory to
#' setwd('Y:/Cris Urcino/CNF experiments/Ep1/')
#' # This means i can access each csv file in the following way
#' tables <- c(
#'     "./190823.M.40.m3.p1 CNF-Ep1/amp.csv",
#'     "./200209.M.31.m3.p1 3uM.Ep1/amp.csv",
#'     "./200209.M.31.m3.p2 3uM.Ep1/amp.csv"
#' )
#' #
#' # Now barPlotter() it,
#' barPlotter(tables)
#' }
#' @export
barPlotter <- function(dat = NULL, cols = 'Dark2', selectCT = T, horiz=T){
    print(deparse(substitute(dat)))
    cat("\nREAD ME\nWelcome to barPlotter to use me,\nbarPlotter(dat=RD.experiment, col = \'YlOrRd\')\n\nCustomize your colors, TRY \n\n\'Set1\', \'Set2\', \'Pastel1\', \'Dark2\', \'PuBu\', \'Reds\'\n\nGo to this webpage to get other names \nhttps://www.datanovia.com/en/wp-content/uploads/dn-tutorials/ggplot2/figures/0101-rcolorbrewer-palette-rcolorbrewer-palettes-1.png")

    cat('\nSelect what you Want to be on your barPlot\n')

    if( class(dat) == 'character' ){
        tableList <- list()
        for( i in 1:length(dat)){
            table <- read.csv(dat[i], row.names = 1)
            tableList[[ row.names(table)[2] ]] <- table
        }
    }else{
        table <- TableBrewer(dat, ,F,F)
        # for each row except the first (this is the number of cells in each cell type)
        tableList <- list()
        for(i in 2:dim(table)[1]){
            tableList[[ row.names(table)[i] ]] <- table[c(1,i),,drop=F]
        }
    }

    # Select the cell type to display
    if(selectCT){
        tableNames <- colnames(table)
        cat("\nSelect Cell Types to display\n")
        tableNames <- select.list(tableNames, multiple=T, title="Select CellTypes")
        table <- table[tableNames]

        tableList <- list()
        for(i in 2:dim(table)[1]){
            tableList[[ row.names(table)[i] ]] <- table[c(1,i),,drop=F]
        }

    }


    #tableList[['ATP']] <- read.csv('./atp.csv', row.names=1)
    #tableList[['Ca-free ATP']] <- read.csv('./caFree.csv', row.names=1)
    #tableList[['MRS-2365']] <- read.csv('./mrs2365.csv', row.names=1)

    tablePercs <- Reduce(rbind,lapply(tableList, function(x) round(x[2,]/x[1,]*100, digits=2)))

    # BARPLOT
    if(!horiz){
        tablePercsMut <- as.matrix(rev(tablePercs[nrow(tablePercs):1,]))
        require(RColorBrewer)
        cols <- brewer.pal(10, cols)[1:length(tableList)]
        graphics.off()
        dev.new(width=5, height=12)
        par(mar=c(5,4,4,7))
        bpDims <- barplot(
            tablePercsMut,
            beside=T,
            horiz=T,
            col=cols,
            yaxt='n',
            xlab='% Cell Class Reponding',
            xlim=c(0,100),
            border=NA)

        # LEGEND
        par(xpd=T)
        responses <- names(tableList)
        legend(
            par('usr')[2]+xinch(.2),
            par('usr')[4]+yinch(.2),
            responses,
            fill=rev(cols),
            border=NA,
            bty='n',
            horiz=F,
            cex=.7)

        # YLAB
        yLocs <- apply(bpDims, 2, mean)
        xLocs <- rep(par('usr')[1]-xinch(.5), length(yLocs))
        tableLabs <- rev(names(tableList[[1]]))
        par(xpd=T)
        text(xLocs, yLocs, tableLabs)

        # BARLABS
        tablePercsMut <- tablePercsMut[nrow(tablePercsMut):1,,drop=F]
        bpDims <- bpDims[nrow(bpDims):1,, drop=F]

        for(i in 1:nrow(bpDims)){
            xLocs <- tablePercsMut[i,] + xinch(.3)
            yLocs <- bpDims[i,] +yinch(.02)
            textToPlace <- paste(rev(tableList[[i]][2,]), '/', rev(tableList[[i]][1,]))
            text(xLocs, yLocs, textToPlace, cex=.6)
        }
        par(xpd=T)
    }else{
        require(RColorBrewer)
        tablePercsMut <- as.matrix(tablePercs)

        #cols <- rev(brewer.pal(10, cols))[1:length(tableList)]
        cols <- brewer.pal(10, cols)[1:length(tableList)]

        graphics.off()
        dev.new(width=12, height=5)
        par(mar=c(5,5,5,2))
        bpDims <- barplot(
            tablePercsMut,
            beside=T,
            horiz=F,
            col=cols,
            xaxt='n',
            ylab='% Cell Class Reponding',
            ylim=c(0,100),
            border=NA)

        # LEGEND
        par(xpd=T)
        responses <- names(tableList)
        legend(
            mean(c(par('usr')[1], par('usr')[2]))-xinch(1.8),
            par('usr')[4]+yinch(.8),
            responses,
            fill=cols,
            border=NA,
            bty='n',
            horiz=T,
            cex=.7)

        # XLAB
        xLocs <- apply(bpDims, 2, mean)
        yLocs <- rep(par('usr')[3]-yinch(.3), length(xLocs))
        tableLabs <- names(tableList[[1]])
        par(xpd=T)
        text(xLocs, yLocs, tableLabs,srt=45)

        # BARLABS
        #tablePercsMut <- tablePercsMut[1:nrow(tablePercsMut),,drop=F]
        #print(tablePercsMut)
        #bpDims <- bpDims[1:nrow(bpDims),, drop=F]

        for(i in 1:nrow(bpDims)){
            yLocs <- tablePercsMut[i,] + yinch(.25)
            xLocs <- bpDims[i,] - xinch(0.02)
            textToPlace <- paste(tableList[[i]][2,], '/', tableList[[i]][1,])
            text(xLocs, yLocs, textToPlace, cex=.6, srt=90)
        }
    }

    dat.name <- deparse(substitute(dat))
    if(any(dat.name %in% c("tmp.rd", "tmpRD","tmp"))){
        dat.name <- ls(pattern = "^RD[.]", envir = .GlobalEnv)
        if(length(dat.name) > 1){
            dat.name <- deparse(substitute(dat))
        }
    }

    text(par('usr')[1]-xinch(.1), par('usr')[3]-yinch(.88), dat.name, cex=.7)

}

#' @export
boxPlotter <- function(mat, xlim_top=NULL, xlim_bottom=NULL, activewindows, controlwindows ){
    tryCatch({
        mainName <- paste(paste(activewindows, collapse = " + "), "/\n", paste(controlwindows, collapse = "+"))
        xLabel <- "Active.Max/Control.Max"
    },error=function(e){
        mainName <<- names(mat)[1]
        xLabel <- NA
    })


    xLabel <- "Active.Max/Control.Max"

    if( is.null(xlim_top) ){
        xlim_top<-max(stat[,1])
    }else{
        xlim_top<-xlim_top
    }

    if( is.null(xlim_bottom) ){
        xlim_bottom<-min(stat[,1])
    }else{
        xlim_bottom<-xlim_bottom
    }

    xlim <- c(xlim_bottom, xlim_top)

    bpStats <- boxplot(
        mat[,1],
        outline=F,
        ylim  = xlim,
        boxfill=rgb(1,1,1,.6,maxColorValue = 1),
        width=5,
        lty=1,
        lwd=3,
        main=mainName,
        xlab=xLabel,
        horizontal=T,
        xaxt='n')
    par(xpd=T)

    tryCatch({
        # Jittered cell names to view ontop of boxplot
        if( length(row.names(mat)) > 500 ){
            abovequartiles <-
                mat[,1] < bpStats$stats[1,1] |
                mat[,1] > bpStats$stats[5,1]

            cellsToView <- row.names(mat)[abovequartiles]
        }else{
            cellsToView <- row.names(mat)
        }
        text(
            jitter(
                rep(
                    1,
                    length(mat[cellsToView, 1])
                ),10
            )~mat[cellsToView,1],
            labels=row.names(mat[cellsToView,]),
            cex=.5,
            col=rgb(1,1,1,2, maxColorValue=10)
        )#,ylim=c(0,2.5), add=T, vertical=T, method="jitter", jitter=.2)
    }, error=function(e) NULL)

}


#' This is a function that returns a regular expression to find the testPulse
#' @export
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


#' This function displayed the newly create stats in an automated way.
#' It depends on the cell types present to display the stat.
#' #' Adding more controls at the begining really improves the view of this
#' @param dat RD.experiment
#' @param stat slection between 'de' which returns the direct effect stat or 'ide' which displays the indirect effect stat
#' @param controlToView a major highlighter of this vixualization is the ability to compare against a control window region
#' @export
ecdfPlotter <- function(dat, cells = NA, controlNames, testNames, legendSep = 0.2, rdName = NA, cell_types = NA, rowLayout = 4){
    # if the rdName is NA
    if(!is.na(rdName)){
        mainName <- rdName
    }else{
        mainName <- deparse(substitute(dat))
    }

    if(is.na(cell_types)){
        cellTypes <- c("L1","L2","L3","L4","L5","L6","G7","G8","G9","G10","R11","R12","R13","N14","N15","N16", "UC")
    }else{
        cellTypes <- cell_types
    }


    if(length(controlNames) > 1 ){
        collumns <- c(controlNames, testNames)

        controlColorFunc <- colorRampPalette(c("black", 'gray50'))
        cols <- c(
            controlColorFunc(length(controlNames)),
            rev(rev(RColorBrewer::brewer.pal(n = length(collumns), 'Dark2')))
        )

        lwds <- c(
            rep(3, length(controlNames)),
            rep(2, length(collumns))
        )
    }else{
        collumns <- c(testNames)
        colorFunc <- colorRampPalette(c("hotpink4", 'plum1'))

        cols <- colorFunc(length(collumns))

        lwds <- rep(2, length(collumns))

    }

    # Now Plot it up!
    mar <- c(0,4,4,0)

    # Creat a layout to fill in
    #rowLayout <- 8
    cellTypeTotal <- length(cell_types)
    # only 6 cell types allowed per collumn
    # calculate the number of collumns
    layoutCollumns <- ceiling(cellTypeTotal / rowLayout)

    ecdfSeq <- c(seq(1, cellTypeTotal, 1) + 1, rep(0, (layoutCollumns * rowLayout) - cellTypeTotal))
    dim(ecdfSeq) <- c(rowLayout, layoutCollumns)

    titlePos <- rep(1, layoutCollumns)
    xaxisPos <- seq(cellTypeTotal+2, length.out = layoutCollumns)
    legendPos <- rep(xaxisPos[length(xaxisPos)]+1, layoutCollumns)

    layoutMat <- rbind(
        titlePos,
        ecdfSeq,
        xaxisPos,
        legendPos
    )

    print(layoutMat)


    #par(mfrow = c( (length(cellTypes) + 4), 1), mar = mar)
    layout(layoutMat)
    plot(0, xlim = c(-1,1), pch = '', bty = 'n', yaxt = 'n', ylab='', xaxt='n', main = mainName)
    mar[3] <- 0

    for(i in 1:length(cellTypes)){
        par(mar = mar, las = 2)
        tryCatch({
            if(!is.na(cells)){
                cellsToView <- intersect(cells, dat$cell_types[[ cellTypes[i] ]])
            }else{
                cellsToView <- dat$cell_types[[ cellTypes[i] ]]
            }
            plot(
                ecdf(dat$scp[cellsToView , collumns[1]]),
                col = cols[1],
                xlim = c(-1,1),
                main = '',
                ylab = cellTypes[i],
                yaxt = 'n',
                xaxt = 'n',
                bty = 'n',
                do.points=F,
                col.01line = NULL,
                lwd = lwds[1]
            )
            for(j in 2:length(collumns)){
                lines(
                    ecdf(dat$scp[cellsToView , collumns[j]]),
                    type = 'l',
                    col = cols[j],
                    do.points=F,
                    col.01line = NULL,
                    lwd = lwds[j]
                )
            }

            legend('topleft', legend = paste("N= ", length(dat$cell_types[[cellTypes[i] ]])), bty = 'n')

        }, error = function(e){
            plot(
                0,0,
                ylim =c(0,1),
                xlim = c(-1,1),
                yaxt = 'n',
                xaxt = 'n',
                bty =
                'n',
                ylab = cellTypes[i],
                main = "", pch =''
            )
            text(0,0.5,"No Values")
            legend(
                'topleft',
                legend = paste("N= ", length(dat$cell_types[[cellTypes[i] ]])),
                bty = 'n'
            )
        })
        abline(v = c(0,-1), h = c(0))
    }

    #Plot of the x axis
    for(i in 1:layoutCollumns){
        mar[1] <- 3
        par(mar = mar)
        plot(0, xlim = c(-1,1), pch = '', bty = 'n', yaxt = 'n', ylab='')
    }

    # plot for legend
    plot(0, xlim = c(-1,1), pch = '', bty = 'n', yaxt = 'n', ylab='', xaxt='n')
    compNames <- Reduce(c, lapply(strsplit(collumns, "_"),function(x){x[1]}))
    par(xpd=T)
    legendVals <- sub("[.]max[.]ide", "", collumns)
    legendVals <- sub("[_]", "\n", legendVals)
    legendVals <- sub("[.]mmnorm", "", legendVals)

    gsub("[.]mmnorm", '', collumns)

    legend('top',
        legendVals[1:4] ,
        fill = cols[1:4],
        bty = 'n',
        border = NA,
        horiz = T,
        text.width = legendSep
    )

    if( length(legendVals) > 4){
        legend('bottom',
            legendVals[5:length(legendVals)] ,
            fill = cols[5:length(legendVals)],
            bty = 'n',
            border = NA,
            horiz = T,
            text.width = legendSep
        )
    }

}
