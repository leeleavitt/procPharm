
#This function plots the stat(data.frame form) of all cells, and individual
#cell type densities.
#dat: RD data
#cells: Total group
#cell_types: How to seperate the groups
#stat: premade statistic in data.frame formate where row names are cell.names
#xlim_top: this is the maximun xlim value to display
#xlim_bottom: this is the minimun xlim value to display
#overlay: will plot the density plot ontop of the cells density plot
#dens_sep: This will plot out the densitys on seperate plots
#plot_new: Will create a new window for this plot
#abline_loc: where to display the added line to help display data better
density_ct_plotter<-function(dat, cells, cell_types,stat=dat$c.dat["area"],xlim_top=NULL, xlim_bottom=NULL,overlay=T,dense_sep=T,plot_new=T,env=NULL,dat.name=NULL, abline_loc=0){	
    par(xpd=F)
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
    require(RColorBrewer)
    color<-brewer.pal(8,"Dark2")
    color<-rep(color,10)
    #color<-sample(rainbow(length(cell_types),start=.2, end=.85))
    
    all.cells.density<-density(stat[,1])
    #Overlay plot used in bp.selector
    if(is.null(cell_types)){
    
        if(!is.null(dat$cell_types)){
            cell_types<-dat$cell_types
        }else{
            overlay=F
            dense_sep=F
        }
        #perform a logical test to determine whether to plot the cells
        #selected_cell_types<-list()
        #for(i in 1:length(cell_types)){
        #	print(paste(names(cell_types)[i],"=",length(cell_types[[i]])))
        #	if(length(cell_types[[i]])>10){
        #		selected_cell_types<-append(selected_cell_types,cell_types[i])
        #	}
        #}
        #cell_types<-selected_cell_types
    }else{
        bringToTop(-1)
        print("which Cell Types would you like to view on the plotter")
        selected_cell_types<-select.list(names(cell_types), multiple=T)
        cell_types<-cell_types[selected_cell_types]
    }

    if(dense_sep==T){
        plot_sep<-ceiling(sqrt(length(cell_types)+1))
        par(mfrow=c(plot_sep,plot_sep),mai=c(.25,.25,.25,.25))
    }
    
    if(is.null(xlim_top)){
        xlim_top<-max(stat[,1])
    }else{xlim_top<-xlim_top}
    
    if(is.null(xlim_bottom)){
        xlim_bottom<-min(stat[,1])
    }else{xlim_bottom<-xlim_bottom}
    density_window<-dev.cur()
    xlim<-c(xlim_bottom,xlim_top)
    dev.set(density_window)
    plot(
        all.cells.density, 
        xlim=xlim, 
        ylim=c(0,max(all.cells.density$y)*1.5), 
        pch="",lwd=3, col="black",
        main=names(stat)
    )
    polygon(all.cells.density,col="red",lwd=1)

    
    #Provide density plots with lines overla
    if(overlay==T){
        for(i in 1:length(cell_types)){
            if(length(cell_types[[i]])>2){
                cell_type_density<-density(stat[cell_types[[i]],])
                lines(cell_type_density, col="black", lwd=5) 
                lines(cell_type_density, col=color[i], lwd=2)  
            }
        }
    legend("topleft",legend=names(cell_types), fill=color, cex=.6,box.col="Black")
    }
    
    if(dense_sep){
    par(xpd=T)
        for(i in 1:length(cell_types)){
            if(length(cell_types[[i]])>2){
                cell_type_density<-density(stat[cell_types[[i]],])

                plot(
                    cell_type_density, col="black", lwd=5,
                    xlim=xlim, 
                    ylim=c(0,max(all.cells.density$y)*1.5),
                    main=paste(names(cell_types[i])," n=",length(cell_types[[i]])),
                    bty="l"
                )   
                abline(v=abline_loc,col="red")
                lines(cell_type_density, col=color[i], lwd=2) 
            }else{
                plot(0,0,pch="",main=paste(names(cell_types[i])," n=",length(cell_types[[i]])),bty="l")
            }
        }
        plot(0,0,main=NULL,xlab=NULL,ylab=NULL,xaxt=NULL,yaxt=NULL,bty="n",pch="")
        #legend("topleft",legend=names(cell_types), fill=color, cex=.8,bg="gray70")
        text(0+xinch(.2),0,dat.name, cex=1.1)

    }
}

# Display census in a barplot
barPlotter <- function(dat = NULL, cols = 'YlOrRd'){
    cat("\nREAD ME\nWelcome to barPlotter to use me,\nbarPlotter(dat=RD.experiment, col = \'YlOrRd\')\nCustomize your colors, go to this webpage to get other names \nhttps://www.datanovia.com/en/wp-content/uploads/dn-tutorials/ggplot2/figures/0101-rcolorbrewer-palette-rcolorbrewer-palettes-1.png")

    cat('\nSelect what you Want to be on your barPlot\n')
    table <- TableBrewer(dat, ,F,F)
    # for each row except the first (this is the number of cells in each cell type)
    tableList <- list()
    for(i in 2:dim(table)[1]){
        tableList[[ row.names(table)[i] ]] <- table[c(1,i),,drop=F]
    }
    #tableList[['ATP']] <- read.csv('./atp.csv', row.names=1)
    #tableList[['Ca-free ATP']] <- read.csv('./caFree.csv', row.names=1)
    #tableList[['MRS-2365']] <- read.csv('./mrs2365.csv', row.names=1)


    tablePercs <- Reduce(rbind,lapply(tableList, function(x) round(x[2,]/x[1,]*100, digits=2)))
    tablePercsMut <- as.matrix(rev(tablePercs[nrow(tablePercs):1,]))

    # BARPLOT
    require(RColorBrewer)
    cols <- rev(brewer.pal(5, cols))[length(tableList):1]
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

}

