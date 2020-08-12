#' New Boxplot Function
#' @param dat : Experiment list RD.
#' @param l.cells: Cells in a list format
#' @param dat.name: Dataframe to pull data from
#' @param col.name: collumn name  to get data for the boxplot
#' @param jitter.f: Factor of jitter to  accomplish
#' @param pts: points to add to boxplot
#' @param notchs: logical (T/F) stand for notch selection 
#' @param col.name c("area","mean.gfp.start","mean.cy5.start")
#' @export 
boxPlotList<-function(dat,l.cells=NULL,dat.name="c.dat",col.name=NULL,jitter.f=.5,pts=T, notchs=F, bplog="y", sort=T){

    #back up operation to fill with cells for the boxplotting
    if(is.null(l.cells)){
        l.cells<-dat$cell_types
    }else{
        l.cells<-l.cells
    }
    l.cells <- l.cells[lengths(l.cells)>0]
    
    if(is.null(dat.name)){
        dat.name<-"c.dat"
    }else{
        dat.name<-dat.name
    }
    if(is.null(col.name)){
        col.name<-select.list(names(dat[[dat.name]]), multiple=T)
    }else{
        col.name<-col.name
    }
    
    #Create a blank list to fill with information
    l.info<-list()
    
    l.cell.types<-names(l.cells)
    l.cell.types<-select.list(l.cell.types,multiple=T)
    
    #First create a boxplot to get the median statistics
    #but first use a for loop to gather the data needed
    for(i in 1:length(l.cell.types)){
        l.info[[ l.cell.types[i] ]]<-dat[[dat.name]][l.cells[[ l.cell.types[i] ]],col.name[1]]
    }

    #open a window
    dev.new()
    bp.stats<-boxplot(l.info)#plot the boxplot and assign it to an object ot gather stats
    colnames(bp.stats$stats)<-bp.stats$names #rename the collumn in the stats portion
    #reorder the data based on the median measure and gather the cell names
    if(sort==T){
        l.cell.types<-colnames(bp.stats$stats)[order(bp.stats$stats[3,], decreasing=T)] 
    }
    dev.off()#tunr off window
    
    #now create an empty list to refill with data
    l.info<-list()
    #once again regather the data
    for(i in l.cell.types){
        l.info[[i]]<-dat[[dat.name]][l.cells[[ i ]],c("id",col.name)]
    }
    
    # Now begin createing a dataframe to creata boxplot that can be intereacted with based on clicking
    l.cell.types<-names(l.info)
    bp.width<-vector()
    for(i in 1:length(l.cell.types)){
        l.info[[i]]["xplot"]<-jitter(rep(i,length(l.cells[[ l.cell.types[i] ]])),jitter.f)
        l.info[[i]]["cell.type"]<-l.cell.types[i]
        l.info[[i]]["cell.type.total"]<-length(l.cells[[ l.cell.types[i] ]])
        l.info[[i]]["cell.type.total.cb"]<-paste(l.cell.types[i],":",length(l.cells[[ l.cell.types[i] ]]),sep=" ")
        bp.width[i]<-length(l.cells[[ l.cell.types[i] ]])
    }
    
    #reduce the list into a dataframe
    bp.df<-Reduce(rbind,l.info)
    #Make the collum of cell types has a levels input for the boxplot below
    #this will allow it to be plotted based on above ordering
    bp.df$cell.type.total.cb<-ordered(bp.df$cell.type.total.cb,levels=unique(as.character(bp.df$cell.type.total.cb)))
    
    #now Boxplot
    dev.new(width=8, height=(3*length(col.name)))
    par(mfrow=c(length(col.name),1), bty="l")
    for(i in 1:length(col.name)){
        boxplot(get(col.name[i])~cell.type.total.cb, data=bp.df, varwidth=T,las=2, lwd=1.5,lty=1, outline=T, log=bplog, notch=notchs,main=tools::toTitleCase(gsub("\\.", " ", col.name[i])))
        
        if(pts){
            text(bp.df[,"xplot"], bp.df[,col.name[i]], "*", cex=.5)
            #text(bp.df[,"xplot"], bp.df[,col.name[i]], bp.df[,"id"], cex=.5)
        }else{}
    }
    
    bp.sel<-select.list(col.name, title="Select a Bp")
    
    tryCatch(windows(width=12, height=6,xpos=0, ypos=10), error=function(e) windows(width=12, height=6))
    bp.win<-dev.cur()
    
    tryCatch(windows(width=14,height=4,xpos=0, ypos=540), error=function(e) windows(width=14,height=4))
    click.window<-dev.cur()
    
    dev.set(bp.win)
    par(mai=c(2,1,1,1), bty="l")
    final.bp<-boxplot(get(bp.sel)~cell.type.total.cb, data=bp.df, varwidth=T,las=2, cex=.8, lwd=1.5,lty=1, outline=T, log=bplog, notch=notchs,main=tools::toTitleCase(gsub("\\.", " ", bp.sel)))
    text(bp.df[,"xplot"], bp.df[,bp.sel], bp.df[,"id"], cex=.5, col=rgb(0,0,0,15,maxColorValue=100))
    
    xreg<-par("usr")[1]
    yreg<-par("usr")[2]
    
    #points(xreg+xinch1)
    
    i<-identify(bp.df[,"xplot"], bp.df[,bp.sel], labels=bp.df[,"id"], n=1)
    ret.list <- NULL
    while(length(i) > 0){
        cell.i<-bp.df[i,"id"]
        dev.set(click.window)
        PeakFunc7(dat,cell.i,t.type="mp.1")
        dev.set(bp.win)
        #i<-identify(bp.df[,"xplot"], bp.df[,bp.sel], labels=bp.df[,"id"], n=1)
        i<-identify(bp.df[,"xplot"], bp.df[,bp.sel],labels="", n=1)
    }
    
    return(list(l.cell.types=l.cell.types,final.bp=final.bp, bp.df=bp.df))
    
    

}

# dat <- tmpRD
# cell <- NULL
# cells <- tmpRD$cell_types$neurons
# groups <- NULL
# dat.name <- NULL
# plot.new=T
# save.bp=F
# view.cells=F
# env=NULL
# statType = "minMax"
#' Interactive statistic maker. This creates a statistic based on peak heights or peak areas.
#' @export
bp.selector<-function(dat, cell=NULL, cells=NULL, groups = NULL, dat.name=NULL,plot.new=T,save.bp=F,view.cells=F, env=NULL, statType = "minMax"){
    #print(environment())
    if(is.null(env)){
        env<-.GlobalEnv
    }else{env<-env}
    
    if(is.null(dat.name)){
        dat.name<-deparse(substitute(dat))
    }else{dat.name<-dat.name}
    
    #grab the RD name from the RD
    if(is.null(dat.name)){
        dat.name<-deparse(substitute(dat))
    }else{dat.name<-dat.name}
    
    #Make sure you have some type of cells
    if(is.null(cells)){
        cells<-dat$c.dat$id
    }else{cells<-cells}
    
    #Choose a cell to display fro selecting stats
    if(is.null(cell)){
        cell<-dat$c.dat[1,'id']
    }else{cell<-cell}
    
    # Group of cells to view on the density cell plotter
    if(!is.null(groups)){
        formals(density_ct_plotter)$overlay <- TRUE
        formals(density_ct_plotter)$cell_types <- groups

    }else{
        formals(density_ct_plotter)$overlay <- FALSE
    }

    ###################################################################
    #This region needs significant work to improve to all data aspects
    ###################################################################
    ## Selcet eith Area or Peak Height
    bringToTop(-1)
    cat("\nSelect the statistic type to make the new statistic.\n")
    sel <- c("Peak Height", "Area", "SNR", "WhereMax")
    type <- sel[menu(sel, title=)]
    if(type == "Peak Height" ){
        type <- ".max"
    }else if(type == "Area"){
        type <- ".tot"
    }else if(type == "SNR"){
        type <- ".snr"
    }else if(type == "WhereMax"){
        type <- ".wm"
    }
    
    #Find the window regions
    levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    #Find the middle region of the windows
    levs.mean<-sort(tapply(dat$t.dat[,"Time"], as.factor(dat$w.dat$wr1), mean))
    #clean up the levs
    levs<-setdiff(names(levs.mean),"")
    #not sre
    levs.mean<-levs.mean[levs]
    #regional asignment for window region labeling
    #ys<-rep(1.05*(max(dat$t.dat[,"X.1"])), length(levs))
    
    #Create a new plot
    if(plot.new){
        dev.new(width=14, height=8)
    }else{}
    #Define the layout of the window region
    layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))
    par(bg="gray90")
    
    peakfunc.window<-dev.cur()
    PeakFunc7(dat,cell, lmain="  ",bcex=1.5, info=F)

    ## Minmax vs Custom
    if(statType == 'custom'){
        formals(density_ct_plotter)$xlim_top <- 3
        formals(density_ct_plotter)$xlim_bottom <- 0
        formals(boxPlotter)$ylim <- c(0,3)
        
        #define the open window
        #plot the trace specified at the beigning
        title(expression("RED"* phantom("/BLUE")), col.main="red")
        title(expression(phantom("RED/")*"BLUE"),col.main="blue")
        title(expression(phantom("RED")*"/"*phantom("BLUE")),col.main="black")

        # add point to the plot to define buttons
        ys <- rep(par("usr")[3],length(levs))
        points(levs.mean, ys, pch=16, cex=2)
        
        ###Selecting Control Windows
        tryCatch(bringToTop(-1), error=function(e)NULL)
        cat("RED:: Choose one or more window regions for the numerator in the equations,\n\nAmplification-or-block = active.window / control.window\n\nCLICK LARGE BLACK DOTS to select\nClick stop in the top left.\n")
        flush.console()

        #Select windows to define numerator
        activewindows <- identify(x=levs.mean,y=ys,labels="X",plot=T, col="red", cex=1.5)
        #collect the names of what you have selected
        activewindows<- levs[activewindows]
        
        ###Selecting Active Windows
        tryCatch(bringToTop(-1), error=function(e)NULL)
        cat("BLUE:: Choose one or more window regions for the denominator in the equations,\n\n Amplification-or-block = active.window / control.window\n\nClick stop in the top left, and then STOP LOCATOR from the drop down\n")
        flush.console()
        
        #change focus back to the peakwindow for active window selection
        dev.set(peakfunc.window)
        controlwindows <- identify(x=levs.mean,y=ys,labels="X",plot=T, col="blue",cex=1.5)
        controlwindows<-levs[controlwindows]
        
        #now if there are multiple control windows selected, 
        if(length(controlwindows)>1){
            #create the name for scp collumn lookup
            controlmax<-paste(controlwindows, type, sep="")
            #add that name to scp, and do a row mean
            controlmaxmean<-data.frame(rowMeans(dat$scp[,controlmax,drop=F]))
        }else{
            controlmax<-paste(controlwindows, type, sep="")
            controlmaxmean<-dat$scp[,controlmax,drop=F]
        }
        #same as above!
        if(length(activewindows)>1){
            activemax<-paste(activewindows, type, sep="")
            activemaxmean<-data.frame(rowMeans(dat$scp[, activemax, drop=F]))
        }else{
            activemax<-paste(activewindows, type, sep="")
            activemaxmean<-dat$scp[,activemax, drop=F]
        }

        max_amp_mean<-activemaxmean/controlmaxmean
        max_amp_mean[,2]<-seq(from=1,to=dim(max_amp_mean)[1],by=1)

        # Calculate percent change and select for cells
        cat("\nWould you like to save this statistic to scp? \n")
        bringToTop(-1)
        sel <- c("yes","no")
        save_stat_op <- sel[menu(sel, title="Save Stat?")]
        if(save_stat_op == 'yes'){
            cat("Enter the name of the statistic to be added to scp \n")
            stat.name<-scan(n=1, what='character')
            dat$scp[stat.name]<-max_amp_mean
            assign(dat.name,dat, envir=env)
        }
    }else if(statType == 'minMax'){
        formals(density_ct_plotter)$xlim_top <- 1
        formals(density_ct_plotter)$xlim_bottom <- -1
        formals(boxPlotter)$ylim <- c(-1,1)

        title("(After-Before)/(After+Before)")
    
        # add point to the plot to define buttons
        ys<-rep(par("usr")[3],length(levs))

        points(levs.mean, ys, pch=16, cex=2)
        #label each point with levs text
        #text(levs.mean,ys,labels=names(levs.mean),pos=c(1,3),cex=1, srt=90)
        
        ###Selecting Control Windows
        tryCatch(bringToTop(-1), error=function(e)NULL)
        cat("\nChoose the Pulse Following the compound of interest.\n\nThis is the AFTER pulse\n\nCLICK LARGE BLACK DOTS to select\nYou Only Get one click.\n\nCLICK ANY KEY TO CONTINUE\n
        ")
        flush.console()

        
        #scan(n=1)
        #Select windows to define numerator
        activewindows <- identify(x=levs.mean,y=ys,labels="X",plot=T, col="red", cex=1.5,n=1)
        #collect the names of what you have selected
        activewindows<- levs[activewindows]
        
        ###Selecting Active Windows
        tryCatch(bringToTop(-1), error=function(e)NULL)
        cat("\n###############################################\nChoose the Pulse Before the compound of interest.\n\nThis is the BEFORE pulse\nYou only get one click.\n\nPress ENTER to continue\n")
        flush.console()

        #scan(n=1)
        #change focus back to the peakwindow for active window selection
        dev.set(peakfunc.window)
        controlwindows <- identify(x=levs.mean,y=ys,labels="X",plot=T, col="blue",cex=1.5,n=1)
        controlwindows<-levs[controlwindows]
        
        #Find the scp collumn to provide the best stat
        aftermax<-paste(activewindows, type, sep="")
        aftermaxmean<-dat$scp[,aftermax, drop=F]
        
        beforemax<-paste(controlwindows, type, sep="")
        beforemaxmean<-dat$scp[,beforemax, drop=F]
        
        max_amp_mean<-(aftermaxmean-beforemaxmean)/(aftermaxmean+beforemaxmean)
        max_amp_mean[,2] <- seq(from=1,to=dim(max_amp_mean)[1],by=1)
        
        # Calculate percent change and select for cells
        # Calculate percent change and select for cells
        cat("\nWould you like to save this statistic to scp? \n")
        bringToTop(-1)
        sel <- c("yes","no")
        save_stat_op<- sel[menu(sel, title="Save Stat?")]
        if(save_stat_op == 'yes'){
            cat("Enter the name of the statistic to be added to scp \n")
            stat.name<-scan(n=1, what='character')
            dat$scp[stat.name]<-max_amp_mean
            assign(dat.name,dat, envir=env)
        }

    }        

    # Visualize the data!
    density_ct_plotter(
        dat, 
        cells,  
        stat = max_amp_mean[,1,drop = F],
        overlay=T,
        dense_sep=F,
        plot_new=F)
    
    par(new = TRUE) 
    
    boxPlotter(
        max_amp_mean[cells,], 
        activewindows = activewindows, 
        controlwindows = controlwindows)

    #170131 adding 2 point localization
    tryCatch({
        bringToTop(-1)
        sel <- c("yes", "no")
        localize_log <- sel[menu(sel, title = "Would you like to localize your boxplot?")]
        if( length(localize_log) == 0 ){
            localize_log <- "F"
        }else{ 
            if(localize_log != "yes"){
                localize_log<-"F"
            }else{
                localize_log <- "T"
            }
        }

        localize <- as.logical(localize_log)

        if(localize){
            cat("\nTo localize the boxplot\n1:: Everything above click will be selected\n2:: Select the bottom range then the top range\n\n")
            bringToTop(-1)
            sel <- c("1", "2")
            selector<- sel[menu(sel)]
            
            if(selector=="1"){loc<-locator(n=1, type="p", pch=15, col="red")}
            if(selector=="2"){loc<-locator(n=2, type="p", pch=15, col="red")}

            abline(v=loc$x,col="red")
            
            if(length(loc$x)==1){
                keepLogic <- max_amp_mean[cells, 1] > loc$x[1]
            }
            if(length(loc$x)==2){
                keepLogic <- max_amp_mean[cells, 1] > loc$x[1] & max_amp_mean[cells, 1] < loc$x[2]
            }
            
            subsetMat <- max_amp_mean[cells, 1, drop=F][keepLogic,,drop=F]
            x.names <- row.names(subsetMat[order(subsetMat[,1], decreasing=T),,drop=F])
        }else{
            print('hihihi')
            x.names <- row.names(max_amp_mean[cells,][order(max_amp_mean[cells,1],decreasing=T),])
        }
        }, error=function(e) {x.names <<- row.names(max_amp_mean[order(max_amp_mean[,1],decreasing=T),])} 
    )
    if(view.cells){
        bringToTop(-1)
        sel <- c("Yes", "No")
        continue<-sel[menu(sel, title="View Selected Cells?")]
    }else{continue<-"No"}
    
    if(continue=="Yes"){
        real.cells <- tcd(dat, x.names,dat.name=dat.name, track = F)
        return(real.cells)
    }else{
        return(x.names)
    }
}

#' Function allows for selection and deselection of cells to build stacked traces
#' @export
XYtrace <- function(dat, cell=NULL, img=NULL, cols=NULL, labs=F, y.var=T){
    graphics.off()
    dat.name<-deparse(substitute(dat))
    x.coor<-grep("\\.x",names(dat$c.dat), value=T, ignore.case=T)
    if(length(x.coor)>1){x.coor<-"center.x"}
    y.coor<-grep("\\.y",names(dat$c.dat), value=T, ignore.case=T)
    if(length(x.coor)>1){y.coor<-"center.y"}
    area<-grep("area",names(dat$c.dat), value=T, ignore.case=T)
    
    lab1<-grep("cgrp",names(dat$c.dat), value=T, ignore.case=T)
    if(length(lab1)==0){lab1<-grep("gfp.1",names(dat$c.dat), value=T, ignore.case=T)}
    
    lab1.1<-grep("cgrp",names(dat$c.dat), value=T, ignore.case=T)
    if(length(lab1.1)==0){lab1.1<-grep("gfp.2",names(dat$c.dat), value=T, ignore.case=T)}
    
    lab2<-grep("ib4",names(dat$c.dat), value=T, ignore.case=T)
    if(length(lab2)==0){lab2<-grep("tritc",names(dat$c.dat), value=T, ignore.case=T)}
    
    if(is.null(cell)){cell<-row.names(dat$c.dat)}
    else{cell<-cell}
    cell.coor<-dat$c.dat[cell,c(x.coor, y.coor)]

    
    
    # select the names of the collumns containing coordinates
    levs <- unique(dat$w.dat[,"wr1"])
    levs<-setdiff(levs, "")
    if(labs==TRUE){
    if(is.null(cols)){cols="orangered1"} else{cols=cols}}
    pch=16
    
    dev.new(height=4,width=12)
    dev.new(width=8, height=8)
    dev.new(height=8,width=12)
    lmain<-"XY ROI"
    
    
    if(is.null(img)){img<-dat$img1}
    img.dim.x<-dim(img)[1]	
    img.dim.y<-dim(img)[2]
    dev.set(dev.list()[2])
    par(mar=c(0,0,0,0))
    plot(0, 0, xlim=c(0,img.dim.x),ylim=c(img.dim.y,0),xaxs="i", yaxs="i",col=cols,pch=".")
    if(!is.null(img)){rasterImage(img, 0, img.dim.y, img.dim.x, 0);points(cell.coor[,1],cell.coor[,2],col=cols,pch=0)}
    else{
    points(cell.coor[,1],cell.coor[,2], col=cols, cex=dat$c.dat[,area]/200)
    points(cell.coor[,1],cell.coor[,2],col=cols, pch=4)}

    i <- identify(cell.coor[,1],cell.coor[,2],n=1,plot=F, col=NA, tolerance=0.05)
    i.names<-row.names(dat$c.dat[cell,])[i]
    while(length(i) > 0)
    {	#selected name of cell
        s.names <- row.names(dat$c.dat[cell,])[i]
        dev.set(dev.list()[1])
        if(y.var){PeakFunc6(dat,s.names, Plotit.both=F)}
        else{PeakFunc5(dat,s.names, Plotit.both=T)}

        dev.set(dev.list()[2])
        # If a cell is selected, that has already been selected, 
        # then remove that cell from the list
        if(length(intersect(i.names,s.names))==1){
            i.names<-setdiff(i.names,s.names)
            points(cell.coor[s.names,1],cell.coor[s.names,2],col="gray70",pch=0,cex=2.4)
            points(cell.coor[i.names,1],cell.coor[i.names,2],col="red",pch=0,cex=2.4)	
        }
        # If it han't been selected, then add it to the list
        else{i.names<-union(i.names,s.names)
        points(cell.coor[i.names,1],cell.coor[i.names,2],col="red",pch=0,cex=2.4)}
        
        if(length(i.names)>=1){
            dev.set(dev.list()[3])
            LinesEvery.5(dat,m.names=i.names, plot.new=F,img="img1", cols="black", dat.n=dat.name)}				
            dev.set(dev.list()[2])
            i <- identify(cell.coor[,1],cell.coor[,2],labels=dat$c.dat[cell,1],n=1,plot=T, pch=0,col="white", tolerance=0.05)
        }
    dev.off()
    graphics.off()
    return(row.names(dat$c.dat[i.names,]))
       
}

#' More advanced image clickers
#' @export 
XYtrace.2<-function(dat, cells=NULL, img=NULL, cols=NULL, zoom=T, labs=T, yvar=F, zf=40, t.type=NULL, sf=1, plot.labs=T){
    dat.name<-deparse(substitute(dat))
    print(class(cells))
    if(is.null(t.type)){t.type<-select.list(names(dat),title="Select a Trace")}
    #setup first windows for analysis and give each of them names
    dev.new(width=8, height=8)
    pic.window<-dev.cur()
        
    #plot image in the window
    if(is.null(cells)){cells<-dat$c.dat$id
    }else{cells<-cells}

    #if(is.null(img)){img<-dat$img1}
    if(is.null(img)){
        img.name<-image.selector(dat)
        img<-dat[[img.name]]
    }
    if(is.null(cols)){cols<-cols}

    img.dim.y<-dim(img)[1]
    img.dim.x<-dim(img)[2]	
    dev.set(which=pic.window)
    par(mar=c(0,0,0,0))
    plot(0, 0, xlim=c(0,img.dim.x),ylim=c(img.dim.y,0),xaxs="i", yaxs="i",col=cols,pch=".")
    rasterImage(img, 0, img.dim.y, img.dim.x, 0)

    if(zoom){
        zoom<-select.list(c("Manual", "Regional"), title="Zoom?  Cancel=NO")
        
        if(zoom=="Manual"){
            #Select regions to zoom on
            print("select X region first, then Y Region")
            x.sel<-locator(n=2, type="p", col="Red")$x
            y.sel<-locator(n=2, type="p", col="Red")$y

            rect(x.sel[1],y.sel[2],x.sel[2],y.sel[1], border="red")

            # before moving on, lets shrink won the image bya factor of 1/2 to have a preview image
            # to refer to
            dev.new(width=4, height=4)
            pic.window.2<-dev.cur()
            par(mar=c(0,0,0,0))
            plot(0, 0, xlim=c(0,img.dim.x),ylim=c(img.dim.y,0),xaxs="i", yaxs="i",col=cols,pch=".")
            if(!is.null(img)){
                rasterImage(img, 0, img.dim.y, img.dim.x, 0)
            }
            rect(x.sel[1],y.sel[2],x.sel[2],y.sel[1], border="red")

            # now i need to clsoe the window and open a new one with the same type of selection
            x.size<-abs(x.sel[1]-x.sel[2])
            y.size<-abs(y.sel[1]-y.sel[2])

            #if you want to mainatin the same aspect ratio
            #width vs height ratio
            x.plot.size<-8*(x.size/img.dim.x)
            y.plot.size<-8*(y.size/img.dim.y)

            #if you want to double the aspect ratio
            #width vs height ratio
            x.plot.size<-16*(x.size/img.dim.x)
            y.plot.size<-16*(y.size/img.dim.y)

            #plot the new image
            dev.off(which=pic.window)
            dev.new(width=x.plot.size, height=y.plot.size)
            pic.window<-dev.cur()

            par(mar=c(0,0,0,0))
            plot(0, 0, xlim=c(x.sel[1],x.sel[2]),ylim=c(y.sel[2],y.sel[1]),xaxs="i", yaxs="i",pch=".")
            rasterImage(img[y.sel[1]:y.sel[2],x.sel[1]:x.sel[2], ], x.sel[1], y.sel[2], x.sel[2], y.sel[1])
        }
        if(zoom=="Regional"){

            rect(0,img.dim.y/2, img.dim.x/2, 0, border="blue",lwd=3)
            rect(img.dim.x/2, img.dim.y/2, img.dim.x, 0, border="red", lwd=3)
            rect(0, img.dim.y, img.dim.x/2, img.dim.y/2, border="green", lwd=3)
            rect(img.dim.x/2, img.dim.y, img.dim.x, img.dim.y/2, border="purple", lwd=3)
            rect(img.dim.x*1/4, img.dim.y*3/4, img.dim.x*3/4, img.dim.y*1/4, border="navy", lwd=3)
            rect(img.dim.x*6/16, img.dim.y*10/16, img.dim.x*10/16, img.dim.y*6/16, border="red", lwd=3)

            text.place.x<-c(.02, .52, .02, .52, .27,.395)
            text.place.x<-text.place.x*img.dim.x
            text.place.y<-c(.02, .02, .52, .52, .27,.395)
            text.place.y<-text.place.y*img.dim.y
            
            #text.y<-img.dim.y*round(text.place$y/img.dim.y, digits=2)
            #text.x<-img.dim.x*round(text.place$x/img.dim.x, digits=2)
            text(text.place.x, text.place.y, c(1,2,3,4,5,6), col=c("blue", "red", "green", "purple","navy","red"), cex=3)
            
            region.selection<-as.numeric(select.list(as.character(c(1,2,3,4,5,6))))

            if(region.selection==1){
                dev.set(which=pic.window)
                par(mar=c(0,0,0,0))
                plot(0, 0, 
                    xlim=c(0, img.dim.x/2),
                    ylim=c(img.dim.y/2,0),xaxs="i", yaxs="i",col=cols,pch="."
                )
                rasterImage(img, 0, img.dim.y, img.dim.x, 0)
            }
            
            if(region.selection==2){
                dev.set(which=pic.window)
                par(mar=c(0,0,0,0))
                plot(0, 0, 
                    xlim=c(img.dim.x/2, img.dim.x),
                    ylim=c(img.dim.y/2,0),xaxs="i", yaxs="i",col=cols,pch="."
                )
                rasterImage(img, 0, img.dim.y, img.dim.x, 0)
            }

            if(region.selection==3){
                dev.set(which=pic.window)
                par(mar=c(0,0,0,0))
                plot(0, 0, 
                    xlim=c(0, img.dim.x/2),
                    ylim=c(img.dim.y/2,img.dim.y),xaxs="i", yaxs="i",col=cols,pch="."
                )
                rasterImage(img, 0, img.dim.y, img.dim.x, 0)
            }
            
            if(region.selection==4){
                dev.set(which=pic.window)
                par(mar=c(0,0,0,0))
                plot(0, 0, 
                    xlim=c(img.dim.x/2, img.dim.x),
                    ylim=c(img.dim.y/2,img.dim.y),xaxs="i", yaxs="i",col=cols,pch="."
                )
                rasterImage(img, 0, img.dim.y, img.dim.x, 0)
                #rasterImage(
                #	img[img.dim.y/2:img.dim.y,img.dim.x/2:img.dim.x,],
                #	img.dim.x/2, img.dim.y, img.dim.x, img.dim.y/2)
            }

            if(region.selection==5){
                dev.set(which=pic.window)
                par(mar=c(0,0,0,0))
                plot(0, 0, 
                    xlim=c(img.dim.x*1/4, img.dim.x*3/4),
                    ylim=c(img.dim.y*3/4,img.dim.y*1/4),xaxs="i", yaxs="i",col=cols,pch="."
                )
                rasterImage(img, 0, img.dim.y, img.dim.x, 0)
            }

            if(region.selection==6){
                dev.set(which=pic.window)
                par(mar=c(0,0,0,0))
                plot(0, 0, 
                    xlim=c(img.dim.x*6/16, img.dim.x*10/16),
                    ylim=c(img.dim.y*10/16,img.dim.y*6/16),xaxs="i", yaxs="i",col=cols,pch="."
                )
                rasterImage(img, 0, img.dim.y, img.dim.x, 0)
            }

        }
    }


    #Define the collumn names
    x.coor<-grep("\\.x",names(dat$c.dat), value=T, ignore.case=T)
    if(length(x.coor)>1){x.coor<-"center.x"}

    y.coor<-grep("\\.y",names(dat$c.dat), value=T, ignore.case=T)
    if(length(y.coor)>1){y.coor<-"center.y"}

    area<-grep("area",names(dat$c.dat), value=T, ignore.case=T)
    if(length(area)>1){area<-"area"}

    #Interactive Plot 
    dev.new(height=4,width=12)
    trace.window<-dev.cur()

    dev.new(height=8,width=12)
    lines.window<-dev.cur()

    cell.coor<-dat$c.dat[cells,c(x.coor, y.coor)]

    dev.set(which=pic.window)
    if(labs){points(cell.coor[,1],cell.coor[,2],col="gold", pch=4, cex=.1)}

    i <- identify(cell.coor[,1],cell.coor[,2],n=1,plot=F, col=NA, tolerance=0.1)
    i.names<-row.names(dat$c.dat[cells,])[i]

    while(length(i) > 0)
    {	#selected name of cell
            s.names <- row.names(dat$c.dat[cells,])[i]
            dev.set(which=trace.window)
            if(yvar){PeakFunc7(dat,s.names, yvar=F, zf=zf, t.type=t.type,dat.n=dat.name)}
            else{PeakFunc7(dat,s.names, yvar=F, zf=zf, t.type=t.type,dat.n=dat.name)}

            dev.set(which=pic.window)
            # If a cell is selected, that has already been selected, 
            # then remove that cell from the list
            if(length(intersect(i.names,s.names))==1){
                i.names<-setdiff(i.names,s.names)
                #points(cell.coor[s.names,1],cell.coor[s.names,2],col="gray70",pch=0,cex=2.4)
                #points(cell.coor[i.names,1],cell.coor[i.names,2],col="red",pch=0,cex=2.4)	
            }
            # If it han't been selected, then add it to the list
            else{i.names<-union(i.names,s.names)
            #points(cell.coor[i.names,1],cell.coor[i.names,2],col="red",pch=0,cex=2.4)
            }
            
            if(length(i.names)>=1){
                dev.set(which=lines.window)
                LinesEvery.5(dat,m.names=i.names, plot.new=F, img=c("img1", "img2", "img6","img7"), cols="black",sf=sf, t.type=t.type)}				
                dev.set(which=pic.window)
                #i <- identify(cell.coor[,1],cell.coor[,2],n=1,plot=F,col="white", tolerance=0.05)
                i <- identify(cell.coor[,1],cell.coor[,2],labels=dat$c.dat[cells,1],n=1,plot=T, pch=0,col="white", tolerance=0.05, cex=.5)
            }
    dev.off()
    graphics.off()
    return(row.names(dat$c.dat[i.names,]))		   
}

