
readkeygraph <- function(prompt){
    getGraphicsEvent(prompt = prompt, 
                 onMouseDown = NULL, onMouseMove = NULL,
                 onMouseUp = NULL, onKeybd = onKeybd,
                 consolePrompt = "uh")
    Sys.sleep(0.01)
    return(keyPressed)
}

onKeybd <- function(key){
    keyPressed <<- key
}

# Fucntion to sample it but maintain current order
mySample <- function(values,N){
        size <- length(values)

        values[sapply(1:size, function(i){
                    select <- as.logical(rbinom(1,1,N/(size+1-i)))
                    if(select) N <<- N - 1
                    select
        })]
}

#' TraceClick.d
#' ev(tcd)
#' 
#' This function is an interactive dashboard for viewing your exeriments.
#'
#' The output of this function returns a list of vectors of cell names
#' There are a few ways to input cell names into this program;
#' 
#' 1)Character; ex. cells=c('X.1','X.2','X.3','X.4')
#' 2)Numeric; ex. cells=c(1,2,3,4)
#' 3)Character Lists; ex. active.cells[[1]]
#' 
#' Character lists/Cell groups, can be handled and displayed in a variety 
#' of ways. Using Keyboard commands (CASE SENSITVE); 
#' 
#' @param up arrow move through list specified in entry
#' @param down arrow Move down through list spcified in entry
#' @param a place cell into the correct group in cell_types
#' @param c add cells to g.names
#' @param r reset g.names
#' @param 1-0  add cells to g.names1 through g.name10
#' @param q Quits the program
#' @param c add cells to g.names
#' @param s stack g.names
#' @param d details for peakfunc	
#' @param D LinesEvery seperation	
#' @param f New trace fitting for pottassium pulses
#' @param F New smoothing factor for fit trace
#' @param i Select image to display on Stacked Traces
#' @param I image for Multiview
#' @param l choose window region to display on stack trace plot
#' @param o order all cells in a new way
#' @param O order cells in Stacked Traces and multiview
#' @param p Toggles points on graph
#' @param P Pick a group/cells to click through
#' @param R reset group specified
#' @param r rename group names
#' @param s stack selected Groups
#' @param t brings up list of RD file. Select Trace (anything starting with t or mp)
#' @param u Underlines the Trace
#' @param v Show where cells are located and give zoomed in view
#' @param V choose cell info to display on traces
#' @param w Change Line Width on plot
#' @param x score this cell as a drop
#' @param X score cell as not dropped.
#' @param y Zoom yaxis automatically
#' @param z image zoom
#' @param F1 advanced stat comparison statistic creator
#' @param F2 advanced stat peak comparison min max norm creator
#' @param F3 Density plot visualization
#' @param F5 does something
#' @export
tcd<-function(dat, cells=NULL,img=dat$img1, l.img=c("img1"), yvar=FALSE, t.type="t.dat", plot.new=F, info=T, pts=T, lns=T, bcex=1, levs=NULL, klevs=NULL, sft=NULL, underline=T, zf=20, lw=2, sf=1, dat.name=NULL, view_func_description=F, save_question = T){
    graphics.off()
    print(environment())
    if(is.null(dat.name)){
        dat.name<-deparse(substitute(dat))
    }else{
        dat.name<-dat.name
    }
    if(view_func_description){
    cat(
    "
    #############################################
    Welcome to Trace.Click.dev
    #############################################

    The output of this function returns a list of vectors of cell names

    There are a few ways to input cell names into this program;

    1)Character; ex. cells=c('X.1','X.2','X.3','X.4')
    2)Numeric; ex. cells=c(1,2,3,4)
    3)Character Lists; ex. active.cells[[1]]

    Character lists/Cell groups, can be handled and displayed in a variety 
    of ways. Using Keyboard commands (CASE SENSITVE); 

    1)s: Stack group of cells 
    2)v: View images of cells
    3)P: Pick a group to scroll through with up and down arrows
        UP ARROW: move through list specified in entry
        DOWN ARROW: Move down through list spcified in entry
        o: reorders traces in the way specified.
    4)r: Rename your group use '.' and a space seperator ex. 'cool.cellz'
    5)R: Empty the specified group of all cells

    UP ARROW: move through list specified in entry
    DOWN ARROW: Move down through list spcified in entry

    ########################
    Stacked Traces Features

    u: Add or remove line under trace
    p: Add or removed points in single trace view
    t: Select the type of trace to display (anythin starting with a t or mp) 
    d: Remove  most information on the single trace view
    D: How much the traces are seperated, Must be greater than 0 ex. 0.2
    i: Image/Images to display on left side of traces
    V: 1.Choose Dataframe 2.Choose Values to display on right side of trace

    ####################
    Viewing cell images

    v: Select the group to view
    I: Change the image

    ##############################
    Making Groups

    1,2,3,4,5,6,7,8,9,0,-,+: add cells to g.names1 through g.name12
    shift+ (above value) removes cell from that group

    To clean up a group press P, select the group of interest
    press 'o' the sort the group in a specified way (ex area) 
    and then use shift + whatever key the cells are stored
    ex('1,2,3,4,5,6,7,8,9,0,-,+')


    q: Quits the program
    c: add cells to g.names
    s: stack g.names


    #d: details for peakfunc	
    #D: LinesEvery seperation	
    #f: New trace fitting for pottassium pulses
    #F: New smoothing factor for fit trace
    #i: Select image to display on Stacked Traces
    #I: image for Multiview
    #l: choose window region to display on stack trace plot
    #o: order all cells in a new way
    #O: order cells in Stacked Traces and multiview
    #p: Toggles points on graph
    #P: Pick a group/cells to click through
    #R: reset group specified
    #r: rename group names
    #s: stack selected Groups
    #t: brings up list of RD file. Select Trace (anything starting with t or mp)
    #u: Underlines the Trace
    #v: Show where cells are located and give zoomed in view
    #V: choose cell info to display on traces
    #w: Change Line Width on plot
    #x: score this cell as a drop
    #X: score cell as not dropped.
    #y: Zoom yaxis automatically
    #z: image zoom
    ")
    }else{}

    dat.tmp<-dat
    if(plot.new){graphics.off()}
    if(is.null(sft)){sft<-7}
    
    tryCatch(windows(width=14,height=4,xpos=0, ypos=50), error=function(e) windows(width=14,height=4))
    click.window<-dev.cur()
    
    # tryCatch(windows(width=10,height=6,xpos=0, ypos=450), error=function(e) windows(width=14,height=4))
    # lines.window<-dev.cur()
    
    # dimx<-dim(img)[2]
    # dimy<-dim(img)[1]
    # haight<-10*dimy/dimx
    # tryCatch(windows(width=haight*dimx/dimy, height=haight,xpos=1130, ypos=200), error=function(e) windows(width=haight*dimx/dimy, height=haight))
    # view.window<-dev.cur()
    
    # tryCatch(windows(width=8, height=8,xpos=1130, ypos=0), error=function(e) windows(width=8, height=8))
    # multipic.window<-dev.cur()
    
    # tryCatch(windows(width=12, height=2,xpos=0, ypos=550), error=function(e) windows(width=12, height=2))
    # traceimpute.window<-dev.cur()
    
    window.flag<-0
    lines.flag <- 0
    cell.i <- 1
    p.names<-NULL
    values<-"area"
    lines.color='black'
    
    #If no cell input collect all cells
    if(is.null(cells)){
    cells<-dat$c.dat$id
        cnames <- names(dat$c.dat$id)
        g.names<-cnames
    }else{}
    
    #If inputing a numeric vector, convert to character by adding a X. to beiging
    if(class(cells)=="numeric"){
        cells<-paste("X.", cells, sep="")
        cnames<-cells
        g.names<-cnames
    }

    #If inputing a list fill in 
    if(class(cells)=="list"){
        #Reduce g.names to combine all cells from the list into g.names
        g.names<-Reduce(union,cells)
        #initialize a list
        gt.names<-list()
        #Now fill in the list
        if( !is.null( names(cells) ) ){
            for(i in 1:length(cells)){
                #Fill in the gt.names with the names of the cells
                gt.names[[ names(cells)[i] ]]<-cells[[i]]
                #assign(names(cells)[i],cells[[i]])
            }
        }else{
            for(i in 1:length(cells)){
                #Fill in the gt.names with the names of the cells
                gt.names[[ paste0("g.names",i) ]]<-cells[[i]]
                #assign(names(cells)[i],cells[[i]])
            }
        }

        #if the length of the cell list is less than 12, fill in the remaining
        #list entries with empty regions
        if(length(gt.names)<12){
            for(i in ( length(gt.names)+1 ):12){
                #fill in with an NA
                gt.names[[paste("g.names",i,sep="")]]<-NA
                #remove the NA to allow for 
                gt.names<-lapply(gt.names, function(x) x[!is.na(x)])
            }
        }
        cells<-dat$c.dat$id
        cnames<-cells
        #gt.names<-list(g.names1=g.names1, g.names2=g.names2, g.names3=g.names3, g.names4=g.names4, g.names5=g.names5, g.names6=g.names6, g.names7=g.names7, g.names8=g.names8,g.names9=g.names9, g.names10=g.names10, g.names11=g.names11, g.names12=g.names12, g.names=g.names)
    }else{
        cnames<-cells
        g.names<-cnames
        g.names1<-NA
        g.names2<-NA
        g.names3<-NA
        g.names4<-NA
        g.names5<-NA
        g.names6<-NA
        g.names7<-NA
        g.names8<-NA
        g.names9<-NA
        g.names10<-NA
        g.names11<-NA
        g.names12<-NA

        gt.names<-list(g.names1=g.names1, g.names2=g.names2, g.names3=g.names3, g.names4=g.names4, g.names5=g.names5, g.names6=g.names6, g.names7=g.names7, g.names8=g.names8,g.names9=g.names9, g.names10=g.names10, g.names11=g.names11, g.names12=g.names12, g.names=g.names)
        gt.names<-lapply(gt.names, function(x) x[!is.na(x)])
        
        cells<-cells
        cnames<-cells
    }
    
    
    keyPressed <- "z"
    #group.names<-NULL
    if(is.null(levs)){
        levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    }else{
        levs<-levs
    }
    
    if(is.null(klevs)){
        klevs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    }else{
        klevs<-levs
    }
    
    while(keyPressed!="q"){
        cell.pick <- cnames[cell.i]
        
        dev.set(which=click.window)
        p1 <- PeakFunc7(dat,cell.pick, t.type=t.type, yvar=yvar, info=info, bcex=bcex, pts=pts, lns=lns, levs=levs, underline=underline, dat.n=dat.name, zf=zf)
        p1.par<-par()
        
        ##LinesEvery
        if(lines.flag == 1){
            if(length(p.names) < 500){
                if(length(p.names)>11){
                    tryCatch(dev.off(which=lines.window), error=function(e)NULL)
                    tryCatch(windows(width=10,height=12,xpos=0, ypos=100), error=function(e)windows(width=10,height=12))
                    lines.window<-dev.cur()
                }else{
                    tryCatch(dev.off(which=lines.window), error=function(e)NULL)
                    tryCatch(windows(width=10,height=7,xpos=0, ypos=250) , error=function(e)windows(width=10,height=7))
                    lines.window<-dev.cur()
                }
                dev.set(which=lines.window)
                tryCatch(LinesEvery.5(dat,p.names,plot.new=F, img=l.img,lmain=paste(gsub("[$]","",p.namez), 'n=',length(p.names)), t.type=t.type, lw=lw, col=lines.color, lns=lns, levs=levs, bcex=1, underline=underline, dat.n=dat.name, zf=zf, sf=sf, values=values),error=function(e) print("You haven't stacked traces yet, yo."))
                lines.flag <- 0
            }
        }
        
        if(lines.flag == 2){
            sample.to.display <- as.numeric(select.list(as.character(c(5,10,20,50,70,100))),title='Sample Number?')
            
            print(p.names)
            print(sample.to.display)
            if(sample.to.display < length(p.names)){
                if(sample.to.display > 20){
                    tryCatch(dev.off(which=lines.window.2), error=function(e)NULL)
                    tryCatch(windows(width=10,height=12,xpos=0, ypos=100), error=function(e)windows(width=10,height=12))
                    lines.window.2<-dev.cur()
                }else{
                    tryCatch(dev.off(which=lines.window.2), error=function(e)NULL)
                    tryCatch(windows(width=10,height=7,xpos=0, ypos=100), error=function(e)windows(width=10,height=12))
                    lines.window.2<-dev.cur()
                }

                tryCatch(
                    LinesEvery.5(
                        dat,
                        mySample(p.names, sample.to.display),
                        plot.new=F,
                        lmain=paste("Sample",sample.to.display,"out of",length(p.names)), 
                        img=l.img, lw=lw, t.type=t.type, col="black", lns=lns, levs=levs, bcex=1, underline=underline, dat.n=dat.name, zf=zf,sf=sf, values=values)
                ,error=function(e) print("You haven't stacked traces yet, yo."))
                lines.flag <- 0
            }else{
                lines.flag <- 1
            }
        }

        ##Pic zoom
        if(window.flag == 1){
            tryCatch(dev.off(which=view.window), error=function(e)NULL)
            tryCatch({
                dimx<-dim(img)[2]
                dimy<-dim(img)[1]
                haight<-10*dimy/dimx
                tryCatch(windows(width=haight*dimx/dimy, height=haight,xpos=1130, ypos=200), error=function(e) windows(width=haight*dimx/dimy, height=haight))
            }, error = function(e)NULL)
            view.window <- dev.cur()

            tryCatch(dev.off(which=mulitpic.window), error=function(e)NULL)
            tryCatch(windows(width=8, height=8,xpos=1130, ypos=0), error=function(e) windows(width=8, height=8))
            multipic.window<-dev.cur()
            
            if(length(p.names) < 500 ){
                dev.set(which=view.window)
                tryCatch(
                    cell.view(dat,
                        cell=p.names, 
                        img=img,
                        cols="Yellow",
                        plot.new=F,
                        cell.name=T, 
                        lmain=paste(gsub("[$]","",p.namez)), 
                        zoom=FALSE)
                    ,error=function(e) print("You haven't collected cells to view")
                )
                
                dev.set(which=multipic.window)
                tryCatch(
                    multi.pic.zoom(
                        dat,
                        p.names,
                        img, 
                        plot.new=F, 
                        zf=zf, 
                        labs=F)
                    ,error=function(e) print("You haven't collected cells to view")
                )
            }else{
                cat("\nThere are too many cells to view try again with less than 500 \n")
            }
            window.flag <- 0
        }
        
        #if(lines.flag==0){
            dev.set(which=click.window)
        #}		
        
        #title(sub=paste("Group ",group.i," n=",g.num," Cell ",cell.i,sep=""))
        ## How many cells are you looking at
        text(par("usr")[1], par("usr")[4]+yinch(.5),paste(cell.i, ":",length(cnames)))
        #click.i <- identify(x=xs,y=ys,n=1,plot=F)
        
        keyPressed <- readkeygraph("[press any key to continue]")
        
        if(keyPressed=="Up")
        {cell.i <- cell.i + 1;if(cell.i>length(cnames)){cell.i<-1};lines.flag<-0}
        
        if(keyPressed=="Down")
        {cell.i <- cell.i - 1;if(cell.i<1){cell.i<-length(cnames)};lines.flag<-0}
        
    #a: Assign this will place the cell of interest into the correct cell class
        if(keyPressed=='a'){
            #need to see if the dat has a cell_types
            cellToReassign <- cnames[cell.i]
            cat('\nReassigning cell', cellToReassign,'\n')
            #cellTypeId <- grep('^cell',names(dat), value=T)
            if(length(cellTypeId) > 0){
                #get only cell types that we want to reassign
                cellTypeNames <- names(dat[[cellTypeId]])
                cellTypesToNotClean <- c('neurons')
                cellTypesToClean <- setdiff(cellTypeNames, cellTypesToNotClean)
                
                #remove it from all groups
                dat[[cellTypeId]][cellTypesToClean] <- lapply(dat[[cellTypeId]][cellTypesToClean], function(X) setdiff(X,cellToReassign))
                
                ###now that we have indicated that we would like to place this cell into a new group
                #First lets find all groups we can assign to\
                tryCatch(bringToTop(-1), error=function(e)NULL)
                cat('\nWhich cell class does this actually belong to?\n')
                correctCellClass <- cellTypesToClean[menu(cellTypesToClean)]
                print(dat[[cellTypeId]][[correctCellClass]])
                dat[[cellTypeId]][[correctCellClass]] <- union(dat[[cellTypeId]][[correctCellClass]], cellToReassign)
                print(dat[[cellTypeId]][correctCellClass])
                assign(dat.name, dat, envir=.GlobalEnv)
            }else{
                cat('\nSorry You haven\'t defined cell types yet. Please do this first!\n')
            }
        }

    #c: add cells to g.names
        if(keyPressed=="c"){
			g.names<-union(g.names,cnames[cell.i]);print(g.names)}
    #C: Remove cells from g.names
        if(keyPressed=="C"){
			g.names<-setdiff(g.names,cnames[cell.i]);print(g.names)}
    #d: details for peakfunc	
        if(keyPressed=="d"){
            if(info){info=F}else{info=T}
            lines.flag<-1
        }
    #D: LinesEvery seperation	
        if(keyPressed=="D"){
            tryCatch(bringToTop(-1), error=function(e)NULL)
            print("change the zoom factor")
            print(paste("This is the current zoom",sf))
            sf<-scan(n=1)
            if(sf==0){sf<-.001}
            lines.flag<-1
        }
    # #f: New trace fitting for pottassium pulses
    #     if(keyPressed=="f"){
    #         lines.flag<-3
    #     }
    #f: Fixes the scoring
    #This function will allow you to fix the region you click on
        # if(keyPressed === 'f'){
            # #first find the middle region of each levs to correct the scoring
            # levsMiddle <- tapply(levsMiddle <- tapply(dat$w.dat[,1], as.factor(dat$w.dat$wr1),mean)[levs]
            # yLocation <- rep(par('usr')[4] + yinch(.5), length(levsMiddle))
            # par(xpd=T)

            # points(levsMiddle, yLocation, pch=7))
            # identify(levsMiddle, yLocation, n=1)



            # par(xpd=F)
            # }


        # }



    #F: New smoothing factor for fit trace
        if(keyPressed=="F"){
            print("Change the loess smoothing factor")
            print(paste("This is the current smoothing",sft))
            sft<-scan(n=1)
            lines.flag<-3
        }
    #h: Change the hue/color of the traces
        if(keyPressed=="h"){
            lines.color<-select.list(c('rainbow','black','brew.pal','topo'))
            if(lines.color==''){
                lines.color<-'black'
            }
            lines.flag<-1
        }
    #i: Select image to display on Stacked Traces
        if(keyPressed=="i"){
            l.img<-image.selector(dat)
            lines.flag<-1
        }
    #I: image for Multiview
        if(keyPressed=="I"){
            img<-dat[[image.selector(dat, multi=F)]]
            #lines.flag<-1
            window.flag<-1
        }
    #l: choose window region to display on stack trace plot
        if(keyPressed=="l"){
            #if(lns){lns<-FALSE}else{lns<-TRUE}
            levs<-select.list(setdiff(unique(as.character(dat$w.dat[,"wr1"])),""), multiple=T)
            if( (levs=="") || identical(levs,character(0)) ){levs<-NULL}#levs<-setdiff(unique(as.character(dat$w.dat$wr1)),"")}
            lines.flag<-1
        }
    #m: Move groups to another group
        if(keyPressed=="m"){
            tryCatch(bringToTop(-1), error=function(e)NULL)
            cat("
            Select the Group you would like to move
            
            ")
            
            gt_to_move<-select.list(names(gt.names), multiple=F)
            print(paste("You Selected Group ",gt_to_move))
            cat("
            Select the Target group to replace
            
            ")
            gt_to_replace<-select.list(names(gt.names), multiple=F)
            print(paste("Group ",gt_to_replace, "was replaced by ", gt_to_move))
            gt.names[[gt_to_replace]]<-gt.names[[gt_to_move]]
        }    

    #o: order all cells in a new way
        if(keyPressed=="o"){
        toMatch<-c("c.dat","bin","scp")
        order_dat<-grep(paste(toMatch,collapse="|"),names(dat),value=TRUE)

        
        datfram<-select.list(order_dat,title="Where is the data?")
        collumn<-select.list(names(dat[[datfram]]),title="Collumn to sort")
        
        tryCatch(cnames<-c.sort.2(dat[[datfram]],cnames,collumn=collumn),error=function(e) print("Something went wrong try again"))
        cell.i<-1
        }	
    #O: order cells in Stacked Traces and multiview
        if(keyPressed=="O"){
            tryCatch({
                    p.names <- c.sort.2(dat, p.names)
                    lines.flag <- 1
                    window.flag <- 1
                },error=function(e) print("You have not stacked traces yet.")
            )
        }	
    #p: Toggles points on graph
        if(keyPressed=="p"){
            if(pts){pts<-FALSE}else{pts<-TRUE}
            lines.flag<-1
        }
    #P: Pick a group/cells to click through
        if(keyPressed=="P"){
            tryCatch(bringToTop(-1), error=function(e)NULL)
            cat("\nPick a Group of cells or a single cell to observe \nIf you Click cancel, all cells will be returned\n")
            selection<-select.list(c("group","cells"))
            if(selection=="group"){
                gt.to.click<-select.list(names(gt.names), multiple=F)
                if( is.null(gt.names[[gt.to.click]]) | is.logical( gt.names[[gt.to.click]]) ){
                    tryCatch(bringToTop(-1), error=function(e)NULL)
                    cat("\nNothing is in this Group\n")
                }else{
                    cell.i<-1
                    print(gt.to.click)
                    cnames <- gt.names[[gt.to.click]]
                    tryCatch({
                        cnames <- c.sort.2(dat[[datfram]],cnames,collumn=collumn)
                    },
                    error=function(e) print("Something went wrong try again") )
                    p.names <- cnames
                    print(cnames)
                }
            }
            if(selection=="cells"){
                cell.i<-1
                cnames<-select.list(as.character(dat$c.dat$id), multiple=T)
                tryCatch({
                    cnames <- c.sort.2(dat[[datfram]],cnames,collumn=collumn)
                },error=function(e) print("Something went wrong try again"))
                p.names <- cnames
            }
            if(selection==""){
                cell.i<-1
                cnames<-dat$c.dat$id
            }
        }

    #R: reset group specified
        if(keyPressed=="R"){
            p.namez<-paste(select.list(names(gt.names)),sep="")
            if(p.namez!=""){
                print(p.namez)
                gt.names[[p.namez]]<-NA
                gt.names[[p.namez]]<-gt.names[[p.namez]][ !is.na(gt.names[[p.namez]]) ]
                #gt.names[[p.namez]]<-lapply(gt.names[[p.namez]], function(x) x[!is.na(x)])
                print(paste("You Emptied", p.namez))
            }else{}
        }
    #r: rename group names
        if(keyPressed=="r"){
            tryCatch(bringToTop(-1), error=function(e)NULL)
            print("Select a group to rename")
            gt.to.rename<-select.list(names(gt.names), multiple=F)
            name.number<-which(names(gt.names)==gt.to.rename,arr.ind=T)
            print("Type in the new name Cannot start with number, no spaces.")
            tryCatch(names(gt.names)[name.number]<-scan(n=1, what='character'),error=function(e) print("You did not enter a name, so this group was not renamed"))
            #assign(names(gt.names)[name.number],gt.names[[name.number]])
            #lines.flag<-1
        }

    #s: stack selected groups
        if(keyPressed=="s"){
            p.namez<-paste(select.list(names(gt.names)),sep="")
            p.names<-gt.names[[p.namez]]
            #p.names<-get(ls(pattern=p.namez))
            lines.flag<-1
        }
        
    #S: Sample selected groups
        if(keyPressed=="S"){
            # p.namez<-paste(select.list(names(gt.names)),sep="")
            # print(p.namez)
            # p.names<-gt.names[[p.namez]]
            # #p.names<-get(ls(pattern=p.namez))
            # print(p.names)
            lines.flag <- 2
        }
    #t: brings up list of RD file. Select Trace (anything starting with t or mp)
        if(keyPressed=="t"){
            toMatch<-c("t[.]","blc","snr","mp")
            trace_dat<-grep(paste(toMatch,collapse="|"),names(dat),value=TRUE)
            t.type1<-t.type
            t.type<-select.list(trace_dat)
            if(t.type==""){t.type<-t.type1}
            lines.flag<-1
        }
    #u: Underlines the Trace
        if(keyPressed=="u"){
            if(underline){underline=F}else{underline=T}
            lines.flag<-1
        }
    #v: Show where cells are located and give zoomed in view
        if(keyPressed=="v"){
            p.namez<-paste(select.list(names(gt.names)),sep="")
            print(p.namez)
            p.names<-gt.names[[p.namez]]
            print(p.names)
            window.flag<-1
        }
    #V: choose cell info to display on traces
        if(keyPressed=="V"){
            #if(lns){lns<-FALSE}else{lns<-TRUE}
            values<-select.list(names(dat$c.dat), multiple=T)
            lines.flag<-1
        }

    #w: Change Line Width on plot
        if(keyPressed=="w"){
            tryCatch(bringToTop(-1), error=function(e)NULL)
            print("change the line width (lw) for LinesEvery")
            print(paste("This is the current lw",lw))
            lw<-scan(n=1)
            lines.flag<-1
        }
    #x: Drop cell
        if(keyPressed=="x")
        {
            print(cnames[cell.i])
            dat$bin[cnames[cell.i], "drop"]<-1
            print(dat$bin[cnames[cell.i], "drop"])
            print(paste("You Dropped Cell",cnames[cell.i]))
            # now that you have dropped a cell, this need to be removed from
            # cell types
            cellTypeId <- grep('^cell',names(dat), value=T)
            if(length(cellTypeId) > 0){
                drops <- row.names(dat$bin[dat$bin$drop==1,])
				print(drops)
                dat[[cellTypeId]] <- lapply(dat[[cellTypeId]], function(X) setdiff(X,drops))
                assign(dat.name,dat, envir=.GlobalEnv)
            }else{assign(dat.name,dat, envir=.GlobalEnv)}
        }
    #X: undrop cell
        if(keyPressed=="X")
        {
            print(cnames[cell.i])
            dat$bin[cnames[cell.i], "drop"]<-0
            print(dat$bin[cnames[cell.i], "drop"])
            print(paste("You Dropped Cell",cnames[cell.i]))
        }
    #y: Zoom yaxis automatically
        if(keyPressed=="y"){
            if(yvar){yvar<-FALSE}else{yvar<-TRUE}
        }
    #z: image zoom
        if(keyPressed=="z"){
            tryCatch(bringToTop(-1), error=function(e)NULL)
            print("change the zoom factor")
            print(paste("This is the current zoom",zf))
            zf<-scan(n=1)
            lines.flag<-1
            window.flag<-1
        }
        
        #if(keyPressed=="k")
        #{
        #	dat$bin[cnames[cell.i], "drop"]<-0
        #	print(paste("You Dropped Cell",cnames[cell.i]))
        #}
        
        #F1: Simple bp.selector. Create the statistic labeled on the plot. The localize question
        #allows you to click the boxplot to select a subset of cells to observe
        if(keyPressed=="F1"){
            tryCatch({
                #first open a new window
                #after undergoing a logical test to see if it exists
                if(length(ls(pattern='bp.selector.window'))==0){
                    dev.new(width=14, height=8)
                    #give this window a name
                    bp.selector.window<-dev.cur()
                }
                #give the focus to the new window
                dev.set(bp.selector.window)
                #empty gt.names[[12]]
                gt.names[[12]]<-NA
                #remove the NA, which will be repalced with a logical(0)
                gt.names[[12]]<-lapply(gt.names[[12]], function(x) x[!is.na(x)])
                #do the function bp.selector to gather data
                tryCatch(bringToTop(-1), error=function(e)NULL)
                cat("##############################################################################\nStat Maker: CUSTOM\n##############################################################################\n\nThis function allows you to create statistics based on the statistic you select.\nThis Function finds a represention of peak amplification and or block \nThis function will take in what ever you are currently scrolling through\n\nYou have the option to localize your boxplot. This means, select cells\nspecifically based on where you click on the boxplot.\n\nTwo clicks means you need\nto specify the lower range followed by the upper range.\nOne click will take everything greater than your click\n\nThe Other option that will arise is, would you like the save the stat.\nIf you do, the console will prompt you to enter a name. Ensure no spaces in the name\nThe next option will be whether you would like to make another statistic.\n")
                dev.set(bp.selector.window)
                gt.names[[12]]<-bp.selector(dat,
                    cnames[cell.i],
                    cnames, 
                    groups = gt.names, 
                    plot.new=F,
                    dat.name=NULL,
                    env=environment(),
                    statType = 'custom')
                #Now fill TCD with the cells just selected.
                cnames<-gt.names[[12]]	
                cell.i<-1
                lines.flag<-1
                windows.flag<-1
            }, error = function(e) cat("\nDid not work. Review documentation\n")
            )
        }
        
        #F2: Advanced Statistic maker This function uses the function (After-Before)/(After+Before)
        #this function allows you to save the stat.  This will be added to the scp dataframe at the bottom.
        #if you have created statistics, be sure to save your RD file before you close
        if(keyPressed=="F2"){
            tryCatch({
                #first open a new window
                #after undergoing a logical test to see if it exists
                if(length(ls(pattern='bp.selector.window'))==0){
                    dev.new(width=14, height=8)
                    #give this window a name
                    bp.selector.window<-dev.cur()
                }
                #give the focus to the new window
                dev.set(bp.selector.window)
                #empty gt.names[[12]]
                gt.names[[12]]<-NA
                #remove the NA, which will be repalced with a logical(0)
                gt.names[[12]]<-lapply(gt.names[[12]], function(x) x[!is.na(x)])
                #do the function bp.selector to gather data
                tryCatch(bringToTop(-1), error=function(e)NULL)
                cat("##############################################################################\nStat Maker: MinMaxnorm\n##############################################################################\n\nThis function allows you to create statistics based on the statistic you select.\nThis Function finds a represention of peak amplification and or block\nThis function will take in what ever you are currently scrolling through\n\nYou have the option to localize your boxplot. This means, select cells\nspecifically based on where you click on the boxplot.\nTwo clicks means you need to specigy the lower range followed by the upper range.\nOne click will take everything greater than your click\nThe Other option that will arise is, 'would you like the save the stat?'\nIf you do, the console will prompt you to enter a name. Ensure no spaces in the name\nThe next option will be whether you would like to make another statistic."
                )
                dev.set(bp.selector.window)
                gt.names[[12]]<-bp.selector(dat,
                    cnames[cell.i],
                    cnames, 
                    groups = gt.names, 
                    plot.new = F,
                    dat.name=NULL,
                    env=environment(),
                    statType = 'minMax')
                #Now fill TCD with the cells just selected.
                cnames<-gt.names[[12]]	
                cell.i<-1
                lines.flag<-1
                windows.flag<-1
            }, error = function(e) cat("\nDid not work. Review documentation\n")
            )
        }
        

        #F3: Plotting the Density plots.  There are many options for this plot
        if(keyPressed=="F3"){
            if(length(ls(pattern="density_win"))==0){
                dev.new(width=10,height=10)
                density_win<-dev.cur()
            }else{
                dev.off(density_win)
                dev.new(width=10,height=10)
                density_win<-dev.cur()
            }
            tryCatch(bringToTop(-1), error=function(e)NULL)
            cat("What dataframe wil contain your stat? \n")
            dense_df_q<-select.list(names(dat))
            cat("What attribute would you like to see the distribution? \n")
            dense_df_att<-menu(names(dat[[dense_df_q]]))
            statz<-dat[[dense_df_q]][dense_df_att]

            #define the top xlim value
            cat("Define Top xlim value \n")
            cat("Enter n to allow default Max value \n")
            xlim_top<-scan(n=1, what = 'raw')
            if(xlim_top == 'n' ){
                xlim_top<-max(dat[[dense_df_q]][dense_df_att])
            }else{
                xlim_top <- as.numeric(xlim_top)
            }
            
            cat("Define bottom xlim value \n")
            cat("Enter n to allow default Max value \n")
            xlim_bottom<-scan(n=1, what = 'raw')
            if(xlim_bottom == 'n'){
                xlim_bottom<-min(dat[[dense_df_q]][dense_df_att])
            }else{
                xlim_bottom <- as.numeric(xlim_bottom)
            }
            
            cat("\nSeperate the density plots?")
            sel <- c('yes', 'no')
            sel <- sel[menu(sel)]

            if(sel == 'yes'){
                formals(density_ct_plotter)$dense_sep <- T
            }else if(sel == 'no'){
                formals(density_ct_plotter)$dense_sep <- F
            }

            dev.set(density_win)
            density_ct_plotter(
                dat,
                g.names,
                cell_types=NULL, 
                stat=statz,
                overlay=T, 
                plot_new=F,
                xlim_top=xlim_top,
                xlim_bottom=xlim_bottom,
                dat.name=dat.name)

            lines.flag<-1
        }
        
        # #F4: Utilizing Topview
        # if(keyPressed=="F4"){
        #     p.namez<-paste(select.list(names(gt.names)),sep="")
        #     p.names<-gt.names[[p.namez]]

        #     aux_var<-c('area')

        #     #What i need to do is selectively import gfp and tritc variables into the 
        #     #topview function
        #     #this means search in the bin data frame for ib4 and gfp
        #     add_vars <- grep('mcherry|cy5|gfp|drop', names(dat$bin),value=T)
        #     aux_var<-c(aux_var, add_vars)
        #     TopView(dat, p.names, 12, 6, dat_name=dat.name, aux.var=aux_var)
        # }
        
        if(keyPressed=="F4"){
            if(length(ls(pattern="density_win"))==0){
                dev.new(width=10,height=10)
                density_win<-dev.cur()
            }else{}
            tryCatch(bringToTop(-1), error=function(e)NULL)
            cat("What dataframe wil contain your stat? \n")
            dense_df_q<-select.list(names(dat))
            cat("What attribute would you like to see the distribution? \n")
            dense_df_att<-menu(names(dat[[dense_df_q]]))
            statz<-dat[[dense_df_q]][dense_df_att]

            #define the top xlim value
            cat("Define Top xlim value \n")
            cat("Enter n to allow default Max value \n")
            xlim_top<-scan(n=1, what = 'raw')
            if(xlim_top == 'n' ){
                xlim_top<-max(dat[[dense_df_q]][dense_df_att])
            }else{
                xlim_top <- as.numeric(xlim_top)
            }
            
            cat("Define bottom xlim value \n")
            cat("Enter n to allow default Max value \n")
            xlim_bottom<-scan(n=1, what = 'raw')
            if(xlim_bottom == 'n'){
                xlim_bottom<-min(dat[[dense_df_q]][dense_df_att])
            }else{
                xlim_bottom <- as.numeric(xlim_bottom)
            }
            
            dev.set(density_win)
            density_ct_plotter(dat, 
                cnames, 
                cell_types = gt.names, 
                stat=statz,
                overlay=T, 
                dense_sep=F,
                plot_new=F,
                xlim_top=xlim_top,
                xlim_bottom=xlim_bottom,
                dat.name=dat.name)

            lines.flag<-1
        }

        #F5: Censusus Viewer
        if(keyPressed=="F5"){
            cat("\nSelect a binary column to add to the 12th group\n")
			cnames_orig <- cnames
			tryCatch({
                cells_to_view <- census_viewer(dat)
                
                if( is.na(cells_to_view) ){
                    cnames <- cnames_orig
                    cat(
                    "\nThere were no cells in that selection\n"
                    )
                }else{ 
                    cell.i<-1
                    cnames <- cells_to_view$cells

                    oldName <- names(gt.names)[12]
                    gt.names[[12]] <- cells_to_view$cells
                    names(gt.names)[12] <- cells_to_view$name
                    p.namez <- cells_to_view$name
                    p.names <- gt.names[[12]]
                    cat("/nThe group ", oldName, ' has been replaced by ', names(gt.names)[12], "\n")
                    lines.flag<-1
                }
            },
                error= function(e){
                    cat('\nYou most likely do not have cell_types made\n')
                }
            )
		}   
        
        #F6: Censusus Viewer
        if(keyPressed=="F6"){
            cat("\n This is a temporary function to view the responses per cell class.\n")

			cnames_orig <- cnames
            cat("Please select the collumn you would like to view\n")
			cells_to_view <- cellzand_tcd(dat$bin)
			if( is.na(cells_to_view) ){
				cnames <- cnames_orig
				cat(
				"There were no cells in that selection"
				)
			}else{ 
				cell.i<-1
				cnames <- cells_to_view$cells
				gt.names[[12]] <- cells_to_view$cells
				names(gt.names)[12] <- cells_to_view$name
				if(length(cells_to_view$cells[1]) > 20 ){
                    p.namez <- cells_to_view$name
				    p.names <- gt.names[[12]]
                    lines.flag<-1
                }
			}
		}  

        #F7: Load cell Types into the groups to pick with 'P'
        if(keyPressed=='F7')  {
           cellTypeId <- grep('^cell',names(dat), value=T)
            if(length(cellTypeId)>0){
				if(length(cellTypeId)>1){
					tryCatch(bringToTop(-1), error=function(e)NULL)
					cat('\n Select the cell type to load in \n')
					cellTypeId <- select.list(cellTypeId, title="Select Cell Type")
				}
                tryCatch(bringToTop(-1), error=function(e)NULL)             
                cat("\nI have filled in your cell_types to choose by pressing \'P\' ENJOY!\n")
                flush.console()
                Sys.sleep(0.5)
                gt.names <- list()
                for(i in 1:length(dat[[cellTypeId]])){
                    #Fill in the gt.names with each cell type
                    gt.names[[ names(dat[[cellTypeId]][i]) ]]<-dat[[cellTypeId]][[i]]
                }
            }else{
                cat('\nSorry you haven\'t defined cell types yet, so i can\'t fill it it for you.\n')
            }
        }

        if(keyPressed=="1")
        {gt.names[[1]]<-union(gt.names[[1]],cnames[cell.i]);print(gt.names[1])}
        if(keyPressed=="!")
        {gt.names[[1]]<-setdiff(gt.names[[1]],cnames[cell.i]);print(gt.names[1])}

        if(keyPressed=="2")
        {gt.names[[2]]<-union(gt.names[[2]],cnames[cell.i]);print(gt.names[2])}
        if(keyPressed=="@")
        {gt.names[[2]]<-setdiff(gt.names[[2]],cnames[cell.i]);print(gt.names[2])}

        if(keyPressed=="3")
        {gt.names[[3]]<-union(gt.names[[3]],cnames[cell.i]);print(gt.names[3])}
        if(keyPressed=="#")
        {gt.names[[3]]<-setdiff(gt.names[[3]],cnames[cell.i]);print(gt.names[3])}

        if(keyPressed=="4")
        {gt.names[[4]]<-union(gt.names[[4]],cnames[cell.i]);print(gt.names[4])}
        if(keyPressed=="$")
        {gt.names[[4]]<-setdiff(gt.names[[4]],cnames[cell.i]);print(gt.names[4])}

        if(keyPressed=="5")
        {gt.names[[5]]<-union(gt.names[[5]],cnames[cell.i]);print(gt.names[5])}
        if(keyPressed=="%")
        {gt.names[[5]]<-setdiff(gt.names[[5]],cnames[cell.i]);print(gt.names[5])}

        if(keyPressed=="6")
        {gt.names[[6]]<-union(gt.names[[6]],cnames[cell.i]);print(gt.names[6])}
        if(keyPressed=="^")
        {gt.names[[6]]<-setdiff(gt.names[[6]],cnames[cell.i]);print(gt.names[6])}

        if(keyPressed=="7")
        {gt.names[[7]]<-union(gt.names[[7]],cnames[cell.i]);print(gt.names[7])}
        if(keyPressed=="&")
        {gt.names[[7]]<-setdiff(gt.names[[7]],cnames[cell.i]);print(gt.names[7])}

        if(keyPressed=="8")
        {gt.names[[8]]<-union(gt.names[[8]],cnames[cell.i]);print(gt.names[8])}
        if(keyPressed=="*")
        {gt.names[[8]]<-setdiff(gt.names[[8]],cnames[cell.i]);print(gt.names[8])}

        if(keyPressed=="9")
        {gt.names[[9]]<-union(gt.names[[9]],cnames[cell.i]);print(gt.names[9])}
        if(keyPressed=="(")
        {gt.names[[9]]<-setdiff(gt.names[[9]],cnames[cell.i]);print(gt.names[9])}
    
        if(keyPressed=="0")
        {gt.names[[10]]<-union(gt.names[[10]],cnames[cell.i]);print(gt.names[10])}
        if(keyPressed==")")
        {gt.names[[10]]<-setdiff(gt.names[[10]],cnames[cell.i]);print(gt.names[10])}

        if(keyPressed=="-")
        {gt.names[[11]]<-union(gt.names[[11]],cnames[cell.i]);print(gt.names[11])}
        if(keyPressed=="_")
        {gt.names[[11]]<-setdiff(gt.names[[11]],cnames[cell.i]);print(gt.names[11])}

        if(keyPressed=="=")
        {gt.names[[12]]<-union(gt.names[[12]],cnames[cell.i]);print(gt.names[12])}
        if(keyPressed=="+")
        {gt.names[[12]]<-setdiff(gt.names[[12]],cnames[cell.i]);print(gt.names[12])}

        BACKUP<<-gt.names 
        if(keyPressed=="q")
        {
            #graphics.off()
            tryCatch({
                dev.off(which=click.window)
                dev.off(which=lines.window)
                dev.off(which=lines.window.2)			
                dev.off(which=view.window)
                dev.off(which=multipic.window)
                dev.off(which=traceimpute.window)
            }, error=function(e) print("this windows hasn't been opened yet"))

        }
    }
    #rd.name <- as.character(substitute(dat))
    #print(rd.name)
    #assign(rd.name, dat, envir=.GlobalEnv)
    #gt.names<-list(g.names1=g.names1, g.names2=g.names2, g.names3=g.names3, g.names4=g.names4, g.names5=g.names5, g.names6=g.names6, g.names7=g.names7, g.names8=g.names8,g.names9=g.names9, g.names10=g.names10, g.names11=g.names11, g.names12=g.names12, g.names=g.names)
    BACKUP<<-gt.names 
    
    assign(dat.name,dat, envir=.GlobalEnv)
    tryCatch(bringToTop(-1), error=function(e)NULL)
    if(save_question){
        print('Would y ou like to save you cell groups?')
        selection<-select.list(c('no','yes'),title='Save Groups?')
        if(selection=='yes'){
            print("Write in your name")
            save.names <- scan(n=1, what='character')
            save_label <- save.names
            assign(save.names, gt.names, envir = .GlobalEnv)  
            assign(save.names , gt.names)
            save(list = save.names ,file=paste(save_label,'.Rdata',sep=''))
            gt.names<<-gt.names
        }else{
            gt.names<<-gt.names
            return(gt.names)
        }
    }else{
        gt.names<<-gt.names
        return(gt.names)
    }
    #print(rd.name)
}

### Function to select rows based on collumn parameters
# dat can be either a raw RD object or an RD dataframe
# ex dat -or- dat$bin

cellzand_tcd<-function(dat,collumn=NULL, parameter=1,cells=NULL){
    
    cells_to_view <- list()
    bob<-list()
    if(is.null(cells)){cells<-dat$c.dat$id}else{cells<-cells}
    if(class(dat)=="list"){
        dat.select<-select.list(names(dat), title="Select DataFrame")
        dat<-dat[[dat.select]]
        if(is.null(cells)){
            cells<-row.names(dat)}else{cells<-cells
        }
    }else{
        dat<-dat
        if(is.null(cells)){cells<-row.names(dat)}else{cells<-cells}
    }
    
    if(is.null(collumn)){
        collumn<-select.list(names(dat), multiple=T, title="Select Collumn")
        cells_to_view$name <- collumn
    }else{collumn<-collumn}
    
    if(is.null(parameter)){
        parameter<-1
    }else{parameter<-parameter}
    
    for(i in collumn){
        bob[[i]]<-row.names(dat)[dat[,i]>=parameter]
    }
    
    bob<-Reduce(union, bob)
    
    bob<-intersect(bob,cells)

    cells_to_view$cells <- bob
    if( length(bob) == 0){
        return(NA)
    }else{
        return(cells_to_view)
    }
}

