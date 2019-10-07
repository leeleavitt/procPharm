graphics.off()
l.col<-"red"

Trace_select_grid<-function(dat, x.names, levs=select.list(names(dat$bin)), t.type="blc", preselect=T, l.col="red", window.w=10, window.h=10, title1="hi"){
    
    x.names<-rev(x.names)

    #Now create 3 extra spaces for buttons
    xn <- length(x.names)
    num.grid <- xn+4
    #This is the number of grids for the rows
    nr <- floor(sqrt(num.grid))
    #this is the number of grids for the rows
    nc <- ceiling((num.grid)/nr)
    #this is the maximun value needed to aquire the matrix of interest
    mtx <- max(nr,nc)
    #this helps to find the center location of each cell
    dx <- seq(0,1,length.out=(mtx+1))[-1]
    #this defines the size between the cells
    sl <- (dx[2]-dx[1])/2
    #This relocates the cells to the far left
    dx <- dx-sl

    all.x <- as.vector(matrix(rep(dx,mtx),byrow=F,ncol=mtx))
    all.y <- as.vector(matrix(rep(dx,mtx),nrow=mtx,byrow=T))

    #Lees trace image plotter
	if(is.null(levs)){
		levs<-setdiff(unique(dat$w.dat$wr1),"")
	}else{levs<-levs}

    levs_min<-min(as.numeric(row.names(which(dat$w.dat["wr1"]==levs,arr.ind=T))))
    levs_max<-max(as.numeric(row.names(which(dat$w.dat["wr1"]==levs,arr.ind=T))))

    levs_min<-which(row.names(dat$blc)==as.character(levs_min))
    levs_max<-which(row.names(dat$blc)==as.character(levs_max))

    peak_min<-min(dat[[t.type]][levs_min:levs_max,dat$c.dat$id])
    peak_max<-max(dat[[t.type]][levs_min:levs_max,dat$c.dat$id])*1.4

	#now loop through the data and create png plots of each region
	png.name<-c()
	start.time<-Sys.time()
    for(i in 1:xn){
		png.name[i]<-paste("tmp_png_",i,".png", sep="")
        png(png.name[i], 30,30, res=10, bg="transparent")
        par(bty="n",mai=c(0,0,0,0))
        plot(dat[[t.type]][ levs_min:levs_max, x.names[i] ],type='l',lwd=3,xaxt='n',yaxt='n',col="white", ylim=c(0,peak_max))
        dev.off()
		#print(i)
	}	
	end_time<-Sys.time()
	
	print(paste("Elapsed time saving:",end_time-start.time))
	
	#now lets open up single view window
    dev.new(width=14,height=4,title="SingleCell")
    trace_view <- dev.cur()

    #Open the grid window
    dev.new(height=window.w,width=window.h,canvas="black",title=title1)
    grid_view <- dev.cur()
    op <- par(mar=c(0,0,0,0))	
    plot(c(0,1),c(0,1),xaxt="n",yaxt="n",type="n",ylab="",xlab="")	
	
	require(png)
	start.time<-Sys.time()
	for(i in 1:xn){
        tmp_img<-readPNG(png.name[i])
        dim(tmp_img)
        
        xl <- all.x[i]-sl*.9
        xr <- all.x[i]+sl*.9
        xt <- all.y[i]-sl*.9
        xb <- all.y[i]+sl*.9
        
        dev.set(grid_view)
        rasterImage(tmp_img,xl,xt,xr,xb)
        unlink(png.name[i])
    }
	end.time<-Sys.time()
	print(paste("Elapsed plot time", end.time-start.time))

    cexr <- sl/.05
	
    text(all.x[xn+1],all.y[xn+1],"Done",col="white",cex= cexr)
    text(all.x[xn+2],all.y[xn+2],"All",col="white",cex= cexr)
    text(all.x[xn+3],all.y[xn+3],"None",col="white",cex= cexr)
    text(all.x[xn+4],all.y[xn+4],"Reset",col="white",cex= cexr)
	
	if(preselect){
		fg <- rep("black",length(all.x))
		all.sel <- dat$bin[x.names,levs]
		names(all.sel) <- x.names	
		fg[1:xn]<-all.sel
		fg[fg=="1"]<-"red"
		fg[fg=="0"]<-"blue"
	}else{
		fg[1:xn]="blue"
	}
	#fg[1:xn] <- "blue"
    symbols(all.x,all.y,squares=rep(sl*1.9,length(all.x)),add=T,inches=F,fg=fg,lwd=3)
	cexd<-4
    #first click defines the split
    #create a named squence, where all are scored as a 0
    #name it
	
	#fg<-all.sel
	#fg[fg==1]="red"
	#fg[fg==0] <- "blue"
	
    #symbols(all.x,all.y,squares=rep(sl*1.9,length(all.x)),add=T,inches=F,fg=fg,lwd=cexr)

    not.done=TRUE
    #Click to define
	if(!preselect){
	    click1 <- locator(n=1)
		#this isnhow kevin find the click location
		dist <- sqrt((click1$x[[1]]-all.x)^2 + (click1$y[[1]]-all.y)^2)
		sel.i <- which.min(dist)
		print(sel.i)

		###Done
		if(sel.i == xn+1){
			not.done=FALSE
			return(all.sel)
		}
		###All
		if(sel.i == xn+2){
			all.sel[1:xn] <- 1
			fg[1:xn] <- l.col
		}
		###None
		if(sel.i == xn+3){
			all.sel[1:xn] <- 0
			fg[1:xn] <- "blue"
		}
		###Reset
		if(sel.i == xn+4){
			#make everything score to a 0
			all.sel[] <- 0
			#now recolor them
			fg[1:xn] <- "blue"
			symbols(all.x,all.y,squares=rep(sl*1.9,length(all.x)),add=T,inches=F,fg=fg,lwd=cexd)
			
			#now click again
			click1 <- locator(n=1)
			#this isnhow kevin find the click location
			dist <- sqrt((click1$x[[1]]-all.x)^2 + (click1$y[[1]]-all.y)^2)
			sel.i <- which.min(dist)

			dev.set(grid_view)
			#now from 1 to the value selected
			pos.i <- 1:max((sel.i-1),1) 
			#make everything above your selection 0
			all.sel[neg.i] <- 0
			#now from selection to the start
			neg.i <- sel.i:xn	
			#score as a 1
			all.sel[pos.i] <- 1
			#define the colors
			fg[neg.i] <- "blue"
			fg[pos.i] <- "red"
			symbols(all.x,all.y,squares=rep(sl*1.9,length(all.x)),add=T,inches=F,fg=fg,lwd=cexd)
		}
		
		if(sel.i <= xn){
			#go to trace view
			dev.set(trace_view)
			#plot the trace
			PeakFunc7(dat,x.names[sel.i], t.type="blc")
			#go back to the grid
			dev.set(grid_view)
			#now from 1 to the value selected
			neg.i <- 1:max((sel.i-1),1) 
			#make everything above your selection 0
			all.sel[neg.i] <- 0
			#now from selection to the start
			pos.i <- sel.i:xn	
			#score as a 1
			all.sel[pos.i] <- 1
			#define the colors
			fg[neg.i] <- "blue"
			fg[pos.i] <- "red"
			symbols(all.x,all.y,squares=rep(sl*1.9,length(all.x)),add=T,inches=F,fg=fg,lwd=cexd)
		}
	}else{}

    while(not.done){
        symbols(all.x,all.y,squares=rep(sl*1.9,length(all.x)),add=T,inches=F,fg=fg,lwd=cexd)
        click1 <- locator(n=1)
        dist <- sqrt((click1$x[[1]]-all.x)^2 + (click1$y[[1]]-all.y)^2)
        sel.i <- which.min(dist)
        
        ###Done
        if(sel.i == xn+1){
            not.done=FALSE
            return(all.sel)
        }
        ###All
        if(sel.i == xn+2){
            all.sel[1:xn] <- 1
            fg[1:xn] <- l.col
        }
        ###None
        if(sel.i == xn+3){
            all.sel[1:xn] <- 0
            fg[1:xn] <- "blue"
        }
        ###Reset
        if(sel.i == xn+4){
            #make everything score to a 0
            all.sel[] <- 0
            #now recolor them
            fg[1:xn] <- "blue"
            symbols(all.x,all.y,squares=rep(sl*1.9,length(all.x)),add=T,inches=F,fg=fg,lwd=cexd)
            
            #now click again
            click1 <- locator(n=1)
            #this isnhow kevin find the click location
            dist <- sqrt((click1$x[[1]]-all.x)^2 + (click1$y[[1]]-all.y)^2)
            sel.i <- which.min(dist)
            print(sel.i)
            dev.set(grid_view)
            #now from 1 to the value selected
            neg.i <- 1:max((sel.i-1),1) 
            #make everything above your selection 0
            all.sel[neg.i] <- 0
            #now from selection to the start
            pos.i <- sel.i:xn	
            #score as a 1
            all.sel[pos.i] <- 1
            #define the colors
            fg[neg.i] <- "blue"
            fg[pos.i] <- "red"
            symbols(all.x,all.y,squares=rep(sl*1.9,length(all.x)),add=T,inches=F,fg=fg,lwd=cexd)
        }

        if(sel.i <= xn){
            #go to trace view
            dev.set(trace_view)
            #plot the trace
            PeakFunc7(dat,x.names[sel.i], t.type="blc")
            #go back to the grid
            dev.set(grid_view)
            
            if(all.sel[sel.i] ==0)
            {
                all.sel[sel.i] <- 1
                fg[sel.i] <- l.col
            }else{
                all.sel[sel.i] <- 0
                fg[sel.i] <- "blue"
            }
        }
    }
}

dice <- function(x, n,min.n=10)
{
	x.lst <- split(x, as.integer((seq_along(x) - 1) / n))
	x.i <- length(x.lst)
	if(length(x.lst[x.i]) < min.n & x.i > 1)
	{
		x.lst[[x.i-1]] <- c(x.lst[[x.i-1]],x.lst[[x.i]])
		x.lst <- x.lst[1:(x.i-1)]
	}
	return(x.lst)
}


#################################
#Welcome to a new method to score cells
#################################
RDView_2<-function(dat,levs=NULL){
	dat.name<-deparse(substitute(dat))
	cat(
"HOWDY partner, we r bout to score some rouwdy responses \n
from your cells. Please selact what we should score \n
and how we should initially sort this data. \n")
	if(is.null(levs)){
		levs<-setdiff(unique(dat$w.dat$wr1),"")
	}else{levs<-levs}
	cat(
"\nWitch window region would you like to score????\n \n What do you say?\n")
	
    lev<-levs[menu(levs)]
	#lev<-levs[26]
	#how would you like to sor this variable?
	cat("#############\nAnd how shall we sort? \n ############### \n")
	sorted.cells<-c.sort.2(dat$scp[grep(lev, names(dat$scp),value=T)],dat$c.dat$id)
	sorted.cells
	subset.list<-dice(sorted.cells, 300, 300/4)
	subset.list


	for(x.names in subset.list){
		graphics.off()
		scored.cells<-Trace_select_grid(dat,x.names, lev, t.type="blc",  preselect=T)
		dat$bin[names(which(scored.cells==1)),lev]=1
		dat$bin[names(which(scored.cells==0)),lev]=0
		
		cat("would you like to continue scoring?")
		choice<-select.list(c("yes","no"))
		if(choice=="yes"){
		}else{
			print("your dun")
			break
		}
	}
	assign(dat.name,dat, envir=.GlobalEnv)
}




















    
