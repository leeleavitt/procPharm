library(utils)

## designate packages to install/load
all_pkgs <-  c('reticulate', 'png', 'RColorBrewer', 'MALDIquant', 'data.table', 'docxtractr', 'xlsx')
## find packages that need to be installed
already_installed <- rownames( installed.packages() )
to_install <- setdiff(all_pkgs, already_installed)
if (length(to_install) > 0) {
    install.packages(to_install, dependencies=TRUE)
}
## now load all packages
sapply(all_pkgs, library, character.only=TRUE, logical.return=TRUE)

##############################################################################################
# Fucntions usied for data Input
##############################################################################################
#Function used in ReadDataDump.
#Converts time to minutes
ConvertTime <- function(x){
    vals <- strsplit(x,":")[[1]]
    retval <- NA
    if(length(vals)==3)
    {
        retval <- as.integer(vals[1])*60+as.integer(vals[2])+as.single(vals[3])/60
    }
    if(length(vals)==2)
    {
        retval <- as.integer(vals[1])+as.single(vals[2])/60
    }
    
    return(retval)
    
}

pharming_harvest <- function(main_dir="Y:/Mario Giacobassi/Manjus experiments", area_conversion=1.625, img_name_vec = NULL, image_question=T){
    cat('#########################################\nPHARM HARVEST\n#########################################\n')
    tryCatch(cat(readLines("Y:/Box Sync/procpharm/farm.txt"),sep='\n'), error=function(e) cat('\nPHARM ON PHARM ANIMAL\n'))
    cat('\nI am your data harvester. In this function I will\nextract data from any video files you have, and package everything together\n')
    #MAIN DIRECTORY SELECTION
    if(is.null(main_dir)){
        cat('\nHOLD UP! Tell me where your Pharm is!\nThis is where each experiment is in a seperate folder.\nExample "Z:/Lee Leavitt"\n')
        main_dir<-scan(n=1,what="character")
    }
    #AREA CONVERSION
    if(is.null(area_conversion)){
        cat('\nYou have NOT told me the area conversion please enter that now,For example\n4x bin = 1.625\n4x nobin = 3.25\n10x bin 3.25\n10x nobin = 6.5\n')
        area_conversion <- scan(n=1, what='numeric')
    }
    #IMPORT RETICULATE
    if( !library(reticulate, logical.return = T) ){
        install.packages('reticulate');library(reticulate)
    }
    if( !library(png, logical.return = T) ){
        install.packages('png');library(png)
    }

    #INITIALIZE EXPERIMENT names
    setwd(main_dir)
    cat('\nI have entered your Pharm,\n',main_dir, '\nselect each experiment I need to harvest.\n')
    exp_dir<-select.list(list.dirs(), multiple=T)
    (exp_dir_1 <- sub("./","",exp_dir))
    (exp_dir_2 <- lapply(strsplit(exp_dir_1, '/'), function(x) x[length(x)] ))
    (exp_dir_3 <- Reduce(c, exp_dir_2))
    (exp_dir_4 <- lapply(strsplit(exp_dir_3,' '), function(x) x[1]))
    (exp_dir_5 <- Reduce(c, exp_dir_4))
    rd.names <- paste0('RD.', exp_dir_5)
    cat('\nThese are the experiments I am going to process for you\n')
    cat(rd.names, sep='\n')
    total_start_time <- Sys.time()
    for( i in 1:length(exp_dir) ){
        setwd(exp_dir[i])

        #Input the names of the files that CP created
        cell_data_name <- list.files(pattern= '^cell.*[.]txt$')
        #cell_data_name<-"celldatacells_filtered.txt"

        ################################################
        #IMAGE IMPORT
        ################################################
        #Names of the files you want loaded into the RD file
        if( is.null(img_name_vec) ){
            img_name_vec<-c(
                "bf.gfp.tritc.start.png",
                "gfp.tritc.start.ci.ltl.rs.png",
                "tritc.start.ci.ltl.rs.png",
                "gfp.start.ci.ltl.rs.png",
                "bf.start.lab.png",
                "fura2.png",
                "fura2.divide.start.png",
                "roi.img.png")
        }
        # Add images
        img_list<-list()
        for( j in 1:length(img_name_vec) ){
            img_list[[ paste0("img",j) ]] <- tryCatch(readPNG(img_name_vec[j]), error=function(e)NULL)
        }
        if(image_question == T){
            cat('\nThese are the images I have attempted to load for you\nIf any are NULL, and want to add different images say yes to the \nnext question. You will be asked to select a png image for each loaction.\n\n')
            cat(img_name_vec, sep='\n')
            cat(str(img_list))

            cat('\nDO YOU WANT DIFFERENT IMAGES[y,n]?\n')
            img_reselect <- scan(n=1,what='character')
            if( img_reselect=='y' ){
                cat("\nAlright buddy I am going to give you options if you don't\nwant any image there just go ahead and press 0\n\n")
                png_imgs <- list.files(pattern='png')
                for( j in 1:8 ){
                    cat("\nWhat do you want for image ", j, '\n')
                    selection <- menu(png_imgs)
                    if(selection==0){
                        img_list[[paste0("img",j)]] <- NULL
                    }else{
                        img_list[[paste0("img",j)]] <- readPNG(png_imgs[selection])
                    }
                    cat('\nI have added ', png_imgs[selection],' to position ',j,'\n')
                }
             }
        }
            
            
            
        ########################################################
        #VIDEO PROCESSING
        ########################################################
        c_dat_make<-file.info(cell_data_name)$mtime
        video_data_name <- list.files(pattern="[vV]ideo.*[.]txt$")
        video_dat_make <- file.info(video_data_name)$mtime

        #If the video was make before the c.dat then you need to make the video_data.txt
        #this means if time_since_make is less than 1 than the 
        time_since_make <- video_dat_make - c_dat_make 

        if(length(video_data_name) < 1 ){
            py_pharm <- import('python_pharmer')
            video <- list.files(pattern="^[Vv]ideo.*nd2$")
            if( length(video) > 1 ){
                cat("\nSelect your video to process\n")
                video_num <- menu(video)
                video <- video[video_num]
            }
            # Now read in video DataS
            py_pharm$video_roi_extractor_faster(video)
        }else{
			if(time_since_make < 0){
				py_pharm <- import('python_pharmer')
				video <- list.files(pattern="^[Vv]ideo.*nd2$")
				if( length(video) > 1 ){
					cat("\nSelect your video to process\n")
					video_num <- menu(video)
					video <- video[video_num]
				}
				# Now read in video Data
				py_pharm$video_roi_extractor_faster(video)
			}
		}

        require(data.table)
        f2_img <- fread("video_data.txt")

        ##DCAST Mean Trace 
        start_time<-Sys.time()
        t.340 <- dcast(data = f2_img, ImageNumber ~ ObjectNumber, value.var = 'Intensity_MeanIntensity_f2_340')
        t.380 <- dcast(data = f2_img, ImageNumber ~ ObjectNumber, value.var = 'Intensity_MeanIntensity_f2_380')
        t.dat <- t.340/t.380
        print(Sys.time()-start_time)
        ###################################################
        #TIME INFO EXTRACTION
        ###################################################
        #Older version of NIS elements is not compatible with the nd2reader python package
        #until then the researcher will need to expor the time information before they start
        #this function. The funciton will sense whether the file is present or not.
		py_pharm <- import('python_pharmer')
		video <- list.files(pattern="^[Vv]ideo.*nd2$")
		
		if( length(list.files(pattern = "^time.*txt$")) < 1 ){
			py_pharm$time_info_gather( video )
		}
		time_info <- read.delim("time.info.txt", sep="\t")
		time_min <- round(time_info[2]/1000/60, digits=3)
        # Create row.names
        cell.names <- paste("X.", colnames(t.340)[-1], sep="")
        traces <- ls(pattern = "^[t.]{2}.*[0-9a-z]{3}")
        for( j in 1:length(traces) ){
            traze <- get( traces[j] )
            traze[,1] <- time_min[,1]
            traze <- as.data.frame(traze)
            colnames(traze) <- c("Time",cell.names)
            row.names(traze)<-time_min[,1]
            assign(traces[j], traze)
        }

        ########################################################
        # wr1 import
        ########################################################
		wrdef <- "wr1.docx"
        wrdef <- list.files(pattern = '^wr1')
		wrdef_logic <- grep(".docx", wrdef, ignore.case=T, value=T)
		if( length(wrdef_logic) == 1 ){
			require(docxtractr)
			wr <- docx.wr1.importer(wrdef)
			w.dat <- MakeWr.docx(t.dat, wr)		

			## Check for duplicated rows	
			if(length(which(duplicated(row.names(t.dat))))>=1){
				dup<-which(duplicated(row.names(t.dat)))
				paste(dup)
				t.dat<-t.dat[-dup,]
				w.dat<-w.dat[-dup,]
			}

		}else{
            tryCatch(
                wr <- ReadResponseWindowFile(wrdef)
                , error=function(e) print("YOU NEED A wr1.docx or wr1.csv")
            )
            Wr<-length(wr[,1])#complete and revise this section

            w.dat <- MakeWr(t.dat,wr)
		}
        ########################################################
        #CELL DATA PROCESSING
        ########################################################\
        c.dat<-read.delim(cell_data_name, header=T, sep="\t")
        #Now use the neewly created roi_checker i built in python
        roiToRemove <- py_pharm$roi_checker()
        if(length(roiToRemove)>0){
            cat("\nremoving ROI:", roiToRemove,'\n')
            c.dat <- c.dat[-roiToRemove,]
        }

        #These are the collumns needed for analysis
        c_dat_names<-c(
            "ObjectNumber", 
            "AreaShape_Area", 
            "AreaShape_Center_X", 
            "AreaShape_Center_Y", 
            "AreaShape_FormFactor", 
            "AreaShape_Perimeter",
            "Intensity_MeanIntensity_CGRP_start_ci",
            "Intensity_MeanIntensity_CGRP_end_ci",
            "Intensity_MeanIntensity_IB4_start_ci",
            "Intensity_MeanIntensity_IB4_end_ci",
			"Intensity_MeanIntensity_TRITC_end_ci",
            "Intensity_MeanIntensity_BF_start",
            "Intensity_MeanIntensity_BF_end",
            "Intensity_MeanIntensity_DAPI_ci",
            "Intensity_MeanIntensity_mcherry_start_ci",
            "Intensity_MeanIntensity_mcherry_end_ci",
            "Location_Center_X",
            "Location_Center_Y")	
        c_dat_rn <- c(
            "id",
            "area",
            "center.x.simplified", 
            "center.y.simplified",
            "circularity",
            "perimeter",
            "mean.gfp.start",
            "mean.gfp.end",
            "mean.cy5.start",
            "mean.cy5.end",
			"mean.tritc.end",
            "mean.bf.start",
            "mean.bf.end",
            "mean.dapi",
            "mean.mcherry.start",
            "mean.mcherry.end",
            "center.x",
            "center.y")
        # Find these collumns in the original c.dat
        (c_dat_names_update <- grep(paste0(c_dat_names, collapse="|"), names(c.dat), value=T, ignore.case=T) )
        # Find which collumn names remain
        cdnp_val<-c()
        for(j in 1:length(c_dat_names_update) ){
            value <- grep(c_dat_names_update[j], c_dat_names, ignore.case=T)
            if(length(value) > 0){
                cdnp_val[j] <- grep(c_dat_names_update[j], c_dat_names, ignore.case=T)
            }else{ cdnp_val[j] <- NA }
        }
        ( cdnp_val <- cdnp_val[ !is.na(cdnp_val) ] )
        # Update the renaming values
        c_dat_rn_update <- c_dat_rn[cdnp_val]
        # Create data subset to rename
        c_dat_1 <- c.dat[ c_dat_names_update ]
        # Rename the collumns in this new data.frame
        names(c_dat_1)<- c_dat_rn_update
        #put this newly renamed data frame at the start of the c.dat
        c.dat<-cbind(c_dat_1, c.dat)	
        #From line 81
        c.dat[,"id"] <- cell.names
        row.names(c.dat) <- cell.names
        #convert area
        c.dat[,"area"]<-c.dat$area*area_conversion

        # Initial and simple Data processing
        tmp.rd <- list(t.dat=t.dat,w.dat=w.dat,c.dat=c.dat)
        levs<-setdiff(unique(as.character(w.dat[,2])),"")
        snr.lim=5; hab.lim=.05; sm=2; ws=3; blc="SNIP"
        pcp <- ProcConstPharm(tmp.rd,sm,ws,blc)
        scp <- ScoreConstPharm(tmp.rd,pcp$blc,pcp$snr,pcp$der,snr.lim,hab.lim,sm)
        bin <- bScore(pcp$blc,pcp$snr,snr.lim,hab.lim,levs,tmp.rd$w.dat[,"wr1"])
        bin <- bin[,levs]
        bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp 

        tmp.rd <- list(t.dat=t.dat,t.340=t.340,t.380=t.380, 
        w.dat=w.dat,c.dat=c.dat, bin=bin, scp=scp, snr=pcp$snr, blc=pcp$blc, der=pcp$der) 
        
        tmp.rd <- TraceBrewer(tmp.rd) 
        tmp.rd <- c(tmp.rd, img_list)
        rd.name <- rd.names[i]
        f.name <- paste(rd.name,".Rdata",sep="")
        assign(rd.name,tmp.rd)
        save(list=rd.name,file=f.name)

        rm(f2.img)
        rm(rd.name)
        rm(tmp.rd)
        rm(c.dat)
        rm(cell.names)
        setwd(main_dir)
        gc()
        alarm()
        cat("\n#####################################################\nYour harvest has gone Successfully. Congratulations.\n#####################################################\n")
    }
    total_end_time <- Sys.time()
    cat("\n#####################################################\nTotal Harvest took. ",total_end_time - total_start_time,"\n#####################################################\n")
}

# readdatadump sam espinosa
ReadDataDump.se <- function(fname=NULL,wrdef=NULL, Wr=NULL, c.dat=NULL,img1=NULL,img2=NULL,img3=NULL,img4=NULL, img5=NULL, img6=NULL, img7=NULL, img8=NULL,rd.name=NULL,sep="\t"){
    require(png)
    require(zoom)
    require(RColorBrewer)
    require(MALDIquant)
        tmp <- read.delim(fname,fileEncoding="UCS-2LE",sep=sep)
        all.names <- names(tmp)
        time.name <- grep("Time",all.names,value=T,ignore=T)[1]
        if(time.name != "Time..ms."){warning(paste(time.name,"assumed to be in ms"))}
        
        id.name <- grep("ID",all.names,value=T,ignore=T)[1]
        if(id.name != "ID"){warning(paste(id.name,"assumed to be it ROI.ID"))}
        
        ratio.name <- grep("Ratio",all.names,value=T,ignore=T)
        if(is.na(ratio.name)){stop("no ratio data")}
        else{if(ratio.name != "Ratio.340.380"){warning(ratio.name,"assumed to be Ratio data")}}
            
        x.names <- unique(tmp[,id.name])
        x.tab <- table(tmp[,id.name])
        if(max(x.tab) != min(x.tab)){error("all ids do not have the same number of data points")}
        x.row <- max(x.tab)
        t.dat <- matrix(tmp[,ratio.name],byrow=FALSE,nrow=x.row)
        time.val <- tmp[tmp[,id.name]==x.names[1],time.name]
        
        if(length(grep(":",time.val[1]))==0)
        {
            x <- as.single(time.val)
            if(max(x) > 1000000)#in ms
            {
                x <- x/60000
            }
            else if(max(x) > 1500) #in seconds
            {
                x <- x/60	
            }		
            time.val <- x
        }
        else{time.val <- sapply(as.character(time.val),ConvertTime)}
        t.dat <- cbind(time.val,t.dat) #note assumption of ms
        t.dat <- as.data.frame(t.dat)
        t.dat<- t.dat[unique(row.names(t.dat)),]
        names(t.dat) <- c("Time",paste("X.",x.names,sep=""))
        
        
        if(!is.null(c.dat)){
        c.dat<-read.delim(file=c.dat,fileEncoding="UCS-2LE", sep=sep)
        c.dat.names<-names(c.dat)
        
        cx.name <- grep("Xpx",c.dat.names,value=T,ignore=T)
        if(is.na(cx.name)){stop("no Center X data")}
        else{if(cx.name != "CentreXpx"){warning(cx.name,"assumed to be Center X data")}}
        
        cy.name <- grep("Ypx",c.dat.names,value=T,ignore=T)
        if(is.na(cy.name)){stop("no Center Y data")}
        else{if(cy.name != "CentreYpx"){warning(cy.name,"assumed to be Center Y data")}}

        area.name <- grep("Area",c.dat.names,value=T,ignore=T)
        if(is.na(area.name)){stop("no Area data")}
        else{if(area.name != "ROIArea"){warning(paste(area.name,"assumed to be Area"))}}

        mean.gfp<-grep("MeanGreen",c.dat.names,value=T,ignore=T)
        if(length(mean.gfp)==0){warning(paste("no gfp.1 data from c.dat"))}
        else{if(mean.gfp!="MeanGFP"){warning(paste(mean.gfp, "assumed to be GFP.1"))}}
        
        mean.tritc<-grep("MeanBlue",c.dat.names,value=T,ignore=T)
        if(length(mean.tritc)==0){warning(paste("no tritc data from c.dat"))}
        else{if(mean.tritc!="MeanTRITC"){warning(paste(mean.tritc, "assumed to be TRITC"))}}
        
        cnames <- c(area.name, cx.name, cy.name, mean.gfp, mean.tritc)
    #	o.names <- setdiff(c.dat.names,c(time.name,id.name,area.name,ratio.name,cx.name,cy.name, mean.gfp, mean.tritc))
    #	if(length(o.names) > 0){warning(paste(o.names,"added to c.dat"));cnames <- c(cnames,o.names)}

        c.dat <- c.dat[,cnames]
        c.dat <- cbind(paste("X.",x.names,sep=""),c.dat)
        c.dat <- data.frame(c.dat)
        colnames(c.dat)[1:4] <- c("id","area","center.x", "center.y")
        
        # If gfp and tritc are not present then evaluate
        # 1st if there is only tritc, name the 6th column mean.tritc
        # 2nd if there is only gfp, name the 6th collumn mean.gfp
        # 3rd if there are both then rename both 6th and 7th collumn
        if(!length(mean.gfp)==0 & !length(mean.tritc)==0){
        if(length(mean.gfp)==0 & length(mean.tritc)==1){colnames(c.dat)[5]<-"mean.tritc"}
        if(length(mean.tritc)==0 & length(mean.gfp)==1){colnames(c.dat)[5]<-c("mean.gfp")}
        if(length(mean.tritc)==1 & length(mean.gfp)==1){colnames(c.dat)[5:6]<-c("mean.gfp","mean.tritc")}
        row.names(c.dat) <- c.dat[,"id"]
        }}
        else{
        area.name <- grep("Area",all.names,value=T,ignore=T)[1]
        if(is.na(area.name)){stop("no ROI.Area data")}
        else{if(area.name != "ROI.Area"){warning(paste(area.name,"assumed to be ROI.Area"))}}
        
        cx.name <- grep("Center.X",all.names,value=T,ignore=T)
        if(is.na(cx.name)){stop("no Center X data")}
        else{if(cx.name != "Center.X"){warning(cx.name,"assumed to be Center X data")}}
        
        cy.name <- grep("Center.Y",all.names,value=T,ignore=T)
        if(is.na(cy.name)){stop("no Center Y data")}
        else{if(cy.name != "Center.Y"){warning(cy.name,"assumed to be Center Y data")}}
        
        cnames <- c(area.name,cx.name,cy.name)
        c.dat <- tmp[match(x.names,tmp[,id.name]),cnames]
        c.dat <- cbind(paste("X.",x.names,sep=""),c.dat)
        c.dat <- data.frame(c.dat)
        names(c.dat)[1:4] <- c("id","area","center.X","center.Y") 
        row.names(c.dat) <- c.dat[,"id"]
    }
    if(!is.null(wrdef))
    {
        wr <- ReadResponseWindowFile(wrdef)
        Wr<-length(wr[,1])#complete and revise this section
        if(length(colnames(wr))<2){w.dat<-WrMultiplex(t.dat,wr,n=Wr)}
        else{w.dat <- MakeWr(t.dat,wr)}
        }
    else
    {
        WrCreate.rdd(t.dat, n=Wr)
        wr <- ReadResponseWindowFile("wr1.csv")
        w.dat <- MakeWr(t.dat,wr)
    }
    
    if(!is.null(img1)){img1<-readPNG(img1)}
    if(!is.null(img2)){img2<-readPNG(img2)}
    if(!is.null(img3)){img3<-readPNG(img3)}
    if(!is.null(img4)){img4<-readPNG(img4)}
    
    if(is.null(rd.name)){rd.name <- paste("RD",make.names(date()),sep="")}
    
    if(length(which(duplicated(row.names(t.dat))))>=1){
    dup<-which(duplicated(row.names(t.dat)))
    paste(dup)
    t.dat<-t.dat[-dup,]
    w.dat<-w.dat[-dup,]
    }
        
    tmp.rd <- list(t.dat=t.dat,w.dat=w.dat,c.dat=c.dat, img1=img1, img2=img2, img3=img3)
    f.name <- paste(rd.name,".Rdata",sep="")
    assign(rd.name,tmp.rd)
    save(list=rd.name,file=f.name)
    return(paste(nrow(tmp.rd$c.dat),"traces read saved to ",f.name))
    #save as RD file
}

# readdatadump Lee Leavitt
#ReadDataDump.lee <- function(fname=NULL,wrdef=NULL, Wr=NULL, c.dat=NULL,img1=NULL,img2=NULL,img3=NULL,img4=NULL,rd.name=NULL,sep="\t")
# fancy added for cell definer
ReadDataDump.lee <- function(rd.name=NULL,img1="bf.f2.png",img2="bf.f2.lab.png",img3="bf.png",img4=NULL,img5=NULL, img6=NULL, img7=NULL, img8=NULL, fancy=F,fname="Data (full).txt",wrdef="wr1.csv", Wr=NULL, c.dat="ROI Data.txt" ,sep="\t"){
        require(png)
        require(zoom)
        require(RColorBrewer)
        require(MALDIquant)

    ##################################################################################
    # Video Data import
    ##################################################################################
        
        if(length(fname)>1){
            tmp1 <- read.delim(fname[1],fileEncoding="UCS-2LE",sep=sep)
            tmp2 <- read.delim(fname[2],fileEncoding="UCS-2LE",sep=sep)
            tmp<-rbind(tmp1, tmp2)
        }else{
            tmp <- read.delim(fname,fileEncoding="UCS-2LE",sep=sep)
        }

        all.names <- names(tmp)
        
        time.name <- grep("Time",all.names,value=T,ignore=T)[1]
        if(time.name != "Time..ms."){warning(paste(time.name,"assumed to be in ms"))}
        
        id.name <- grep("ID",all.names,value=T,ignore=T)[1]
        if(id.name != "ID"){warning(paste(id.name,"assumed to be it ROI.ID"))}
        
        ratio.name <- grep("Ratio",all.names,value=T,ignore=T)
        if(is.na(ratio.name)){stop("no ratio data")}
        else{if(ratio.name != "Ratio.340.380"){warning(ratio.name,"assumed to be Ratio data")}}
            
        x.names <- unique(tmp[,id.name])
        x.tab <- table(tmp[,id.name])
        if(max(x.tab) != min(x.tab)){warning("all ids do not have the same number of data points")}
        x.row <- max(x.tab)
        t.dat <- matrix(tmp[,ratio.name],byrow=FALSE,nrow=x.row)
        time.val <- tmp[tmp[,id.name]==x.names[1],time.name]
        
        if(length(grep(":",time.val[1]))==0)
        {
            x <- as.single(time.val)
            if(max(x) > 1000000)#in ms
            {
                x <- x/60000
            }
            else if(max(x) > 1500) #in seconds
            {
                x <- x/60	
            }		
            time.val <- x
        }
        else{time.val <- sapply(as.character(time.val),ConvertTime)}
        t.dat <- cbind(time.val,t.dat) #note assumption of ms
        t.dat <- as.data.frame(t.dat)
        t.dat<- t.dat[unique(row.names(t.dat)),]
        names(t.dat) <- c("Time",paste("X.",x.names,sep=""))
        
    ##################################################################################
    # Cell Data import
    ##################################################################################

    if(!is.null(c.dat)){
        c.dat<-read.delim(file=c.dat,fileEncoding="UCS-2LE", sep=sep)
        c.dat.names<-names(c.dat)
        
        id.name <- grep("id",c.dat.names,value=T,ignore=T)
        if(is.na(id.name)){stop("no ID data")}
        else{if(id.name != "RoiID"){warning(cx.name,"assumed to be ID data")}}

        cx.name <- grep("Xpx",c.dat.names,value=T,ignore=T)
        if(is.na(cx.name)){stop("no Center X data")}
        else{if(cx.name != "CentreXpx"){warning(cx.name,"assumed to be Center X data")}}
        
        cy.name <- grep("Ypx",c.dat.names,value=T,ignore=T)
        if(is.na(cy.name)){stop("no Center Y data")}
        else{if(cy.name != "CentreYpx"){warning(cy.name,"assumed to be Center Y data")}}

        perimeter.name<-grep("perimeter", c.dat.names, value=T, ignore=T)
        if(is.na(perimeter.name)){stop("no Perimeter data")}
        else{if(perimeter.name != "Perimeter"){warning(paste(perimeter.name,"assumed to be Perimeter"))}}
        
        area.name <- grep("Area",c.dat.names,value=T,ignore=T)
        if(is.na(area.name)){stop("no Area data")}
        else{if(area.name != "ROIArea"){warning(paste(area.name,"assumed to be Area"))}}

        
        #mean.gfp<-grep("gfp.1",c.dat.names,value=T,ignore=T)
        mean.gfp<-grep("GFP",c.dat.names,value=T,ignore=F)
        if(length(mean.gfp)==0){mean.gfp<-grep("gfp",c.dat.names,value=T,ignore=T);warning(paste("no gfp.1 data from c.dat"))}
        else{if(mean.gfp!="MeanGFP"){warning(paste(mean.gfp, "assumed to be GFP.1"))}}
        
        mean.gfp.2<-grep("gfp.2",c.dat.names,value=T,ignore=T)
        if(length(mean.gfp.2)==0){warning(paste("no gfp.2 data from c.dat"))}
        else{if(mean.gfp.2!="MeanGFP"){warning(paste(mean.gfp.2, "assumed to be GFP.2"))}}
        
        mean.tritc<-grep("TRITC",c.dat.names,value=T,ignore=F)
        if(length(mean.tritc)==0){warning(paste("no tritc data from c.dat"))}
        else{if(mean.tritc!="MeanTRITC"){warning(paste(mean.tritc, "assumed to be TRITC"))}}
        
        mean.dapi<-grep("DAPI",c.dat.names,value=T,ignore=F)
        if(length(mean.dapi)==0){warning(paste("no dapi data from c.dat"))}
        else{if(mean.dapi!="MeanDAPI"){warning(paste(mean.dapi, "assumed to be DAPI"))}}

        cnames <- c(id.name,area.name, perimeter.name, cx.name, cy.name, mean.gfp, mean.gfp.2, mean.tritc, mean.dapi)
    #	o.names <- setdiff(c.dat.names,c(time.name,id.name,area.name,ratio.name,cx.name,cy.name, mean.gfp, mean.tritc))
    #	if(length(o.names) > 0){warning(paste(o.names,"added to c.dat"));cnames <- c(cnames,o.names)}
        
        c.dat<-c.dat[cnames]#create c.dat with specified collumns from cnames
        c.dat <- c.dat[order(c.dat[,id.name]),] # order rows by ROIid
        c.dat[,id.name] <- paste("X.",c.dat[,id.name],sep="")#rename ROIid with a X.cell#
        row.names(c.dat)<-c.dat[,id.name]# assign row.names the ROIid name
        c.dat <- data.frame(c.dat)#convert to data frame
        colnames(c.dat)[1:5] <- c("id","area","perimeter","center.x", "center.y")#rename collumns these names
        c.dat["circularity"]<-((c.dat$perimeter^2)/(4*pi*c.dat$area)) # create a circularity measurement

        ## If the class of the collumn is a factor, then the collumn is filled with "N/A"
        # therefore make the NULL/ remove it.  If not, then perform an unecessarily complex 
        # set of selection to rename the collumn what you want.
        if(class(c.dat[,mean.gfp])=="factor"){c.dat[,mean.gfp]<-NULL
        }else{
        colnames(c.dat)[which(colnames(c.dat)==mean.gfp)]<-"mean.gfp"}
        
        if(class(c.dat[,mean.gfp.2])=="factor"){c.dat[,mean.gfp.2]<-NULL
        }else{colnames(c.dat)[which(colnames(c.dat)==mean.gfp.2)]<-"mean.gfp.2"}
        
        if(class(c.dat[,mean.tritc])=="factor"){c.dat[,mean.tritc]<-NULL
        }else{colnames(c.dat)[which(colnames(c.dat)==mean.tritc)]<-"mean.tritc"}
        
        if(class(c.dat[,mean.dapi])=="factor"){c.dat[,mean.dapi]<-NULL
        }else{colnames(c.dat)[which(colnames(c.dat)==mean.dapi)]<-"mean.dapi"}

        }
        else{
        area.name <- grep("Area",all.names,value=T,ignore=T)[1]
        if(is.na(area.name)){stop("no ROI.Area data")}
        else{if(area.name != "ROI.Area"){warning(paste(area.name,"assumed to be ROI.Area"))}}
        
        cx.name <- grep("Center.X",all.names,value=T,ignore=T)
        if(is.na(cx.name)){stop("no Center X data")}
        else{if(cx.name != "Center.X"){warning(cx.name,"assumed to be Center X data")}}
        
        cy.name <- grep("Center.Y",all.names,value=T,ignore=T)
        if(is.na(cy.name)){stop("no Center Y data")}
        else{if(cy.name != "Center.Y"){warning(cy.name,"assumed to be Center Y data")}}
        
        cnames <- c(area.name,cx.name,cy.name)
        c.dat <- tmp[match(x.names,tmp[,id.name]),cnames]
        c.dat <- cbind(paste("X.",x.names,sep=""),c.dat)
        c.dat <- data.frame(c.dat)
        names(c.dat)[1:4] <- c("id","area","center.x","center.y") 
        row.names(c.dat) <- c.dat[,"id"]
    }
    #####################################################
    # Window Region Definition
    #####################################################

    if(!is.null(wrdef))
        {
            wr <- ReadResponseWindowFile(wrdef)
            Wr<-length(wr[,1])#complete and revise this section
            if(length(colnames(wr))<2){w.dat<-WrMultiplex(t.dat,wr,n=Wr)}
            else{w.dat <- MakeWr(t.dat,wr)}
            }
        else
        {
            WrCreate.rdd(t.dat, n=Wr)
            wr <- ReadResponseWindowFile("wr1.csv")
            w.dat <- MakeWr(t.dat,wr)
        }
        tmp.rd <- list(t.dat=t.dat,w.dat=w.dat,c.dat=c.dat)
        #####################################################
        #Create Despiked data
        #####################################################
        wts <- tmp.rd$t.dat
        for(i in 1:5) #run the despike 5 times.
        {
            wt.mn3 <- Mean3(wts)
            wts <- SpikeTrim2(wts,1,-1)
            print(sum(is.na(wts))) #this prints out the number of points removed should be close to 0 after 5 loops.
            wts[is.na(wts)] <- wt.mn3[is.na(wts)]
        }
        tmp.rd$mp <- wts

        # Initial Data processing
        levs<-setdiff(unique(as.character(w.dat[,2])),"")
        snr.lim=4;hab.lim=.05;sm=3;ws=30;blc="SNIP"
        
        pcp <- ProcConstPharm(tmp.rd$mp,sm,ws,blc)
        scp <- ScoreConstPharm(tmp.rd,pcp$blc,pcp$snr,pcp$der,snr.lim,hab.lim,sm)
        bin <- bScore(pcp$blc,pcp$snr,snr.lim,hab.lim,levs,tmp.rd$w.dat[,"wr1"])
        bin <- bin[,levs]
        bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
        bin<-pf.function(bin,levs)
        
        tmp.rd$t.dat<-t.dat
        tmp.rd$w.dat<-w.dat
        tmp.rd$c.dat<-c.dat
        tmp.rd$bin<-bin
        tmp.rd$scp<-scp
        tmp.rd$snr<-pcp$snr
        tmp.rd$blc<-pcp$blc
        tmp.rd$der<-pcp$der
        # Add images
        if(!is.null(img1)){tmp.rd$img1<-readPNG(img1)}
        if(!is.null(img2)){tmp.rd$img2<-readPNG(img2)}
        if(!is.null(img3)){tmp.rd$img3<-readPNG(img3)}
        if(!is.null(img4)){tmp.rd$img4<-readPNG(img4)}
        if(!is.null(img5)){tmp.rd$img5<-readPNG(img5)}
        if(!is.null(img6)){tmp.rd$img6<-readPNG(img6)}
        if(!is.null(img7)){tmp.rd$img7<-readPNG(img7)}
        if(!is.null(img8)){tmp.rd$img8<-readPNG(img8)}


    #####################################################
    # Cell Label Scoring	
    #####################################################

        if(fancy==TRUE){tmp.rd<-cell.creator(tmp.rd)}		# Create list of binary  labeled neurons}
        else{tmp.rd$cells<-NULL}
        
        if(is.null(rd.name)){rd.name <- paste("RD",make.names(date()),sep="")}
        
        if(length(which(duplicated(row.names(t.dat))))>=1){
        dup<-which(duplicated(row.names(t.dat)))
        paste(dup)
        t.dat<-t.dat[-dup,]
        w.dat<-w.dat[-dup,]
        }
        
        
        f.name <- paste(rd.name,".Rdata",sep="")
        assign(rd.name,tmp.rd)
        save(list=rd.name,file=f.name)
        return(paste(nrow(tmp.rd$c.dat),"traces read saved to ",f.name))
        #save as RD file
    }

    #HAS TRACEBREWER
    ReadDataDump.lee.2 <- function(rd.name=NULL,img1="bf.f2.png",img2="bf.f2.lab.png",img3=NULL,img4=NULL,img5=NULL, img6=NULL, img7=NULL, img8=NULL, fancy=F,fname="Data (full).txt",wrdef="wr1.docx", Wr=NULL, c.dat="ROI Data.txt" ,sep="\t")
    {
    require(png)
    require(MALDIquant)

    ##################################################################################
    # Video Data import
    ##################################################################################
        
        if(length(fname)>1){
            tmp1 <- read.delim(fname[1],fileEncoding="UCS-2LE",sep=sep)
            tmp2 <- read.delim(fname[2],fileEncoding="UCS-2LE",sep=sep)
            tmp<-rbind(tmp1, tmp2)
        }else{
            tmp <- read.delim(fname,fileEncoding="UCS-2LE",sep=sep)
        }

        all.names <- names(tmp)
        
        time.name <- grep("Time",all.names,value=T,ignore=T)[1]
        if(time.name != "Time..ms."){warning(paste(time.name,"assumed to be in ms"))}
        
        id.name <- grep("ID",all.names,value=T,ignore=T)[1]
        if(id.name != "ID"){warning(paste(id.name,"assumed to be it ROI.ID"))}
        
        ratio.name <- grep("Ratio",all.names,value=T,ignore=T)
        if(is.na(ratio.name)){stop("no ratio data")}else{if(ratio.name != "Ratio.340.380"){warning(ratio.name,"assumed to be Ratio data")}}
            
        x.names <- unique(tmp[,id.name])
        x.tab <- table(tmp[,id.name])
        if(max(x.tab) != min(x.tab)){warning("all ids do not have the same number of data points")}
        x.row <- max(x.tab)
        t.dat <- matrix(tmp[,ratio.name],byrow=FALSE,nrow=x.row)
        time.val <- tmp[tmp[,id.name]==x.names[1],time.name]
        
        if(length(grep(":",time.val[1]))==0)
        {
            x <- as.single(time.val)
            if(max(x) > 1000000)#in ms
            {
                x <- x/60000
            }
            else if(max(x) > 1500) #in seconds
            {
                x <- x/60	
            }		
            time.val <- x
        }else{time.val <- sapply(as.character(time.val),ConvertTime)}
        t.dat <- cbind(time.val,t.dat) #note assumption of ms
        t.dat <- as.data.frame(t.dat)
        t.dat<- t.dat[unique(row.names(t.dat)),]
        names(t.dat) <- c("Time",paste("X.",x.names,sep=""))
        
    ##################################################################################
    # Cell Data import
    ##################################################################################

    if(!is.null(c.dat)){
        c.dat<-read.delim(file=c.dat,fileEncoding="UCS-2LE", sep=sep)
        c.dat.names<-names(c.dat)
        
        id.name <- grep("id",c.dat.names,value=T,ignore=T)
        if(is.na(id.name)){stop("no ID data")
        }else{if(id.name != "RoiID"){warning(id.name,"assumed to be ID data")}}

        cx.name <- grep("Xpx",c.dat.names,value=T,ignore=T)
        if(is.na(cx.name)){stop("no Center X data")}else{if(cx.name != "CentreXpx"){warning(cx.name,"assumed to be Center X data")}}
        
        cy.name <- grep("Ypx",c.dat.names,value=T,ignore=T)
        if(is.na(cy.name)){stop("no Center Y data")}else{if(cy.name != "CentreYpx"){warning(cy.name,"assumed to be Center Y data")}}

        perimeter.name<-grep("perimeter", c.dat.names, value=T, ignore=T)
        if(is.na(perimeter.name)){stop("no Perimeter data")}else{if(perimeter.name != "Perimeter"){warning(paste(perimeter.name,"assumed to be Perimeter"))}}
        
        area.name <- grep("Area",c.dat.names,value=T,ignore=T)
        if(is.na(area.name)){stop("no Area data")}else{if(area.name != "ROIArea"){warning(paste(area.name,"assumed to be Area"))}}

        
        #mean.gfp<-grep("gfp.1",c.dat.names,value=T,ignore=T)
        mean.gfp<-grep("GFP",c.dat.names,value=T,ignore=F)
        if(length(mean.gfp)==0){mean.gfp<-grep("gfp",c.dat.names,value=T,ignore=T);warning(paste("no gfp.1 data from c.dat"))}else{if(mean.gfp!="MeanGFP"){warning(paste(mean.gfp, "assumed to be GFP.1"))}}
        
        mean.gfp.2<-grep("gfp.2",c.dat.names,value=T,ignore=T)
        if(length(mean.gfp.2)==0){warning(paste("no gfp.2 data from c.dat"))}else{if(mean.gfp.2!="MeanGFP"){warning(paste(mean.gfp.2, "assumed to be GFP.2"))}}
        
        mean.tritc<-grep("TRITC",c.dat.names,value=T,ignore=F)
        if(length(mean.tritc)==0){warning(paste("no tritc data from c.dat"))}else{if(mean.tritc!="MeanTRITC"){warning(paste(mean.tritc, "assumed to be TRITC"))}}
        
        mean.cy5<-grep("TRITC",c.dat.names,value=T,ignore=F)
        if(length(mean.cy5)==0){warning(paste("no cy5 data from c.dat"))}else{if(mean.cy5!="MeanTRITC"){warning(paste(mean.cy5, "assumed to be TRITC"))}}

        mean.dapi<-grep("DAPI",c.dat.names,value=T,ignore=F)
        if(length(mean.dapi)==0){warning(paste("no dapi data from c.dat"))}
        else{if(mean.dapi!="MeanDAPI"){warning(paste(mean.dapi, "assumed to be DAPI"))}}

        cnames <- c(id.name,area.name, perimeter.name, cx.name, cy.name, mean.gfp, mean.gfp.2, mean.tritc,mean.cy5, mean.dapi)
    #	o.names <- setdiff(c.dat.names,c(time.name,id.name,area.name,ratio.name,cx.name,cy.name, mean.gfp, mean.tritc))
    #	if(length(o.names) > 0){warning(paste(o.names,"added to c.dat"));cnames <- c(cnames,o.names)}
        
        c.dat<-c.dat[cnames]#create c.dat with specified collumns from cnames
        c.dat <- c.dat[order(c.dat[,id.name]),] # order rows by ROIid
        c.dat[,id.name] <- paste("X.",c.dat[,id.name],sep="")#rename ROIid with a X.cell#
        row.names(c.dat)<-c.dat[,"RoiID"]# assign row.names the ROIid name
        c.dat <- data.frame(c.dat)#convert to data frame
        colnames(c.dat)[1:5] <- c("id","area","perimeter","center.x", "center.y")#rename collumns these names
        c.dat["circularity"]<-((c.dat$perimeter^2)/(4*pi*c.dat$area)) # create a circularity measurement

        ## If the class of the collumn is a factor, then the collumn is filled with "N/A"
        # therefore make the NULL/ remove it.  If not, then perform an unecessarily complex 
        # set of selection to rename the collumn what you want.
        if(class(c.dat[,mean.gfp])=="factor"){c.dat[,mean.gfp]<-NULL
        }else{
        colnames(c.dat)[which(colnames(c.dat)==mean.gfp)]<-"mean.gfp"}
        
        if(class(c.dat[,mean.gfp.2])=="factor"){c.dat[,mean.gfp.2]<-NULL
        }else{colnames(c.dat)[which(colnames(c.dat)==mean.gfp.2)]<-"mean.gfp.2"}
        
        if(class(c.dat[,mean.tritc])=="factor"){c.dat[,mean.tritc]<-NULL
        }else{colnames(c.dat)[which(colnames(c.dat)==mean.tritc)]<-"mean.tritc"}
        
        if(class(c.dat[,mean.cy5])=="factor"){c.dat[,mean.cy5]<-NULL
        }else{colnames(c.dat)[which(colnames(c.dat)==mean.cy5)]<-"mean.cy5"}

        if(class(c.dat[,mean.dapi])=="factor"){c.dat[,mean.dapi]<-NULL
        }else{colnames(c.dat)[which(colnames(c.dat)==mean.dapi)]<-"mean.dapi"}

        }
        else{
        area.name <- grep("Area",all.names,value=T,ignore=T)[1]
        if(is.na(area.name)){stop("no ROI.Area data")}
        else{if(area.name != "ROI.Area"){warning(paste(area.name,"assumed to be ROI.Area"))}}
        
        cx.name <- grep("Center.X",all.names,value=T,ignore=T)
        if(is.na(cx.name)){stop("no Center X data")}
        else{if(cx.name != "Center.X"){warning(cx.name,"assumed to be Center X data")}}
        
        cy.name <- grep("Center.Y",all.names,value=T,ignore=T)
        if(is.na(cy.name)){stop("no Center Y data")}
        else{if(cy.name != "Center.Y"){warning(cy.name,"assumed to be Center Y data")}}
        
        cnames <- c(area.name,cx.name,cy.name)
        c.dat <- tmp[match(x.names,tmp[,id.name]),cnames]
        c.dat <- cbind(paste("X.",x.names,sep=""),c.dat)
        c.dat <- data.frame(c.dat)
        names(c.dat)[1:4] <- c("id","area","center.x","center.y") 
        row.names(c.dat) <- c.dat[,"id"]
    }
    #####################################################
    # Window Region Definition
    #####################################################
    ############ wr1 import
    wrdef<-"wr1.docx"
    require(docxtractr)
    if(!is.null(wrdef)){
            wr<-docx.wr1.importer(wrdef)
            w.dat<-MakeWr.docx(t.dat,wr)
    }


    #if(!is.null(wrdef))
    #	{
    ##		wr <- ReadResponseWindowFile(wrdef)
    #		Wr<-length(wr[,1])#complete and revise this section
    #		if(length(colnames(wr))<2){w.dat<-WrMultiplex(t.dat,wr,n=Wr)}
    #		else{w.dat <- MakeWr(t.dat,wr)}
    #		}
    #	else
    #	{
    #		WrCreate.rdd(t.dat, n=Wr)
    #		wr <- ReadResponseWindowFile("wr1.csv")
    #		w.dat <- MakeWr(t.dat,wr)
    #	}
        tmp.rd <- list(t.dat=t.dat,w.dat=w.dat,c.dat=c.dat)
        #####################################################
        #Create Despiked data
        #####################################################
        wts <- tmp.rd$t.dat
        for(i in 1:5) #run the despike 5 times.
        {
            wt.mn3 <- Mean3(wts)
            wts <- SpikeTrim2(wts,1,-1)
            print(sum(is.na(wts))) #this prints out the number of points removed should be close to 0 after 5 loops.
            wts[is.na(wts)] <- wt.mn3[is.na(wts)]
        }
        tmp.rd$mp <- wts

        # Initial Data processing
        levs<-setdiff(unique(as.character(w.dat[,2])),"")
        snr.lim=4;hab.lim=.05;sm=2;ws=20;blc="SNIP"
        
        pcp <- ProcConstPharm(tmp.rd$mp,sm,ws,blc)
        scp <- ScoreConstPharm(tmp.rd,pcp$blc,pcp$snr,pcp$der,snr.lim,hab.lim,sm)
        bin <- bScore(pcp$blc,pcp$snr,snr.lim,hab.lim,levs,tmp.rd$w.dat[,"wr1"])
        bin <- bin[,levs]
        bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
        bin<-pf.function(bin,levs)
        
        tmp.rd$t.dat<-t.dat
        tmp.rd$w.dat<-w.dat
        tmp.rd$c.dat<-c.dat
        tmp.rd$bin<-bin
        tmp.rd$scp<-scp
        tmp.rd$snr<-pcp$snr
        tmp.rd$blc<-pcp$blc
        tmp.rd$der<-pcp$der
        
        tmp.rd<-TraceBrewer(tmp.rd)

        # Add images
        if(!is.null(img1)){tmp.rd$img1<-readPNG(img1)}
        if(!is.null(img2)){tmp.rd$img2<-readPNG(img2)}
        if(!is.null(img3)){tmp.rd$img3<-readPNG(img3)}
        if(!is.null(img4)){tmp.rd$img4<-readPNG(img4)}
        if(!is.null(img5)){tmp.rd$img5<-readPNG(img5)}
        if(!is.null(img6)){tmp.rd$img6<-readPNG(img6)}
        if(!is.null(img7)){tmp.rd$img7<-readPNG(img7)}
        if(!is.null(img8)){tmp.rd$img8<-readPNG(img8)}


    #####################################################
    # Cell Label Scoring	
    #####################################################

    if(fancy==TRUE){tmp.rd<-cell.creator(tmp.rd)}		# Create list of binary  labeled neurons}
    else{tmp.rd$cells<-NULL}
    
    if(is.null(rd.name)){rd.name <- paste("RD",make.names(date()),sep="")}
    
    if(length(which(duplicated(row.names(t.dat))))>=1){
    dup<-which(duplicated(row.names(t.dat)))
    paste(dup)
    t.dat<-t.dat[-dup,]
    w.dat<-w.dat[-dup,]
    }
    
    
    f.name <- paste(rd.name,".Rdata",sep="")
    assign(rd.name,tmp.rd)
    save(list=rd.name,file=f.name)
    return(paste(nrow(tmp.rd$c.dat),"traces read saved to ",f.name))
    #save as RD file
}
# readdatadump Lee Leavitt 170209

# readdatadump Lee Leavitt 170209
#ReadDataDump.lee <- function(fname=NULL,wrdef=NULL, Wr=NULL, c.dat=NULL,img1=NULL,img2=NULL,img3=NULL,img4=NULL,rd.name=NULL,sep="\t")
# fancy added for cell definer
#this import now has a 340 and 380 ch
ReadDataDump.microglia <- function(rd.name=NULL,img1="bf.f2.png",img2="bf.f2.lab.png",img3="bf.png",img4=NULL,img5=NULL, img6=NULL, img7=NULL, img8=NULL, fancy=F,fname="Data (full).txt",wrdef="wr1.docx", Wr=NULL, c.dat="ROI Data.txt" ,sep="\t"){
    require(png)

    require(RColorBrewer)
    require(MALDIquant)

    ##################################################################################
    # Video Data import
    ##################################################################################
        
        if(length(fname)>1){
            tmp1 <- read.delim(fname[1],fileEncoding="UCS-2LE",sep=sep)
            tmp2 <- read.delim(fname[2],fileEncoding="UCS-2LE",sep=sep)
            tmp<-rbind(tmp1, tmp2)
        }else{
            tmp <- read.delim(fname,fileEncoding="UCS-2LE",sep=sep)
        }

        all.names <- names(tmp)
        
        time.name <- grep("Time",all.names,value=T,ignore=T)[1]
        if(time.name != "Time..ms."){warning(paste(time.name,"assumed to be in ms"))}
        
        id.name <- grep("ID",all.names,value=T,ignore=T)[1]
        if(id.name != "ID"){warning(paste(id.name,"assumed to be it ROI.ID"))}
        
        ratio.name <- grep("Ratio",all.names,value=T,ignore=T)
        if(is.na(ratio.name)){stop("no ratio data")}
        else{if(ratio.name != "Ratio.340.380"){warning(ratio.name,"assumed to be Ratio data")}}
            
        x.names <- unique(tmp[,id.name])
        x.tab <- table(tmp[,id.name])
        if(max(x.tab) != min(x.tab)){warning("all ids do not have the same number of data points")}
        x.row <- max(x.tab)
        t.dat <- matrix(tmp[,ratio.name],byrow=FALSE,nrow=x.row)
        time.val <- tmp[tmp[,id.name]==x.names[1],time.name]
        
        if(length(grep(":",time.val[1]))==0)
        {
            x <- as.single(time.val)
            if(max(x) > 1000000)#in ms
            {
                x <- x/60000
            }
            else if(max(x) > 1500) #in seconds
            {
                x <- x/60	
            }		
            time.val <- x
        }
        else{time.val <- sapply(as.character(time.val),ConvertTime)}
        t.dat <- cbind(time.val,t.dat) #note assumption of ms
        t.dat <- as.data.frame(t.dat)
        t.dat<- t.dat[unique(row.names(t.dat)),]
        names(t.dat) <- c("Time",paste("X.",x.names,sep=""))
        
    ##################################################################################
    # Cell Data import
    ##################################################################################

    if(!is.null(c.dat)){
        c.dat<-read.delim(file=c.dat,fileEncoding="UCS-2LE", sep=sep)
        c.dat.names<-names(c.dat)
        
        id.name <- grep("id",c.dat.names,value=T,ignore=T)
        if(is.na(id.name)){stop("no ID data")}
        else{if(id.name != "RoiID"){warning(id.name,"assumed to be ID data")}}

        cx.name <- grep("Xpx",c.dat.names,value=T,ignore=T)
        if(is.na(cx.name)){stop("no Center X data")}
        else{if(cx.name != "CentreXpx"){warning(cx.name,"assumed to be Center X data")}}
        
        cy.name <- grep("Ypx",c.dat.names,value=T,ignore=T)
        if(is.na(cy.name)){stop("no Center Y data")}
        else{if(cy.name != "CentreYpx"){warning(cy.name,"assumed to be Center Y data")}}

        perimeter.name<-grep("Perimeter", c.dat.names, value=T, ignore=T)
        if(is.na(perimeter.name)){stop("no Perimeter data")}
        else{if(perimeter.name != "Perimeter"){warning(paste(perimeter.name,"assumed to be Perimeter"))}}
        
        area.name <- grep("ROIArea",c.dat.names,value=T,ignore=T)
        if(is.na(area.name)){stop("no Area data")}
        else{if(area.name != "ROIArea"){warning(paste(area.name,"assumed to be Area"))}}

        
        #mean.gfp<-grep("gfp.1",c.dat.names,value=T,ignore=T)
        mean.gfp.start<-grep("MeanGFP.start",c.dat.names,value=T,ignore=F)
        if(length(mean.gfp.start)==0){mean.gfp.start<-grep("gfp",c.dat.names,value=T,ignore=T);warning(paste("no gfp.1 data from c.dat"))}
        else{if(mean.gfp.start!="MeanGFP"){warning(paste(mean.gfp.start, "assumed to be GFP.1"))}}
        
        mean.gfp.end<-grep("MeanGFP.end",c.dat.names,value=T,ignore=T)
        if(length(mean.gfp.end)==0){warning(paste("no gfp.2 data from c.dat"))}
        else{if(mean.gfp.end!="MeanGFP"){warning(paste(mean.gfp.end, "assumed to be GFP.2"))}}
        
        mean.tritc.start<-grep("MeanTRITC.start",c.dat.names,value=T,ignore=F)
        if(length(mean.tritc.start)==0){warning(paste("no tritc data from c.dat"))}
        else{if(mean.tritc.start!="MeanTRITC"){warning(paste(mean.tritc.start, "assumed to be TRITC"))}}
        
        mean.tritc.end<-grep("MeanTRITC.end",c.dat.names,value=T,ignore=F)
        if(length(mean.tritc.end)==0){warning(paste("no tritc data from c.dat"))}
        else{if(mean.tritc.end!="MeanTRITC"){warning(paste(mean.tritc.end, "assumed to be TRITC"))}}

        mean.dapi<-grep("DAPI",c.dat.names,value=T,ignore=F)
        if(length(mean.dapi)==0){warning(paste("no dapi data from c.dat"))}
        else{if(mean.dapi!="MeanDAPI"){warning(paste(mean.dapi, "assumed to be DAPI"))}}

        cnames <- c(id.name,area.name, perimeter.name, cx.name, cy.name, mean.gfp.start, mean.gfp.end, mean.tritc.start, mean.tritc.end, mean.dapi)
    #	o.names <- setdiff(c.dat.names,c(time.name,id.name,area.name,ratio.name,cx.name,cy.name, mean.gfp, mean.tritc))
    #	if(length(o.names) > 0){warning(paste(o.names,"added to c.dat"));cnames <- c(cnames,o.names)}
        
        c.dat<-c.dat[cnames]#create c.dat with specified collumns from cnames
        c.dat <- c.dat[order(c.dat[,id.name]),] # order rows by ROIid
        c.dat[,id.name] <- paste("X.",c.dat[,id.name],sep="")#rename ROIid with a X.cell#
        row.names(c.dat)<-c.dat[,id.name]# assign row.names the ROIid name
        c.dat <- data.frame(c.dat)#convert to data frame
        colnames(c.dat)[1:5] <- c("id","area","perimeter","center.x", "center.y")#rename collumns these names
        c.dat["circularity"]<-((c.dat$perimeter^2)/(4*pi*c.dat$area)) # create a circularity measurement

        ## If the class of the collumn is a factor, then the collumn is filled with "N/A"
        # therefore make the NULL/ remove it.  If not, then perform an unecessarily complex 
        # set of selection to rename the collumn what you want.
        if(class(c.dat[,mean.gfp.start])=="factor"){c.dat[,mean.gfp.start]<-NULL
        }else{
        colnames(c.dat)[which(colnames(c.dat)==mean.gfp.start)]<-"mean.gfp.start"}
        
        if(class(c.dat[,mean.gfp.end])=="factor"){c.dat[,mean.gfp.end]<-NULL
        }else{colnames(c.dat)[which(colnames(c.dat)==mean.gfp.end)]<-"mean.gfp.end"}
        
        if(class(c.dat[,mean.tritc.start])=="factor"){c.dat[,mean.tritc.start]<-NULL
        }else{colnames(c.dat)[which(colnames(c.dat)==mean.tritc.start)]<-"mean.tritc.start"}
        
        if(class(c.dat[,mean.tritc.end])=="factor"){c.dat[,mean.tritc.end]<-NULL
        }else{colnames(c.dat)[which(colnames(c.dat)==mean.tritc.end)]<-"mean.tritc.end"}


        if(class(c.dat[,mean.dapi])=="factor"){c.dat[,mean.dapi]<-NULL
        }else{colnames(c.dat)[which(colnames(c.dat)==mean.dapi)]<-"mean.dapi"}

        }
        else{
        area.name <- grep("Area",all.names,value=T,ignore=T)[1]
        if(is.na(area.name)){stop("no ROI.Area data")}
        else{if(area.name != "ROI.Area"){warning(paste(area.name,"assumed to be ROI.Area"))}}
        
        cx.name <- grep("Center.X",all.names,value=T,ignore=T)
        if(is.na(cx.name)){stop("no Center X data")}
        else{if(cx.name != "Center.X"){warning(cx.name,"assumed to be Center X data")}}
        
        cy.name <- grep("Center.Y",all.names,value=T,ignore=T)
        if(is.na(cy.name)){stop("no Center Y data")}
        else{if(cy.name != "Center.Y"){warning(cy.name,"assumed to be Center Y data")}}
        
        cnames <- c(area.name,cx.name,cy.name)
        c.dat <- tmp[match(x.names,tmp[,id.name]),cnames]
        c.dat <- cbind(paste("X.",x.names,sep=""),c.dat)
        c.dat <- data.frame(c.dat)
        names(c.dat)[1:4] <- c("id","area","center.x","center.y") 
        row.names(c.dat) <- c.dat[,"id"]
    }
    #####################################################
    # Window Region Definition
    #####################################################
        ########################################################
        # wr1 import
        ########################################################
        
		
		wrdef <- "wr1.docx"
        wrdef <- list.files(pattern = '^wr1')
		wrdef_logic <- grep(".docx", wrdef, ignore.case=T, value=T)
		if( length(wrdef_logic) >= 1 ){
			require(docxtractr)
			wr <- docx.wr1.importer(wrdef)
			w.dat <- MakeWr.docx(t.dat, wr)		

			## Check for duplicated rows	
			if(length(which(duplicated(row.names(t.dat))))>=1){
				dup<-which(duplicated(row.names(t.dat)))
				paste(dup)
				t.dat<-t.dat[-dup,]
				w.dat<-w.dat[-dup,]
			}

		}else{
            wr <- ReadResponseWindowFile(wrdef)
            Wr<-length(wr[,1])#complete and revise this section
			w.dat <- MakeWr(t.dat,wr)
		}
        #####################################################
        #Create Despiked data
        #####################################################
        wts <- tmp.rd$t.dat
        for(i in 1:5) #run the despike 5 times.
        {
            wt.mn3 <- Mean3(wts)
            wts <- SpikeTrim2(wts,1,-1)
            print(sum(is.na(wts))) #this prints out the number of points removed should be close to 0 after 5 loops.
            wts[is.na(wts)] <- wt.mn3[is.na(wts)]
        }
        tmp.rd$mp <- wts
        
        #170127
        # Take the despiked data, subtract the minimum value from the trace, then divide by the maximun value
        # to create traces that are all on the same 0 to 1 scale
        tmp.dat<-tmp.rd$mp

        for(k in 1:length(colnames(tmp.rd$mp))){

            tmp.dat[,k]<-tmp.rd$mp[,k]-min(tmp.rd$mp[,k])
            tmp.dat[,k]<-tmp.dat[,k]/max(tmp.dat[,k])

        }
        tmp.dat[,1]<-tmp.rd$t.dat[,1]

        tmp.rd$mp.1<-tmp.dat



        # Initial Data processing
        levs<-setdiff(unique(as.character(w.dat[,2])),"")
        snr.lim=4;hab.lim=.05;sm=3;ws=30;blc="SNIP"
        
        pcp <- ProcConstPharm(tmp.rd$mp,sm,ws,blc)
        scp <- ScoreConstPharm(tmp.rd,pcp$blc,pcp$snr,pcp$der,snr.lim,hab.lim,sm)
        bin <- bScore(pcp$blc,pcp$snr,snr.lim,hab.lim,levs,tmp.rd$w.dat[,"wr1"])
        bin <- bin[,levs]
        bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
        bin<-pf.function(bin,levs)
        
        tmp.rd$t.dat<-t.dat
        tmp.rd$w.dat<-w.dat
        tmp.rd$c.dat<-c.dat
        tmp.rd$bin<-bin
        tmp.rd$scp<-scp
        tmp.rd$snr<-pcp$snr
        tmp.rd$blc<-pcp$blc
        tmp.rd$der<-pcp$der
        # Add images
        if(!is.null(img1)){tmp.rd$img1<-readPNG(img1)}
        if(!is.null(img2)){tmp.rd$img2<-readPNG(img2)}
        if(!is.null(img3)){tmp.rd$img3<-readPNG(img3)}
        if(!is.null(img4)){tmp.rd$img4<-readPNG(img4)}
        if(!is.null(img5)){tmp.rd$img5<-readPNG(img5)}
        if(!is.null(img6)){tmp.rd$img6<-readPNG(img6)}
        if(!is.null(img7)){tmp.rd$img7<-readPNG(img7)}
        if(!is.null(img8)){tmp.rd$img8<-readPNG(img8)}


    #####################################################
    # Cell Label Scoring	
    #####################################################

        if(fancy==TRUE){tmp.rd<-cell.creator(tmp.rd)}		# Create list of binary  labeled neurons}
        else{tmp.rd$cells<-NULL}
        
        if(is.null(rd.name)){rd.name <- paste("RD",make.names(date()),sep="")}
        
        if(length(which(duplicated(row.names(t.dat))))>=1){
        dup<-which(duplicated(row.names(t.dat)))
        paste(dup)
        t.dat<-t.dat[-dup,]
        w.dat<-w.dat[-dup,]
        }
        
        
        f.name <- paste(rd.name,".Rdata",sep="")
        assign(rd.name,tmp.rd)
        save(list=rd.name,file=f.name)
        return(paste(nrow(tmp.rd$c.dat),"traces read saved to ",f.name))
        #save as RD file
}

# creates despiked trace from t.dat trace
# While using function, save back to RD, or a tmp.rd object
mp.brewer<-function(tmp.rd){

    dat.name<-deparse(substitute(tmp.rd))

    #####################################################
    #Create Despiked data
    #####################################################
    wts <- tmp.rd$t.dat
    for(i in 1:5) #run the despike 5 times.
    {
        wt.mn3 <- Mean3(wts)
        wts <- SpikeTrim2(wts,1,-1)
        print(sum(is.na(wts))) #this prints out the number of points removed should be close to 0 after 5 loops.
        wts[is.na(wts)] <- wt.mn3[is.na(wts)]
    }
    tmp.rd$mp <- wts
    
    #170127
    # Take the despiked data, subtract the minimum value from the trace, then divide by the maximun value
    # to create traces that are all on the same 0 to 1 scale
    tmp.dat<-tmp.rd$mp

    for(k in 1:length(colnames(tmp.rd$mp))){

        tmp.dat[,k]<-tmp.rd$mp[,k]-min(tmp.rd$mp[,k])
        tmp.dat[,k]<-tmp.dat[,k]/max(tmp.dat[,k])

    }
    tmp.dat[,1]<-tmp.rd$t.dat[,1]

    tmp.rd$mp.1<-tmp.dat
    
    #return(tmp.rd)
    assign(dat.name, tmp.rd, envir=.GlobalEnv)
}

#develope cellular binary score and place the binary label of the cells into a cell list called cells
cell.creator<-function(dat, score=F, subset.n=250){		
        if(is.null(subset.n)){subset.n<-250
        }else{subset.n<-subset.n}
        
        if(score){dat<-ROIreview(dat, subset.n=subset.n, pad=5)
        }else{dat<-dat}
        
        levs<-setdiff(unique(as.character(dat$w.dat$wr1)),"")
        cells<-list()
        neuron.response<-select.list(levs, title="What defines Neurons?", multiple=T)
        neurons<-cellzand(dat$bin,neuron.response, 1)
        drop<-cellzand(dat$bin, "drop", 1)
        neurons<-setdiff(neurons,drop)
        pf<-apply(dat$bin[,c("gfp.bin", "tritc.bin")],1,paste, collapse="")
        dat$bin["lab.pf"]<-as.factor(pf)
        #lab.groups<-unique(dat$bin$lab.pf)[-grep(pattern="NA",unique(dat$bin$lab.pf))]
        lab.groups<-as.character(unique(dat$bin$lab.pf))
        cells<-list()
        for(i in lab.groups){
            x.names<-row.names(dat$bin[which(dat$bin[,"lab.pf"]==i, arr.ind=T),])
            cells[[i]]<-x.names
        }
        
        glia.response<-select.list(c(levs, "none"), title="What defines glia?", multiple=T)
        if(glia.response!="none"){
            drop<-cellzand(dat$bin, "drop", 1)
            glia<-cellzand(dat$bin,glia.response, 1)
            glia<-setdiff(glia,drop)
            cells[["000"]]<-setdiff(glia, neurons)
        } 
        else {cells[["000"]]<-setdiff(row.names(dat$c.dat), neurons)}
        dat$cells<-cells
        return(dat)
        
}

ReadResponseWindowFile <- function(fname){
    dat <- read.csv(fname)
    return(dat)
}

#wr file should be
#NEW format for wr file (three column table, treatment, at, duration)

#IMPORT DOCX tables
docx.wr1.importer<-function(file.name='wr1.docx'){
    if( !library(docxtractr, logical.return=T) ){install.packages(docxtractr)}else{}
    #read in docx
    wr1<-read_docx(file.name)
    #Extract each table
    wr1<-docx_extract_all_tbls(wr1, guess_header=F)
    #out table is the third one
    wr1<-Reduce(c,wr1[[3]])
    #split up each vaue based on a single space
    wr1<-strsplit(wr1, ' ')
    
    #Now perform a test and provide a wait if there is an error where the window 
    #region has to little information
    error<-0
    print("There is an >1 <2 info error at")
    for(i in 1:length(wr1)){
        if(length(wr1[[i]])>1 & length(wr1[[i]])<3){
            print(wr1[[i]])
            error<-error+1
        }
    }
    
    print("There is an >3 info error at")
    for(i in 1:length(wr1)){
        if(length(wr1[[i]])>3){
            print(wr1[[i]])
            error<-error+1
        }
    }

    
    print(paste("You have a total of",error,"errors")) 
    if(error>0){
        print("Fix these Errors")
        print("PRESS ANY KEY TO CONTINUE")
        scan(n=1)
        cat("These are your window region definitions
        If you would like to make anymore changes do so now
        ")
        
        for(i in 1:length(wr1)){
            if(length(wr1[[i]])==3){
                print(wr1[[i]])
            }
        }
        print("PRESS ANY KEY TO CONTINUE")
        scan(n=1)
        
        wr1<-read_docx(file.name)
        #Extract each table
        wr1<-docx_extract_all_tbls(wr1, guess_header=F)
        #out table is the third one
        wr1<-Reduce(c,wr1[[3]])
        #split up each vaue based on a single space
        wr1<-strsplit(wr1, ' ')
    }else{}

    wr1.logic<-unlist(lapply(wr1, function(x) length(x)>1 ))
    wr1.locations<-which(wr1.logic, arr.ind=T)
    wr1<-wr1[wr1.locations]
    #wr1<-Reduce(rbind,wr1)
    wr1<-do.call(cbind, lapply(wr1, data.frame, stringsAsFactors=F))
    row.names(wr1)<-c('at','treatment','duration')
    wr1['at',]<-as.numeric(wr1['at',])-(10/60)
    return(wr1)
}

MakeWr.docx <- function(t.dat,wr1,padL=0,padR=0)
{
    w.dat <- t.dat[,1:2]
    names(w.dat)[2] <- "wr1"
    w.dat["wr1"] <- ""
    wr1["treatment",] <- make.names(wr1["treatment",],unique=T)
    for(i in 1:ncol(wr1))
    {
        x1 <- which.min(abs(as.numeric(wr1["at",i])-t.dat[,"Time"]))
        x2 <- which.min(abs(( as.numeric( wr1["at",i] ) + as.numeric( wr1["duration",i] ) )-t.dat[,"Time"]))
        w.dat[max((x1-padL),1):min((x2+padR),nrow(t.dat)),"wr1"] <- as.character(wr1["treatment",i])
    }
    return(w.dat)
}

GetWr <- function(fname)
{
    wr1 <- read.csv(fname)
    return(wr1)
}

MakeWr <- function(t.dat,wr1,padL=0,padR=0)
{
    w.dat <- t.dat[,1:2]
    names(w.dat)[2] <- "wr1"
    w.dat["wr1"] <- ""
    wr1["treatment"] <- make.names(wr1[,"treatment"],unique=T)
    for(i in 1:nrow(wr1))
    {
        x1 <- which.min(abs(wr1[i,"at"]-t.dat[,"Time"]))
        x2 <- which.min(abs((wr1[i,"at"]+wr1[i,"duration"])-t.dat[,"Time"]))
        w.dat[max((x1-padL),1):min((x2+padR),nrow(t.dat)),"wr1"] <- as.character(wr1[i,"treatment"])
    }
    return(w.dat)
}

#fill forward for flevs in the window region.
FillWR <- function(wr1,flevs)
{
    u.names <- unique(wr1)
    wr2 <- NumBlanks(wr1)

    u2.names <- unique(wr2)
    b.names <- grep("blank",u2.names,value=T)
    for(i in flevs)
    {
        for(j in 1:(length(u2.names)-1))
        {
            if(u2.names[j]==i & is.element(u2.names[j+1],b.names) )
            {
                wr1[wr2==u2.names[j+1]] <- i
            }
        }
    }
    return(wr1)
}

#adjust the windows to maximize shift regions and peak regions
#try to minimize the false positive rates but growing/shinking windows
#works reasonably well, but is only counting peaks.  It is not accountinf
#for shape aspects of the trace.
WrAdjust <- function(dat,pcp=NULL,wr=NULL,wr.levs=NULL,snrT=4,minT=10)
{
    gtrfunc <- function(x,a){sum(x>a)}
    if(is.null(wr)){wr <- dat$w.dat[,"wr1"]}
    wr.new <- wr
    wrb <- NumBlanks(wr)
    wi <- 1:(length(wrb))
    x.names <- names(dat$t.dat[,-1])
    if(is.element("bin",names(dat)))
        if(is.element("drop",names(dat$bin)))
        {
            x.names <- row.names(dat$bin[dat$bin[,"drop"]==0,])
        }
    if(is.null(wr.levs))
    {
        wr.levs <- unique(wr)
        wr.levs <- wr.levs[wr.levs != ""]
    }
    if(is.null(pcp))
    {
        pcp <- ProcConstPharm(dat)
    }
    #OK expand/contract each window to give best false positive ratio.
    #keep a min width.
    hits <- apply(pcp$snr[,x.names],1,gtrfunc,a=snrT)
    wrb.levs <- unique(wrb)
    b.levs <- grep("blank",wrb.levs,value=T)
    for(i in wr.levs[wr.levs != wrb.levs[length(wrb.levs)]])
    {
        i1 <- match(i,wrb.levs)
        if(is.element(wrb.levs[i1+1],b.levs))
        {
            targs <- hits[wrb==i | wrb==wrb.levs[i1+1]]
            tval <- NULL
            endT <- length(targs)
            lp <- 0
            for(j in minT:(endT-1))
            {
                lp <- lp+1
                #tval[lp] <- mean(targs[1:j])/((sum(targs[(j+1):endT])+1)/length(targs[(j+1):endT]))
                tval[lp] <- 1/((sum(targs[(j+1):endT])+1)/length(targs[(j+1):endT]))				
            }			
            iopt <- match(i,wr)+which.max(tval)+(minT-1)
        }
        else
        {iopt <- max(wi[wr==i])}
        wr.new[wr==i] <- ""
        wr.new[match(i,wr):iopt] <- i
    }	
    return(wr.new)
}

WrCreate.rdd<-function(t.dat, n=NULL){
    
    window.dat<-data.frame()
    #dev.new(width=10,height=6) 
    x.names<-names(t.dat)[-1]
    LinesSome(t.dat,m.names=x.names,lmain="",subset.n=15)
    ## Plot the total sum of all peaks
    #t.sum<-apply(t.dat[-1], 1, sum)
    #plot(t.dat[,1], t.sum, type="l", lwd=2)
    
    i<-1
    for(i in i:n){
    dose<-locator(n=2, type="o", pch=15, col="red")
    abline(v=c(dose$x[1],dose$x[2]), col="red", lwd=1)
    dose.type<-scan(file="", what="character", n=1, quiet=T)
    duration<-dose$x[2]-dose$x[1]
    window.dat[i,1]<-dose.type
    window.dat[i,2]<-dose$x[1]
    window.dat[i,3]<-duration
    window.dat<-print(window.dat)
    names(window.dat)<-c("treatment", "at", "duration")
    }
    graphics.off()
    write.csv(window.dat, file="wr1.csv", row.names=F)
}

# General Read data dump for an already created RD file without window data
WrCreate.1<-function(dat, n=14, cell=NULL){
    window.dat<-data.frame()
    if(is.null(cell)){cell<-"X.1"}
    else(cell<-cell)
    t.sum<-apply(dat$t.dat[-1], 1, sum)
    dev.new(width=14,height=4) 
    ymax<-max(dat$t.dat[,cell])*1.05
    ymin<-min(dat$t.dat[,cell])*.95
    yrange<-ymax-ymin

    ylim <- c(ymin,ymax)
    xlim <- range(dat$t.dat[,1]) # use same xlim on all plots for better comparison
    
    par(mar=c(6,4.5,3.5,11))
    plot(dat$t.dat[,cell]~dat$t.dat[,1], main=cell,xlim=xlim,ylim=ylim,xlab="", ylab="",pch=16, lwd=1, cex=.5)
    #axis(1, at=seq(0, length(dat$t.dat[,1]), 5),tick=TRUE )  

    

    for(i in 1:n){
        dose<-locator(n=2, type="o", pch=15, col="red")
        abline(v=c(dose$x[1],dose$x[2]), col="red", lwd=1)
        dose.type<-scan(file="", what="character", n=1, quiet=T)
        duration<-dose$x[2]-dose$x[1]
        window.dat[i,1]<-dose.type
        window.dat[i,2]<-dose$x[1]
        window.dat[i,3]<-duration
        window.dat<-print(window.dat)
        names(window.dat)<-c("treatment", "at", "duration")
        wr1<-window.dat
    }
    t.dat<-return(MakeWr(dat$t.dat,wr1,padL=0,padR=0))
}

WrMultiplex<-function(t.dat, wr, n=NULL){
    w.dat<-t.dat[,1:2]
    names(w.dat)[2]<-"wr1"
    w.dat["wr1"]<-""
    if(is.null(n)){n=length(wr[,1])}
    library(cluster)
    pamk<-pam(w.dat[,1], k=n)
    wr[1] <- make.names(wr[,1],unique=T)
    levs<-wr[,1]
    w.dat[,"wr1"]<-levs[pamk$clustering]
    return(w.dat)}


#Made a mistake in you window regions?
#If it is only the time sequence, but all other info is corrected
#then make complete F.  This will allow you to select the windows that need reapri
#if the naming is off, then make complete=F.  You will need to do a complete reapri
# you will lose all information from RDView
WindowRepair<-function(dat, complete=T){

    tmp<-dat #first create a tmp to repair

    tmp.rd<-dat # then create a tmp.rd to completely screwup for repairs

    #Now do all of this to tmp.rd
    wrdef<-"wr1.csv"
    t.dat<-tmp.rd$t.dat
    if(!is.null(wrdef))
        {
            wr <- ReadResponseWindowFile(wrdef)
            Wr<-length(wr[,1])#complete and revise this section
            if(length(colnames(wr))<2){w.dat<-WrMultiplex(t.dat,wr,n=Wr)}
            else{w.dat <- MakeWr(t.dat,wr)}
            tmp.rd$w.dat<-w.dat
        }
    levs<-setdiff(unique(as.character(tmp.rd$w.dat$wr1)),"")

    #5 set the thresholds for scoring and run the automatic scoring
    sm <- 2 #smooth window size set
    ws <- 30 #window peak size
    snr.lim <- 4 #signal to noise threshold
    hab.lim <- .05 #height above baseline threshold
    blc="SNIP"

    pcp <- ProcConstPharm(tmp.rd,sm,ws,blc)
    scp <- ScoreConstPharm(tmp.rd,pcp$blc,pcp$snr,pcp$der,snr.lim,hab.lim,sm)
    bin <- bScore(pcp$blc,pcp$snr,snr.lim,hab.lim,levs,tmp.rd$w.dat[,"wr1"])
    bin <- bin[,levs]
    bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
    bin<-pf.function(bin,levs)

    tmp.rd$bin<-bin
    tmp.rd$blc<-pcp$blc
    tmp.rd$snr<-pcp$snr
    tmp.rd$der<-pcp$der
    tmp.rd$bin<-bin
    tmp.rd$scp<-scp

    #Now surgically add the selected corrected data from tmp.rd to the tmp
    #starting with the window region
    tmp$w.dat$wr1<-tmp.rd$w.dat$wr1

    #select the window region you want to repair
    if(complete==T){
        tmp$scp<-tmp.rd$scp
        tmp$bin<-tmp.rd$bin
    }else{
        print("Select windows to repair")
        windows.tp<-select.list(names(tmp.rd$bin), multiple=T)
    
    if(length(windows.tp)>1){
        #next add the binary information
        for(i in windows.tp){
            tmp$bin[i]<-tmp.rd$bin[i]
            win.stats<-grep(i,names(tmp$scp), value=T)
            tmp$scp[win.stats]<-tmp.rd$scp[win.stats]
        }
    }else{}
    }
    #now save back to the RD object
    dat<-tmp
    return(dat)
}

#Made a mistake in you window regions?
#If it is only the time sequence, but all other info is corrected
#then make complete F.  This will allow you to select the windows that need reapri
#if the naming is off, then make complete=F.  You will need to do a complete reapri
# you will lose all information from RDView

WindowRepair_docx<-function(dat, complete=F, trace_brew=F){
    require(docxtractr)
    tmp<-dat #first create a tmp to repair
    tmp.rd<-dat # then create a tmp.rd to completely screwup for repairs

    #Now do all of this to tmp.rd
    wrdef<-"wr1.docx"
    t.dat<-tmp.rd$t.dat
    if(!is.null(wrdef)){
        wr <- docx.wr1.importer(wrdef)
        Wr<-length(wr[,1])#complete and revise this section
        if(length(colnames(wr))<2){
            w.dat<-WrMultiplex(t.dat,wr,n=Wr)
        }else{
            w.dat <- MakeWr.docx(t.dat,wr)
        }
        tmp.rd$w.dat<-w.dat
    }
    levs<-setdiff(unique(as.character(tmp.rd$w.dat$wr1)),"")

    #5 set the thresholds for scoring and run the automatic scoring
    sm <- 2 #smooth window size set
    ws <- 30 #window peak size
    snr.lim <- 4 #signal to noise threshold
    hab.lim <- .05 #height above baseline threshold
    blc="SNIP"
    
    if(length(summary(names(tmp.rd)=='t.norm')) == 2 | trace_brew ){
        tmp.rd <- TraceBrewer(tmp.rd)
    }
    pcp <- ProcConstPharm(tmp.rd,sm,ws,blc)
    scp <- ScoreConstPharm(tmp.rd,pcp$blc,pcp$snr,pcp$der,snr.lim,hab.lim,sm)
    bin <- bScore(pcp$blc,pcp$snr,snr.lim,hab.lim,levs,tmp.rd$w.dat[,"wr1"])
    bin <- bin[,levs]
    bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
    bin<-pf.function(bin,levs)

    tmp.rd$bin<-bin
    tmp.rd$blc<-pcp$blc
    tmp.rd$snr<-pcp$snr
    tmp.rd$der<-pcp$der
    tmp.rd$bin<-bin
    tmp.rd$scp<-scp

    #Now surgically add the selected corrected data from tmp.rd to the tmp
    #starting with the window region
    tmp$w.dat$wr1<-tmp.rd$w.dat$wr1

    #select the window region you want to repair
    if(complete==T){
        tmp$scp<-tmp.rd$scp
        tmp$bin<-tmp.rd$bin
    }else{
        print("Select windows to repair")
        windows.tp<-select.list(names(tmp.rd$bin), multiple=T)
    
        if(length(windows.tp)>1){
            #next add the binary information
            for(i in windows.tp){
                tmp$bin[i]<-tmp.rd$bin[i]
                win.stats<-grep(i,names(tmp$scp), value=T)
                tmp$scp[win.stats]<-tmp.rd$scp[win.stats]
            }
        }else{}
    }
    #now save back to the RD object
    dat<-tmp
    return(dat)
}

#Something happened within you experiment and you need to delete a region of row fro the traces
#to recover the data and clean up the display.
#dat: This is the RD object
#cell: Cell the display on the plot, if left empty, X.1 will be used
#complete: logical, If you say true, all window regions will be reassessed. If False window region can be selected.

RegionDeleter<-function(dat, cell=NULL, complete=T){

    if(is.null(cell)){cell<-"X.1"}
    tmp<-dat #first create a tmp to repair

    tmp.rd<-dat # then create a tmp.rd to completely screwup for repairs

    #Now do all of this to tmp.rd
    dev.new(width=16, height=4)
    PeakFunc7(tmp.rd, cell)
    x.del<-locator(n=2, type="o", col="red", pch=4, lwd=2)$x
    
    rows.to.remove<-which(tmp.rd$t.dat[,1]>=x.del[1] & tmp.rd$t.dat[,1]<=x.del[2],arr.ind=T)
    
    tmp.rd$t.dat<-tmp.rd$t.dat<-dat$t.dat[-rows.to.remove,]
    tmp.rd$w.dat<-tmp.rd$w.dat<-dat$w.dat[-rows.to.remove,]

    tmp.rd<-mp.brewer(tmp.rd)
    levs<-setdiff(unique(as.character(tmp.rd$w.dat$wr1)),"")

    #5 set the thresholds for scoring and run the automatic scoring
    sm <- 2 #smooth window size set
    ws <- 30 #window peak size
    snr.lim <- 4 #signal to noise threshold
    hab.lim <- .05 #height above baseline threshold
    blc="SNIP"

    pcp <- ProcConstPharm(tmp.rd,sm,ws,blc)
    scp <- ScoreConstPharm(tmp.rd,pcp$blc,pcp$snr,pcp$der,snr.lim,hab.lim,sm)
    bin <- bScore(pcp$blc,pcp$snr,snr.lim,hab.lim,levs,tmp.rd$w.dat[,"wr1"])
    bin <- bin[,levs]
    bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
    bin<-pf.function(bin,levs)

    tmp.rd$bin<-bin
    tmp.rd$blc<-pcp$blc
    tmp.rd$snr<-pcp$snr
    tmp.rd$der<-pcp$der
    tmp.rd$bin<-bin
    tmp.rd$scp<-scp

    #Now surgically add the selected corrected data from tmp.rd to the tmp
    #starting with the window region
    #tmp$w.dat$wr1<-tmp.rd$w.dat$wr1

    #select the window region you want to repair
    if(complete==T){
        tmp<-tmp.rd
    }else{
        tmp$t.dat<-tmp.rd$t.dat
        tmp$w.dat<-tmp.rd$w.dat
        tmp$blc<-tmp.rd$blc
        tmp$snr<-tmp.rd$snr
        tmp$mp<-tmp.rd$mp
        tmp$mp.1<-tmp.rd$mp.1
        
        print("Select windows to repair")
        windows.tp<-select.list(names(tmp.rd$bin), multiple=T)
    
        if(length(windows.tp)>1){
            #next add the binary information
            for(i in windows.tp){

            tmp$bin[i]<-tmp.rd$bin[i]
            win.stats<-grep(i,names(tmp$scp), value=T)
            tmp$scp[win.stats]<-tmp.rd$scp[win.stats]
        }
    }else{}
    }
    #now save back to the RD object
    dat<-tmp
    return(dat)
}
    
#This is a function to rename mislabeled pulses
WindowRenamer<-function(dat){
    dat.name<-deparse(substitute(dat))

    pulsenames<-setdiff(unique(as.character(dat$w.dat$wr1)),"")
    pulse<-select.list(pulsenames,title="Pulse To Rename", multiple=T)
    pulze<-pulse
    
    bringToTop(-1)
    print("#############These are your pulses###############")
    print(pulsenames)
    print("#############This is the pulse to rename:")
    print(pulse)
    print("#############Enter the new name")
    pulse.rn<-scan(n=length(pulse),what="character")
    pulse.rn<-make.names(pulse.rn, unique=T)
    print(pulse.rn)
    
    for(i in length(pulse):1){
        #Rename Bin dataframe
        colnames(dat$bin)[colnames(dat$bin)==pulse[i]] <- pulse.rn[i]
        #Rename w.dat
        dat$w.dat[which(dat$w.dat[,"wr1"]==pulse[i], arr.ind=T),"wr1"]<-pulse.rn[i]
        #Rename SCP
        pulze[i]<-paste(pulse[i],".", sep="")
        scp.col.to.rn<-grep(pulze[i],colnames(dat$scp),fixed=T)
        names(dat$scp)<-sub(pulze[i],paste(pulse.rn[i],".", sep=""), names(dat$scp))
    }
    #assign(dat.name,dat, envir=.GlobalEnv)
    return(dat)
}

##############################################################################################
##############################################################################################
#take a defined window vector and
#number of the contiguos blank regions ("")
NumBlanks <- function(x){
    nw <-  as.character(x)
    xlen <- length(x)
    bl.cnt <- 1
    mi <- match("",nw)
    while(!is.na(mi) & bl.cnt < 20)
    {
        mi2 <- mi+1
        while((x[mi2] == "") & mi2 <= xlen){mi2 <- mi2+1}
        nw[mi:(mi2-1)] <- paste("blank",bl.cnt,sep="")
        bl.cnt <- bl.cnt+1
        mi <- match("",nw)
    }
    return(as.factor(nw))
}

NumAny <- function(x,targ){
    nw <-  as.character(x)
    xlen <- length(x)
    bl.cnt <- 1
    mi <- match(targ,nw)
    while(!is.na(mi) & bl.cnt < 20)
    {
        mi2 <- mi+1
        while((x[mi2] == targ) & mi2 <= xlen){mi2 <- mi2+1}
        nw[mi:(mi2-1)] <- paste(targ,bl.cnt,sep="")
        bl.cnt <- bl.cnt+1
        mi <- match(targ,nw)
    }
    return(as.factor(nw))	
}

#measure the variance about the smooth
SmoothVar <- function(x,shws=2){
    s1 <- createMassSpectrum(1:length(x),x)
    s3 <- smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=shws)
    return(sd(x-intensity(s3)))
}

#########################################################################
#Trace Spike removal, smoothing, normalizing
##########################################################################
#adds an mp data.frame to the tmp RD object.
#RD object must have t.dat and w.dat[,"wr1"] 
#this will have the data for the Despiked and smoothed traces.
#ulim and dlim define the limits of two point rise and fall.
#if they are not set the function will try to estimate them from 
#the blank regions of the traces, RD must have blank regions (e.g. "")
#currently only the negative spikes are eliminated.  That is points
#that have large increases followed by large decreases (or vice versa)
#smoothing is done using 13 points.
smoothfunc <- function(y,ts,pts=min(13,sum(!is.na(y)))){
    yloe <- loess(y ~ ts,span=pts/length(y))
    ymp <- predict(yloe,newdata=data.frame(ts=ts))
    return(ymp)
}

DespikeSmooth <- function(tmp,ulim=NULL,dlim= NULL){
    print("Starting Despike-smooth timing")
    start_time <- Sys.time()
    wt <- tmp$t.dat
    ts <- tmp$t.dat[,1]
    wtd <- wt[-1,] - wt[-nrow(wt),]
    print(paste("step 1 time",Sys.time()-start_time))
        wtd <- sweep(wtd[,-1],1,wtd[,1],'/')
    print(paste("step 2 time", Sys.time()-start_time))
        wtm <- wtd[-1,]*wtd[-nrow(wtd),]
    print(paste("step 3 time", Sys.time()-start_time))
        wrb <- NumBlanks(tmp$w.dat[,"wr1"])
    print(paste("step 4 time", Sys.time()-start_time)	)
        if(is.null(ulim) | is.null(dlim))
        {
            qvals <- quantile(
                        as.vector(
                            as.matrix(
                                wtm[grep("blank",wrb[1:nrow(wtm)]),]
                            )
                        )
                    ,probs=c(0,.001,.5,.999,1))
        }
    print(paste("step 5 time",Sys.time()-start_time))
        if(is.null(dlim)){dlim <- qvals[2]}
        if(is.null(ulim)){ulim <- qvals[4]}	
        wtm.z <- wtm[1,]
        wtm.z[,] <- 0
        wtm <- rbind(wtm.z,wtm,wtm.z)
        wtrm <- wtm < dlim
        x <- wt[,1]
        wt <- wt[,-1]
        wt[wtrm] <- NA
        mp <- sapply(wt, smoothfunc, ts=x)
    print(paste("step 6 time",Sys.time()-start_time))
        tmp$mp <- tmp$t.dat
        tmp$mp[,-1] <- mp
        return(tmp)
}

LinearBLextend <- function(x,y,intvl=3,plotit=F){
    require(MASS)
    #break into intervals of 3 minute.
    time.tot <- max(x)-min(x)
    gaps <- ceiling(time.tot/intvl)
    gap.fac <- sort(rep(seq(1,gaps),length.out=length(x)))
    gap.i <- tapply(y,gap.fac,which.min)
    gap.i <- gap.i + match(unique(gap.fac),gap.fac)-1
    mindat <- data.frame(x=x[gap.i],y=y[gap.i])
    rlt <- rlm(y ~ x,data=mindat)	
    if(plotit)
    {
        plot(x,y)
        points(x[gap.i],y[gap.i],pch=16,col="red")
        abline(rlt,lty=2,col="blue",lwd=2)
        #lines(x[gap.i],predict(xloe))
        #points(x1,predict(rlt,newdata=data.frame(x=x1)),pch=16,col="blue",cex=2)
    }
    return(coefficients(rlt))
}

#add points to the end of the t.dat
PadTdat <- function(tmp,n=5){
    xdat <- tmp$t.dat
    r1 <- sapply(tmp$t.dat[,-1],LinearBLextend,x=tmp$t.dat[,1],plotit=F)
    r1 <- data.frame(t(r1))
    tmp$scp[row.names(r1),"rlm.a"] <- r1[,1]
    tmp$scp[row.names(r1),"rlm.b"] <- r1[,2]
    x <- xdat[,1]
    dx <- x[-1]-x[-length(x)]
    tmax <- max(x)
    tseq <- seq(1,n)*median(dx)+tmax
    r2 <- data.frame(t(r1[,"x"] %*% t(tseq) + r1[,1]))
    names(r2) <- row.names(r1)
    r2[,"Time"] <- tseq
    xdat <- rbind(xdat,r2[,names(xdat)])
    w1 <- tmp$w.dat[1:n,]
    w1[,"Time"] <- tseq	
    for(i in names(w1)[sapply(w1,is.character)]){w1[,i] <- "epad"}
    tmp$w.dat <- rbind(tmp$w.dat,w1)
    row.names(tmp$w.dat)<-tmp$w.dat$Time #Lee's additions
    tmp$t.dat <- xdat
    row.names(tmp$t.dat)<-tmp$t.dat$Time #Lee's additions
    tmp$bin['epad']<-0 #Lee's additions
    return(tmp)
}

TraceNormal<-function(dat, t.type='blc'){
    tmp.dat<-dat[[t.type]]

    for(k in 1:length(colnames(dat$mp))){

        tmp.dat[,k]<-dat$mp[,k]-min(dat$mp[,k])
        tmp.dat[,k]<-tmp.dat[,k]/max(tmp.dat[,k])

    }
    tmp.dat[,1]<-dat$t.dat[,1]

    dat$t.norm<-tmp.dat
    return(dat)
}

TraceBrewer<-function(dat){
    cat('
    #The current flow of our trace cleaning protocol is as follows, and this is
    #what the function automatically fills in for the RD list

    #t.dat>                 Raw data                   t.dat
    #t.dat.pad>             3 End points added at end  NA
    #t.dat.pad.ds.s>        despike and smooth         mp
    #t.dat.pad.ds.s.n>      normalize 0 to 1           t.norm
    #t.dat.pad.ds.s.n.blc>  Baseline Corrected         blc
    ')
    tmp.rd<-dat
    start.time<-proc.time()
    # Kevin has created a new way of creating cleaned up traces.
    # Add a 3 point padding to the end of the experiment which return the trace back to baseline
    # This helps preserve the response shape
    tmp.rd<-PadTdat(tmp.rd)
    print(paste("Completed Padding at:",(proc.time()-start.time))[3])
    # Kevin now uses this to despike and smooth the data
    tmp.rd<-DespikeSmooth(tmp.rd)
    print(paste("Completed Despiking at:",(proc.time()-start.time)[3]))
    # Now what we need to do is provide the analysis with some type of normalized trace
    tmp.rd<-TraceNormal(tmp.rd,'mp')
    print(paste("Completed Normalizing at:",(proc.time()-start.time)[3]))


    # Now do a baseline correction based on the trace padding, and the despike smooth function, and
    # the normalized trace 
    pcp.tmp<-ProcConstPharm(tmp.rd$t.norm)
    print(paste("Completed Baseline Correction at:",(proc.time()-start.time)[3]))
    tmp.rd$blc<-pcp.tmp$blc
    tmp.rd$snr<-pcp.tmp$snr
    # Now perform trace statistics on the specified trace
    tmp.rd$scp<-ScoreConstPharm.2(tmp.rd,'blc')
    print(paste("Completed Window Statistics at:",(proc.time()-start.time)[3]))
    return(tmp.rd)
}
##############################################################################################
# Cornerstones of trace washing, peak detection, and binary scoring
##############################################################################################

#the first argument is the raw data
#the second argument is the halfwindow size for smoothing (shws)
#the third argument is the peak detection halfwindow size (phws)
#the last argument is the baseline correction method (TopHat = blue line SNIP = red line)
#Note that you should use the RoughReview function to determine the best values for
#arguments 2,3 and 4.

#returns a list with two dataframes: snr and blc.
#snr has the peaks detected for all cells, blc has the baseline corrected data for all cells. 
SpikeTrim2 <- function(wt,ulim=NULL,dlim= NULL){
    
    wtd <- wt[-1,]-wt[-nrow(wt),]
    wtd <- sweep(wtd[,-1],1,wtd[,1],'/')
    if(is.null(ulim) | is.null(dlim))
    {
        qvals <- quantile(as.vector(as.matrix(wtd)),probs=c(0,.01,.5,.99,1))
    }
    if(is.null(dlim)){dlim <- qvals[2]}
    if(is.null(ulim)){ulim <- qvals[4]}	
    wt.up <- wtd > ulim
    wt.dn <- wtd < dlim
    wt.ud <- wt.up[-nrow(wt.up),] + wt.dn[-1,]
    wt.du <- wt.up[-1,] + wt.dn[-nrow(wt.dn),]
    wt.na <- wt[2:(nrow(wt)-1),-1]
    wt.na[wt.ud==2] <- NA
    wt.na[wt.du==2] <- NA	
    sum(is.na(wt.na))
    wt[2:(nrow(wt)-1),-1] <- wt.na

    #impute missing using mean of flanking.
    #consider replicating first and last columns and doing this all as a vector
    
    return(wt)
}
#each point is replaced with the mean of the two neighboring points
Mean3 <- function(wt){
    wt.mn <- (wt[-c(1,2),]+wt[-c(nrow(wt),(nrow(wt)-1)),])/2
    wt[2:(nrow(wt)-1),] <- wt.mn
    return(wt)
}

ProcConstPharm <- function(dat,shws=2,phws=20,bl.meth="SNIP"){
    if(class(dat)=="data.frame"){(dat1<-dat)}else{dat1 <- dat$t.dat}
    t.names <- names(dat1)[-1]#Time in first column
    dat1.snr <- dat1 #peak calls stored as SNR
    dat1.snr[,t.names] <- 0
    dat1.bc <- dat1.snr #baseline corrected data

    for(i in t.names)
    {
        p1 <- PeakFunc2(dat1,i,shws=shws,phws=phws,Plotit=F,bl.meth=bl.meth)
        dat1.snr[match(mass(p1$peaks),dat1[,1]),i] <- snr(p1$peaks)
        dat1.bc[i] <- intensity(p1$dat)
    }
    dat1.der<-dat1.bc[-1,]-dat1.bc[-nrow(dat1.bc),]
    dat1.der <- sweep(dat1.der[,-1],1,dat1.der[,1],'/')

    #    dat1.crr <- allCRR(dat1,t.names,Plotit=F) #leave off advanced processing for now
    return(list(snr=dat1.snr,blc=dat1.bc, der=dat1.der))
}
#binary score for all cells for the regions of interest bScore
#argument 1 is the baseline corrected data
#argument 2 is the snr peak data
#argument 3 is the threshold for significance on the peaks
#argument 4 is the intensity above baseline theshold
#argument 5 indicates the regions of interest. (e.g. the response windows for which the cells will be scored)
#argument 6 indicates the response windows. 
#argument 7 indicates the cells to score (if null all cells will be scored)
#returns the scoring for all cells subject to the above parameters.
#as well as the sum for the snr scores and the sd for the snr scores.
bScore <- function(blc,snr,snr.lim,blc.lim,levs,wr,cnames=NULL){
    notzero <- function(x){as.integer(sum(x) > 0)}
    if(is.null(cnames)){cnames <- names(blc)[-1]}
    wr2 <- wr[is.element(wr,levs)]
    b.snr <- snr[is.element(wr,levs),cnames]
    b.blc <- blc[is.element(wr,levs),cnames]
    b.call <- b.blc
    b.call[,] <- 0
    b.call[b.snr > snr.lim & b.blc > blc.lim] <- 1
    b.score <- data.frame(tot=apply(b.snr,2,sum))
    b.score["sd"] <- apply(b.snr,2,sd)
    for(i in levs)
    {
        b.score[i] <- apply(b.call[wr2==i,],2,notzero)
    }
    return(b.score)
}

# Binary scoring dependent upon score const pharm talbe values
# Best way to determine parameters is to look through trace click before hand
# snr.min = minimun signal to noise value
# max.min= minimun above baseline threshold
# tot.min= area minimun to consider
# wm.min= which max, Where within the window region does the maximun value occur
# wm.max= where to stop looking for the maximun value
bscore2<-function(dat, levs.1=NULL, snr.min=2.8, max.min=.03, wm.min=0, wm.max=600){
    scp<-dat$scp
    levs<-setdiff(unique(as.character(dat$w.dat[,2])),"")
    if(is.null(levs.1)){levs.1<-levs}
    else{levs.1<-levs.1}
    #dat2<-matrix(0, nrow=length(dat$c.dat[,1]), ncol=length(levs))
    dat2<-dat$bin[levs]
    #row.names(dat2)<-dat$c.dat[,1]
    #colnames(dat2)<-levs
    x.names<-dat$c.dat[,1]
    for(j in x.names){	
        for(i in levs.1){
            snr.name<-grep(paste(i,".snr", sep=""), names(dat$scp), value=T)
            tot.name<-grep(paste(i,".tot", sep=""), names(dat$scp), value=T)
            max.name<-grep(paste(i,".max", sep=""), names(dat$scp), value=T)
            wm.name<-grep(paste(i,".wm", sep=""), names(dat$scp), value=T)
            
            if(dat$scp[j,snr.name]>=snr.min &
                dat$scp[j,max.name]>=max.min &
                dat$scp[j,wm.name]>=wm.min &
                dat$scp[j,wm.name]<=wm.max)
            {dat2[j,i]<-1}
            else{dat2[j,i]<-0}
            }
            }
            return(dat2)}

    # calculate a table of cell characteristics globally and 
    # within specific windows
    # these specifics should include
    # mean and sd, sum of in window peaks, sum of out of window peaks
    # 1)  some measure of dead cell
    # 2)  yes/no peak response for each window
    # 3) peak height
    # 4) max peak SNR
    # 5) peak timing in window
    # 6)
    # variance of smoothed - raw in window
    # define and number blank windows.
    ScoreConstPharm <- function(dat,blc=NULL, snr=NULL, der=NULL, snr.lim=3,blc.lim=.03,shws=2)
    {
    t.dat<-dat$t.dat
    if(is.null(blc)){blc<-dat$blc
    }else{blc<-blc}
    if(is.null(snr)){snr<-dat$snr
    }else{snr<-snr}
    if(is.null(der)){der<-dat$der
    }else{der<-der}


    wr<-dat$w.dat$wr1

        gtfunc <- function(x,alph){sum(x > alph,na.rm=T)}
        
    lt5func <- function(x,y)
    {
        ltfunc <- function(i){summary(lm(y[i:(i+5)] ~ x[i:(i+5)]))$coefficients[2,3]}
        iseq <- 1:(length(x)-5)
        res <- sapply(iseq,ltfunc)
        return(range(res))
    }

        levs <- setdiff(unique(wr),"")
        cnames <- names(t.dat)[-1]
        res.tab <- data.frame(mean=apply(blc[,cnames],2,mean))
        res.tab["sd"] <- apply(blc[,cnames],2,sd)
        res.tab["snr.iws"] <- apply(snr[is.element(wr,levs),cnames],2,sum)
        res.tab["snr.ows"] <- apply(snr[!is.element(wr,levs),cnames],2,sum)
        res.tab["snr.iwc"] <- apply(snr[is.element(wr,levs),cnames],2,gtfunc,alph=snr.lim)
        res.tab["snr.owc"] <- apply(snr[!is.element(wr,levs),cnames],2,gtfunc,alph=snr.lim)

        dat.der<-der
        
        for(i in cnames)
        {
            s1 <- createMassSpectrum(t.dat[,"Time"],t.dat[,i])
            s3 <- smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=shws)
            bl.th <- estimateBaseline(s3, method="TopHat")[,"intensity"]
            bl.snp <- estimateBaseline(s3, method="SNIP")[,"intensity"]
            eseq <- 1:ceiling((nrow(t.dat)/2))
            lseq <- max(eseq):nrow(t.dat)
            res.tab[i,"bl.diff"] <- mean(bl.th-bl.snp)
            res.tab[i,"earl.bl.diff"] <- mean(bl.th[eseq]-bl.snp[eseq])
            res.tab[i,"late.bl.diff"] <- mean(bl.th[lseq]-bl.snp[lseq])        
        }
        for(i in levs)
        {
            res.tab[paste(i,".snr",sep="")] <- apply(snr[wr==i,cnames],2,max)
            res.tab[paste(i,".tot",sep="")] <- apply(blc[wr==i,cnames],2,sum)
            res.tab[paste(i,".max",sep="")] <- apply(blc[wr==i,cnames],2,max)
            res.tab[paste(i,".ph.a.r",sep="")] <-res.tab[paste(i,".tot",sep="")]/res.tab[paste(i,".max",sep="")]

            res.tab[paste(i,".wm",sep="")] <- apply(blc[wr==i,cnames],2,which.max)
            
            ## Derviative measures
            #res.tab[paste(i,".der.tot",sep="")] <- apply(dat.der[wr==i,cnames],2,sum)
            res.tab[paste(i,".der.tot",sep="")] <- apply(dat.der[wr==i,cnames],2,sum)
            #res.tab[paste(i,".der.tot",sep="")] <- apply(na.omit(dat.der[wr==i,cnames]),2,function(x){sum(x[x>0])})
            res.tab[paste(i,".der.max",sep="")] <- apply(na.omit(dat.der[wr==i,cnames]),2,max)
            res.tab[paste(i,".der.min",sep="")] <- apply(na.omit(dat.der[wr==i,cnames]),2,min)
            res.tab[paste(i,".der.wmax",sep="")] <- apply(na.omit(dat.der[wr==i,cnames]),2,which.max)#function(x){which.max(x[5:length(row.names(x))])})
            res.tab[paste(i,".der.wmin",sep="")] <- apply(na.omit(dat.der[wr==i,cnames]),2,which.min)


            
    #        res.tab[c(paste(i,".dn5",sep=""),paste(i,".up5",sep=""))] <- t(apply(t.dat[wr==i,cnames],2,lt5func,x=t.dat[wr==i,1]))
    #        res.tab[paste(i,".dn5",sep="")] <- apply(blc[wr==i,cnames],2,dn5func)                
        }
        return(res.tab)
}

ScoreConstPharm.2 <- function(dat,t.type=NULL, snr=NULL, der=NULL, snr.lim=3,blc.lim=.03,shws=2){
    require(MALDIquant)
    t.dat<-dat$t.dat
    if(is.null(t.type)){t.type<-'blc'
    }else{t.type<-t.type}
    if(is.null(snr)){snr<-dat$snr
    }else{snr<-snr}
    if(is.null(der)){der<-dat$der
    }else{der<-der}


    wr<-dat$w.dat$wr1

        gtfunc <- function(x,alph){sum(x > alph,na.rm=T)}
        
    lt5func <- function(x,y)
    {
        ltfunc <- function(i){summary(lm(y[i:(i+5)] ~ x[i:(i+5)]))$coefficients[2,3]}
        iseq <- 1:(length(x)-5)
        res <- sapply(iseq,ltfunc)
        return(range(res))
    }

    levs <- setdiff(unique(wr),"")
    cnames <- names(t.dat)[-1]
    res.tab <- data.frame(mean=apply(dat[['blc']][,cnames],2,mean))
    res.tab["sd"] <- apply(dat$blc[,cnames],2,sd)
    res.tab["snr.iws"] <- apply(snr[is.element(wr,levs),cnames],2,sum)
    res.tab["snr.ows"] <- apply(snr[!is.element(wr,levs),cnames],2,sum)
    res.tab["snr.iwc"] <- apply(snr[is.element(wr,levs),cnames],2,gtfunc,alph=snr.lim)
    res.tab["snr.owc"] <- apply(snr[!is.element(wr,levs),cnames],2,gtfunc,alph=snr.lim)

    dat.der<-der
    
    for(i in cnames)
    {
        s1 <- createMassSpectrum(t.dat[,"Time"],t.dat[,i])
        s3 <- smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=shws)
        bl.th <- estimateBaseline(s3, method="TopHat")[,"intensity"]
        bl.snp <- estimateBaseline(s3, method="SNIP")[,"intensity"]
        eseq <- 1:ceiling((nrow(t.dat)/2))
        lseq <- max(eseq):nrow(t.dat)
        res.tab[i,"bl.diff"] <- mean(bl.th-bl.snp)
        res.tab[i,"earl.bl.diff"] <- mean(bl.th[eseq]-bl.snp[eseq])
        res.tab[i,"late.bl.diff"] <- mean(bl.th[lseq]-bl.snp[lseq])        
    }
    for(i in levs)
    {
        res.tab[paste(i,".snr",sep="")] <- apply(snr[wr==i,cnames],2,max)
        res.tab[paste(i,".tot",sep="")] <- apply(dat[[t.type]][wr==i,cnames],2,sum)
        res.tab[paste(i,".max",sep="")] <- apply(dat[[t.type]][wr==i,cnames],2,max)
        res.tab[paste(i,".ph.a.r",sep="")] <-res.tab[paste(i,".tot",sep="")]/res.tab[paste(i,".max",sep="")]

        res.tab[paste(i,".wm",sep="")] <- apply(dat[[t.type]][wr==i,cnames],2,which.max)
    }
    return(res.tab)
}

##############################################################################################
##############################################################################################


##############################################################################################
# Response Scoring
##############################################################################################
#should probably break this into ScoreMulti and ReviewMulti
#Score all RD...Rdata files in a given directory with review
#check for an existing bin file and just review that.
#add a "drop" column to the bin file
# Needs work on drop cells
ScoreMulti <- function(dir.name=NULL,snr.lim=4,hab.lim=.05,sm=3,ws=30,review=T){
    if(is.null(dir.name)){dir.name <- getwd()}
    setwd(dir.name)
    f.names <- list.files(pattern="RD.*\\.Rdata$")
    if(length(f.names) == 0){stop("no RD...Rdata files in given directory")}
    rd.list <- sub("\\.Rdata*","",f.names)
    RD.names <- rd.list #paste(rd.list,".b",sep="")
    RD.f.names <- paste(RD.names,".Rdata",sep="")
    sel.i <- menu(rd.list,title="Select Data to review")			
    while(sel.i != 0)
    {

        j <- sel.i
        load(f.names[j])
        i <- rd.list[j]
        tmp <- get(i)
        tlevs <- c(as.character(unique(tmp$w.dat[,"wr1"])[-1]),"drop")
        if(is.null(tmp$bin))
        {
        tmp.pcp <- ProcConstPharm(tmp,sm,ws,"TopHat")
        tmp.scp <- ScoreConstPharm(tmp$t.dat,tmp.pcp$blc,tmp.pcp$snr,snr.lim,hab.lim,tmp$w.dat[,"wr1"],sm)
        tmp.bin <- bScore(tmp.pcp$blc,tmp.pcp$snr,snr.lim,hab.lim,tlevs,tmp$w.dat[,"wr1"])
        tmp.bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
        }
        else
        {
            tmp.bin <- tmp$bin
            tmp.scp <- tmp$scp
            tmp.blc <- tmp$blc
        }
        if(review)
        {
        tmp.bin <- ScoreReview1(tmp$t.dat,tmp.bin[,tlevs],tmp$w.dat[,"wr1"])
        tmp.bin <- ScoreReview0(tmp$t.dat,tmp.bin[,tlevs],tmp$w.dat[,"wr1"])
        }
        
        tmp$bin <- tmp.bin[,tlevs]
        pf<-apply(tmp$bin[,tlevs],1,paste,collapse="")	
        pf.sum<-summary(as.factor(pf),maxsum=500)
        pf.sum<-pf.sum[order(pf.sum,decreasing=T)]
        pf.ord<-pf.sum
        pf.ord[]<-seq(1,length(pf.sum))
        tmp$c.dat["pf"]<-as.factor(pf)
        tmp$c.dat["pf.sum"]<-pf.sum[pf]
        tmp$c.dat["pf.ord"]<-pf.ord[pf]
        tmp$c.dat<-cbind(tmp$c.dat, tmp$bin)
        
        
        tmp$scp <- tmp.scp
        tmp$snr<-tmp.pcp$snr
        tmp$blc <- tmp.pcp$blc
        assign(RD.names[j],tmp)
        save(list=RD.names[j],file=RD.f.names[j])
        print(paste("DONE REVIEWING ",RD.names[j]," CHANGES SAVED TO FILE.",sep=""))
        sel.i <- menu(rd.list,title="Select Data to review")			
    }
    return(RD.f.names)		
}

ScoreSelect <- function(t.dat,snr=NULL,m.names,wr,levs=NULL,lmain=""){
    sf <- .8
    library(RColorBrewer)
    m.names <- intersect(m.names,names(t.dat))
    lwds <- 3
    if(length(m.names) == 0)
    {stop("no named traces exist in trace dataframe.")}
    
    xseq <- t.dat[,1]
    cols <-brewer.pal(8,"Dark2")
    cols <- rep(cols,ceiling(length(m.names)/length(cols)))
    cols <- cols[1:length(m.names)]
    dev.new(width=14,height=8)
    m.pca <- prcomp(t(t.dat[,m.names]),scale=F,center=T)
    m.names <- m.names[order(m.pca$x[,1],decreasing=sum(m.pca$rot[,1]) < 0)]
    hbc <- length(m.names)*sf+min(2,max(t.dat[,m.names]))
    hb <- ceiling(hbc)
    
    plot(xseq,t.dat[,m.names[1]],ylim=c(-sf,hbc),xlab="Time (min)",ylab="Ratio with shift",main=lmain,type="n", xaxt="n")
    axis(1, at=seq(0, length(t.dat[,1]), 5))

    if(length(wr) > 0)
    {
        if(is.null(levs)){levs <- setdiff(unique(wr),"")}
        x1s <- tapply(xseq,as.factor(wr),min)[levs]
        x2s <- tapply(xseq,as.factor(wr),max)[levs]
        y1s <- rep(-.3,length(x1s))
        y2s <- rep(hbc+.2,length(x1s))
        rect(x1s,y1s,x2s,y2s,col="lightgrey")
        text(xseq[match(levs,wr)],rep(-.1,length(levs)),levs,pos=4,offset=0,cex=1)
    }
    x.sel <- NULL
    xs <-rep(0,(length(m.names)+4))
    ys <- seq(1,length(m.names))*sf+t.dat[1,m.names]
    ys <- as.vector(c(ys,c(2*sf,sf,0,-sf)))
    #    xs[(length(xs)-2):length(xs)] <- c(0,5,10)
    p.names <- c(m.names,"ALL","NONE","FINISH","DROP")
    drop.i <- length(p.names)
    done.n <- drop.i-1
    none.i <- drop.i-2
    all.i <- drop.i-3

    p.cols <- c(cols,c("black","black","black","black"))
    for(i in 1:length(m.names))
    {
        lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lwds)
        if(!is.null(snr))
        {
        pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
        pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
        points(xseq[pp1],t.dat[pp1,m.names[i]]+i*sf,pch=1,col=cols[i])
        points(xseq[pp2],t.dat[pp2,m.names[i]]+i*sf,pch=0,col=cols[i])
        }
    }
    text(x=xs,y=ys,labels=p.names,pos=2,cex=.7,col=p.cols)
    points(x=xs,y=ys,pch=16,col=p.cols)
    click.i <- 1    
    while(click.i < done.n)
    {
        click.i <- identify(xs,ys,n=1,plot=F)
        if(click.i < (length(m.names)+1) & click.i > 0)
        {
            i <- click.i
            if(is.element(i,x.sel))
            {
                lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lwds)
                x.sel <- setdiff(x.sel,i)
            }
                else
                {
                lines(xseq,t.dat[,m.names[i]]+i*sf,col="black",lwd=lwds)
                #lines(xseq,t.dat[,m.names[i]]+i*sf,col="white",lwd=2,lty=2)
                x.sel <- union(x.sel,i)
            }
        }
        if(click.i == none.i)
        {
            x.sel <- NULL
            for(i in 1:length(m.names))
            {
                lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lwds)
            }
        }
        if(click.i == all.i)	
        {
            x.sel <- seq(1,length(m.names))
            for(i in 1:length(m.names))
            {
                lines(xseq,t.dat[,m.names[i]]+i*sf,col="black",lwd=lwds)
            }
            
        }
    }
    return(list(cells=m.names[x.sel],click = p.names[click.i]))
}

##review binary scoring file and toggle 1/0
##names of binary scoring bin must be in wr
##NO NAs
ScoreReview1 <- function(tdat,bin,wr,maxt=20){
    subD <- function(xdat)#trace dat with names NO TIME COL
    {
        s.x <- apply(xdat,2,sum)
        s.names <- names(xdat)[order(s.x)]
        sub.list <- list()
        sub.i <- seq(1,ncol(xdat),by=(maxt+1))
        if(length(sub.i) > 1)
        {
        for(i in 1:(length(sub.i)-1))
        {
            sub.list[[i]] <- s.names[sub.i[i]:(sub.i[i]+maxt)]
        }
        }
        i <- length(sub.i)
        sub.list[[i]] <- s.names[sub.i[i]:(ncol(xdat))]
        return(sub.list)
    }
    
    b.levs <- names(bin)[names(bin) != "drop"]
    drop <- rep(0,nrow(bin))
    if(is.element("drop",names(bin))){drop <- bin[,"drop"]}
    names(drop) <- row.names(bin)
    for(i in b.levs)
    {
        lmain <- paste("Scored as 1 for ",i,sep="")
        b.1 <- row.names(bin)[bin[,i]==1 & drop==0]
        if(length(b.1) > 0)
        {
        if(length(b.1) < maxt){sub1 <- list(b.1)}else{sub1 <- subD(tdat[wr==i,b.1])}
        for(x.names in sub1)
        {
            no.names <- NULL
            dropit <- TRUE
            while(dropit==TRUE & (length(x.names) > 0))
            {
                inp <- ScoreSelect(tdat,,x.names,wr,i,lmain)
                no.names <- inp[["cells"]]
                dropit <- (inp[["click"]]=="DROP")
                if(dropit){drop[no.names] <- 1;x.names <- setdiff(x.names,no.names)}
                dev.off()
            }
            if(length(no.names) > 0)
            {
                bin[no.names,i] <- 0
            }
        }
        }
    }
    bin["drop"] <- drop
    return(bin)
        
}

ScoreReview0 <- function(tdat,bin,wr,maxt=20){
    subD <- function(xdat)#trace dat with names NO TIME COL
    {
        s.x <- apply(xdat,2,sum)
        s.names <- names(xdat)[order(s.x)]
        sub.list <- list()
        sub.i <- seq(1,ncol(xdat),by=(maxt+1))
        if(length(sub.i) > 1)
        {
        for(i in 1:(length(sub.i)-1))
        {
            sub.list[[i]] <- s.names[sub.i[i]:(sub.i[i]+maxt)]
        }
        }
        i <- length(sub.i)
        sub.list[[i]] <- s.names[sub.i[i]:(ncol(xdat))]
        return(sub.list)
    }
    
    b.levs <- names(bin)[names(bin) != "drop"]
    drop <- rep(0,nrow(bin))
    if(is.element("drop",names(bin))){drop <- bin[,"drop"]}
    names(drop) <- row.names(bin)
    for(i in b.levs)
    {
        lmain <- paste("Scored as 0 for ",i,sep="")
        b.1 <- row.names(bin)[bin[,i]==0 & drop==0]
        if(length(b.1) > 0)
        {
        if(length(b.1) < maxt){sub1 <- list(b.1)}else{sub1 <- subD(tdat[wr==i,b.1])}			
        for(x.names in sub1)
        {
            no.names <- NULL
            dropit <- TRUE
            while(dropit==TRUE & (length(x.names)>0))
            {
                inp <- ScoreSelect(tdat,,x.names,wr,i,lmain)
                no.names <- inp[["cells"]]
                dropit <- (inp[["click"]]=="DROP")
                if(dropit){drop[no.names] <- 1;x.names <- setdiff(x.names,no.names)}
                dev.off()
            }
            if(length(no.names) > 0)
            {
                bin[no.names,i] <- 1
            }
        }
        }
    }
    bin["drop"] <- drop
    return(bin)
}

#Now lets create the function to put premade groups into a binary collumns
col_binner<-function(dat,cells){
    cell.names<-select.list(names(cells),multiple=T)
    cells<-cells[cell.names]
    for(i in 1:length(cell.names)){
        dat$bin[cell.names[i]]<-0
        dat$bin[ cells[[i]], cell.names[i] ]<-1
    }
    return(dat)
}

# Create Binary Classes of cells
#AdvANCED, ANY WORK DOWNE IS AUTOMATICALLY SAVED TO THE dAT
#<- IS NO LOINGER REUQIRED
#dAT IS THE rd. INPUT
#LEVS, IS AN OPTIONAL ARGUEMENT, IF LEFT BLACK THE FUNCTION WILL LOOK IN THE BIN DATA.FRAME COLLUMN ANMES
# TO ALLOW TO MANUALLY SECT THE COLLUMNS YOU WANT TO COMBINE
combinner<-function(dat, levs=NULL, bin_it=T){     
    tmp<-dat$bin
    if(is.null(levs)){
        levs<-select.list(names(dat$bin), multiple=T)
    }else{}
    
    newcolnames<-paste(levs,collapse="___")
    pf<-apply(tmp[,levs],1,paste, collapse="")
    pf.sum<-summary(as.factor(pf), maxsum=1500)
    pf.sum<-pf.sum[order(pf.sum, decreasing=T)]
    pf.ord<-pf.sum
    pf.ord[]<-seq(1,length(pf.sum))
    
    dat$scp[newcolnames]<-as.factor(pf)
    
    dat$bin<-tmp
    cat("We have added this barcode to the scp dataframes","\n")
    cat(newcolnames,"\n")
    cat(sort(summary(dat$scp[newcolnames], maxsum=1500), T), sep="\n" )

    if(bin_it){
        dat <- pf_summary(dat,,ncol(dat$scp))
    }

    return(dat)
}

# Create Binary Classes of cells
pf.function<-function(dat, levs){
    tmp<-dat
    pf<-apply(tmp[,levs],1,paste, collapse="")
    pf.sum<-summary(as.factor(pf), maxsum=1500)
    pf.sum<-pf.sum[order(pf.sum, decreasing=T)]
    pf.ord<-pf.sum
    pf.ord[]<-seq(1,length(pf.sum))
    tmp["pf"]<-as.factor(pf)
    tmp["pf.sum"]<-pf.sum[pf]
    tmp["pf.ord"]<-pf.ord[pf]
    return(tmp)
 }
 
bin_to_group<-function(dat){
    bin<-dat$bin
    cat("
    Select the collumns you would like to collect the rows that are scored as 1's.\n")

    cols_sel<-select.list(names(dat$bin), multiple=T)
    
    cell_group<-list()
    for(i in 1:length(cols_sel)){
        cell_group[[ cols_sel[i] ]]<-row.names(which(dat$bin[ cols_sel[i] ]==1,arr.ind=T))
    }
    
    return(cell_group)
}

#This takes a pf and allows you to create a binarry table based on the barcode
#Created in pf.function
pf_summary<-function(dat, response_classes = NULL, pf_col = NULL){
    if(is.null(pf_col)){
        pf_col <- menu( colnames(dat$scp) )
    }else{ pf_col <- pf_col }

    if(is.null(response_classes)){
        response_classes <- unique(dat$scp[,pf_col])
    }else{}
    
    for(i in 1:length(response_classes)){
        response.types<-row.names(
            which(
                dat$scp[pf_col] == as.character(response_classes[i])
            , arr.ind=T)
        )
        dat$bin[ as.character(response_classes[i]) ]<-0
        dat$bin[ response.types, as.character(response_classes[i]) ]<-1
    }
    cat("I Have added new rows to your bin dataframe based off of this \nresponse combination","\n\n")
    cat(colnames(dat$scp)[pf_col], sep="\n")
    cat(as.character(response_classes), sep='\n')
    return(dat)
}

census_viewer<-function(dat, census){
    cat("This is a function where you select your census of interest and then a cell type of interest, use the table /n to reference your choices here")
    sel.i<-1
    while(sel.i!=0){
        cells_to_view<-select.list(names(census))
        cell_type_to_reference<-select.list(names(dat$cell_types))
    
        cells_of_interest<-intersect(dat$cell_types[[cell_type_to_reference]],census[[cells_to_view]])
        
        if(length(cells_of_interest)>1){
            tcd(dat, cells_of_interest)
        }else{}
        
        cat("Would you like to look at t another cell in your census? enter 1, if not enter 0 \n")
        sel.i<-scan(n=1)
    }
}

##############################################################################################
##############################################################################################
#tmp is an RD object, x.names are the cell ids to investiage
#pad is the extra amount of image to select around the cell e.g. 1 = at cell bondaries 1.05 = 5% extra
#stain.name is the stain to display ("tritc","gfp","dapi") anything else defaults to yellow ROI boundaries
#title1 will be the title of the grid selection window.
SelectGrid <- function(tmp,x.names,pad=1.05,stain.name="area",title1="SelectRed",window.h=7,window.w=7,l.col="red",roi.img=NULL){
        imgs <- grep("img",names(tmp),value=T)	
        imgs.yes <- rep(F,length(imgs))
        for(i in 1:length(imgs)){imgs.yes[i] <- length(dim(tmp[[imgs[i]]])) == 3}
        imgs <- imgs[imgs.yes]
        if(length(imgs) < 1){stop("no image data")}	
        imgs.yes <- rep(F,length(imgs))
        for(i in 1:length(imgs)){imgs.yes[i] <- dim(tmp[[imgs[i]]])[3] == 3}
        imgs <- imgs[imgs.yes]
        
        if(length(imgs) < 1){stop("no image data")}
        img.rgb <- data.frame(name=imgs)
        img.rgb["r"] <- 0
        img.rgb["g"] <- 0
        img.rgb["b"] <- 0
        
        for(j in 1:nrow(img.rgb))
        {
            img.rgb[j,"r"] <- as.numeric(mean(tmp[[imgs[j]]][,,1]))
            img.rgb[j,"g"] <- as.numeric(mean(tmp[[imgs[j]]][,,2]))
            img.rgb[j,"b"] <- as.numeric(mean(tmp[[imgs[j]]][,,3]))
        }
        img.rgb["rgb"]<-rowSums(img.rgb[,2:4])
        #set the channel to use and subtract the others. red=1, green=2, blue=3
        #also select the best image.
        
        
        img.red <- imgs[which.max((img.rgb[,"r"]-img.rgb[,"g"]-img.rgb[,"b"])/img.rgb[,"rgb"])]
        
        #This is the old way for finding the green image.
        img.green <- imgs[which.max((img.rgb[,"g"]-(img.rgb[,"r"]-img.rgb[,"b"]))/img.rgb[,"rgb"])]
        #Find Green Image by finding the red neagtive images first
        #red.negative<-which(img.rgb['r']==min(img.rgb['r']))
        #then find the row in a red.neagitive matrix where the green is maximized
        #img.green<-img.rgb[which.max(img.rgb[red.negative,'g']), 1]
        
        img.blue <- imgs[which.max((img.rgb[,"b"]-(img.rgb[,"r"]-img.rgb[,"g"]))/img.rgb[,"rgb"])]
        #img.yellow <- imgs[which.max(img.rgb[,"r"]+img.rgb[,"g"]-img.rgb[,"b"])]
        if(is.null(roi.img)){img.yellow<-"img7"}
        else(img.yellow<-roi.img)	

        if(is.element(stain.name,c("tritc","gfp","dapi","mcherry","cy5","tritc.immuno")))
        {
            sn <- grep(stain.name,names(tmp$c.dat),ignore.case=T,value=T)[1]
            print(sn)
            if(is.null(sn)){stop("no stain value data")}
            x.names <- x.names[order(tmp$c.dat[x.names,sn])]
            if(stain.name=="tritc")
            {
                img.name <- imgs[which.max((img.rgb[,"r"]-img.rgb[,"g"]-img.rgb[,"b"])/img.rgb[,"rgb"])]
                chn <- 1
            }
            if(stain.name=="mcherry")
            {
                img.name <- imgs[which.max((img.rgb[,"r"]-img.rgb[,"g"]-img.rgb[,"b"])/img.rgb[,"rgb"])]
                chn <- 1
            }
            if(stain.name=="cy5")
            {
                img.name <- imgs[which.max((img.rgb[,"r"]-img.rgb[,"g"]-img.rgb[,"b"])/img.rgb[,"rgb"])]
                chn <- 1
            }
            if(stain.name=="gfp")
            {
                img.name <- imgs[which.max((img.rgb[,"g"]-(img.rgb[,"r"]+img.rgb[,"g"]))/img.rgb[,"rgb"])]
                chn <- 2		
            }
            if(stain.name=="dapi")
            {
                img.name <- imgs[which.max((img.rgb[,"b"]-img.rgb[,"r"]-img.rgb[,"g"])/img.rgb[,"rgb"])]
                chn <- 3
            }
            if(stain.name=="tritc.immuno")
            {
                img.name <- imgs[which.max((img.rgb[,"b"]-img.rgb[,"r"]-img.rgb[,"g"])/img.rgb[,"rgb"])]
                chn <- 3
            }
            
            img <- tmp[[img.name]]
            img.dat <- img[,,chn]
            for(i in setdiff(c(1,2,3),chn)){gt.mat <- img.dat < img[,,i];img.dat[gt.mat] <- 0} 
        }else{
            img.name <- img.yellow
            if(is.null(img.name)){img.name <- imgs[which.max(img.rgb[,"b"]+img.rgb[,"r"]-img.rgb[,"g"])]}
            
            sn <- intersect(c("area","circularity"),names(tmp$c.dat))[1]
            x.names <- x.names[order(tmp$c.dat[x.names,sn])]
            img <- tmp[[img.name]]
            img.dat <- (img[,,1]+img[,,2])/2
            med.r <- .99
            med.b <- .99
            if(sum(as.vector(img[,,1]) > med.r)==0){med.r <- quantile(as.vector(img[,,1]),probs=c(.95))[1]}
            if(sum(as.vector(img[,,2]) > med.b)==0){med.b <- quantile(as.vector(img[,,2]),probs=c(.95))[1]}
            img.dat[img[,,1] < med.r] <- 0
            img.dat[img[,,2] < med.b] <- 0		

            #single.img <- tmp$img4
        }
        
        #set up two devices
        graphics.off()
        dev.new(height=window.h,width=window.w,canvas="black",title="SingleCell")
        dev.single <- dev.cur()
        op <- par(mar=c(0,0,0,0))	
        plot(c(0,1),c(0,1),xaxt="n",yaxt="n",type="n",ylab="",xlab="")	
        
        dev.new(height=window.w,width=window.h,canvas="black",title=title1)
        dev.grid <- dev.cur()
        op <- par(mar=c(0,0,0,0))	
        plot(c(0,1),c(0,1),xaxt="n",yaxt="n",type="n",ylab="",xlab="")	
        xn <- length(x.names)
        num.grid <- xn+3
        nr <- floor(sqrt(num.grid))
        nc <- ceiling((num.grid)/nr)
        mtx <- max(nr,nc)
        dx <- seq(0,1,length.out=(mtx+1))[-1]
        sl <- (dx[2]-dx[1])/2
        dx <- dx-sl
        all.x <- as.vector(matrix(rep(dx,mtx),byrow=F,ncol=mtx))
        all.y <- as.vector(matrix(rep(dx,mtx),nrow=mtx,byrow=T))
        
        zf<-(sqrt(tmp$c.dat[x.names,"area"])/pi)*pad
        x <- tmp$c.dat[x.names,"center.x"]
        y <- tmp$c.dat[x.names,"center.y"]
        img.dimx<-dim(tmp$img1)[2]
        img.dimy<-dim(tmp$img1)[1]

        zf[zf > x] <- x[zf > x]
        zf[zf > y] <- y[zf > y]
        zf[x+zf > img.dimx] <- img.dimx-x[x+zf > img.dimx]
        zf[y+zf > img.dimy] <- img.dimy-y[y+zf > img.dimy]
        
        img.left<- x-zf
        img.left[img.left < 1] <- 1
        img.right<- x+zf
        img.right[img.right > img.dimx] <- img.dimx
        img.top<- y-zf
        img.top[img.top < 1] <- 1
        img.bottom<-y+zf
        img.bottom[img.bottom > img.dimy] <- img.dimy

        img.bottom[img.top>=img.bottom & img.top<img.dimy] <- img.top[img.top>=img.bottom] + 1
        img.right[img.left>=img.right & img.left<img.dimx] <- img.left[img.left>=img.right] + 1

        img.top[img.top == img.dimy] <- img.dimy-1
        img.left[img.left == img.dimx] <- img.dimx-1
            
        for(i in 1:xn)
        {
            xl <- all.x[i]-sl*.9
            xr <- all.x[i]+sl*.9
            xt <- all.y[i]-sl*.9
            xb <- all.y[i]+sl*.9
            #rasterImage(tmp$img1[img.bottom[i]:img.top[i],img.left[i]:img.right[i],],xl,xb,xr,xt)
            rasterImage(img.dat[img.bottom[i]:img.top[i],img.left[i]:img.right[i]],xl,xb,xr,xt)
        }
        fg <- rep("black",length(all.x))
        fg[1:xn] <- "grey"
        cexr <- sl/.04
        symbols(all.x,all.y,squares=rep(sl*1.9,length(all.x)),add=T,inches=F,fg=fg,lwd=cexr)
        text(all.x[xn+1],all.y[xn+1],"Done",col="white",cex= cexr)
        text(all.x[xn+2],all.y[xn+2],"All",col="white",cex= cexr)
        text(all.x[xn+3],all.y[xn+3],"None",col="white",cex= cexr)

        #first click defines the split
        all.sel <- rep(0,xn)
        names(all.sel) <- x.names	
        not.done=TRUE
        click1 <- locator(n=1)
        dist <- sqrt((click1$x[[1]]-all.x)^2 + (click1$y[[1]]-all.y)^2)
        sel.i <- which.min(dist)
        if(sel.i == xn+1){not.done=FALSE;return(all.sel)}
        if(sel.i == xn+2){all.sel[1:xn] <- 1;fg[1:xn] <- l.col}
        if(sel.i == xn+3){all.sel[1:xn] <- 0;fg[1:xn] <- "grey"}
        if(sel.i <= xn)
        {
        dev.set(which=dev.single)
        rasterImage(tmp[[img.red]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],0,0,.5,.5,interpolate=F)
        rasterImage(tmp[[img.green]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],0,.5,.5,1,interpolate=F)
        rasterImage(tmp[[img.blue]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],.5,0,1,.5,interpolate=F)
        rasterImage(tmp[[img.yellow]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],.5,.5,1,1,interpolate=F)
        abline(h=.5,col="grey")
        abline(v=.5,col="grey")
        
        dev.set(which=dev.grid)	
        neg.i <- 1:max((sel.i-1),1) 
        all.sel[neg.i] <- 0
        pos.i <- sel.i:xn	
        all.sel[pos.i] <- 1
        fg[neg.i] <- "grey"
        fg[pos.i] <- l.col
        }
        while(not.done)
        {
            symbols(all.x,all.y,squares=rep(sl*1.9,length(all.x)),add=T,inches=F,fg=fg,lwd=cexr)
            click1 <- locator(n=1)
            dist <- sqrt((click1$x[[1]]-all.x)^2 + (click1$y[[1]]-all.y)^2)
            sel.i <- which.min(dist)
            if(sel.i == xn+1){not.done=FALSE;return(all.sel)}
            if(sel.i == xn+2){all.sel[1:xn] <- 1;fg[1:xn] <- l.col}
            if(sel.i == xn+3){all.sel[1:xn] <- 0;fg[1:xn] <- "grey"}
            if(sel.i <= xn)
            {
            dev.set(which=dev.single)
        rasterImage(tmp[[img.red]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],0,0,.5,.5,interpolate=F)
        rasterImage(tmp[[img.green]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],0,.5,.5,1,interpolate=F)
        rasterImage(tmp[[img.blue]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],.5,0,1,.5,interpolate=F)
        rasterImage(tmp[[img.yellow]][img.bottom[sel.i]:img.top[sel.i],img.left[sel.i]:img.right[sel.i],],.5,.5,1,1,interpolate=F)
        abline(h=.5,col="grey")
        abline(v=.5,col="grey")

        dev.set(which=dev.grid)	
        if(all.sel[sel.i] ==0)
        {
            all.sel[sel.i] <- 1
            fg[sel.i] <- l.col
        }
        else
        {
            all.sel[sel.i] <- 0
            fg[sel.i] <- "grey"
        }
        }
    }		
    
}

#three tests Drop (confirm), Red (confirm) and Green (confirm)
#return and RD object with the changes made to c.dat and bin
#tmp is an RD object with images, "tritc.mean" and "gfp.mean" in c.dat
#x.names is a list of specific cells to review
#pad is the expansion factor about the center of the cell.
#subset.n is number of cells to review at once instead of all at once.
ROIreview <- function(tmp,x.names=NULL,pad=2,wh=7,hh=7,subset.n=500, roi.img=NULL){	
    print(names(tmp$c.dat)[1:20])
    choices<-select.list(
        title="Score what?", 
        choices=c("CGRP.GFP", "IB4.TRITC", "IB4.CY5", "NF200.TRITC", "MCHERRY", "Drops"), 
        multiple=T)
    print("how to display ROI")
    if(is.null(roi.img)){roi.img<-image.selector(tmp)}else{roi.img<-roi.img}

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

    if(is.null(x.names)){x.names <- row.names(tmp$c.dat)}
    x.names <- x.names[tmp$bin[x.names,"drop"]==0]
    if(is.na(subset.n) | subset.n > length(x.names)){subset.n=length(x.names)}
    subset.list <- dice(x.names,subset.n,subset.n/4)
    for(x.names in subset.list)
    {

        #drop cells
        if(length(grep("TRUE",choices=="Drops"))>0){
            d.names <- SelectGrid(tmp,x.names,pad,"area","SelectDrops",window.h=hh,window.w=wh,roi.img=roi.img)
            d1.names <- names(d.names[d.names==1])
            if(length(d1.names) > 5)
            {
                d1.names <- SelectGrid(tmp,d1.names,pad,"area","ConfirmDrops",window.h=hh,window.w=wh,roi.img=roi.img) 
                d1.names <- names(d1.names)[d1.names==1]
                if(length(d1.names) > 0){tmp$bin[d1.names,"drop"] <- 1;x.names <- setdiff(x.names,d1.names)}
            }
        }else{}

        #Red Cells
        if(length(grep("TRUE",choices=="IB4.TRITC"))>0){
            r.names <- SelectGrid(tmp,x.names,pad,"tritc","SelectRed",window.h=hh,window.w=wh,roi.img=roi.img)
            r1.names <- names(r.names[r.names==1])
            q1 <- 1:floor(length(r1.names)*.25)
            r2.names <- r1.names[q1]
            if(length(r2.names) > 5)
            {
                r2.names <- SelectGrid(tmp,r2.names,pad*2,"tritc","ConfirmRed",window.h=hh,window.w=wh,roi.img=roi.img)
                r.names[names(r2.names)] <- r2.names
            }
            tmp$bin[names(r.names),"tritc.bin"] <- r.names
        }else{}
        
        #Red Cells
        if(length(grep("TRUE",choices=="IB4.CY5"))>0){
            r.names <- SelectGrid(tmp,x.names,pad,"cy5","SelectRed",window.h=hh,window.w=wh,roi.img=roi.img)
            r1.names <- names(r.names[r.names==1])
            q1 <- 1:floor(length(r1.names)*.25)
            r2.names <- r1.names[q1]
            if(length(r2.names) > 5)
            {
                r2.names <- SelectGrid(tmp,r2.names,pad*2,"cy5","ConfirmRed",window.h=hh,window.w=wh,roi.img=roi.img)
                r.names[names(r2.names)] <- r2.names
            }
            tmp$bin[names(r.names),"cy5.bin"] <- r.names
        }else{}

        #Green Cells
        if(length(grep("TRUE",choices=="CGRP.GFP"))>0){
            r.names <- SelectGrid(tmp,x.names,pad,"gfp","SelectGreen",window.h=hh,window.w=wh,l.col="green",roi.img=roi.img)
            r1.names <- names(r.names[r.names==1])
            q1 <- 1:floor(length(r1.names)*.25)
            r2.names <- r1.names[q1]
            if(length(r2.names) > 5)
            {
                r2.names <- SelectGrid(tmp,r2.names,pad*2,"gfp","ConfirmGreen",window.h=hh,window.w=wh,l.col="green",roi.img=roi.img)
                r.names[names(r2.names)] <- r2.names
            }
            tmp$bin[names(r.names),"gfp.bin"] <- r.names
        }else{}
        
        #NF200
        if(length(grep("TRUE",choices=="NF200.TRITC"))>0){
            r.names <- SelectGrid(tmp,x.names,pad,"tritc.immuno","SelectBlue",window.h=hh,window.w=wh,l.col="blue",roi.img=roi.img)
            r1.names <- names(r.names[r.names==1])
            q1 <- 1:floor(length(r1.names)*.25)
            r2.names <- r1.names[q1]
            if(length(r2.names) > 5)
            {
                r2.names <- SelectGrid(tmp,r2.names,pad*2,"tritc.immuno","ConfirmBlue",window.h=hh,window.w=wh,l.col="blue",roi.img=roi.img)
                r.names[names(r2.names)] <- r2.names
            }
            tmp$bin[names(r.names),"tritc.bin"] <- r.names
        }else{}

        #MCHERRY
        if(length(grep("TRUE",choices=="MCHERRY"))>0){
            r.names <- SelectGrid(tmp,x.names,pad,"mcherry","SelectRed",window.h=hh,window.w=wh,l.col="red",roi.img=roi.img)
            r1.names <- names(r.names[r.names==1])
            q1 <- 1:floor(length(r1.names)*.25)
            r2.names <- r1.names[q1]
            if(length(r2.names) > 5)
            {
                r2.names <- SelectGrid(tmp,r2.names,pad*2,"mcherry","ConfirmRed",window.h=hh,window.w=wh,l.col="red",roi.img=roi.img)
                r.names[names(r2.names)] <- r2.names
            }
            tmp$bin[names(r.names),"mcherry.bin"] <- r.names
        }else{}

        
    }

    graphics.off()
    return(tmp)			
}

#three tests Drop (confirm), Red (confirm) and Green (confirm)
#return and RD object with the changes made to c.dat and bin
#tmp is an RD object with images, "tritc.mean" and "gfp.mean" in c.dat
#x.names is a list of specific cells to review
#pad is the expansion factor about the center of the cell.
#subset.n is number of cells to review at once instead of all at once.
ROIreview2 <- function(tmp,x.names=NULL,pad=2,wh=7,hh=7,subset.n=NA){
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

    if(is.null(x.names)){x.names <- row.names(tmp$c.dat)}
    x.names <- x.names[tmp$bin[x.names,"drop"]==0]
    if(is.na(subset.n) | subset.n > length(x.names)){subset.n=length(x.names)}
    subset.list <- dice(x.names,subset.n,subset.n/4)
    for(x.names in subset.list)
    {
        #drop cells
        d.names <- SelectGrid(tmp,x.names,pad,"area","SelectDrops",window.h=hh,window.w=wh)
        d1.names <- names(d.names[d.names==1])
        if(length(d1.names) > 5)
        {
            d1.names <- SelectGrid(tmp,d1.names,pad,"area","ConfirmDrops",window.h=hh,window.w=wh) 
            d1.names <- names(d1.names)[d1.names==1]
            if(length(d1.names) > 0){tmp$bin[d1.names,"drop"] <- 1;x.names <- setdiff(x.names,d1.names)}
        }
        r.names <- SelectGrid(tmp,x.names,pad,"tritc","SelectRed",window.h=hh,window.w=wh)
        r1.names <- names(r.names[r.names==1])
        q1 <- 1:floor(length(r1.names)*.25)
        r2.names <- r1.names[q1]
        if(length(r2.names) > 5)
        {
            r2.names <- SelectGrid(tmp,r2.names,pad*2,"tritc","ConfirmRed",window.h=hh,window.w=wh)
            r.names[names(r2.names)] <- r2.names
        }
        tmp$bin[names(r.names),"tritc.bin"] <- r.names
    
        #r.names <- SelectGrid(tmp,x.names,pad,"gfp","SelectGreen",window.h=hh,window.w=wh,l.col="green")
        #r1.names <- names(r.names[r.names==1])
        #q1 <- 1:floor(length(r1.names)*.25)
        #r2.names <- r1.names[q1]
        #if(length(r2.names) > 5)
        #{
        #	r2.names <- SelectGrid(tmp,r2.names,pad*2,"gfp","ConfirmGreen",window.h=hh,window.w=wh,l.col="green")
        #	r.names[names(r2.names)] <- r2.names
        #}
        #tmp$bin[names(r.names),"gfp.bin"] <- r.names
        }
    return(tmp)			
}

##############################################################################################
# Drop Scoring
##############################################################################################
# Functions to allow for dropping of cells.  Main function is DropTestMulti
# Drops based on spikey traces, out of window peaks, and baselineshifts

SpikyNorm <- function(xdat){
        shapfunc <- function(x){shapiro.test(x)$p.value}
        i1 <- seq(1,nrow(xdat))
        s1 <- xdat[c(1,i1[-length(i1)]),] #shift 1 time interval forward
        s2 <- xdat[c(i1[-1],i1[length(i1)]),] #shift 1 time interval back
        s3 <- xdat-((s1+s2)/2)
        s.x <- apply(abs(s3),2,shapfunc)
    return(s.x)	
}

DropPick <- function(tdat,bin,wr,maxt=10,s.x=NULL,lmain="Select Cells to drop"){
    #order traces by spikey trait.
    #allow drop selection until 0 selected.
    #spikes are defined as single point deviations from previous and next.
    subD <- function(s.x)#trace dat with names NO TIME COL
    {
        s.names <- names(s.x)[order(s.x)]
        sub.list <- list()
        sub.i <- seq(1,length(s.x),by=(maxt+1))
        if(length(sub.i) > 1)
        {
        for(i in 1:(length(sub.i)-1))
        {
            sub.list[[i]] <- s.names[sub.i[i]:(sub.i[i]+maxt)]
        }
        }
        i <- length(sub.i)
        sub.list[[i]] <- s.names[sub.i[i]:(length(s.x))]
        return(sub.list)
    }
    
    b.levs <- c("drop") #names(bin)[names(bin) != "drop"]
    drop <- rep(0,nrow(bin))
    if(is.element("drop",names(bin))){drop <- bin[,"drop"]}
    names(drop) <- row.names(bin)
    for(i in b.levs)
    {

        b.1 <- row.names(bin)[bin[,i]==0 & drop==0]
        if(is.null(s.x)){s.x <- SpikyNorm(tdat[,-1])}

        if(length(b.1) > 0)
        {
            s.x <-s.x[b.1]
        if(length(b.1) < maxt){sub1 <- list(b.1)}else{sub1 <- subD(s.x)}

        for(x.names in sub1)
        {
            no.names <- NULL
            dropit <- TRUE
            nd <- 0			
            while(dropit==TRUE & (length(x.names)>0))
            {

                inp <- ScoreSelect(tdat,,x.names,wr,,lmain)
                no.names <- inp[["cells"]]
                dropit <- (inp[["click"]]=="DROP")
                if(dropit){drop[no.names] <- 1;x.names <- setdiff(x.names,no.names);nd=1}
                dev.off()
            }
            if(length(no.names) > 0)
            {
                drop[no.names] <- 1
            }
            if(length(no.names)==0 & nd==0)
            {break}
        }
        }
    }
    return(drop)
}

DropTestList <- function(tmp){
        #tmp <- get(rd.name)
        x1 <- DropPick(tmp$t.dat,tmp$bin,tmp$w.dat[,"wr1"],lmain="Select spikey traces to Drop") #defaults to spiky test
        tmp$bin[,"drop"] <- x1
        x1 <- DropPick(tmp$t.dat,tmp$bin,tmp$w.dat[,"wr1"],s.x= -apply(tmp$scp[,"snr.owc",drop=F],1,mean),lmain="Select out of window peaks to Drop")
        tmp$bin[,"drop"] <- x1		
        x1 <- DropPick(tmp$t.dat,tmp$bin,tmp$w.dat[,"wr1"],s.x= -apply(tmp$scp[,"bl.diff",drop=F],1,mean),lmain="Select Baseline Drops")
        tmp$bin[,"drop"] <- x1		
        if(sum(x1 > 0)) #check highest correlations with dropped cells.
        {
            d.names <- names(x1[x1==1])
            ct <- cor(tmp$t.dat[,-1])
            mn <- -apply(ct[,d.names],1,max)
            x1 <- DropPick(tmp$t.dat,tmp$bin,tmp$w.dat[,"wr1"],s.x= mn,lmain="Correlated with other drops")
            tmp$bin[,"drop"] <- x1		
        }
        return(tmp)
}

DropTestMulti <- function(dir.name=NULL,snr.lim=4,hab.lim=.05,sm=3,ws=30,review=F){
    if(is.null(dir.name)){dir.name <- getwd()}
    setwd(dir.name)
    f.names <- list.files(pattern="RD.*\\.Rdata$")
    if(length(f.names) == 0){stop("no RD...Rdata files in given directory")}
    rd.list <- sub("\\.Rdata*","",f.names)
    RD.names <- rd.list #paste(rd.list,".b",sep="")
    RD.f.names <- paste(RD.names,".Rdata",sep="")
    sel.i <- menu(rd.list,title="Select Data to review")
    while(sel.i != 0)
    {
        j <- sel.i
        load(f.names[j])
        i <- rd.list[j]
        tmp <- get(i)
        tlevs <- c(as.character(unique(tmp$w.dat[,"wr1"])[-1]),"drop")
        if(is.null(tmp$bin))
        {
        tmp.pcp <- ProcConstPharm(tmp,sm,ws,"TopHat")
        tmp.scp <- ScoreConstPharm(tmp$t.dat,tmp.pcp$blc,tmp.pcp$snr,snr.lim,hab.lim,tmp$w.dat[,"wr1"],sm)
        tmp.bin <- bScore(tmp.pcp$blc,tmp.pcp$snr,snr.lim,hab.lim,tlevs,tmp$w.dat[,"wr1"])
        tmp.bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
        }
        else
        {
            tmp.pcp <- ProcConstPharm(tmp,sm,ws,"TopHat")
            tmp.scp <- ScoreConstPharm(tmp$t.dat,tmp.pcp$blc,tmp.pcp$snr,snr.lim,hab.lim,tmp$w.dat[,"wr1"],sm)
            tmp.bin <- tmp$bin
            tmp.scp <- tmp$scp
            #tmp.blc <- tmp$blc
        }

        tmp$bin <- tmp.bin[,tlevs]
        tmp$scp <- tmp.scp
        #tmp$blc <- tmp.blc
        
        tmp <- DropTestList(tmp)
        if(review)
        {
              tmp.bin <- ScoreReview1(tmp$t.dat,tmp.bin[,tlevs],tmp$w.dat[,"wr1"])
              tmp.bin <- ScoreReview0(tmp$t.dat,tmp.bin[,tlevs],tmp$w.dat[,"wr1"])
            tmp$bin <- tmp.bin[,tlevs]
        }
        pf<-apply(tmp$bin[,tlevs],1,paste,collapse="")
        pf.sum<-summary(as.factor(pf),maxsum=500)
        pf.sum<-pf.sum[order(pf.sum,decreasing=T)]
        pf.ord<-pf.sum
        pf.ord[]<-seq(1,length(pf.sum))
        tmp$c.dat["pf"]<-as.factor(pf)
        tmp$c.dat["pf.sum"]<-pf.sum[pf]
        tmp$c.dat["pf.ord"]<-pf.ord[pf]
        
        
        tmp$scp <- tmp.scp
        tmp$snr<-tmp.pcp$snr
        tmp$blc <- tmp.pcp$blc
    
        assign(RD.names[j],tmp)		
        save(list=RD.names[j],file=RD.f.names[j])
        print(paste("DONE REVIEWING ",RD.names[j]," CHANGES SAVED TO FILE.",sep=""))
        print(paste("Dropped Cells:", table(tmp$bin[,"drop"])[2]))
        sel.i <- menu(rd.list,title="Select Data to review")			
    }
    return(RD.f.names)		
}

##############################################################################################
##############################################################################################

##############################################################################################
# No Scoring, only processing
##############################################################################################

Trace.prep<-function(dir.name=NULL,snr.lim=4,hab.lim=.05,sm=3,ws=30,blc="SNIP"){
    if(is.null(dir.name)){dir.name <- getwd()}
    setwd(dir.name)
    f.names <- list.files(pattern="RD.*\\.Rdata$")
    if(length(f.names) == 0){stop("no RD...Rdata files in given directory")}
    rd.list <- sub("\\.Rdata*","",f.names)
    RD.names <- rd.list #paste(rd.list,".b",sep="")
    RD.f.names <- paste(RD.names,".Rdata",sep="")
    sel.i <- menu(rd.list,title="Select Data to review")
    while(sel.i != 0)
    {
        j <- sel.i
        load(f.names[j])
        i <- rd.list[j]
        tmp <- get(i)
        tlevs<-c(setdiff(unique(as.character(tmp$w.dat[,2])),""),"drop")
        
        
        tmp.pcp <- ProcConstPharm(tmp,sm,ws,blc)
        tmp.scp <- ScoreConstPharm(tmp,tmp.pcp$blc,tmp.pcp$snr, tmp.pcp$der,snr.lim,hab.lim,sm)
        tmp.bin <- bScore(tmp.pcp$blc,tmp.pcp$snr,snr.lim,hab.lim,tlevs,tmp$w.dat[,"wr1"])
        tmp.bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
        
        pf<-apply(tmp.bin[,tlevs],1,paste,collapse="")
        pf.sum<-summary(as.factor(pf),maxsum=500)
        pf.sum<-pf.sum[order(pf.sum,decreasing=T)]
        pf.ord<-pf.sum
        pf.ord[]<-seq(1,length(pf.sum))
        tmp$c.dat["pf"]<-as.factor(pf)
        tmp$c.dat["pf.sum"]<-pf.sum[pf]
        tmp$c.dat["pf.ord"]<-pf.ord[pf]
        
        tmp$bin<-tmp.bin
        tmp$scp <- tmp.scp
        tmp$snr<-tmp.pcp$snr
        tmp$blc <- tmp.pcp$blc
        tmp$der<-tmp.pcp$der
        assign(RD.names[j],tmp)		
        save(list=RD.names[j],file=RD.f.names[j])
        print(paste("DONE REVIEWING ",RD.names[j]," CHANGES SAVED TO FILE.",sep=""))
        sel.i <- menu(rd.list,title="Select Data to review")			
    }
    return(RD.f.names)	
}

#this is not complete
#condi is the indicator for the conditional frequency table
#this is bad
#####add selection section of selection of experiments to include/exclude
#####conditional expresion tables.
SummarizeMulti <- function(dir.name=NULL,condi=1,recur=F){
    if(is.null(dir.name)){stop("not a directory")}
    setwd(dir.name)
    f.names <- list.files(pattern=".*RD.*\\.Rdata$",recursive=recur,full.names=T)
    f.names <- select.list(f.names,multiple=T,title="Select Experiments For Analysis")
    if(length(f.names) == 0){stop("no RD...Rdata files in given directory")}
    for(i in f.names){load(i)}
    rd.list <- sub("\\.Rdata*","",basename(f.names))
    RD.names <- ls(pat="^RD")
    RD.names <- intersect(rd.list,RD.names)
    if(!setequal(RD.names,rd.list)){stop("dataframes loaded do not match files listed in directory")}
    RD.f.names <- paste(RD.names,".Rdata",sep="")
    
    i <- rd.list[1]
    tmp <- get(i)
    if(sum(is.element(c("bin","scp"),names(tmp))) < 2){stop("Data frame has not been scored")}

    if(names(tmp$bin)[c(1,2)]==c("tot","sd"))
    {tmp$bin <- tmp$bin[,-c(1,2)]}
    freq.tab <- data.frame(mean=apply(tmp$bin[tmp$bin[,"drop"]==0,],2,mean))
    kfreq.tab <- data.frame(mean=apply(tmp$bin[tmp$bin[,"drop"]==0 & tmp$bin[,condi]==1,],2,mean))

    b.names <- row.names(freq.tab)[row.names(freq.tab) != "drop"]
    q.names <- paste(b.names,".max",sep="")
    resp.tab <- data.frame(mean=apply(tmp$scp[tmp$bin[,"drop"]==0,q.names],2,mean))
    for(rn in row.names(resp.tab)){resp.tab[rn,"mean"] <- mean(tmp$scp[tmp$bin[,"drop"]==0 & tmp$bin[,sub("\\.max$","",rn)]==1,rn],na.rm=T)}
    pf.tot <- data.frame(str = apply(tmp$bin[tmp$bin[,"drop"]==0,names(tmp$bin)!="drop"],1,paste,collapse=""))
    pf.tot["exp"] <- i
    for(j in 2:length(RD.names))
    {
        i <- rd.list[j]
        tmp <- get(i)
        if(names(tmp$bin)[c(1,2)]==c("tot","sd"))
        {tmp$bin <- tmp$bin[,-c(1,2)]}
        
        m1 <- apply(tmp$bin[tmp$bin[,"drop"]==0,],2,mean)
        freq.tab[i] <- m1[row.names(freq.tab)]
        m2 <- apply(tmp$bin[tmp$bin[,"drop"]==0 & tmp$bin[,condi]==1,],2,mean)
        kfreq.tab[i] <- m2[row.names(kfreq.tab)]
        resp.tab[i] <- NA
        for(rn in intersect(row.names(resp.tab),names(tmp$scp))){resp.tab[rn,i] <- mean(tmp$scp[tmp$bin[,"drop"]==0 & tmp$bin[,sub("\\.max$","",rn)]==1,rn],na.rm=T)}
        pf.tmp <- data.frame(str = apply(tmp$bin[tmp$bin[,"drop"]==0,names(tmp$bin)!="drop"],1,paste,collapse=""))		
        pf.tmp["exp"] <- i
        pf.tot <- rbind(pf.tot,pf.tmp)
    }
    names(freq.tab)[1] <- rd.list[1]
    names(kfreq.tab)[1] <- rd.list[1]
    names(resp.tab)[1] <- rd.list[1]
    pf.tab <- table(pf.tot[,1],pf.tot[,2])
    return(list(freq.tab=freq.tab,kfreq.tab=kfreq.tab,resp.tab=resp.tab,pf.tab=pf.tab))
}

##############################################################################################
# Stacked traces Plotting
##############################################################################################
LinesSome <- function(t.dat,snr=NULL,m.names,wr=NULL,levs=NULL,lmain="",pdf.name=NULL,morder=NULL,subset.n=5,sf=.25,lw=2,bcex=.6){
    library(cluster)
    if(length(m.names) < subset.n)
    {stop("group size lower than subset size")}
    pam5 <- pam(t(t.dat[,m.names]),k=subset.n)
    s.names <- row.names(pam5$medoids)
    if(!is.null(morder))
    {
        names(morder) <- m.names
        morder <- morder[s.names]
        }
    pam5.tab <- table(pam5$clustering)
    tags <- paste(paste("#",names(pam5.tab),sep=""),as.vector(pam5.tab),sep=":")
    LinesEvery(t.dat,snr,s.names,wr,levs,lmain,pdf.name,morder,rtag=tags,sf,lw,bcex)
    return(pam5$clustering)
}

LinesEvery <- function(t.dat,snr=NULL,m.names,wr,levs=NULL,lmain="",pdf.name=NULL,morder=NULL,rtag=NULL,sf=.7,lw=3,bcex=1,p.ht=7,p.wd=10){
    m.names <- intersect(m.names,names(t.dat))
    xseq <- t.dat[,1]
    library(RColorBrewer)

    if(length(m.names) > 0)
    {
        if(is.null(pdf.name))
        {dev.new(width=14,height=8)}
        else
        {if(length(grep("\\.pdf",pdf.name))>0){pdf(pdf.name,width=p.wd,height=p.ht)}else{png(pdf.name,width=1200,height=600)}}#pdf(pdf.name,width=28,height=16)}
        if(is.null(morder))
        {
            m.pca <- prcomp(t(t.dat[,m.names]),scale=F,center=T)
            morder <- m.pca$x[,1] * c(1,-1)[(sum(m.pca$rot[,1]) < 0)+1]
            #m.names <- m.names[order(m.pca$x[,1],decreasing=sum(m.pca$rot[,1]) < 0)]
        }
        m.names <- m.names[order(morder)]
        
        hbc <- length(m.names)*sf+max(t.dat[,m.names])
        hb <- ceiling(hbc)
        #cols <- rainbow(length(m.names),start=.55)
        cols <-brewer.pal(8,"Dark2")
        cols <- rep(cols,ceiling(length(m.names)/length(cols)))
        cols <- cols[1:length(m.names)]
        par(mar=c(4,1,4,1))
        plot(xseq,t.dat[,m.names[1]],ylim=c(0,hbc),xlab="Time (min)",main=lmain,type="n", xaxt="n",yaxt="n",xlim=c(min(xseq)-1.5,max(xseq)+1.5))#-sf
        axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
        if(!is.null(wr))
        {
            if(!is.null(levs))
            {
            #levs <- setdiff(unique(wr),"")
            x1s <- tapply(xseq,as.factor(wr),min)[levs]
            x2s <- tapply(xseq,as.factor(wr),max)[levs]
            y1s <- rep(-.3,length(x1s))
            y2s <- rep(hbc+.2,length(x1s))
            rect(x1s,y1s,x2s,y2s,col=NA,border="darkgrey")
            cpx <- xseq[match(levs,wr)+round(table(wr)[levs]/2,0)]
            offs <- nchar(levs)*.5
            text(cpx,rep(c(sf/2,sf),length=length(levs)),levs,pos=1,cex=bcex)#,offset=-offs
            }
        }
        for(i in 1:length(m.names))
        {
            lines(xseq,t.dat[,m.names[i]]+i*sf, cex=.5,col=cols[i],lty=1, lwd=lw)
            points(xseq,t.dat[,m.names[i]]+i*sf,pch=15, cex=.5,col=cols[i])
            if(!is.null(snr))
            {
            pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
            pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
                                        #                pp3 <- dat$crr[,m.names[i]] > 0
            points(xseq[pp1],t.dat[pp1,m.names[i]]+i/10,pch=1,col=cols[i])
            points(xseq[pp2],t.dat[pp2,m.names[i]]+i/10,pch=0,col=cols[i])
                                        #                points(xseq[pp3],t.dat[pp3,m.names[i]]+i/10,pch=2,col=cols[i],cex=.5)
                                        }    
        }
        text(rep(0,length(m.names)),seq(1,length(m.names))*sf+t.dat[1,m.names],m.names,cex=.8*bcex,col=cols,pos=2)
        if(!is.null(rtag))
        {
            rtag <- rtag[order(morder)]
            text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag,cex=.8*bcex,col=cols,pos=4)
        }
      if(!is.null(pdf.name))
        {dev.off()}
    }
    
}

#Simplified LinesEvery which only needs 2 entries; RD and m.names.
LinesEvery.2 <- function(dat,m.names, blc=FALSE, snr=NULL,lmain="",cols=NULL, levs=NULL,m.order=NULL,rtag=NULL,rtag2=NULL,rtag3=NULL, plot.new=TRUE,sf=.7,lw=.9,bcex=.8,p.ht=7,p.wd=10){
    if(blc){t.dat<-dat$blc}
    else{t.dat<-dat$t.dat}
    wr<-dat$w.dat[,2]
    if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")}
    else{levs<-levs}
    m.names <- intersect(m.names,names(t.dat))
    xseq <- t.dat[,1]
    if(plot.new){dev.new(width=10,height=6)}
    library(RColorBrewer)
    
    ## Tool for Sorting cells based on c.dat collumn name
    if(length(m.names) > 0)
    {        
        if(!is.null(m.order)){	
            tmp<-dat$c.dat[m.names,]
            n.order<-tmp[order(tmp[,m.order]),]
            m.names <- row.names(n.order)
        }
        else{
            m.pca <- prcomp(t(t.dat[,m.names]),scale=F,center=T)
            morder <- m.pca$x[,1] * c(1,-1)[(sum(m.pca$rot[,1]) < 0)+1]
            m.names <- m.names[order(m.pca$x[,1],decreasing=sum(m.pca$rot[,1]) < 0)]
            m.names <- m.names[order(morder)]
        }
        
    ## Tool for color labeleing
        if(is.null(cols)){
            #cols <- rainbow(length(m.names),start=.55)
            cols <-brewer.pal(8,"Dark2")
            cols <- rep(cols,ceiling(length(m.names)/length(cols)))
            cols <- cols[1:length(m.names)]
        } 
    ## Tool for single color labeling
        else {cols<-cols
            cols <- rep(cols,ceiling(length(m.names)/length(cols)))
            cols <- cols[1:length(m.names)]
        }
        
        hbc <- length(m.names)*sf+max(t.dat[,m.names])
        hb <- ceiling(hbc)
        #par(xpd=TRUE)
        par(mar=c(4,2,4,3))
        plot(xseq,t.dat[,m.names[1]],ylim=c(0,hbc),xlab="Time (min)",main=lmain,type="n", xaxt="n",yaxt="n",xlim=c(min(xseq)-1.5,max(xseq)))#-sf
        axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
        axis(2, 1.4, )
        text(rep(0,length(m.names)),seq(1,length(m.names))*sf+t.dat[1,m.names],m.names,cex=.8*bcex,col=cols,pos=2)

    ## Tool for adding window region labeling
        if(length(wr) > 0){
            #levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
            x1s <- tapply(xseq,as.factor(wr),min)[levs]
            x2s <- tapply(xseq,as.factor(wr),max)[levs]
            y1s <- rep(-.3,length(x1s))
            y2s <- rep(hbc+.2,length(x1s))
            rect(x1s,y1s,x2s,y2s,col="grey90",border="black")
            cpx <- xseq[match(levs,wr)+round(table(wr)[levs]/2,0)]
            offs <- nchar(levs)*.5
            text(dat$t.dat[match(levs,wr),"Time"],rep(c(sf/2,sf),length=length(levs)),levs,pos=4,offset=0,cex=bcex)#,offset=-offs}
    ## Tool for adding line and point plot for graph
            for(i in 1:length(m.names)){
                lines(xseq,t.dat[,m.names[i]]+i*sf, lty=1,col=cols[i],lwd=lw)
                points(xseq,t.dat[,m.names[i]]+i*sf,pch=16,col=cols[i],cex=.3)

                if(!is.null(snr)){
                    pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
                    pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
                    points(xseq[pp1],t.dat[pp1,m.names[i]]+i/10,pch=1,col=cols[i])
                    points(xseq[pp2],t.dat[pp2,m.names[i]]+i/10,pch=0,col=cols[i])
                }    
            }
        }
    ## Tool for adding cell data labeling to end of graph
            if(!is.null(dat$c.dat[m.names, "area"])){rtag<-"area";rtag <- round(dat$c.dat[m.names,rtag], digits=0)}
            else{rtag<-NULL}
            if(!is.null(dat$c.dat[m.names, "CGRP"])){rtag2<-"CGRP";rtag2 <- round(dat$c.dat[m.names,rtag2], digits=0)}
            else{rtag2<-NULL}
            #if(!is.null(dat$c.dat[m.names, "mean.gfp"])){rtag2<-"mean.gfp";rtag2 <- round(dat$c.dat[m.names,rtag2], digits=0)}
            #else{rtag2<-NULL}
            if(!is.null(dat$c.dat[m.names, "mean.gfp"])){rtag2<-"mean.gfp";rtag2 <- round(dat$c.dat[m.names,rtag2], digits=0)}
            else{rtag2<-NULL}
            if(!is.null(dat$c.dat[m.names, "IB4"])){rtag3<-"IB4";rtag3 <- round(dat$c.dat[m.names,rtag3], digits=0)}
            else{rtag3<-NULL}
            if(!is.null(dat$c.dat[m.names, "mean.tritc"])){rtag3<-"mean.tritc";rtag3 <- round(dat$c.dat[m.names,rtag3], digits=0)}
            else{rtag3<-NULL}
            if(!is.null(dat$c.dat[m.names, "mean.gfp.2"])){rtag4<-"mean.gfp.2";rtag4 <- round(dat$c.dat[m.names,rtag4], digits=0)}
            else{rtag4<-NULL}
            text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag,cex=.9*bcex,col=cols,pos=4)
            text(rep(max(xseq)*1.04,length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag2,cex=.9*bcex,col="darkgreen",pos=4)
            text(rep(max(xseq)*1.08,length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag3,cex=.9*bcex,col="red",pos=4)


         
    }
}

# pic.plot=T plots images next to trace, unles more than 10 traces
# XY.plot, shows cells in image
LinesEvery.3 <- function(dat,m.names, img=NULL,pic.plot=TRUE, XY.plot=TRUE, blc=T, snr=NULL,lmain="",cols=NULL, levs=NULL, levs.cols="grey90",m.order=NULL,rtag=NULL,rtag2=NULL,rtag3=NULL, plot.new=TRUE,sf=.7,lw=.9,bcex=.6,p.ht=7,p.wd=10){
    if(blc){t.dat<-dat$blc}
    else{t.dat<-dat$t.dat}
    wr<-dat$w.dat[,2]
    if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")}
    else{levs<-levs}
    m.names <- intersect(m.names,names(t.dat))
    xseq <- t.dat[,1]
    hbc <- length(m.names)*sf+max(t.dat[,m.names])
    hb <- ceiling(hbc)
    library(RColorBrewer)
    
    ## Tool for Sorting cells based on c.dat collumn name
        if(length(m.names) > 0)
        {

            if(!is.null(m.order)){	
                tmp<-dat$c.dat[m.names,]
                n.order<-tmp[order(tmp[,m.order]),]
                m.names <- row.names(n.order)
            }
            else{
                m.pca <- prcomp(t(t.dat[,m.names]),scale=F,center=T)
                morder <- m.pca$x[,1] * c(1,-1)[(sum(m.pca$rot[,1]) < 0)+1]
                m.names <- m.names[order(m.pca$x[,1],decreasing=sum(m.pca$rot[,1]) < 0)]
                m.names <- m.names[order(morder)]
            }
    ### Picture Plotting!		
        if(XY.plot==T){cell.zoom.2048(dat, cell=m.names,img=img, cols="white",zoom=F, plot.new=T)}
        ## Tool for color labeleing
            if(is.null(cols)){
                #cols <- rainbow(length(m.names),start=.55)
                cols <-brewer.pal(8,"Dark2")
                cols <- rep(cols,ceiling(length(m.names)/length(cols)))
                cols <- cols[1:length(m.names)]
            } 
        ## Tool for single color labeling
            else {cols<-cols
                cols <- rep(cols,ceiling(length(m.names)/length(cols)))
                cols <- cols[1:length(m.names)]
            }
            
            if(plot.new){dev.new(width=10,height=6)}
            par(xpd=FALSE)
            par(mar=c(4,2,4,5))
            plot(xseq,t.dat[,m.names[1]],ylim=c(0,hbc),xlab="Time (min)",main=lmain,type="n", xaxt="n",yaxt="n",xlim=c(min(xseq)-1.5,max(xseq)))#-sf
            axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
            axis(2, 1.4, )
            text(rep(0,length(m.names)),seq(1,length(m.names))*sf+t.dat[1,m.names],m.names,cex=.8*bcex,col=cols,pos=2)

        ## Tool for adding window region labeling
            if(length(wr) > 0){
                #levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
                x1s <- tapply(xseq,as.factor(wr),min)[levs]
                x2s <- tapply(xseq,as.factor(wr),max)[levs]
                y1s <- rep(-.3,length(x1s))
                y2s <- rep(hbc+.2,length(x1s))
                rect(x1s,y1s,x2s,y2s,col=levs.cols,border="black")
                cpx <- xseq[match(levs,wr)+round(table(wr)[levs]/2,0)]
                offs <- nchar(levs)*.5
                par(xpd=TRUE)
                text(dat$t.dat[match(levs,wr),"Time"],rep(c((sf*.7)/5,(sf*.7)),length=length(levs)),levs,pos=4,offset=0,cex=bcex)#,offset=-offs}
                par(xpd=FALSE)
            }
        
        ## Tool for adding line, point and picture to the plot
            for(i in 1:length(m.names)){
                ypos<-t.dat[,m.names[i]]+i*sf
                lines(xseq,ypos, lty=1,col=cols[i],lwd=lw)
                points(xseq,ypos,pch=16,col=cols[i],cex=.3)
                if(!is.null(snr)){
                    pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
                    pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
                    points(xseq[pp1],t.dat[pp1,m.names[i]]+i/10,pch=1,col=cols[i])
                    points(xseq[pp2],t.dat[pp2,m.names[i]]+i/10,pch=0,col=cols[i])
                }
            }
            par(xpd=TRUE)
            if(!is.null(dat$c.dat[m.names, "area"])){rtag<-"area";rtag <- round(dat$c.dat[m.names,rtag], digits=0)
            text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],paste(rtag),cex=.9*bcex,col=cols,pos=4)}

            if(!is.null(dat$c.dat[m.names, "mean.gfp"])){rtag2<-"mean.gfp";rtag2 <- round(dat$c.dat[m.names,rtag2], digits=0)
            text(rep(max(xseq)*1.04,length(m.names)),seq(1,length(m.names))*sf+(t.dat[nrow(t.dat),m.names]),paste(rtag2),cex=.9*bcex,col="springgreen3",pos=4)}

            if(!is.null(dat$c.dat[m.names, "mean.gfp.1"])){rtag2<-"mean.gfp.1";rtag2 <- round(dat$c.dat[m.names,rtag2], digits=0)
            text(rep(max(xseq)*1.04,length(m.names)),seq(1,length(m.names))*sf+(t.dat[nrow(t.dat),m.names]),paste(rtag2),cex=.9*bcex,col="springgreen3",pos=4)}

            if(!is.null(dat$c.dat[m.names, "mean.tritc"])){rtag3<-"mean.tritc";rtag3 <- round(dat$c.dat[m.names,rtag3], digits=0)
            text(rep(max(xseq)*1.08,length(m.names)),seq(1,length(m.names))*sf+(t.dat[nrow(t.dat),m.names]),paste(rtag3),cex=.9*bcex,col="red1",pos=4)}
            
        if(is.null(img)){img<-dat$img.gtd}
        if(pic.plot==TRUE & length(m.names)<5){
            pic.pos<-list()
            for(i in 1:length(m.names)){
                ypos<-t.dat[,m.names[i]]+i*sf
                pic.pos[[i]]<-mean(ypos)}

                for(i in 1:length(m.names)){		
                    zf<-20
                    x<-dat$c.dat[m.names[i],"center.x"]
                    left<-x-zf
                    if(left<=0){left=0; right=2*zf}
                    right<-x+zf
                    if(right>=2048){left=2048-(2*zf);right=2048}
                    
                    y<-dat$c.dat[m.names[i],"center.y"]
                    top<-y-zf
                    if(top<=0){top=0; bottom=2*zf}
                    bottom<-y+zf
                    if(bottom>=2048){top=2048-(2*zf);bottom=2048}
                    
                    par(xpd=TRUE)
                    xleft<-max(dat$t.dat[,1])*1.05
                    xright<-max(dat$t.dat[,1])*1.13
                    ytop<-pic.pos[[i]]+(.06*hb)
                    ybottom<-pic.pos[[i]]-(.06*hb)
                    rasterImage(img[top:bottom,left:right,],xleft,ytop,xright,ybottom)
                }
            }
        else{multi.pic.zoom(dat, m.names,img=img, plot.new=T)}
        }
    #return(pic.pos)
}

# LinesEvery With all inputs into a single window, excpet XY plot
LinesEvery.4 <- function(dat,m.names, img=NULL,pic.plot=TRUE, zf=NULL, t.type=FALSE, snr=NULL,lmain="",cols=NULL, levs=NULL, levs.cols="grey90",m.order=NULL,rtag=NULL,rtag2=NULL,rtag3=NULL,plot.new=T,sf=.7,lw=.9,bcex=.6,p.ht=7,p.wd=10){
    require(png)
    #if(blc){t.dat<-dat$blc}
    if(t.type){t.type<-menu(names(dat));t.dat<-dat[[t.type]]}# if trace type is empty select the data, you would like your trace to be
    else{t.dat<-dat$blc}
    wr<-dat$w.dat[,2]
    if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")}
    else{levs<-levs}
    m.names <- intersect(m.names,names(t.dat))
    xseq <- t.dat[,1]
    hbc <- length(m.names)*sf+max(t.dat[,m.names])
    hb <- ceiling(hbc)
    library(RColorBrewer)
    
    ## Tool for Sorting cells based on c.dat collumn name
        if(length(m.names) > 0)
        {

            if(!is.null(m.order)){	
                tmp<-dat$c.dat[m.names,]
                n.order<-tmp[order(tmp[,m.order]),]
                m.names <- row.names(n.order)
            }
            else{
                #m.pca <- prcomp(t(t.dat[,m.names]),scale=F,center=T)
                #morder <- m.pca$x[,1] * c(1,-1)[(sum(m.pca$rot[,1]) < 0)+1]
                #m.names <- m.names[order(m.pca$x[,1],decreasing=sum(m.pca$rot[,1]) < 0)]
                #m.names <- m.names[order(morder)]
                m.names<-m.names
            }
    ### Picture Plotting!		
        #if(XY.plot==T){cell.zoom.2048(dat, cell=m.names,img=img, cols="white",zoom=F, plot.new=T)}
        ## Tool for color labeleing
            if(is.null(cols)){
                #cols <- rainbow(length(m.names),start=.55)
                cols <-brewer.pal(8,"Dark2")
                cols <- rep(cols,ceiling(length(m.names)/length(cols)))
                cols <- cols[1:length(m.names)]
            } 
        ## Tool for single color labeling
            else {cols<-cols
                cols <- rep(cols,ceiling(length(m.names)/length(cols)))
                cols <- cols[1:length(m.names)]
            }
            
            if(plot.new){
                if(length(m.names)>5){dev.new(width=16,height=6);layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(10,6), heights=c(6,6))}
                else(dev.new(width=10,height=6))
            }
            else{
                if(length(m.names)>5){layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(10,6), heights=c(6,6))}
                }
            par(xpd=FALSE,mar=c(4,2,4,5), bty="l")
            plot(xseq,t.dat[,m.names[1]],ylim=c(0,hbc),xlab="Time (min)",main=lmain,type="n", xaxt="n",yaxt="n",xlim=c(min(xseq)-1.5,max(xseq)))#-sf
            bob<-dev.cur()
            axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
            axis(2, 1.4, )
            text(rep(0,length(m.names)),seq(1,length(m.names))*sf+t.dat[1,m.names],m.names,cex=.8*bcex,col=cols,pos=2)

        ## Tool for adding window region labeling
            if(length(wr) > 0){
                #levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
                x1s <- tapply(xseq,as.factor(wr),min)[levs]
                x2s <- tapply(xseq,as.factor(wr),max)[levs]
                y1s <- rep(-.3,length(x1s))
                y2s <- rep(hbc+.2,length(x1s))
                rect(x1s,y1s,x2s,y2s,col=levs.cols,border="black")
                cpx <- xseq[match(levs,wr)+round(table(wr)[levs]/2,0)]
                offs <- nchar(levs)*.5
                par(xpd=TRUE)
                text(dat$t.dat[match(levs,wr),"Time"],rep(c((sf*.7)/5,(sf*.7)),length=length(levs)),levs,pos=4,offset=0,cex=bcex)#,offset=-offs}
                par(xpd=FALSE)
            }
        
        ## Tool for adding line, point and picture to the plot
            for(i in 1:length(m.names)){
                ypos<-t.dat[,m.names[i]]+i*sf
                lines(xseq,ypos, lty=1,col=cols[i],lwd=lw)
                points(xseq,ypos,pch=16,col=cols[i],cex=.3)
                if(!is.null(snr)){
                    pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
                    pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
                    points(xseq[pp1],t.dat[pp1,m.names[i]]+i/10,pch=1,col=cols[i])
                    points(xseq[pp2],t.dat[pp2,m.names[i]]+i/10,pch=0,col=cols[i])
                }
            }
            par(xpd=TRUE)
            if(!is.null(dat$c.dat[m.names, "area"])){rtag<-"area";rtag <- round(dat$c.dat[m.names,rtag], digits=0)
            text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],paste(rtag),cex=.9*bcex,col=cols,pos=4)}

            if(!is.null(dat$c.dat[m.names, "mean.gfp"])){rtag2<-"mean.gfp.bin";rtag2 <- round(dat$bin[m.names,rtag2], digits=0)
            text(rep(max(xseq)*1.04,length(m.names)),seq(1,length(m.names))*sf+(t.dat[nrow(t.dat),m.names]),paste(rtag2),cex=.9*bcex,col="springgreen3",pos=4)}

            if(!is.null(dat$c.dat[m.names, "mean.gfp.1"])){rtag2<-"mean.gfp.1";rtag2 <- round(dat$c.dat[m.names,rtag2], digits=0)
            text(rep(max(xseq)*1.04,length(m.names)),seq(1,length(m.names))*sf+(t.dat[nrow(t.dat),m.names]),paste(rtag2),cex=.9*bcex,col="springgreen3",pos=4)}

            if(!is.null(dat$c.dat[m.names, "mean.tritc"])){rtag3<-"mean.tritc.bin";rtag3 <- round(dat$bin[m.names,rtag3], digits=0)
            text(rep(max(xseq)*1.08,length(m.names)),seq(1,length(m.names))*sf+(t.dat[nrow(t.dat),m.names]),paste(rtag3),cex=.9*bcex,col="red1",pos=4)}
            
        if(is.null(img)){img<-dat[[select.list(grep("img",names(dat), value=T))]]}
        if(pic.plot==TRUE & length(m.names)<=5){
            pic.pos<-list()
            for(i in 1:length(m.names)){
                ypos<-t.dat[,m.names[i]]+i*sf
                pic.pos[[i]]<-mean(ypos)}

                for(i in 1:length(m.names)){		
                    #if(dat$bin[m.names[1],"mean.gfp.bin"]!=1 & dat$bin[m.names[1],"mean.tritc.bin"]!=1){img.p<-dat$img.gtd #if the cell is neither red or green, then make the img to plot img.gtd
                    #}else{img.p<-img} 
                    img.p<-img
                    
                    if(is.null(zf)){zf<-20}else{zf<-zf}
                    x<-dat$c.dat[m.names[i],"center.x"]
                    left<-x-zf
                    if(left<=0){left=0; right=2*zf}
                    right<-x+zf
                    if(right>=2048){left=2048-(2*zf);right=2048}
                    
                    y<-dat$c.dat[m.names[i],"center.y"]
                    top<-y-zf
                    if(top<=0){top=0; bottom=2*zf}
                    bottom<-y+zf
                    if(bottom>=2048){top=2048-(2*zf);bottom=2048}
                    
                    par(xpd=TRUE)
                    xleft<-max(dat$t.dat[,1])*1.05
                    xright<-max(dat$t.dat[,1])*1.13
                    ytop<-pic.pos[[i]]+(.06*hb)
                    ybottom<-pic.pos[[i]]-(.06*hb)
                    if(length(dim(img))>2){rasterImage(img.p[top:bottom,left:right,],xleft,ytop,xright,ybottom)
                    }else{rasterImage(img.p[top:bottom,left:right],xleft,ytop,xright,ybottom)}
                }
            }
        else{
            par(mar=c(0,0,0,0))
            plot(0,0,xlim=c(0,6), ylim=c(0,6), xaxs="i",yaxs="i", xaxt='n', yaxt='n')
            tmp.img<-multi.pic.zoom.2(dat, m.names,img=img)
            dev.set(bob) # FUCK THIS!
            rasterImage(tmp.img, 0,0,6,6)
        }
        }
    #return(pic.pos)
}

# LinesEvery same as .4 but has image at begining of trace and moves to pic plot at >10 
LinesEvery.5.1 <- function(dat,m.names, img=dat$img1,pic.plot=TRUE, multi.pic=T,zf=NULL, t.type="mp", snr=NULL,lmain="",cols=NULL, levs=NULL, levs.cols="grey90",  m.order=NULL,rtag=NULL,rtag2=NULL,rtag3=NULL,plot.new=T,sf=1,lw=2,bcex=.6,p.ht=7,p.wd=10, lns=T, pts=F){
    require(png)
    #if(blc){t.dat<-dat$blc}
    if(class(t.type)=="character"){t.dat<-dat[[t.type]]}# if trace type is empty select the data, you would like your trace to be
    else{t.type<-menu(names(dat));t.dat<-dat[[t.type]]}
    wr<-dat$w.dat[,2]
    if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")}
    else{levs<-levs}
    m.names <- intersect(m.names,names(t.dat))
    xseq <- t.dat[,1]
    hbc <- length(m.names)*sf+max(t.dat[,m.names])
    hb <- ceiling(hbc)
    library(RColorBrewer)
    
    ## Tool for Sorting cells based on c.dat collumn name
        if(length(m.names) > 0)
        {
            #if(is.null(pdf.name))
            #	{dev.new(width=14,height=8)}
            #else
                #{if(length(grep("\\.pdf",pdf.name))>0){pdf(pdf.name,width=p.wd,height=p.ht)}else{png(pdf.name,width=1200,height=600)}}

            if(!is.null(m.order)){	
                tmp<-dat$c.dat[m.names,]
                n.order<-tmp[order(tmp[,m.order]),]
                m.names <- row.names(n.order)
            }
            else{
                #m.pca <- prcomp(t(t.dat[,m.names]),scale=F,center=T)
                #morder <- m.pca$x[,1] * c(1,-1)[(sum(m.pca$rot[,1]) < 0)+1]
                #m.names <- m.names[order(m.pca$x[,1],decreasing=sum(m.pca$rot[,1]) < 0)]
                #m.names <- m.names[order(morder)]
                m.names<-m.names
            }
    ### Picture Plotting!		
        #if(XY.plot==T){cell.zoom.2048(dat, cell=m.names,img=img, cols="white",zoom=F, plot.new=T)}
        ## Tool for color labeleing
            if(is.null(cols)){
                #cols <- rainbow(length(m.names),start=.55)
                cols <-brewer.pal(8,"Dark2")
                cols <- rep(cols,ceiling(length(m.names)/length(cols)))
                cols <- cols[1:length(m.names)]
            } 
        ## Tool for single color labeling
            else {cols<-cols
                cols <- rep(cols,ceiling(length(m.names)/length(cols)))
                cols <- cols[1:length(m.names)]
            }
            
            if(multi.pic){
                if(plot.new){
                    if(length(m.names)>10){dev.new(width=16,height=6);layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(10,6), heights=c(6,6))}
                    else(dev.new(width=10,height=6))
                }
                else{
                    if(length(m.names)>10){layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(10,6), heights=c(6,6))}
                    }
            }else{dev.new(width=10,height=6)}
            par(xpd=FALSE,mar=c(4,3,4,10), bty="l")
            plot(xseq,t.dat[,m.names[1]],ylim=c(0,hbc),xlab="Time (min)",main=lmain,type="n", xaxt="n",yaxt="n",xlim=c(min(xseq)-1.5,max(xseq)))#-sf
            bob<-dev.cur()
            axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
            axis(2, 1.4, )
            text(rep(0,length(m.names)),seq(1,length(m.names))*sf+t.dat[1,m.names],m.names, cex=.5,col=cols,pos=2)

        ## Tool for adding window region labeling
            if(length(wr) > 0){
                #levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
                x1s <- tapply(dat$w.dat[,1],as.factor(wr),min)[levs]
                x2s <- tapply(dat$w.dat[,1],as.factor(wr),max)[levs]
                y1s <- rep(-.3,length(x1s))
                y2s <- rep(hbc+.2,length(x1s))
                rect(x1s,y1s,x2s,y2s,col=levs.cols,border="black")
                cpx <- xseq[match(levs,wr)+round(table(wr)[levs]/2,0)]
                offs <- nchar(levs)*.5
                par(xpd=TRUE)
                text(dat$t.dat[match(levs,wr),"Time"],rep(c((sf*.7)*.5,(sf*.7),(sf*.7)/5),length=length(levs)),levs,pos=4,offset=0,cex=bcex*.8)#,offset=-offs}
                par(xpd=FALSE)
            }
        
        ## Tool for adding line, point and picture to the plot
            for(i in 1:length(m.names)){
                ypos<-t.dat[,m.names[i]]+i*sf
                if(lns){lines(xseq,ypos, lty=1,col=cols[i],lwd=lw)}
                if(pts){points(xseq,ypos,pch=16,col=cols[i],cex=.3)}
                if(!is.null(snr)){
                    pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
                    pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
                    points(xseq[pp1],t.dat[pp1,m.names[i]]+i/10,pch=1,col=cols[i])
                    points(xseq[pp2],t.dat[pp2,m.names[i]]+i/10,pch=0,col=cols[i])
                }
            }
            par(xpd=TRUE)
            if(!is.null(dat$c.dat[m.names, "area"])){rtag<-"area";rtag <- round(dat$c.dat[m.names,rtag], digits=0)
            text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],paste(rtag),cex=.9*bcex,col=cols,pos=4)}

            if(!is.null(dat$c.dat[m.names, "mean.gfp.start"])){rtag2<-"mean.gfp.start";rtag2 <- round(dat$c.dat[m.names,rtag2], digits=6)
            text(rep(max(xseq)+xinch(.5),length(m.names)),seq(1,length(m.names))*sf+(t.dat[nrow(t.dat),m.names]),paste(rtag2),cex=.6*bcex,col="darkgreen",pos=4)}

            #if(!is.null(dat$c.dat[m.names, "mean.gfp.end"])){rtag2<-"mean.gfp.end";rtag2 <- round(dat$c.dat[m.names,rtag2], digits=6)
            #text(rep(max(xseq)+xinch(1),length(m.names)),seq(1,length(m.names))*sf+(t.dat[nrow(t.dat),m.names]),paste(rtag2),cex=.6*bcex,col="darkgreen",pos=4)}

            if(!is.null(dat$c.dat[m.names, "mean.tritc.start"])){rtag3<-"mean.tritc.start";rtag3 <- round(dat$c.dat[m.names,rtag3], digits=6)
            text(rep(max(xseq)+xinch(1),length(m.names)),seq(1,length(m.names))*sf+(t.dat[nrow(t.dat),m.names]),paste(rtag3),cex=.6*bcex,col=,pos=4)}
            
            if(!is.null(dat$c.dat[m.names, "mean.tritc.end"])){rtag3<-"mean.tritc.end";rtag3 <- round(dat$c.dat[m.names,rtag3], digits=6)
            text(rep(max(xseq)+xinch(1.5),length(m.names)),seq(1,length(m.names))*sf+(t.dat[nrow(t.dat),m.names]),paste(rtag3),cex=.6*bcex,col=,pos=4)}

        if(is.null(img)){
        img.p<-dat[[select.list(grep("img",names(dat), value=T))]]
        if(is.null(img.p)){img.p<-dat$img1}
        }else{img.p<-img}
        if(is.null(zf)){zf<-20}else{zf<-zf}
        
        #if(pic.plot==TRUE & length(m.names)<=10){
        if(pic.plot==TRUE){
        if(length(m.names)<=10){

            pic.pos<-list()
            for(i in 1:length(m.names)){
                ypos<-t.dat[1,m.names[i]]+i*sf
                pic.pos[[i]]<-ypos}

                for(i in 1:length(m.names)){		
                    #if(dat$bin[m.names[1],"mean.gfp.bin"]!=1 & dat$bin[m.names[1],"mean.tritc.bin"]!=1){img.p<-dat$img.gtd #if the cell is neither red or green, then make the img to plot img.gtd
                    #}else{img.p<-img} 
                    #img.p<-img
                    
                    img.dim<-dim(dat$img1)[1]

                    x<-dat$c.dat[m.names[i],"center.x"]
                    left<-x-zf
                    if(left<=0){left=0; right=2*zf}
                    right<-x+zf
                    if(right>=img.dim){left=img.dim-(2*zf);right=img.dim}
                    
                    y<-dat$c.dat[m.names[i],"center.y"]
                    top<-y-zf
                    if(top<=0){top=0; bottom=2*zf}
                    bottom<-y+zf
                    if(bottom>=img.dim){top=img.dim-(2*zf);bottom=img.dim}
                    
                    par(xpd=TRUE)
                    xleft<-min(dat$t.dat[,1])-xinch(1)
                    xright<-min(dat$t.dat[,1])-xinch(.5)
                    ytop<-pic.pos[[i]]+yinch(.25)
                    ybottom<-pic.pos[[i]]-yinch(.25)
                    
                    tryCatch(rasterImage(img.p[top:bottom,left:right,],xleft,ybottom,xright,ytop),error=function(e) rasterImage(img.p[top:bottom,left:right],xleft,ybottom,xright,ytop))
                }
            }
        else{
            par(mar=c(0,0,0,0))
            plot(0,0,xlim=c(0,6), ylim=c(0,6), xaxs="i",yaxs="i", xaxt='n', yaxt='n')
            tmp.img<-multi.pic.zoom.2(dat, m.names,img=img.p, labs=T, zf=zf, cols=cols)
            dev.set(bob) # FUCK THIS!
            rasterImage(tmp.img, 0,0,6,6)
        }
        }
        }
            #if(!is.null(pdf.name))
            #{dev.off()}

    #return(pic.pos)
}

# LinesEvery same as .4 but has image at begining of trace and moves to pic plot at >10 
#multipi does not work on this.  Instead, if greater than 10, the traces are plotted as a portrait orientation
# Also window labels are rotated on axis and place on the bottom of the plot
# I am also adding two more images to the left side of the plot
#171009 added underline.  Helps to show irreversibility
#171031 added dat.n for the name fotthe experiment
LinesEvery.5 <- function(dat,m.names, img="img1",channel=NULL,pic.plot=TRUE,zf=NULL, t.type="mp.1", snr=NULL,lmain="",cols="black", levs=NULL, levs.cols="grey90", values=NULL,plot.new=T,sf=1,lw=1,bcex=1,p.ht=7,p.wd=10, lns=T, pts=F, underline=T,dat.n=NULL){
    #require(RColorBrewer)
    dat.name<-deparse(substitute(dat))
    if(dat.name=="dat" | dat.name == "tmp.rd" | dat.name == "tmp_rd"){
        dat.name<-dat.n
    }else{
        dat.name<-dat.name
    }
    
    #Trace Selector if t.type is empty.  t.type must be character input
    if(class(t.type)=="character"){
        t.dat<-dat[[t.type]]# if trace type is empty select the data, you would like your trace to be
    }else{
        t.type<-menu(names(dat));t.dat<-dat[[t.type]]
    }
    
    m.names <- intersect(m.names,names(t.dat))
    xseq <- t.dat[,1]
    #upper ylimit 
    hbc <- length(m.names)*sf+max(t.dat[,m.names])

    #Selecting multiple images
    if(is.null(img)){
        img.l<-select.list(grep("img",names(dat), value=T), multiple=T)
    }else{
        img.l<-img
    }

    if(length(m.names) > 0){
        #For pdf output
            #if(is.null(pdf.name))
            #	{dev.new(width=14,height=8)}
            #else
                #{if(length(grep("\\.pdf",pdf.name))>0){pdf(pdf.name,width=p.wd,height=p.ht)}else{png(pdf.name,width=1200,height=600)}}
        ## Tool for addind value tags displayed on the right side of trace
        #See line 3016 for where values come into play
        #values<-c("area", "mean.gfp.start", "mean.gfp.end" "mean.tritc.start", "mean.tritc.end")
            if(is.null(values)){
                values<-c("area")
            }else{values<-values}

        ## Tool for color labeleing
        ## Tool for single color labeling
            if(cols=="brew.pal"){
                #cols <- rainbow(length(m.names),start=.55)
                require(RColorBrewer)
                cols <-brewer.pal(8,"Dark2")
                cols <- rep(cols,ceiling(length(m.names)/length(cols)))
                cols <- cols[1:length(m.names)]
            }
            if(cols=="rainbow"){
                cols<-rainbow(length(m.names),start=.7,end=.1)
            }
            if(cols=="topo"){
                cols<-topo.colors(length(m.names))
            }else{		
                cols<-cols
                cols <- rep(cols,ceiling(length(m.names)/length(cols)))
                cols <- cols[1:length(m.names)]
            }
            if(plot.new){
                if(length(m.names)>10){dev.new(width=10+length(img)+(length(values)*.6),height=12)}
                else(dev.new(width=10+length(img)+(length(values)*.6),height=8))
            }
 
            xinch(length(img))
            par(xpd=FALSE,mai=c(2,.5+(.5*length(img.l)), 1, 0.6*length(values)), bty="l")
            plot(xseq,t.dat[,m.names[1]],ylim=c(0,hbc),xlab="",main=lmain,type="n",yaxt="n",xlim=c(min(xseq)-1.5,max(xseq)), ylab="")#-sf
            #axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
            #axis(2, 1.4, )
            #Label cell names
            text(rep(0,length(m.names))-xinch(.1),seq(1,length(m.names))*sf+t.dat[1,m.names],m.names, cex=.5,col=cols,pos=3)

        ## Tool for adding window region labeling
            if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
            }else{levs<-levs}
            wr<-dat$w.dat$wr1
            if(length(wr) > 0){
                x1s <- tapply(dat$w.dat[,1],as.factor(wr),min)[levs]
                x2s <- tapply(dat$w.dat[,1],as.factor(wr),max)[levs]
                y1s <- rep(par("usr")[3],length(x1s))
                y2s <- rep(par("usr")[4],length(x1s))
                rect(x1s,y1s,x2s,y2s,col=levs.cols,border="black")
                par(xpd=TRUE)
                #text(x1s-xinch(.1),par("usr")[3]-yinch(1),levs,cex=.8*bcex, srt=90)	
                #dat$t.dat[match(levs,wr),"Time"]
                levs.loc<-tapply(dat$w.dat[,"Time"],as.factor(wr),mean)[levs]
				levs_cex <- nchar(levs)
				levs_cex[ levs_cex <= 12*1.3 ] <- 1
				levs_cex[ levs_cex > (12*1.3) ] <- 12/levs_cex[ levs_cex>(12*1.3) ]*1.3
                text(levs.loc,par("usr")[3],levs,pos=3,offset=-4.3,cex=levs_cex, srt=90)	
                par(xpd=FALSE)
            }
        
        ## Tool for adding line, point and picture to the plot
            for(i in 1:length(m.names)){
                ypos<-t.dat[,m.names[i]]+i*sf
                if(lns){lines(xseq,ypos, lty=1,col=cols[i],lwd=lw)}
                if(pts){points(xseq,ypos,pch=16,col=cols[i],cex=.3)}
                if(!is.null(snr)){
                    pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
                    pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
                    points(xseq[pp1],t.dat[pp1,m.names[i]]+i/10,pch=1,col=cols[i])
                    points(xseq[pp2],t.dat[pp2,m.names[i]]+i/10,pch=0,col=cols[i])
                }
                if(underline){abline(h=min(ypos), col="black")}else{}
            }
            par(xpd=TRUE)
            
        ## Tool for adding Value info on right side of trace
            placement<-seq(0,length(values),.5)
            digits<-c(0,rep(4,length(values)))
            text(max(xseq)+xinch(placement[1:length(values)]), par("usr")[4]+yinch(.2), pos=4, values,cex=bcex*.75, srt=30)
            
            for(i in 1:length(values)){
                if(!is.null(dat$c.dat[m.names, values[i]])){
                rtag<-values[i]
                rtag <- round(dat$c.dat[m.names,rtag], digits=digits[i])
                text(
                    rep(max(xseq)+xinch(placement[i]),length(m.names)),
                    seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],
                    paste(rtag),
                    cex=.65*bcex,
                    col=cols,
                    pos=4)
                }
            }
            
        ##Tool for adding images to the left side of the plot
        if(is.null(zf)){
            zf<-20
        }else{zf<-zf}
        
        pic.pos<-list()
        for(i in 1:length(m.names)){
            ypos<-t.dat[1,m.names[i]]+i*sf
            pic.pos[[i]]<-ypos
        }

        xinchseq1<-seq(1,5,.5)
        xinchseq2<-seq(.5,5,.5)
        
        if(is.null(channel)){channel<-rep(list(c(1:3)),length(img.l))
        }else{channel<-channel}

        for(j in 1:length(img.l)){
            for(i in 1:length(m.names)){		
                img.dim<-dim(dat$img1)[1]
                x<-dat$c.dat[m.names[i],"center.x"]
                left<-x-zf
                if(left<=0){
                    left=0
                    right=2*zf
                }
                
                right<-x+zf
                if(right>=img.dim){
                    left=img.dim-(2*zf)
                    right=img.dim
                }
                    
                y<-dat$c.dat[m.names[i],"center.y"]
                top<-y-zf
                if(top<=0){
                    top=0
                    bottom=2*zf
                }
                bottom<-y+zf
                
                if(bottom>=img.dim){
                    top=img.dim-(2*zf)
                    bottom=img.dim
                }
                
                par(xpd=TRUE)
                xleft<-min(dat$t.dat[,1])-xinch(xinchseq1[j])
                xright<-min(dat$t.dat[,1])-xinch(xinchseq2[j])
                ytop<-pic.pos[[i]]+yinch(.25)
                ybottom<-pic.pos[[i]]-yinch(.25)
                
                tryCatch(
                    rasterImage(dat[[img.l[j]]][top:bottom,left:right,channel[[j]]],xleft,ybottom,xright,ytop),
                    error=function(e) rasterImage(dat[[img.l[j]]][top:bottom,left:right],xleft,ybottom,xright,ytop)
                )
            }
        }
    }

        tryCatch(
			legend(x=par("usr")[2]-xinch(1.2), y=par("usr")[3]-yinch(1.6), xpd=TRUE, inset=c(0,-.14), bty="n", cex=.7, legend=dat.name),
		error=function(e) NULL)


        #if(!is.null(pdf.name))
        #{dev.off()}

    #return(pic.pos)
}

#How to display single or multiple window regions as specified by you
#performs pam analysis around as many mediods as you want
#displays the information as a heat map with red represeenting most populace group
# and white as least populace
#legend
#xlim= Logical added to have option to group traces around window regions
#medios= Logical.  If true then groupw ill be split into subset.n groups
LevsViewer <- function(dat,m.names=NULL, ylim=c(0,1.4), xlim=F, mediods=T, min.group=0,linez=F,subset.n=5,img=NULL, pic.plot=FALSE, t.type="mp.1",lmain="",cols=NULL, levs=NULL, levs.cols="grey90",plot.new=T,lw=.9,bcex=.8,opacity=3){
    # if trace type is empty select the data, you would like your trace to be
        if(is.null(t.type)){
            t.type<-"t.dat"
            t.dat<-dat[[t.type]]
        }else{
            t.dat<-dat[[t.type]]
        }
        
    ## select the region to plot from window region
        levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    ## This is how i determine the x region to plot
        if(xlim){
            plot.region<-select.list(levs, multiple=T)
        }else{
            plot.region<-levs
        }
        
    ## What cells would you like to plot?
    # if you do not have cells entered, then all cells will be selected.
        if(is.null(m.names)){
            m.names<-dat$c.dat$id
        }else{
            m.names<-m.names}
    ## Create new plot if you have nothing selected.
        if(plot.new){dev.new(height=5,width=5*length(plot.region))}

    ## If you have decided to select window regions based on levs, then this code will select
    # the minimun x value from all selected levs regions, and the ultimate maximun region from the selected levs region
    ## If not, then the maximun and minimun x regions will be the maximun and minimun region of the entire time frame
        if(xlim){
            x.min<-which(t.dat$Time==
                min(tapply(t.dat[,"Time"], as.factor(dat$w.dat$wr1), min)[plot.region]))
            x.max<-which(t.dat[,"Time"]==
                max(tapply(t.dat[,"Time"], as.factor(dat$w.dat$wr1), max)[plot.region]))
        }else{
            x.min<-min(t.dat[,"Time"])
            x.min<-which(t.dat[,"Time"]==x.min)
            x.max<-max(t.dat[,"Time"])
            x.max<-which(t.dat[,"Time"]==x.max)
        }

    ## Y limits will be automatically decided if this input is NULL
    # if not, the y limits can be specified via c()
        if(is.null(ylim)){
            ylim<-c(min(t.dat[,m.names]),max(t.dat[,m.names]))
        }else{
            ylim<-ylim
        }
        
    ## xseq is an arguement input for the plotting function
        xseq<-t.dat[x.min:x.max,"Time"]
    ## I'm not sure why i added this to the code.  Perhaps I have had issues in the past
        m.names <- intersect(m.names,names(t.dat))
    ## Make sure this is loaded for pretty colors, although we are currently using heat.colors	
        library(RColorBrewer)

    ###PLOTTING####
    ##parameters for plot
        par(xpd=FALSE,mar=c(4,2,4,5), bty="l", bg="gray70")
        plot(xseq,t.dat[x.min:x.max,m.names[1]],ylim=ylim,xlab="Time (min)",main=lmain,type="n")#-sf
        axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), .5))
        axis(2, tick=T, )
        abline(h=seq(0,max(ylim),.2), lty=3, lwd=.1	)
        #abline(v=seq(
        #	floor(t.dat$Time[x.min]),
        #	ceiling(t.dat$Time[x.max]),.5), 
        #	lty=3, lwd=.1
        #)

    ## Add all lines to the plot with light opacity.
    # in the future this may need a type of calculation to allow for a realtionship between 
    # opacity and number of lines
        if(linez){
            for(i in 1:length(m.names)){
                ypos<-t.dat[x.min:x.max,m.names[i]]
                color<-rgb(1,1,1, opacity, maxColorValue=10)
                lines(xseq,ypos, lty=1,col=color,lwd=lw)
                #points(xseq,ypos,pch=16,col=color,cex=.3)
            }
        }
        
    ##Tool for adding trace difference averages using Partitioning about mediods
        if(mediods){
            library(cluster)
            if(subset.n>=(length(m.names)/2)){
                subset.n<-floor(length(m.names)/4)
            }else{subset.n<-subset.n}
            
            pam5 <- pam(t(t.dat[x.min:x.max,m.names]),k=subset.n)
            s.names <- row.names(pam5$medoids)
            pam5.tab <- table(pam5$clustering)
            #Create Labels for the PAM groupings
            #tags <- paste(paste("#",names(pam5.tab),sep=""),as.vector(pam5.tab),sep=":")
            group.means<-list()
            group.names<-list()
            
        ## Create naming for PAM clusters
            for(i in 1:subset.n){
                x.names<-names(which(pam5$clustering==i, arr.ind=T))
                #group.info<-paste(i,":",length(x.names), sep="")
                group.info<-length(x.names)
                group.names[[i]]<-x.names
                names(group.names)[i]<-group.info
            }

        #only select groups that have more than 2 traces
        #this takes some magic to happen
            bg<-summary(group.names)
            bg[,1]<-as.numeric(bg[,1])
            
            if(subset.n<=length(m.names)){
                bg.big<-which(as.numeric(bg[,1])>min.group)
            }else{
                bg.big<-which(as.numeric(bg[,1])>min.group)
            }
            
            if(length(bg.big)>1){
                bg.2<-bg[bg.big,]
                bg.2<-bg.2[order(as.numeric(bg.2[,1]),decreasing=T),]
                bg.2.names<-row.names(bg.2)	
                group.names<-group.names[bg.2.names]
            }else{
                bg.2<-bg
                bg.2<-bg.2[order(as.numeric(bg.2[,1]),decreasing=T),]
                bg.2.names<-row.names(bg.2)	
                group.names<-group.names[bg.2.names]
            }
            
            cols <-brewer.pal(8,"Dark2")
            cols<-rainbow(length(group.names))
            #cols<-rep(cols,length(group.names))
            #cols <- rep(cols,ceiling(length(s.names)/length(cols)))
            #cols <- cols[1:length(s.names)]
            #cols<-heat.colors(length(group.names))
            #cols<-rev(cols)
            
            for(i in 1:length(group.names)){
                if(length(group.names[[i]])>1){
                    lines(xinch(i*0)+xseq, yinch(i*.15)+apply(t.dat[x.min:x.max,group.names[[i]]],1,mean), col=cols[i], lwd=2)
                }else{
                    lines(xinch(i*0)+xseq, yinch(i*.15)+t.dat[x.min:x.max,group.names[[i]]], col=cols[i], lwd=2)
                }
            }
            
            legend("topright",legend=rev(names(group.names)),title="Cell Total", cex=.8,lty=1,lwd=2, bty="", col=rev(cols))
        }
    ## Tool for adding window region labeling
    wr<-dat$w.dat[x.min:x.max,2]
    if(length(wr) > 0){
        #levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
        x1s <- tapply(xseq,as.factor(wr),min)[plot.region]
        #abline(v=x1s, lwd=1.2)
        par(xpd=TRUE)
        text(x1s,rep(.9,length=length(plot.region)),plot.region,pos=4,offset=0,cex=bcex)#,offset=-offs}
        par(xpd=FALSE)
    }
    

    if(is.null(img)){img<-dat$img2}
    if( (pic.plot) & (length(m.names)>=5)){
        dev.new()
        par(mar=c(0,0,0,0))
        plot(0,0,xlim=c(0,6), ylim=c(0,6), xaxs="i",yaxs="i", xaxt='n', yaxt='n')
        multi.pic.zoom(dat, m.names,img=img)
    }
    return(group.names)
}

LinesSome.2 <- function(dat,m.names,snr=NULL,lmain="",pdf.name=NULL,morder=NULL,subset.n=5,sf=1,lw=3,bcex=1)
{
    library(cluster)
    t.dat<-dat$t.dat
    wr<-dat$w.dat[,2]
    levs<-unique(as.character(wr))[-1]
    
    if(length(m.names) < subset.n)
    {stop("group size lower than subset size")}
    pam5 <- pam(t(t.dat[,m.names]),k=subset.n)
    s.names <- row.names(pam5$medoids)
    if(!is.null(morder))
    {
        names(morder) <- m.names
        morder <- morder[s.names]
        }
    pam5.tab <- table(pam5$clustering)
    tags <- paste(paste("#",names(pam5.tab),sep=""),as.vector(pam5.tab),sep=":")
    LinesEvery(t.dat,snr,s.names,wr,levs,lmain,pdf.name,morder,rtag=tags,sf,lw,bcex)
    return(pam5$clustering)
}

TraceSelect <- function(dat,m.names,blc=NULL,snr=NULL,wr=NULL,levs=NULL,lmain="",m.order=NULL,rtag=NULL,rtag2=NULL,rtag3=NULL){
    if(!is.null(blc)){t.dat<-dat$blc}
    else{t.dat<-dat$t.dat}
    
    if(is.null(wr)){wr<-dat$w.dat[,2]}
    sf <- .2
    bcex<-1
    library(RColorBrewer)
    m.names <- intersect(m.names,names(t.dat))
    lwds <- 3
    if(length(m.names) > 0)
    {
    
    xseq <- t.dat[,1]
    cols <-brewer.pal(8,"Dark2")
    cols <- rep(cols,ceiling(length(m.names)/length(cols)))
    cols <- cols[1:length(m.names)]
    dev.new(width=14,height=8)
    
    if(!is.null(m.order)){
        (tmp<-dat$c.dat[m.names,])
        (n.order<-tmp[order(tmp[,m.order]),])
        (m.names <- row.names(n.order))
        }
    else{
        m.pca <- prcomp(t(t.dat[,m.names]),scale=F,center=T)
        m.names <- m.names[order(m.pca$x[,1],decreasing=sum(m.pca$rot[,1]) < 0)]
        }
    
    
    hbc <- length(m.names)*sf+max(t.dat[,m.names])
    hb <- ceiling(hbc)
    
    plot(xseq,t.dat[,m.names[1]],ylim=c(-sf,hbc),xlab="Time (min)",ylab="Ratio with shift",main=lmain,type="n", xaxt="n")
    axis(1, at=seq(0, length(t.dat[,1]), 5))

    if(length(wr) > 0)
    {
        if(is.null(levs)){levs <- setdiff(unique(wr),"")}
        x1s <- tapply(xseq,as.factor(wr),min)[levs]
        x2s <- tapply(xseq,as.factor(wr),max)[levs]
        y1s <- rep(-.3,length(x1s))
        y2s <- rep(hbc+.2,length(x1s))
        rect(x1s,y1s,x2s,y2s,col="lightgrey")
        text(xseq[match(levs,wr)],rep(-.1,length(levs)),levs,pos=4,offset=0,cex=1)
    }
    x.sel <- NULL
    xs <-c(rep(0,length(m.names)),c(.1,.1,.1))
    ys <- seq(1,length(m.names))*sf+t.dat[1,m.names]
    ys <- as.vector(c(ys,c(sf,0,-sf)))
    #    xs[(length(xs)-2):length(xs)] <- c(0,5,10)
    p.names <- c(m.names,"ALL","NONE","FINISH")
    done.n <- length(p.names)
    none.i <- done.n-1
    all.i <- none.i-1
    p.cols <- c(cols,c("black","black","black"))
    for(i in 1:length(m.names))
    {
        lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lwds)
        if(!is.null(snr))
        {
        pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
        pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
        points(xseq[pp1],t.dat[pp1,m.names[i]]+i*sf,pch=1,col=cols[i])
        points(xseq[pp2],t.dat[pp2,m.names[i]]+i*sf,pch=0,col=cols[i])
        }
    }
    text(x=xs,y=ys,labels=p.names,pos=2,cex=.7,col=p.cols)
    points(x=xs,y=ys,pch=16,col=p.cols)
        
        if(is.null(rtag)){
        if(!is.null(m.order)){
            rtag <- dat$c.dat[m.names,m.order]
            text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag,cex=.8*bcex,col=cols,pos=4)
        }}
        else{
            rtag <- round(dat$c.dat[m.names,rtag], digits=0)
            text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag,cex=.8*bcex,col=cols,pos=4)
         }

        if(!is.null(rtag2)){
            (rtag2 <- round(dat$c.dat[m.names,rtag2], digits=0))
            text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag2,cex=.8*bcex,col=cols,pos=3)
        }
        if(!is.null(rtag3)){
            rtag3 <- round(dat$c.dat[m.names,rtag3], digits=0)
            text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag3,cex=.8*bcex,col=cols,pos=1)
        }
    
    click.i <- 1    
    while(click.i != done.n)
    {
        click.i <- identify(xs,ys,n=1,plot=F)
        if(click.i < (length(m.names)+1) & click.i > 0)
        {
            i <- click.i
            if(is.element(i,x.sel))
            {
                lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lwds)
                x.sel <- setdiff(x.sel,i)
            }
                else
                {
                lines(xseq,t.dat[,m.names[i]]+i*sf,col="black",lwd=lwds)
                #lines(xseq,t.dat[,m.names[i]]+i*sf,col="white",lwd=2,lty=2)
                x.sel <- union(x.sel,i)
            }
        }
        if(click.i == none.i)
        {
            x.sel <- NULL
            for(i in 1:length(m.names))
            {
                lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lwds)
            }
        }
        if(click.i == all.i)	
        {
            x.sel <- seq(1,length(m.names))
            for(i in 1:length(m.names))
            {
                lines(xseq,t.dat[,m.names[i]]+i*sf,col="black",lwd=lwds)
            }
            
        }
    }
    dev.off()
    return(m.names[x.sel])
    }
}

#TraceSelectLarge takes a large list of traces
#subsets it and passes each on to Trace select
TraceSelectLarge <- function(t.dat,snr=NULL,m.names,wr,levs=NULL,lmain="",subset.n=10,rtag=NULL){
    sel.names <- NULL
    s <- ceiling(length(m.names)/subset.n)
    for(i in 1:s)
    {
        x1 <- (i-1)*subset.n+1
        x2 <- min(length(m.names),x1+subset.n)
        x.sel <- TraceSelect(t.dat,snr,m.names[x1:x2],wr,levs,lmain,rtag[x1:x2])
        sel.names <- union(sel.names,x.sel)
    }		
    return(sel.names)	
}

LinesStack <- function(dat,m.names,lmain="",levs=NULL, plot.new=TRUE,bcex=.7, sf=.2, subset.n=5){
    if(plot.new){dev.new(width=10,height=6)}
    if(length(m.names)>subset.n){
        t.dat<-dat$t.dat
        wr<-dat$w.dat[,2]
        if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")}
        else{levs<-levs}
        m.names <- intersect(m.names,names(t.dat))
        hbc <- subset.n*sf+max(t.dat[,m.names])
        xseq <- t.dat[,1]
        library(RColorBrewer)
        par(mar=c(4,2,4,4))
        hbc <- (subset.n*(.8*sf))+max(t.dat[,m.names])
        #ylim <- c(-.1,2.5)
        ylim<-c(-.1,hbc)
        plot(xseq,t.dat[,m.names[1]],ylim=ylim,xlab="Time (min)",main=lmain,type="n", xaxt="n",xlim=c(min(xseq)-1.5,max(xseq)+25))#-sf
        axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
        ## Tool for adding window region labeling
        if(length(wr) > 0){
            #levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
            x1s <- tapply(xseq,as.factor(wr),min)[levs]
            x2s <- tapply(xseq,as.factor(wr),max)[levs]
            y1s <- rep(min(ylim)-.2,length(x1s))
            y2s <- rep(max(ylim)+.2,length(x1s))
            rect(x1s,y1s,x2s,y2s,col="grey90",border="black")
            text(dat$t.dat[match(levs,wr),"Time"],rep(c(-.05, abs(min(ylim))),length=length(levs)),levs,cex=bcex,offset=0, pos=4)#,offset=-offs}
        }
        blc<-dat$blc
        ## Tool for adding line and point plot for all lines
            #matlines(xseq, blc[,m.names], col=rgb(0,0,0,3, maxColorValue=100), lwd=.01)
            #matpoints(xseq, blc[,m.names], col=rgb(0,0,0,3, maxColorValue=100), pch=16, cex=.03)
        
        #cols <- rainbow(length(m.names),start=.55)
     
        library(cluster)
        blc<-dat$blc
        pam5 <- pam(t(blc[,m.names]),k=subset.n)
        s.names <- row.names(pam5$medoids)
        pam5.tab <- table(pam5$clustering)
        #tags <- paste(paste("#",names(pam5.tab),sep=""),as.vector(pam5.tab),sep=":")
        group.means<-list()
        group.names<-list()
        for(i in 1:subset.n){
            x.names<-names(which(pam5$clustering==i, arr.ind=T))
            group.names[[i]]<-x.names
            group.means[i]<-paste(
            round(mean(dat$c.dat[x.names, "area"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "area"]), digits=0)," : ",
            round(mean(dat$c.dat[x.names, "mean.gfp"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.gfp"]), digits=0)," : ",	
            round(mean(dat$c.dat[x.names, "mean.tritc"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.tritc"]), digits=0), sep="")
            # adding standard deviation,"\u00b1",round(sd(dat$c.dat[x.names, "area"]), digits=0)," : ", 
        }
        
        
        tags <- paste(as.vector(pam5.tab),":",group.means)
        info<-pam5$clustering
        
        ## Tool For adding color to selected Traces
        cols <-brewer.pal(8,"Dark2")
        cols <- rep(cols,ceiling(length(s.names)/length(cols)))
        cols <- cols[1:length(s.names)]

        ## Tool for adding labeling for single line within stacked traces
        for(i in 1:length(s.names)){
            lines(xseq, blc[,s.names[i]]+i*sf, col=cols[i], lwd=.2)
            points(xseq, blc[,s.names[i]]+i*sf, col=cols[i], pch=16, cex=.02)
                    matlines(xseq, blc[,names(which(info==i, arr.ind=T))]+i*sf, col=rgb(0,0,0,50, maxColorValue=100), lwd=.01)
            text(x=min(blc[,1]), y=blc[nrow(t.dat),s.names[i]]+i*sf, labels=s.names[i], col=cols[i], pos=2, cex=bcex)
            text(x=max(blc[,1]), y=blc[nrow(t.dat),s.names[i]]+i*sf, labels=tags[i], col=cols[i], pos=4, cex=bcex)
        }
    #return(pam5$clustering)	
    return(group.names)
    }
} 
### Best Linesstack Yet
# stack traces accorinding to # of groups defined
# uses pam clustering
# and bp.func2
LinesStack.2.1 <- function(dat,m.names=NULL,lmain="",levs=NULL, plot.new=TRUE,bcex=.7, sf=.1, subset.n=5, img=NULL, cols=NULL,bp.param=NULL){
    graphics.off()
    if(is.null(img)){img<-dat$img1}
    if(is.null(m.names)){m.names<-dat$c.dat$id}
    if(plot.new){dev.new(width=10,height=6)}
    if(length(m.names)>subset.n){
        
        t.dat<-dat$t.dat
        wr<-dat$w.dat[,2]
        if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")}
        else{levs<-levs}
        m.names <- intersect(m.names,names(t.dat))
        hbc <-max(t.dat[,m.names])*subset.n *.643
        xseq <- t.dat[,1]
        library(RColorBrewer)
        par(mar=c(4,2,4,4),bty="l")
        #hbc <- (subset.n*(.8*sf))+max(t.dat[,m.names])
        #ylim <- c(-.1,2.5)
        #ylim<-c(-.1,hbc)
        ylim<-c(0,subset.n+subset.n*sf)
        plot(xseq,t.dat[,m.names[1]],ylim=ylim,xlab="Time (min)",main=lmain,type="n", xaxt="n",xlim=c(min(xseq)-1.5,max(xseq)+25))#-sf
        axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
        
        ## Tool for adding window region labeling
        if(length(wr) > 0){
            levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
            x1s <- tapply(xseq,as.factor(wr),min)[levs]
            #x2s <- tapply(xseq,as.factor(wr),max)[levs]
            #y1s <- rep(min(ylim)-.2,length(x1s))
            #y2s <- rep(max(ylim)+.2,length(x1s))
            #rect(x1s,y1s,x2s,y2s,col="grey90",border="black")
        }
            abline(v=x1s)
        text(dat$t.dat[match(levs,wr),"Time"]+.5,rep(c(-.05, abs(min(ylim)),abs(min(ylim))+.1),length=length(levs)),levs,cex=bcex,offset=0, pos=4)#,offset=-offs}

        blc<-dat$blc
        ## Tool for adding line and point plot for all lines
            #matlines(xseq, blc[,m.names], col=rgb(0,0,0,3, maxColorValue=100), lwd=.01)
            #matpoints(xseq, blc[,m.names], col=rgb(0,0,0,3, maxColorValue=100), pch=16, cex=.03)
        
        #cols <- rainbow(length(m.names),start=.55)
     
        library(cluster)
        blc<-dat$blc
        pam5 <- pam(t(blc[,m.names]),k=subset.n)
        s.names <- row.names(pam5$medoids)
        pam5.tab <- table(pam5$clustering)
        #tags <- paste(paste("#",names(pam5.tab),sep=""),as.vector(pam5.tab),sep=":")
        group.means<-list()
        group.names<-list()
        for(i in 1:subset.n){
            x.names<-names(which(pam5$clustering==i, arr.ind=T))
            group.names[[i]]<-x.names
            group.means[i]<-paste(
            round(mean(dat$c.dat[x.names, "area"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "area"]), digits=0)," : ",
            round(mean(dat$c.dat[x.names, "mean.gfp"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.gfp"]), digits=0)," : ",	
            round(mean(dat$c.dat[x.names, "mean.tritc"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.tritc"]), digits=0), sep="")
            # adding standard deviation,"\u00b1",round(sd(dat$c.dat[x.names, "area"]), digits=0)," : ", 
        }
        
        
        tags <- paste(as.vector(pam5.tab),":",group.means)
        info<-pam5$clustering
        
        ## Tool For adding color to selected Traces
        if(is.null(cols)){
            cols <-brewer.pal(8,"Dark2")
            cols <- rep(cols,ceiling(length(s.names)/length(cols)))
            cols <- cols[1:length(s.names)]
        }else{cols<-rep(cols, length(m.names))}

        ## Tool for adding labeling for single line within stacked traces
        par(xpd=T)
        for(i in 1:length(s.names)){
            if(length(group.names[[i]])>=2){
                matlines(xseq, blc[,group.names[[i]]]+i+sf, col=rgb(0,0,0,20, maxColorValue=100), lwd=.01)
                lines(xseq, apply(blc[,group.names[[i]]],1,mean)+i+sf, col=cols[i], lwd=.2)
                points(xseq, apply(blc[,group.names[[i]]],1,mean)+i+sf, col=cols[i], pch=16, cex=.02)
                text(x=min(blc[,1]), y=blc[nrow(t.dat),s.names[i]]+i+sf, labels=i, col=cols[i], pos=2, cex=bcex)
                text(x=max(blc[,1]), y=blc[nrow(t.dat),s.names[i]]+i+sf, labels=tags[i], col=cols[i], pos=4, cex=bcex)
            }else{
                lines(xseq, blc[,group.names[[i]]]+i+sf, col=cols[i], lwd=.2)
                points(xseq, blc[,group.names[[i]]]+i+sf, col=cols[i], pch=16, cex=.02)
                text(x=min(blc[,1]), y=blc[nrow(t.dat),s.names[i]]+i+sf, labels=i, col=cols[i], pos=2, cex=bcex)
                text(x=max(blc[,1]), y=blc[nrow(t.dat),s.names[i]]+i+sf, labels=tags[i], col=cols[i], pos=4, cex=bcex)
                }
        }
        par(xpd=F)
    

    # Tool for adding boxplot
        par(xpd=T)
        dev.current<-dev.cur()
        if(is.null(bp.param)){
            dat.select<-"c.dat"
            bp.param<-c(
            grep("area",names(dat$c.dat),value=T),
            #tryCatch(grep("mean.gfp",names(dat$c.dat)),error=function(e) NULL),
            grep("mean.gfp",names(dat$c.dat),value=T),
            grep("mean.tritc",names(dat$c.dat),value=T))
            
            cols<-c("blue", "darkgreen","red")
        #}else{
        #	dat.select<-select.list(names(dat))
        #	bp.param<-as.character(select.list(names(dat[[dat.select]]), multiple=T))
        #	cols<-NULL
        #	}
        }else{
        dat.select<-"c.dat"
        bp.param<-bp.param}
        
        
        for(i in 1:length(s.names)){
            xleft<-max(blc[,1])+xinch(.3)
            xright<-xleft+xinch(1)*length(bp.param)
            y<-(blc[nrow(t.dat),group.names[[i]]]+i+sf)
            ybottom<- y-yinch(.5)
            ytop<-y+yinch(.5)

            bp.img<-bpfunc.3(dat,group.names[[i]],dat.select, bp.param, print.out=T, cols=cols, bcex=bcex)
            dev.set(dev.current)
            rasterImage(bp.img,xleft, ybottom, xright, ytop)
        }
        
    continue<-select.list(c("yes", "no"))
    if(continue=="yes"){
        while(i!=length(s.names)+1){
            i<-scan(n=1)
            if(i>length(s.names)| i==0){i<-length(s.names)+1
            }else{LinesEvery.5(dat,sample(names(which(info==i, arr.ind=T)))[1:15], img=NULL, pic.plot=T, sf=.3, lmain=i,m.order="area")}
            #multi.pic.zoom(dat, names(which(info==i, arr.ind=T)), img, plot.new=T)
        }
    }	
    }
    else{LinesEvery.5(dat, m.names,img)}
    #return(pam5$clustering)	
    return(group.names)
        
} 

##170109
#intereact: LOGICAL; 
#TRUE select cell groups to work though and return list of groups of cells
#FALSE only plot out the groups, and dont return group of cells

##region.group: Select a region to group the cells around.  Brings up option to select region to group around
#170403 bp logical: lets you choose whether to boxplot

#170508 Allows to select the trace you would like to use for grouping with option:
#t.type:input character

#170605:  Adding a drop function to this.  It will automatically update the RD.file.  I need something to drop cells much faster
#
LinesStack.2<- function(dat,m.names=NULL,t.type=NULL,lmain="", interact=T, region.group=T,levs=NULL, plot.new=TRUE,bcex=.7, sf=1.1, subset.n=NULL, img=NULL,bp.param=NULL, bp=F, bp.pts=F){
    #graphics.off()
    if(is.null(img)){img<-dat$img1}
    if(is.null(m.names)){
        dropped.cells<-cellzand(dat$bin, "drop",1)
        m.names<-setdiff(dat$c.dat$id, dropped.cells)
    }else{
        dropped.cells<-cellzand(dat$bin, "drop",1)
        m.names<-setdiff(m.names, dropped.cells)
        }
        
    if(is.null(subset.n)){subset.n<-as.numeric(select.list(as.character(c(5,10,15,20,25,30))))}
    if(plot.new){
    
        if(subset.n>=10){
            dev.new(width=14,height=10)
            }
        else{dev.new(width=14,height=10)}
        
        linesstack.win<-dev.cur()
    }
    if(length(m.names)>subset.n){
        
        if(is.null(t.type)){t.dat<-dat$t.dat}
        else{t.dat<-dat[[t.type]]}
        
        wr<-dat$w.dat[,2]
        if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")}
        else{levs<-levs}
        m.names <- intersect(m.names,names(t.dat))
        hbc <-(max(t.dat[,m.names])+subset.n)*sf
        #hbc <- (subset.n*(.8*sf))+max(t.dat[,m.names])
        
        xseq <- t.dat[,1]
        library(RColorBrewer)
        par(mai=c(2,1,1,1))
        
        ylim<-c(-.1,hbc)
        #ylim <- c(-.1,2.5)
        
        plot(xseq,t.dat[,m.names[1]],ylim=ylim,xlab="",main=lmain,type="n",xlim=c(min(xseq)-1.5,max(xseq)+10))#-sf
        #axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
        
        ## Tool for adding window region labeling
        if(length(wr) > 0){
            #levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
            x1s <- tapply(xseq,as.factor(wr),min)[levs]
            x2s <- tapply(xseq,as.factor(wr),max)[levs]
            y1s <- rep(min(ylim)-.2,length(x1s))
            y2s <- rep(max(ylim)+.2,length(x1s))
            rect(x1s,y1s,x2s,y2s,col="grey90",border="black")
            levs.loc<-tapply(dat$w.dat[,"Time"],as.factor(wr),mean)[levs]
            par(xpd=T)
            text(levs.loc,par("usr")[3],levs,pos=3,offset=-4.3,cex=bcex, srt=90)	
            #text(levs.loc,par("usr")[3],levs,pos=3,offset=-4.2,cex=bcex, srt=90)	
            #text(t.dat[match(levs,wr),"Time"], par('usr')[4]-yinch(.5),levs, cex=bcex, offset=0, pos=4 )
            #text(t.dat[match(levs,wr),"Time"],rep(c(-.05, abs(min(ylim)),abs(min(ylim))+.1),length=length(levs)),levs,cex=bcex,offset=0, pos=4)#,offset=-offs}
        }
        ## Tool for adding line and point plot for all lines
            #matlines(xseq, blc[,m.names], col=rgb(0,0,0,3, maxColorValue=100), lwd=.01)
            #matpoints(xseq, blc[,m.names], col=rgb(0,0,0,3, maxColorValue=100), pch=16, cex=.03)
        
        #cols <- rainbow(length(m.names),start=.55)
     
        library(cluster)	
        ## To select data within the experiment to group around
        if(region.group){
            dev.new(width=15, height=10)
            LinesEvery.5(dat, sample(row.names(dat$c.dat)[1:5]), plot.new=F, lmain="Click to Select region to Groups Cells", t.type="t.dat", img="img1")
            #LinesEvery.4(dat, sample(row.names(dat$c.dat)[1:5]), plot.new=F, lmain="Click to Select region to Groups Cells", img=dat$img1)
            b.xseq<-locator(n=2, type="o", pch=15, col="red")$x
            dev.off()
            x.min<-which(abs(t.dat$Time-b.xseq[1])==min(abs(t.dat$Time-b.xseq[1])))
            x.max<-which(abs(t.dat$Time-b.xseq[2])==min(abs(t.dat$Time-b.xseq[2])))
            
            pam5 <- pam(t(t.dat[x.min:x.max,m.names]),k=subset.n)
            s.names <- row.names(pam5$medoids)
            pam5.tab <- table(pam5$clustering)
            #tags <- paste(paste("#",names(pam5.tab),sep=""),as.vector(pam5.tab),sep=":")
            group.means<-list()
            group.names<-list()
            for(i in 1:subset.n){
                x.names<-names(which(pam5$clustering==i, arr.ind=T))
                group.names[[i]]<-x.names
                group.means[i]<-paste(
                tryCatch(round(mean(dat$c.dat[x.names, "area"]), digits=0),error=function(e) NULL),
                "\u00b1",
                tryCatch(round(sd(dat$c.dat[x.names, "area"]), digits=1),error=function(e) NULL))#," : ",
                #tryCatch(round(mean(dat$c.dat[x.names, "mean.gfp"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.gfp"]), digits=0),error=function(e) NULL)," : ",	
                #tryCatch(round(mean(dat$c.dat[x.names, "mean.tritc"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.tritc"]), digits=0),error=function(e) NULL), sep="")
                #tryCatch(round(sd(dat$c.dat[x.names, "area"]), digits=0),"\u00b1",error=function(e) NULL)
            }
    
        }else{
            library(cluster)
            pam5 <- pam(t(t.dat[,m.names]),k=subset.n)
            s.names <- row.names(pam5$medoids)
            pam5.tab <- table(pam5$clustering)
            #tags <- paste(paste("#",names(pam5.tab),sep=""),as.vector(pam5.tab),sep=":")
            group.means<-list()
            group.names<-list()
            for(i in 1:subset.n){
                x.names<-names(which(pam5$clustering==i, arr.ind=T))
                group.names[[i]]<-x.names
                group.means[i]<-paste(
                tryCatch(round(mean(dat$c.dat[x.names, "area"]), digits=0),error=function(e) NULL),
                "\u00b1",
                tryCatch(round(sd(dat$c.dat[x.names, "area"]), digits=1),error=function(e) NULL))				
                #round(mean(dat$c.dat[x.names, "mean.gfp"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.gfp"]), digits=0)," : ",	
                #round(mean(dat$c.dat[x.names, "mean.tritc"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.tritc"]), digits=0)
                #adding standard deviation
                #"\u00b1",round(sd(dat$c.dat[x.names, "area"]), digits=0), sep="")
            }
        }			
        tags <- paste(as.vector(pam5.tab),":",group.means)
        info<-pam5$clustering
        
        ## Tool For adding color to selected Traces
        cols <-brewer.pal(8,"Dark2")
        cols <- rep(cols,ceiling(length(s.names)/length(cols)))
        cols <- cols[1:length(s.names)]

        ## Tool for adding labeling for single line within stacked traces
        par(xpd=T)
        dev.set(which=linesstack.win)
        for(i in 1:length(s.names)){
            if(length(group.names[[i]])>=2){
                matlines(xseq, (t.dat[,group.names[[i]]]+i)*sf, col=rgb(0,0,0,20, maxColorValue=100), lwd=.01)
                lines(xseq, apply(t.dat[,group.names[[i]]],1,mean)+i*sf, col=cols[i], lwd=.2)
                points(xseq, apply(t.dat[,group.names[[i]]],1,mean)+i*sf, col=cols[i], pch=16, cex=.02)
                text(x=min(t.dat[,1]), y=t.dat[1,s.names[i]]+i*sf, labels=i, col=cols[i], pos=2, cex=bcex)
                text(x=max(t.dat[,1]), y=t.dat[nrow(dat$t.dat),s.names[i]]+i*sf, labels=tags[i], col=cols[i], pos=4, cex=bcex)
            }else{
                lines(xseq, t.dat[,group.names[[i]]]+i*sf, col=cols[i], lwd=.2)
                points(xseq, t.dat[,group.names[[i]]]+i*sf, col=cols[i], pch=16, cex=.02)
                text(x=min(t.dat[,1]), y=t.dat[1,s.names[i]]+i*sf, labels=i, col=cols[i], pos=2, cex=bcex)
                text(x=max(t.dat[,1]), y=t.dat[nrow(dat$t.dat),s.names[i]]+i*sf, labels=tags[i], col=cols[i], pos=4, cex=bcex)
                }
        }
        if(region.group){
            points(b.xseq,rep(min(ylim),2),pch=15, col="blue", cex=.5)
            abline(v=b.xseq, col="blue")
        }else{}

        par(xpd=F)

    #### Tool for adding boxplot
        if(bp){
                par(xpd=T)
                dev.current<-dev.cur()
                if(is.null(bp.param)){
                    #dat.select<-"c.dat"
                    #bp.param<-c(
                    #grep("area",names(dat$c.dat),value=T),
                    ##tryCatch(grep("mean.gfp",names(dat$c.dat)),error=function(e) NULL),
                    #grep("mean.gfp",names(dat$c.dat),value=T),
                    #grep("mean.tritc",names(dat$c.dat),value=T))
                    
                    #cols<-c("blue", "darkgreen","red")
                #}else{
                    dat.select<-select.list(names(dat))
                    bp.param<-as.character(select.list(names(dat[[dat.select]]), multiple=T))
                    cols<-NULL
                }else{
                dat.select<-"c.dat"
                bp.param<-bp.param
                }
                
                #for(i in 1:length(group.names)){
                    #if(length(group.names[[i]])>5){
                    #	xleft<-max(t.dat[,1])+xinch(.3)
                    #	xright<-xleft+xinch(1)*length(bp.param)
                    #	y<-(apply(t.dat[nrow(t.dat),group.names[[i]]],1,mean)+i)*sf
                    #	ybottom<- y-yinch(.5)
                    #	ytop<-y+yinch(.5)

                    #	bp.img<-bpfunc.3(dat,group.names[[i]],dat.select, bp.param, print.out=T, cols=cols, bcex=bcex)
                    #	dev.set(dev.current)
                    #	rasterImage(bp.img,xleft, ybottom, xright, ytop)
                    #}else{}
                    
                    #170509 How to create a new window with these boxplots
                dev.new(width=length(bp.param), height=subset.n)
                bp.win<-dev.cur()
                par(mfrow=c(subset.n,1))
                group.names.rev<-rev(group.names)
                    
                for(i in 1:length(group.names.rev)){
                        par(mar=c(0,0,0,0))
                        plot(0,0)
                        dim<-par("usr")
                        xleft<-par("usr")[1]
                        xright<-par("usr")[2]
                        ybottom<- par("usr")[3]
                        ytop<-par("usr")[4]
                        bp.img<-bpfunc.3(dat,group.names.rev[[i]],dat.select, bp.param, print.out=T, cols=cols, bcex=bcex, bp.pts=bp.pts)
                        dev.set(bp.win)
                        rasterImage(bp.img,xleft, ybottom, xright, ytop)
                        text(xleft+xinch(.1), 0, subset.n-i+1, cex=2)
                    }
            
            }
        }
        
            if(interact){
                continue<-select.list(c("yes", "no"))
                if(continue=="yes"){
                    while(i!=length(s.names)+1){
                        i<-scan(n=1)
                        if(i>length(s.names)| i==0){i<-length(s.names)+1}
                        else{
                            assesment.selection<-select.list(c("Trace.Click","LinesEvery","LinesStack", "drop"))
                            
                            if(assesment.selection=="Trace.Click"){
                                Trace.Click.dev(dat,names(which(info==i, arr.ind=T)))
                            }
                            
                            if(assesment.selection=="LinesEvery"){
                                number.to.display<-as.numeric(select.list(as.character(c(5,10,20))))
                                LinesEvery.5(dat,sample(names(which(info==i, arr.ind=T)))[1:number.to.display], img, pic.plot=T, lmain=i,m.order="area", plot.new=T, col="black")
                            }
                            
                            if(assesment.selection=="LinesStack"){
                            
                                LinesStack.2(dat,names(which(info==i, arr.ind=T)),bp=F,lmain=i, interact=T, region.group=T,levs=NULL, plot.new=TRUE,bcex=.7, img=dat$img1, t.type="mp.1")
                            }
                            if(assesment.selection=="drop"){
                                rd.namels2 <- as.character(substitute(dat))
                                dat$bin[names(which(info==i, arr.ind=T)), "drop"]<-1
                                assign(rd.namels2, dat, envir=.GlobalEnv)
                                print(paste("You Dropped Group",i))
                            }
                        }
                    }
                }
                #return(pam5$clustering)	

            }
        
    #dev.off(which=linesstack.win)
    return(group.names)
}

# Stacked Traces, 
# Input is a list of cells
# Currently Created for the 5 cell classes of,
# +ib4+cgrp, +IB4, +CGRP, -/-, glia
LinesStack.3 <- function(dat,cells=NULL,lmain="",levs=NULL, plot.new=TRUE,bcex=.7, sf=.9, img=NULL, sample.num=NULL){
    graphics.off()
    if(is.null(img)){img<-dat$img1}
    if(is.null(sample.num)){sample.num<-10}
    
    if(is.null(cells)){
        cells<-dat$cells
        cells<-cells[c('000','00','01','10','11')]
    }else{
        cells.main<-dat$cells
        cells.main<-cells.main[c('000','00','01','10','11')]
        bob<-list()
        
        for(i in 1:length(cells.main)){
            x.names<-intersect(cells,cells.main[[i]])
            bob[[i]]<-x.names
        }

        cells<-bob
        names(cells)<-c('000','00','01','10','11')
    }
    #cells<-cells[c('000','00','01','10','11')]
    
    if(plot.new){dev.new(width=10,height=6)}
        t.dat<-dat$t.dat
        wr<-dat$w.dat[,2]
        if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")}
        else{levs<-levs}
        #m.names <- intersect(m.names,names(t.dat))
        xseq <- t.dat[,1]
        library(RColorBrewer)
        par(mar=c(4,2,4,4), bty="L")
        #hbc <- (5*(.8*sf))+max(t.dat[,Reduce(c,stack(cells)[1])])
        #hbc <- 5*sf+max(t.dat[,Reduce(c,stack(cells)[1])])

        ylim <- c(.5,5.2)
        #ylim<-c(-.1,hbc)
        plot(xseq,t.dat[,cells[[1]][1]],ylim=ylim,xlab="Time (min)",main=lmain,type="n", xaxt="n",xlim=c(min(xseq), max(xseq)*1.5))#-sf
        axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
        
        ## Tool for adding window region labeling
        if(length(wr) > 0){
            #levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
            x1s <- tapply(xseq,as.factor(wr),min)[levs]
            x2s <- tapply(xseq,as.factor(wr),max)[levs]
            y1s <- rep(min(ylim),length(x1s))*1.03-rep(min(ylim),length(x1s))
            y2s <- rep(max(ylim),length(x1s))*1.03
            rect(x1s,y1s,x2s,y2s,col="grey90",border="black")
            text(dat$w.dat[match(levs,wr),"Time"],rep(c(.5,.6,.7),length=length(levs)),levs,cex=.6,offset=0, pos=4)#,offset=-offs}
        }
        blc<-dat$blc	 
     
        ##  Tool for creating mean and st.dev calculation
        library(cluster)
        blc<-dat$blc
        group.means<-list()
        group.names<-list()
        
        for(i in 1:length(cells)){
            if(length(cells[[i]])>1){
                x.names<-cells[[i]]
                group.names[[i]]<-names(cells[i])
                group.means[i]<-paste(
                length(cells[[i]]),":",
                round(mean(dat$c.dat[x.names, "area"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "area"]), digits=0),"   :   ",
                round(mean(dat$c.dat[x.names, "mean.gfp"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.gfp"]), digits=0),"   :   ",	
                round(mean(dat$c.dat[x.names, "mean.tritc"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.tritc"]), digits=0), sep="")
                #adding standard deviation,"\u00b1",round(sd(dat$c.dat[x.names, "area"]), digits=0)," : ", 
            }
            else{
                x.names<-cells[[i]]
                group.names[[i]]<-names(cells[i])
                group.means[i]<-paste(
                length(cells[[i]]))
            }}
        
        
        
        ## Tool For adding color to selected Traces
        cols <-brewer.pal(8,"Dark2")
        cols <- rep(cols,5)
        cols <- cols[1:5]
        
        cols<-c("mediumpurple1","goldenrod1", "firebrick1", "limegreen", "steelblue3")

    
        ## Tool for adding labeling for single line within stacked traces
        for(i in 1:length(cells)){
            if(length(cells[[i]])>1){
                matlines(xseq, blc[,cells[[i]]]+i*sf, col=rgb(0,0,0,20, maxColorValue=100), lwd=.3)
                lines(xseq, apply(blc[,cells[[i]]],1,mean)+i*sf, col=cols[i], lwd=1.2)
                text(x=min(blc[,1]), y=blc[nrow(t.dat),cells[[i]]]+i*sf, labels=group.names[i], col=cols[i], pos=2, cex=bcex)
                text(x=max(blc[,1]), y=blc[nrow(t.dat),cells[[i]]]+i*sf, labels=group.means[i], col="black", pos=4, cex=bcex)
                }
            else{
                lines(xseq, blc[,cells[[i]]]+i*sf, col=rgb(0,0,0,20, maxColorValue=100), lwd=.3)
                text(x=min(blc[,1]), y=blc[nrow(t.dat),cells[[i]]]+i*sf, labels=group.names[i], col=cols[i], pos=2, cex=bcex)
            }}

                
        ## Tool for adding boxplot to plot
        dev.current<-dev.cur()

        for(i in 1:length(cells)){
            xleft<-max(blc[,1])*1.05
            xright<-xleft+xinch(2.74)
            y<-(blc[nrow(t.dat),cells[[i]]]+i*sf)
            ybottom<- y-.55
            ytop<-ybottom+yinch(.85)

            #dev.set(dev.list()[length(dev.list())])
            bp.img<-bpfunc.2(dat,cells[[i]])
            dev.set(dev.current)
            rasterImage(bp.img,xleft, ybottom, xright, ytop)
        }

    continue<-select.list(c("yes", "no"))
    if(continue=="yes"){
        i<-1
        while(i!=00){
            i<-scan(n=1)
            cells.tp<-cells[[i]]
            LinesEvery.4(dat,sample(cells.tp)[1:15], img, pic.plot=T, sf=.3, lmain=i,m.order="area")
            #multi.pic.zoom(dat, names(which(info==i, arr.ind=T)), img, plot.new=T)
        }


    #for(i in 1:length(cells)){
    #		if(length(cells[[i]])<20){
        #		LinesEvery.4(dat,cells[[i]], img, pic.plot=T, lmain=names(cells[i]), m.order="area", levs=levs, sf=.6)
    #		}
    #		else{
    #		# select the range of
    #			sample.num<-ceiling(sample.num/2)
    #			cells.n<-sort(c(ceiling(seq(1,length(cells[[i]]), length.out=5)),ceiling(seq(1,length(cells[[i]]), length.out=5))+1))
    #			cells.rs<-c.sort(dat$c.dat[cells[[i]],], "area")
    #			LinesEvery.4(dat, cells.rs[cells.n],img, lmain=names(cells[i]), m.order="area",levs=levs, sf=.4)}
    #		#multi.pic.zoom(dat, names(which(info==i, arr.ind=T)), img, plot.new=T)
        #}
    }

    #else{LinesEvery.4(dat, m.names,img)}
    #return(pam5$clustering)	
    #return(group.names)
    #print(group.means)
}
 
LinesStack.4 <- function(dat,cells=NULL,lmain="",levs=NULL, plot.new=TRUE,bcex=.7, sf=.9, img=NULL, sample.num=NULL){
    if(is.null(img)){img<-dat$img.gtd}
    if(is.null(sample.num)){sample.num<-10}
    if(is.null(cells)){cells<-dat$cells}
    else{
        cells.main<-dat$cells
        cells.main<-cells.main[c('000','00','01','10','11')]
        bob<-list()
        for(i in 1:length(cells.main)){
            x.names<-intersect(cells,cells.main[[i]])
            bob[[i]]<-x.names
            }

        cells<-bob
        names(cells)<-c('000','00','01','10','11')
    }
    #cells<-cells[c('000','00','01','10','11')]
    
    if(plot.new){dev.new(width=10,height=6)}
        t.dat<-dat$t.dat
        wr<-dat$w.dat[,2]
        if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")}
        else{levs<-levs}
        #m.names <- intersect(m.names,names(t.dat))
        xseq <- t.dat[,1]
        library(RColorBrewer)
        par(mar=c(4,2,4,4), bty="L")
        #hbc <- (5*(.8*sf))+max(t.dat[,Reduce(c,stack(cells)[1])])
        #hbc <- 5*sf+max(t.dat[,Reduce(c,stack(cells)[1])])

        ylim <- c(.5,5.2)
        #ylim<-c(-.1,hbc)
        plot(xseq,t.dat[,cells[[1]][1]],ylim=ylim,xlab="Time (min)",main=lmain,type="n", xaxt="n",xlim=c(min(xseq), max(xseq)*1.5))#-sf
        axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
        
        ## Tool for adding window region labeling
        if(length(wr) > 0){
            #levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
            x1s <- tapply(xseq,as.factor(wr),min)[levs]
            x2s <- tapply(xseq,as.factor(wr),max)[levs]
            y1s <- rep(min(ylim),length(x1s))*1.03-rep(min(ylim),length(x1s))
            y2s <- rep(max(ylim),length(x1s))*1.03
            rect(x1s,y1s,x2s,y2s,col="grey90",border="black")
            text(dat$t.dat[match(levs,wr),"Time"],rep(c(.5,.6,.7),length=length(levs)),levs,cex=.6,offset=0, pos=4)#,offset=-offs}
        }
        blc<-dat$blc	 
     
        ##  Tool for creating mean and st.dev calculation
        library(cluster)
        blc<-dat$blc
        group.means<-list()
        group.names<-list()
        
        for(i in 1:length(cells)){
            if(length(cells[[i]])>1){
                x.names<-cells[[i]]
                group.names[[i]]<-names(cells[i])
                group.means[i]<-paste(
                length(cells[[i]]),":",
                round(mean(dat$c.dat[x.names, "area"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "area"]), digits=0),"   :   ",
                round(mean(dat$c.dat[x.names, "mean.gfp"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.gfp"]), digits=0),"   :   ",	
                round(mean(dat$c.dat[x.names, "mean.tritc"]), digits=0),"\u00b1",round(sd(dat$c.dat[x.names, "mean.tritc"]), digits=0), sep="")
                #adding standard deviation,"\u00b1",round(sd(dat$c.dat[x.names, "area"]), digits=0)," : ", 
            }
            else{
                x.names<-cells[[i]]
                group.names[[i]]<-names(cells[i])
                group.means[i]<-paste(
                length(cells[[i]]))
            }}
        
        
        
        ## Tool For adding color to selected Traces
        cols <-brewer.pal(8,"Dark2")
        cols <- rep(cols,5)
        cols <- cols[1:5]
        
        cols<-c("mediumpurple1","goldenrod1", "firebrick1", "limegreen", "steelblue3")

    
        ## Tool for adding labeling for single line within stacked traces
        for(i in 1:length(cells)){
            if(length(cells[[i]])>1){
                matlines(xseq, blc[,cells[[i]]]+i*sf, col=rgb(0,0,0,10, maxColorValue=100), lwd=.3)
                lines(xseq, apply(blc[,cells[[i]]],1,mean)+i*sf, col=cols[i], lwd=1.2)
                text(x=min(blc[,1]), y=blc[nrow(t.dat),cells[[i]]]+i*sf, labels=group.names[i], col=cols[i], pos=2, cex=bcex)
                text(x=max(blc[,1]), y=blc[nrow(t.dat),cells[[i]]]+i*sf, labels=group.means[i], col="black", pos=4, cex=bcex)
                }
            else{
                lines(xseq, blc[,cells[[i]]]+i*sf, col=rgb(0,0,0,80, maxColorValue=100), lwd=.3)
                text(x=min(blc[,1]), y=blc[nrow(t.dat),cells[[i]]]+i*sf, labels=group.names[i], col=cols[i], pos=2, cex=bcex)
            }}

                
        ## Tool for adding boxplot to plot
        for(i in 1:length(cells)){
            xleft<-max(blc[,1])*1.05
            xright<-xleft+xinch(2.74)
            y<-(blc[nrow(t.dat),cells[[i]]]+i*sf)
            ybottom<- y-.55
            ytop<-ybottom+yinch(.85)

            #dev.set(dev.list()[length(dev.list())])
            rasterImage(bpfunc.2(dat,cells[[i]]),xleft, ybottom, xright, ytop)
        }
    
    
    continue<-select.list(c("yes", "no"))
    if(continue=="yes"){
    for(i in 1:length(cells)){
        
            if(length(cells[[i]])<20){
                LinesEvery.4(dat,cells[[i]], img, pic.plot=T, lmain=names(cells[i]), m.order="area", levs=levs, sf=.6)
            }
            else{
            # select the range of
                sample.num<-ceiling(sample.num/2)
                cells.n<-sort(c(ceiling(seq(1,length(cells[[i]]), length.out=5)),ceiling(seq(1,length(cells[[i]]), length.out=5))+1))
                cells.rs<-c.sort(dat$c.dat[cells[[i]],], "area")
                LinesEvery.4(dat, cells.rs[cells.n],img, lmain=names(cells[i]), m.order="area",levs=levs, sf=.4)}
            #multi.pic.zoom(dat, names(which(info==i, arr.ind=T)), img, plot.new=T)
        }
    }

    #else{LinesEvery.4(dat, m.names,img)}
    #return(pam5$clustering)	
    return(group.names)
    print(group.means)
    
    
    
} 

bpfunc<-function(dat,n.names){
        if(length(n.names)>4){
        #par(width=12, height=4.5)
        par(mfrow=c(2,3))
        par(mar=c(2.5,2.5,2.5,2.5))
        par(cex=.8)
        dat.names<-names(dat$c.dat)
        #lab.1<-grep("gfp.1",dat.names,ignore.case=T, value=T)
        #lab.2<-grep("gfp.2",dat.names, ignore.case=T, value=T)
        #lab.3<-grep("tritc",dat.names, ignore.case=T, value=T)
        #lab.4<-grep("area",dat.names, ignore.case=T, value=T)
        
        #if(dat$c.dat[n.names, lab.1]!="N/A"){lab.1<-lab.1}
        #else{rm(lab.1)}
        #if(dat$c.dat[n.names, lab.2]!="N/A"){}
        #if(dat$c.dat[n.names, lab.3]!="N/A"){ }
        #if(dat$c.dat[n.names, lab.4]!="N/A"){}
    
        
        ##Color intensity 1
        boxplot(dat$c.dat[n.names,"mean.gfp"],main="GFP",bty="n",ylim=c(0,max(dat$c.dat["mean.gfp"])), col="springgreen4", outline=F)
        text(x=jitter(rep(1, length(dat$c.dat[n.names,"mean.gfp"])), factor=10),
        y=dat$c.dat[n.names,"mean.gfp"], 
        labels=as.character(dat$c.dat[n.names,"id"]))
        
        #Color Intensity 2
        boxplot(dat$c.dat[n.names,"mean.tritc"],main="IB4",ylim=c(0,max(dat$c.dat["mean.tritc"])), col="firebrick4", outline=F)
        text(x=jitter(rep(1, length(dat$c.dat[n.names,"mean.tritc"])), factor=10),
        y=dat$c.dat[n.names,"mean.tritc"], 
        labels=as.character(dat$c.dat[n.names,"id"]))
    
        # area
        boxplot(dat$c.dat[n.names,"area"],main="Area",ylim=c(0,max(dat$c.dat["area"])), col="lightslateblue", outline=F)
        text(x=jitter(rep(1, length(dat$c.dat[n.names,"area"])), factor=10),
        y=dat$c.dat[n.names,"area"], 
        labels=as.character(dat$c.dat[n.names,"id"]))
        
        ##Color intensity 1 log
        boxplot(1+dat$c.dat[n.names,"mean.gfp"],main="GFP",bty="n", col="springgreen4", outline=T, log="y")
        text(x=jitter(rep(1, length(dat$c.dat[n.names,"mean.gfp"])), factor=10),
        y=1+dat$c.dat[n.names,"mean.gfp"], 
        labels=as.character(dat$c.dat[n.names,"id"]))
        
        #Color Intensity 2 log
        boxplot(1+dat$c.dat[n.names,"mean.tritc"],main="IB4", col="firebrick4", outline=T, log="y")
        text(x=jitter(rep(1, length(dat$c.dat[n.names,"mean.tritc"])), factor=10),
        y=1+dat$c.dat[n.names,"mean.tritc"], 
        labels=as.character(dat$c.dat[n.names,"id"]))
        
        # area log
         boxplot(1+dat$c.dat[n.names,"area"],main="Area", col="lightslateblue", outline=T, log="y")
        text(x=jitter(rep(1, length(dat$c.dat[n.names,"area"])), factor=10),
        y=1+dat$c.dat[n.names,"area"], 
        labels=as.character(dat$c.dat[n.names,"id"]))
        
        dev.set(dev.list()[1])}
        else{
        par(mfrow=c(1,3))
        par(mar=c(2,2,2,2))
        
        stripchart(dat$c.dat[n.names,"mean.gfp"],main="GFP",ylim=c(0,max(dat$c.dat["mean.gfp"])),cex=2, col=c("green4"), outline=T, vertical=T, pch=".")
        text(x=1,
        y=dat$c.dat[n.names,"mean.gfp"], 
        labels=as.character(dat$c.dat[n.names,"id"]), col="green4")
        
        stripchart(dat$c.dat[n.names,"mean.tritc"],main="IB4",ylim=c(0,max(dat$c.dat["mean.tritc"])), ,cex=2,col="red", outline=F, vertical=T, pch=".")
        text(x=1,
        y=dat$c.dat[n.names,"mean.tritc"], 
        labels=as.character(dat$c.dat[n.names,"id"]), col="red")
        
        stripchart(dat$c.dat[n.names,"area"],main="Area",ylim=c(0,max(dat$c.dat["area"])), ,cex=2,col="lightslateblue", outline=F, vertical=T, pch=".")
        text(x=1,
        y=dat$c.dat[n.names,"area"], 
        labels=as.character(dat$c.dat[n.names,"id"]), col="lightslateblue")
        
        dev.set(dev.list()[5])}
}

bpfunc.2<-function(dat,n.names, bp.pts=T){
    require(png)
    
    png('tmp.png', width=2.74, height=.85, units="in", res=200)
    #dev.new(width=2.74, height=1)
    if(length(n.names)>4){

        
        par(mfrow=c(1,3),mar=c(1,3,2,0), bty="n",lwd=1, lty=1, cex.axis=.8, cex=.6)
        dat.names<-names(dat$c.dat)
        #lab.1<-grep("gfp.1",dat.names,ignore.case=T, value=T)
        #lab.2<-grep("gfp.2",dat.names, ignore.case=T, value=T)
        #lab.3<-grep("tritc",dat.names, ignore.case=T, value=T)
        #lab.4<-grep("area",dat.names, ignore.case=T, value=T)
        
        #if(dat$c.dat[n.names, lab.1]!="N/A"){lab.1<-lab.1}
        #else{rm(lab.1)}
        #if(dat$c.dat[n.names, lab.2]!="N/A"){}
        #if(dat$c.dat[n.names, lab.3]!="N/A"){ }
        #if(dat$c.dat[n.names, lab.4]!="N/A"){}
    
        
        ##Color intensity 1
        boxplot(dat$c.dat[n.names,"mean.gfp"],main="GFP",
        ylim=c(min(dat$c.dat["mean.gfp"]),max(dat$c.dat["mean.gfp"])), col="springgreen4", outline=F,yaxt="n", boxwex=.8, medlwd=.4,whisklty=1)
        if(bp.pts==T){stripchart(dat$c.dat[n.names,"mean.gfp"], add=T, method="jitter", vertical=T, jitter=.2, pch=18, cex=.7)}
        mtext(paste(round(mean(dat$c.dat[n.names, "mean.gfp"]), digits=3),"\u00b1",round(sd(dat$c.dat[n.names, "mean.gfp"]), digits=3)),1, cex=.5)
        #text(x=jitter(rep(1, length(dat$c.dat[n.names,"mean.gfp"])), factor=10),
        #y=dat$c.dat[n.names,"mean.gfp"], 
        #labels=as.character(dat$c.dat[n.names,"id"]), cex=.4)
        axis(2, at=c(round(min(dat$c.dat["mean.gfp"]), digits=3),round(max(dat$c.dat["mean.gfp"]), digits=3)))#,labels=x, col.axis="red", las=2)
        box("figure")
        
        #Color Intensity 2
        boxplot(dat$c.dat[n.names,"mean.tritc"],main="IB4",
        ylim=c(min(dat$c.dat["mean.tritc"]),max(dat$c.dat["mean.tritc"])), col="red", outline=F, boxwex=.8, yaxt="n", medlwd=.4,whisklty=1)
        if(bp.pts==T){stripchart(dat$c.dat[n.names,"mean.tritc"], add=T, method="jitter", vertical=T, jitter=.2, pch=18, cex=.7)}
        #text(x=jitter(rep(1, length(dat$c.dat[n.names,"mean.tritc"])), factor=10),
        #y=dat$c.dat[n.names,"mean.tritc"], 
        #labels=as.character(dat$c.dat[n.names,"id"]), cex=.4)
        mtext(paste(round(mean(dat$c.dat[n.names, "mean.tritc"]), digits=3),"\u00b1",round(sd(dat$c.dat[n.names, "mean.tritc"]), digits=3)),1, cex=.5)
        axis(2, at=c(round(min(dat$c.dat["mean.tritc"]), digits=3),round(max(dat$c.dat["mean.tritc"]), digits=3)))#,labels=x, col.axis="red", las=2)
        box("figure")
        
        # area
        boxplot(dat$c.dat[n.names,"area"],main="Area",
        ylim=c(min(dat$c.dat["area"]),max(dat$c.dat["area"])), col="lightslateblue", outline=F, boxwex=.8, yaxt="n", medlwd=.4,whisklty=1)
        if(bp.pts==T){stripchart(dat$c.dat[n.names,"area"], add=T, method="jitter", vertical=T, jitter=.2, pch=18, cex=.7)}
        #text(x=jitter(rep(1, length(dat$c.dat[n.names,"mean.tritc"])), factor=10),
        #y=dat$c.dat[n.names,"mean.tritc"], 
        #labels=as.character(dat$c.dat[n.names,"id"]), cex=.4)
        mtext(paste(round(mean(dat$c.dat[n.names, "area"]), digits=0),"\u00b1",round(sd(dat$c.dat[n.names, "mean.tritc"]), digits=0)),1, cex=.5)
        axis(2, at=c(round(min(dat$c.dat["area"]), digits=0),round(max(dat$c.dat["area"]), digits=0)))#,labels=x, col.axis="red", las=2)
        box("figure")

    

    }

    else{
        par(mfrow=c(1,3))
        par(mar=c(2,2,2,2))
        
        stripchart(dat$c.dat[n.names,"mean.gfp"],main="GFP",ylim=c(0,max(dat$c.dat["mean.gfp"])),cex=1, col=c("green4"), outline=T, vertical=T, pch=".")
        text(x=1,
        y=dat$c.dat[n.names,"mean.gfp"], 
        labels=as.character(dat$c.dat[n.names,"id"]), col="green4", cex=.8)
        box("figure")
        
        stripchart(dat$c.dat[n.names,"mean.tritc"],main="IB4",ylim=c(0,max(dat$c.dat["mean.tritc"])),cex=1,col="red", outline=T, vertical=T, pch=".")
        text(x=1,
        y=dat$c.dat[n.names,"mean.tritc"], 
        labels=as.character(dat$c.dat[n.names,"id"]), col="red", cex=.8)
        box("figure")
        
        stripchart(dat$c.dat[n.names,"area"],main="Area",ylim=c(0,max(dat$c.dat["area"])),cex=1,col="lightslateblue", outline=T, vertical=T, pch=".")
        text(x=1,
        y=dat$c.dat[n.names,"area"], 
        labels=as.character(dat$c.dat[n.names,"id"]), col="lightslateblue", cex=.8)
        box("figure")
    }
    
    
    dev.off()
    tmp.png <- readPNG("tmp.png")
    dim(tmp.png)
    unlink("tmp.png")
    return(tmp.png)			

}

bpfunc.3<-function(dat,n.names=NULL, dat.select=NULL, parameters=NULL,bp.pts=F, print.out=F, plot.new=T,bcex=NULL, ylim=NULL, cols=NULL){
    if(class(dat)=="list"){
        if(is.null(dat.select)){dat.select<-select.list(names(dat))}else{dat.select<-dat.select}
        dat<-dat[[dat.select]]
    }else{dat<-dat}
    
    require(png)
    #dev.new(width=2.74, height=1)
    if(is.null(n.names)){n.names<-dat$c.dat$id}else{n.names<-n.names}
    
    if(is.null(parameters)){parameters<-select.list(names(dat), multiple=T)
    }else{parameters<-parameters}
    
    if(length(parameters)>6){width=ceiling(sqrt(length(parameters)));height=ceiling(sqrt(length(parameters)))
    }else{width=length(parameters);height=1}
    
    if(print.out){png('tmp.png', width=width*1.5, height=height*1.5, units="in", res=200,type="cairo")
    }else{if(plot.new){dev.new(width=width*1.5, height=height*1.5)}}
    
    if(is.null(bcex)){bcex<-.8}else{bcex<-bcex}
    if(is.null(cols)){cols<-"blue";cols<-(rep(cols, length(parameters)))}else{cols<-cols}
    

    
    par(mfrow=c(height,width),mar=c(1,3,3,0), bty="n",lwd=1, lty=1, cex.axis=.8, cex=.6)


    #loop through selected parameters, and applies it to selected dataframe for dat list	
    if(length(n.names)>4){
        for(i in 1:length(parameters)){
                if(is.null(ylim)){ylim.i<-c(min(dat[,parameters[i]]),max(dat[,parameters[i]]))
                }else{ylim.i<-ylim}
        
                main.name<-strsplit(parameters[i], "[.]")[[1]]
                boxplot(dat[n.names,parameters[i]],main=paste(main.name, collapse=" "),
                ylim=ylim.i, col=cols[i], outline=F,yaxt="n", boxwex=.8, medlwd=.4,whisklty=1, cex=bcex)
                
                if(bp.pts==T){stripchart(dat[n.names,parameters[i]], add=T, method="jitter", vertical=T, jitter=.2, pch=18, cex=.5)
                }else{
                text(x=jitter(rep(1, length(dat[n.names,parameters[i]])), factor=10),
                    y=dat[n.names,parameters[i]], 
                    labels=as.character(row.names(dat[n.names,])), cex=bcex,
                    col=rgb(0,0,0,4,maxColorValue=10))
                    }
                mtext(paste(round(mean(dat[n.names, parameters[i]]), digits=3),"\u00b1",round(sd(dat[n.names, parameters[i]]), digits=3)),1, cex=bcex)
                    axis(2, at=c(min(ylim.i),max(ylim.i)), cex=bcex)#,labels=x, col.axis="red", las=2)
                    box("figure")
        }
    }else{
        for(i in 1:length(parameters)){
            if(is.null(ylim)){ylim.i<-c(min(dat[,parameters[i]]),max(dat[,parameters[i]]))
            }else{ylim.i<-ylim}
                    main.name<-strsplit(parameters[i], "[.]")[[1]]
                    stripchart(dat[n.names,parameters[i]],main=paste(main.name, collapse=" "),ylim=ylim.i,cex=1, col=c("green4"), outline=T, vertical=T, pch=".")
                    text(x=1,
                    y=dat[n.names,parameters[i]], 
                    labels=as.character(dat[n.names,"id"]), col=cols[i], cex=2)
                    box("figure")
        }
    }
    if(print.out){
        dev.off()
        tmp.png <- readPNG("tmp.png")
        dim(tmp.png)
        unlink("tmp.png")
        return(tmp.png)		
    }	
}

# New Boxplot Function
# dat : Experiment list RD.
# l.cells: Cells in a list format
# dat.name: Dataframe to pull data from
# col.name: collumn name  to get data for the boxplot
# jitter.f: Factor of jitter to  accomplish
# pts: points to add to boxplot
# notchs: logical (T/F) stand for notch selection 
#c("area","mean.gfp.start","mean.cy5.start")

boxplotlist<-function(dat,l.cells=NULL,dat.name="c.dat",col.name=NULL,jitter.f=.5,pts=T, notchs=F, bplog="y", sort=T){

    #back up operation to fill with cells for the boxplotting
    if(is.null(l.cells)){
        l.cells<-dat$cell.types
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
    
    windows(width=12, height=6,xpos=0, ypos=10)
    bp.win<-dev.cur()
    
    windows(width=14,height=4,xpos=0, ypos=540)
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

LinesStack.select <- function(dat,m.names,lmain="",levs=NULL, plot.new=TRUE,bcex=.8, sf=.2, subset.n=5){
    t.dat<-dat$t.dat
    wr<-dat$w.dat[,2]
    if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")}
    else{levs<-levs}
    m.names <- intersect(m.names,names(t.dat))
    hbc <- subset.n*sf+max(t.dat[,m.names])
    xseq <- t.dat[,1]
    if(plot.new){dev.new(width=10,height=6)}
    library(RColorBrewer)
     par(mar=c(4,2,4,4))
    #ylim <- c(-.1,1.4)
    ylim<-c(-.1,hbc)
    plot(xseq,t.dat[,m.names[1]],ylim=ylim,xlab="Time (min)",main=lmain,type="n", xaxt="n",xlim=c(min(xseq)-1.5,max(xseq)+1.5))#-sf
    axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
    ## Tool for adding window region labeling
    if(length(wr) > 0){
        #levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
        x1s <- tapply(xseq,as.factor(wr),min)[levs]
        x2s <- tapply(xseq,as.factor(wr),max)[levs]
        y1s <- rep(min(ylim)-.2,length(x1s))
        y2s <- rep(max(ylim)+.2,length(x1s))
        rect(x1s,y1s,x2s,y2s,col="grey90",border="black")
        text(dat$t.dat[match(levs,wr),"Time"],rep(c(abs(min(ylim)), abs(min(ylim*1.5))),length=length(levs)),levs,cex=bcex,offset=0, pos=4)#,offset=-offs}
    }
    blc<-dat$blc
    ## Tool for adding line and point plot for all lines
        #matlines(xseq, blc[,m.names], col=rgb(0,0,0,3, maxColorValue=100), lwd=.01)
        #matpoints(xseq, blc[,m.names], col=rgb(0,0,0,3, maxColorValue=100), pch=16, cex=.03)
    
    #cols <- rainbow(length(m.names),start=.55)
 
    library(cluster)
    blc<-dat$blc
    pam5 <- pam(t(blc[,m.names]),k=subset.n)
    s.names <- row.names(pam5$medoids)
    pam5.tab <- table(pam5$clustering)
    tags <- paste(paste("#",names(pam5.tab),sep=""),as.vector(pam5.tab),sep=":")
    info<-pam5$clustering
    
    ## Tool For adding color to selected Traces
    cols <-brewer.pal(8,"Dark2")
    cols <- rep(cols,ceiling(length(s.names)/length(cols)))
    cols <- cols[1:length(s.names)]

    ## Tool for adding labeling for single line within stacked traces
    for(i in 1:length(s.names)){
        matlines(xseq, blc[,names(which(info==i, arr.ind=T))]+i*sf, col=rgb(0,0,0,10, maxColorValue=100), lwd=.01)
        lines(xseq, blc[,s.names[i]]+i*sf, col=cols[i], lwd=.5)
        points(xseq, blc[,s.names[i]]+i*sf, col=cols[i], pch=16, cex=.03)
        text(x=min(blc[,1]), y=blc[nrow(t.dat),s.names[i]]+i*sf, labels=s.names[i], col=cols[i], pos=2, cex=bcex)
        text(x=max(blc[,1]), y=blc[nrow(t.dat),s.names[i]]+i*sf, labels=tags[i], col=cols[i], pos=4, cex=bcex)
    }
    
    return(pam5$clustering)	
} 
 
Lines.Multi<-function(dat,n.names){
    dev.new(width=2, height=2)
    par(mar=c(0,0,0,0))
    plot(0,0, pch=NA, xlim=c(0,2), ylim=c(0,2))
    points(x=c(1,1), y=c(1.5,1), pch=15)
    text(x=c(1,1), y=c(1.5,1), c("next", "off"), pos=2)

    dev.new()
    click.i<-0
    i<-1

    while(click.i!=2){
        dev.set(dev.list()[2])
        LinesEvery.2(dat,n.names[i:(10+i)], m.order="area", plot.new=F)
        dev.set(dev.list()[1])
        click.i<-identify(x=c(1,1), y=c(1.5,1), n=1)
        if(click.i==1){i<-i+10}
    }
    graphics.off()
}

linesmean<-function(dat, x.names,t.type=NULL, ylim=NULL, bcex=NULL, cols=NULL,lmain=NULL, lines.all=T, pic.plot=F){
    if(is.null(ylim)){ylim<-c(0,1.5)}else{ylim<-ylim}
    if(is.null(bcex)){bcex<-.9}else{bcex<-bcex}
    if(is.null(cols)){cols<-"red"}else{cols<-cols}

    if(is.null(t.type)){t.type<-select.list(names(dat))
    }else{t.type<-t.type}

    dat.t<-dat[[t.type]]

    dev.new(width=10,height=4)

    x.mean<-apply(dat.t[,x.names],1,mean)
    xseq<-dat$blc[,1]

    plot(xseq, x.mean, col="white", lwd=.2, ylim=ylim,main=lmain)

    levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    wr<-dat$w.dat$wr1
    x1s <- tapply(xseq,as.factor(wr),min)[levs]
    x2s <- tapply(xseq,as.factor(wr),max)[levs]
    y1s <- rep(par("usr")[3],length(x1s))
    y2s <- rep(par("usr")[4],length(x1s))
    rect(x1s,y1s,x2s,y2s,col="grey90",border="black")
    par(xpd=T)
    text(dat$w.dat[match(levs,wr),"Time"],rep(c(par("usr")[3]-yinch(.4),par("usr")[3]-yinch(.65)),length=length(levs)),levs,cex=bcex,offset=0, pos=4)#,offset=-offs}

    if(lines.all){
        matlines(xseq, dat.t[,x.names], col=rgb(0,0,0,20, maxColorValue=100), lwd=.01)
        }
    lines(xseq, x.mean, col=cols, lwd=.2)
    points(xseq, x.mean, col=cols, pch=16, cex=.02)

    if(pic.plot){
    cell.view()}
}

PulseViewer<-function(dat, cell, window.min=NULL, select.trace="t.dat"){
    
    if(class(select.trace)=="character"){
        dat.select<-select.trace
        dat.t<-dat[[dat.select]]
    }
    else{
        dat.select<-menu(names(dat))
        dat.t<-dat[[dat.select]]
    }


    window.region<-select.list(setdiff(unique(dat$w.dat$wr1),""))
    
    if(is.null(window.min)){window.min=7
    }else{window.min<-window.min}
    
    #What is the time frame from w.dat?
    window.time<-dat$w.dat[which(dat$w.dat$wr1==window.region,arr.ind=T,useNames=T),"Time"]
    #What is the maximun window defined
    window.max<-min(window.time)+window.min
    #what is the actual value
    window.region<-row.names(dat$w.dat[which(dat$w.dat$Time>=window.min & dat$w.dat$Time<=window.max, useNames=T),"Time"])
    
    plot(dat.t[window.region,"Time"], dat.t[window.region, cell])
    dat.t[,c]~dat.t[,1]
    
}

#Display the analysis of a single trace 
#dat is the trace dataframe with "Time" in the first column and cell trace intensities in subsequent columns
#i is the index column to be analyzed and displayed.
#shws is the smoothing half window size
#Plotit is a flag indicating that the results should be ploted or not.
#wr is the response window factor 
#SNR.lim is the signal to noise ratio limit for peak detection
#bl.meth is the method for baseline correction.
PeakFunc2 <- function(dat,i,shws=2,phws=20,Plotit=F,wr=NULL,SNR.lim=2,bl.meth="TopHat",lmain=NULL){
    library("MALDIquant")
    s1 <- createMassSpectrum(dat[,"Time"],dat[,i])
    if(shws > 1)
        s3 <- smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=shws)
    else
        s3 <- s1
    if(Plotit)
    {
        bSnip <- estimateBaseline(s3, method="SNIP")
        bTopHat <- estimateBaseline(s3, method="TopHat")
    }
    s4 <- removeBaseline(s3, method=bl.meth)
    Baseline <- estimateBaseline(s3, method=bl.meth)
    p <- detectPeaks(s4, method="MAD", halfWindowSize=phws, SNR=SNR.lim)
    if(Plotit)
    {
        xlim <- range(mass(s1)) # use same xlim on all plots for better comparison
        ylim <- c(-.1,1.4)
    #        ylim <- range(intensity(s1))
        plot(s1, main=paste(lmain,i),xlim=xlim,ylim=ylim,xlab="Time (min)", xaxt="n")
        axis(1, at=seq(0, length(dat[,1]), 5))  
        if(length(wr) > 0)
        {
            levs <- setdiff(unique(wr),"")
            levs <- setdiff(levs,grep("blank",levs,value=T))
            x1s <- tapply(dat[,"Time"],as.factor(wr),min)[levs]
            x2s <- tapply(dat[,"Time"],as.factor(wr),max)[levs]
            y1s <- rep(min(ylim)-.2,length(x1s))
            y2s <- rep(max(ylim)+.2,length(x1s))
    #            cols <- rainbow(length(x1s))
            rect(x1s,y1s,x2s,y2s,col="lightgrey")
    #            points(dat[,"Time"],as.integer(wr=="")*-1,pch=15,cex=.6)
            ## for(j in levs)
            ## {
            ##     x1 <- mass(s3)[min(grep(j,wr))]
            ##     x2 <- mass(s3)[max(grep(j,wr))]
            ##     y1 <- min(ylim)-.2
            ##     y2 <- max(ylim)+.2
            ##     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col="lightgrey",lwd=.1)
            ## }
            text(dat[match(levs,wr),"Time"],rep(-.1,length(levs)),levs,pos=4,offset=0,cex=.5)
        }
        
        lines(s3,lwd=3,col="cyan")
        lines(s1)
        lines(bSnip, lwd=2, col="red")
        lines(bTopHat, lwd=2, col="blue")
        lines(s4,lwd=2)
    }
    if((length(p) > 0)&Plotit)
    {
        points(p)
        ## label top 40 peaks
        top40 <- intensity(p) %in% sort(intensity(p), decreasing=TRUE)[1:40]
        labelPeaks(p, index=top40, underline=TRUE,labels=round(snr(p)[top40],2))
    }
    return(list(peaks=p,baseline=Baseline,dat=s4))
}

PeakFunc3 <- function(dat,n.names,shws=2,phws=20,wr=NULL,SNR.lim=2,bl.meth="TopHat",lmain=NULL){

    xlim <- range(dat$t.dat[,1]) # use same xlim on all plots for better comparison
    ylim <- c(-.1,1.4)
    #   ylim <- range(intensity(s1))
    levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    plot(dat$t.dat[,n.names],dat$t.dat[,1], main=paste(lmain,n.names),xlim=xlim,ylim=ylim,xlab="", xaxt="n",pch=16, lwd=1, cex=.5)
    axis(1, at=seq(0, length(dat$t.dat[,1]), 5))  
    lines(dat$t.dat[,n.names]~dat$t.dat[,1])
    points(dat$t.dat[,n.names]~dat$t.dat[,1], pch=16, cex=.4)
    
    lines(dat$blc[,n.names]~dat$t.dat[,1], lwd=1, cex=.5)
    points(dat$blc[,n.names]~dat$t.dat[,1], pch=16, cex=.4)
    
    # Tool for labeling window regions
    if(is.null(wr)){
        wr<-dat$w.dat[,"wr1"]
        levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
        x1s <- tapply(dat$t.dat[,"Time"],as.factor(wr),min)[levs]
        x2s <- tapply(dat$t.dat[,"Time"],as.factor(wr),max)[levs]
        y1s <- rep(min(ylim)-.2,length(x1s))
        y2s <- rep(max(ylim)+.2,length(x1s))
        rect(x1s,y1s,x2s,y2s,col="grey95")
        text(dat$t.dat[match(levs,wr),"Time"],rep(-.1,length(levs)),levs,pos=4,offset=0,cex=.5)
    }
    
    # Tool for labeling the binary score
    if(length(levs)>0){
        levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
        z<-t(dat$bin[n.names,levs])
        zz<-z==1
        zi<-attributes(zz)
        zzz<-which(zz, arr.ind=T)
        #levs<-zi$dimnames[[2]][zzz[,2]]
        levs<-unique(as.character(row.names(zzz)))
        x1s <- tapply(dat$t.dat[,"Time"],as.factor(wr),min)[levs]
        x2s <- tapply(dat$t.dat[,"Time"],as.factor(wr),max)[levs]
        y1s <- rep(min(ylim)-.2,length(x1s))
        y2s <- rep(max(ylim)+.2,length(x1s))
        rect(x1s,y1s,x2s,y2s,col="grey69")
        levs <- setdiff(unique(wr),"")
        text(dat$t.dat[match(levs,wr),"Time"],rep(-.1,length(levs)),levs,pos=4,offset=0,cex=.5)
    }
    
    # Tool for labeling cellular aspects, gfp.1, gfp.2, tritc, area
    legend("topright", xpd=TRUE, inset=c(0,-.14), legend=c(
        if(!is.null(dat$c.dat[n.names, "CGRP"])){paste("CGRP","",round(dat$c.dat[n.names,"CGRP"],digits=0))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp"])){paste("GFP","",round(dat$c.dat[n.names,"mean.gfp"],digits=0))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.2"])){paste("GFP.2","",round(dat$c.dat[n.names,"mean.gfp.2"],digits=0))},
        if(!is.null(dat$c.dat[n.names, "IB4"])){paste("IB4","",round(dat$c.dat[n.names,"IB4"],digits=0))},
        if(!is.null(dat$c.dat[n.names, "mean.tritc"])){paste("IB4","",round(dat$c.dat[n.names, "mean.tritc"], digits=0))}, 
        if(!is.null(dat$c.dat[n.names, "area"])){paste("area","", round(dat$c.dat[n.names, "area"], digits=0))})
    ,bty="n", cex=.8)

    # Tool for lableing window region information
    x.name<-n.names
    levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])), "")
    levs.loc<-tapply(dat$t.dat[,"Time"],as.factor(wr),mean)[levs]
    mtext(c("snr", "tot", "max", "wm"), side=1, at=-1, line=c(1.4, 2.1, 2.8, 3.5), cex=.6)
    for(i in levs){
        snr.name<-grep(paste(i,".snr", sep=""), names(dat$scp), value=T)
        tot.name<-grep(paste(i,".tot", sep=""), names(dat$scp), value=T)
        max.name<-grep(paste(i,".max", sep=""), names(dat$scp), value=T)
        wm.name<-grep(paste(i,".wm", sep=""), names(dat$scp), value=T)
        snr.val<-round(dat$scp[x.name, snr.name], digits=1)
        tot.val<-round(dat$scp[x.name, tot.name], digits=2)
        max.val<-round(dat$scp[x.name, max.name], digits=2)
        wm.val<-round(dat$scp[x.name, wm.name], digits=1)
        mtext(snr.val, side=1, at=levs.loc[i], line=1.4, cex=.6)
        mtext(tot.val, side=1, at=levs.loc[i], line=2.1, cex=.6)
        mtext(max.val, side=1, at=levs.loc[i], line=2.8, cex=.6)
        mtext(wm.val, side=1, at=levs.loc[i], line=3.5, cex=.6)
    }
}

PeakFunc4 <- function(dat,n.names,Plotit.maldi=T,Plotit.der=T,lmain=NULL){
    
    par(mfrow=c(2,1))
    
    if(Plotit.der)
    {	
    ylim<-c(-1, 2)
    plot(dat$der[,n.names]~dat$t.dat[-1,1], ylim=ylim,type="l",ylab=expression(paste(Delta," (340/380)/time")),xlab="",main=paste("Derivative",n.names), xaxt="n",pch=16, lwd=1, cex=.5)	
    
    # Tool for labeling window regions
    wr<-dat$w.dat[,"wr1"]
    levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    x1s <- tapply(dat$t.dat[,"Time"],as.factor(wr),min)[levs]
    x2s <- tapply(dat$t.dat[,"Time"],as.factor(wr),max)[levs]
    y1s <- rep(min(ylim)-.2,length(x1s))
    y2s <- rep(max(ylim)+.2,length(x1s))
    rect(x1s,y1s,x2s,y2s,col="grey95")
    text(dat$t.dat[match(levs,wr),"Time"],rep(-1,length(levs)),levs,pos=4,offset=0,cex=.5)
    
    # Tool for labeling the binary score
    if(length(levs)>0){
        levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
        z<-t(dat$bin[n.names,levs])
        zz<-z==1
        zi<-attributes(zz)
        zzz<-which(zz, arr.ind=T)
        #levs<-zi$dimnames[[2]][zzz[,2]]
        levs<-unique(as.character(row.names(zzz)))
        x1s <- tapply(dat$t.dat[,"Time"],as.factor(wr),min)[levs]
        x2s <- tapply(dat$t.dat[,"Time"],as.factor(wr),max)[levs]
        y1s <- rep(min(ylim)-.2,length(x1s))
        y2s <- rep(max(ylim)+.2,length(x1s))
        rect(x1s,y1s,x2s,y2s,col="grey69")
        levs <- setdiff(unique(wr),"")
        text(dat$t.dat[match(levs,wr),"Time"],rep(-1,length(levs)),levs,pos=4,offset=0,cex=.5)
    }
    
    # Tool for lableing window region information
    x.name<-n.names
    levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])), "")
    levs.loc<-tapply(dat$t.dat[,"Time"],as.factor(wr),mean)[levs]
    mtext(c("tot", "max", "min", "wmax", "wmin"), side=1, at=-1, line=c(0.7,1.4, 2.1, 2.8, 3.5), cex=.6)
    for(i in levs){
        tot.name<-grep(paste(i,".der.tot", sep=""), names(dat$scp), value=T)
        max.name<-grep(paste(i,".der.max", sep=""), names(dat$scp), value=T)
        min.name<-grep(paste(i,".der.min", sep=""), names(dat$scp), value=T)
        wmax.name<-grep(paste(i,".der.wmax", sep=""), names(dat$scp), value=T)
        wmin.name<-grep(paste(i,".der.wmin", sep=""), names(dat$scp), value=T)
        
        tot.val<-round(dat$scp[x.name, tot.name], digits=2)
        max.val<-round(dat$scp[x.name, max.name], digits=2)
        min.val<-round(dat$scp[x.name, min.name], digits=2)
        wmax.val<-round(dat$scp[x.name, wmax.name], digits=2)
        wmin.val<-round(dat$scp[x.name, wmin.name], digits=2)

        mtext(tot.val, side=1, at=levs.loc[i], line=0.7, cex=.6)
        mtext(max.val, side=1, at=levs.loc[i], line=1.4, cex=.6)
        mtext(min.val, side=1, at=levs.loc[i], line=2.1, cex=.6)
        mtext(wmax.val, side=1, at=levs.loc[i], line=2.8, cex=.6)
        mtext(wmin.val, side=1, at=levs.loc[i], line=3.5, cex=.6)
    }
    lines(dat$der[,n.names]~dat$t.dat[-1,1], lwd=.01, col="black")
    abline(h=0.5)

    #axis(1, at=seq(0, length(dat$t.dat[,1]), 5))  
    }

    
    if(Plotit.maldi)
    {
        xlim <- range(dat$t.dat[,1]) # use same xlim on all plots for better comparison
        ylim <- c(0,1.4)
        #   ylim <- range(intensity(s1))
        levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
        plot(dat$t.dat[,n.names]~dat$t.dat[,1], main=paste(lmain,n.names),xlim=xlim,ylim=ylim,xlab="", ylab="(340/380)", xaxt="n",pch=16, lwd=1, cex=.5)
        axis(1, at=seq(0, length(dat$t.dat[,1]), 5))  

        
        # Tool for labeling window regions
        wr<-dat$w.dat[,"wr1"]
        levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
        x1s <- tapply(dat$t.dat[,"Time"],as.factor(wr),min)[levs]
        x2s <- tapply(dat$t.dat[,"Time"],as.factor(wr),max)[levs]
        y1s <- rep(min(ylim)-.2,length(x1s))
        y2s <- rep(max(ylim)+.2,length(x1s))
        rect(x1s,y1s,x2s,y2s,col="grey95")
        #text(dat$t.dat[match(levs,wr),"Time"],rep(-.1,length(levs)),levs,pos=4,offset=0,cex=.5)
        
        # Tool for labeling the binary score
        levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
        z<-t(dat$bin[n.names,levs])
        zz<-z==1
        zi<-attributes(zz)
        zzz<-which(zz, arr.ind=T)
        #levs<-zi$dimnames[[2]][zzz[,2]]
        levs<-unique(as.character(row.names(zzz)))
        x1s <- tapply(dat$t.dat[,"Time"],as.factor(wr),min)[levs]
        x2s <- tapply(dat$t.dat[,"Time"],as.factor(wr),max)[levs]
        y1s <- rep(min(ylim)-.2,length(x1s))
        y2s <- rep(max(ylim)+.2,length(x1s))
        rect(x1s,y1s,x2s,y2s,col="grey69")
        levs <- setdiff(unique(wr),"")
        text(dat$t.dat[match(levs,wr),"Time"],rep(-1,length(levs)),levs,pos=4,offset=0,cex=.5)
        
        # Tool for labeling cellular aspects, gfp.1, gfp.2, tritc, area
        legend("topright", xpd=TRUE, inset=c(0,-.14), legend=c(
            if(!is.null(dat$c.dat[n.names, "CGRP"])){paste("CGRP","",round(dat$c.dat[n.names,"CGRP"],digits=0))},
            if(!is.null(dat$c.dat[n.names, "mean.gfp"])){paste("GFP","",round(dat$c.dat[n.names,"mean.gfp"],digits=0))},
            if(!is.null(dat$c.dat[n.names, "mean.gfp.1"])){paste("GFP.1","",round(dat$c.dat[n.names,"mean.gfp.1"],digits=0))},
            if(!is.null(dat$c.dat[n.names, "mean.gfp.2"])){paste("GFP.2","",round(dat$c.dat[n.names,"mean.gfp.2"],digits=0))},
            if(!is.null(dat$c.dat[n.names, "IB4"])){paste("IB4","",round(dat$c.dat[n.names,"IB4"],digits=0))},
            if(!is.null(dat$c.dat[n.names, "mean.tritc"])){paste("IB4","",round(dat$c.dat[n.names, "mean.tritc"], digits=0))}, 
            if(!is.null(dat$c.dat[n.names, "area"])){paste("area","", round(dat$c.dat[n.names, "area"], digits=0))})
        ,bty="n", cex=.8)

        # Tool for lableing window region information
        x.name<-n.names
        levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])), "")
        levs.loc<-tapply(dat$t.dat[,"Time"],as.factor(wr),mean)[levs]
        mtext(c("snr", "tot", "max", "wm"), side=1, at=-1, line=c(1.4, 2.1, 2.8, 3.5), cex=.6)
        for(i in levs){
            snr.name<-grep(paste(i,".snr", sep=""), names(dat$scp), value=T)
            tot.name<-grep(paste(i,".tot", sep=""), names(dat$scp), value=T)
            max.name<-grep(paste(i,".max", sep=""), names(dat$scp), value=T)
            wm.name<-grep(paste(i,".wm", sep=""), names(dat$scp), value=T)
            snr.val<-round(dat$scp[x.name, snr.name], digits=1)
            tot.val<-round(dat$scp[x.name, tot.name], digits=2)
            max.val<-round(dat$scp[x.name, max.name], digits=2)
            wm.val<-round(dat$scp[x.name, wm.name], digits=1)
            mtext(snr.val, side=1, at=levs.loc[i], line=1.4, cex=.6)
            mtext(tot.val, side=1, at=levs.loc[i], line=2.1, cex=.6)
            mtext(max.val, side=1, at=levs.loc[i], line=2.8, cex=.6)
            mtext(wm.val, side=1, at=levs.loc[i], line=3.5, cex=.6)
        }
        
        lines(dat$t.dat[,n.names]~dat$t.dat[,1])
        points(dat$t.dat[,n.names]~dat$t.dat[,1], pch=16, cex=.4)
        
        lines(dat$blc[,n.names]~dat$t.dat[,1], lwd=1, cex=.5)
        points(dat$blc[,n.names]~dat$t.dat[,1], pch=16, cex=.4)
        
        
        #abline(h=.5)	

    }	
    
   # return(list(peaks=p,baseline=Baseline,dat=s4))
}

# Fixed y axis
# Photo addition
# Derivative plot
# win
PeakFunc5 <- function(dat,n.names,select.trace=F,Plotit.trace=T,Plotit.both=F, info=T,lmain=NULL, bcex=.7, ylim.max=1.6){
    if(is.null(ylim.max)){ylim.max<-1.4}else{ylim.max<-ylim.max}
    if(Plotit.trace){ylim <- c(-.1,ylim.max)}
    if(Plotit.both){ylim <- c(-.5,ylim.max)}
    par(xpd=FALSE)
    if(select.trace==TRUE){
        dat.select<-menu(names(dat))
        dat.t<-dat[[dat.select]]
        }
    else(dat.t<-dat$t.dat)
    xlim <- range(dat.t[,1]) # use same xlim on all plots for better comparison
    
    #   ylim <- range(intensity(s1))
    levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    par(mar=c(6,4.5,3.5,11))
    plot(dat.t[,n.names]~dat.t[,1], main=paste(lmain,n.names),xlim=xlim,ylim=ylim,xlab="", ylab="",pch=16, lwd=1, cex=.5)
    #axis(1, at=seq(0, length(dat.t[,1]), 5),tick=TRUE )  
    
    # Tool for labeling window regions
    wr<-dat$w.dat[,"wr1"]
    levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs]
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs]
    y1s <- rep(par("usr")[4],length(x1s))
    y2s <- rep(par("usr")[3],length(x1s))
    rect(x1s,y1s,x2s,y2s,col="grey95")
    #text(dat.t[match(levs,wr),"Time"],rep(-.1,length(levs)),levs,pos=4,offset=0,cex=.5)
    
    # Tool for labeling the binary score
    levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    z<-t(dat$bin[n.names,levs])
    zz<-z==1
    zi<-attributes(zz)
    zzz<-which(zz, arr.ind=T)
    #levs<-zi$dimnames[[2]][zzz[,2]]
    levs<-unique(as.character(row.names(zzz)))
    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs]
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs]
    y1s <- rep(par("usr")[4],length(x1s))
    y2s <- rep(par("usr")[3],length(x1s))
    rect(x1s,y1s,x2s,y2s,col="grey69")
    levs <- setdiff(unique(wr),"")
    text(dat.t[match(levs,wr),"Time"],c(min(ylim), .1),levs,pos=4,offset=0,cex=bcex)
    
    # Tool for labeling cellular aspects, gfp.1, gfp.2, tritc, area
    
    legend(x=par("usr")[2]-xinch(.4), y=par("usr")[4]+yinch(.5), xpd=TRUE, inset=c(0,-.14), legend=c(
        if(!is.null(dat$c.dat[n.names, "CGRP"])){paste("CGRP","",round(dat$c.dat[n.names,"CGRP"],digits=0))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp"])){paste("GFP","",round(dat$c.dat[n.names,"mean.gfp"],digits=0))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.1"])){paste("GFP.1","",round(dat$c.dat[n.names,"mean.gfp.1"],digits=0))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.2"])){paste("GFP.2","",round(dat$c.dat[n.names,"mean.gfp.2"],digits=0))},
        if(!is.null(dat$c.dat[n.names, "mean.dapi"])){paste("DAPI","",round(dat$c.dat[n.names,"mean.dapi"],digits=0))},
        if(!is.null(dat$c.dat[n.names, "IB4"])){paste("IB4","",round(dat$c.dat[n.names,"IB4"],digits=0))},
        if(!is.null(dat$c.dat[n.names, "mean.tritc"])){paste("IB4","",round(dat$c.dat[n.names, "mean.tritc"], digits=0))}, 
        if(!is.null(dat$c.dat[n.names, "area"])){paste("area","", round(dat$c.dat[n.names, "area"], digits=0))},
        if(!is.null(dat$c.dat[n.names, "ROI.Area"])){paste("area","", round(dat$c.dat[n.names, "ROI.Area"], digits=0))},
        #if(!is.null(dat$c.dat[n.names, "perimeter"])){paste("perimeter","", round(dat$c.dat[n.names, "perimeter"], digits=0))},
        if(!is.null(dat$c.dat[n.names, "circularity"])){paste("circularity","", round(dat$c.dat[n.names, "circularity"], digits=3))}
        )
    ,bty="n", cex=.7)

    #Adding binary scoring for labeling to plot
    par(xpd=TRUE)
    if(!is.null(dat$bin[n.names, "gfp.bin"])){text(y=1.9, x=max(dat.t[,1])*1.09, paste("mean.gfp :",dat$bin[n.names,"gfp.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "tritc.bin"])){text(y=1.9, x=max(dat.t[,1])*1.19, paste("IB4 :",dat$bin[n.names,"tritc.bin"]), cex=.7)}

    
    # Tool for lableing window region information
    if(info){
        x.name<-n.names
        levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])), "")
        levs.loc<-tapply(dat$w.dat[,"Time"],as.factor(wr),mean)[levs]
        mtext(c("max","tot","snr"), side=1, at=-max(dat.t[,1])*.05, line=c(1.4, 2.1, 2.8), cex=.6)
        for(i in levs){
            max.name<-paste(i,".max", sep="")
            max.val<-round(dat$scp[x.name, max.name], digits=3)
            mtext(max.val, side=1, at=levs.loc[i], line=1.4, cex=.6)
            
            tot.name<-paste(i,".tot", sep="")
            tot.val<-round(dat$scp[x.name, tot.name], digits=3)
            mtext(tot.val, side=1, at=levs.loc[i], line=2.1, cex=.6)
            
            snr.name<-paste(i,".snr", sep="")
            snr.val<-round(dat$scp[x.name, snr.name], digits=3)
            mtext(snr.val, side=1, at=levs.loc[i], line=2.8, cex=.6)

        }
    }
    
    
    par(xpd=FALSE)
    if(Plotit.both){
        if(!is.null(dat$der)){lines(dat$der[,n.names]~dat.t[-1,1], lwd=.01, col="paleturquoise4")}
        abline(h=0)
        lines(dat.t[,n.names]~dat.t[,1])
        points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
    }
    
    if(Plotit.trace){
        lines(dat.t[,n.names]~dat.t[,1])
        points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
    }
    
    ## Tool for adding rasterImages to plot
    img.dim<-dim(dat$img1)[1]
    zf<-20
    x<-dat$c.dat[n.names,"center.x"]
    left<-x-zf
    if(left<=0){left=0; right=2*zf}
    right<-x+zf
    if(right>=img.dim){left=img.dim-(2*zf);right=img.dim}
    
    y<-dat$c.dat[n.names,"center.y"]
    top<-y-zf
    if(top<=0){top=0; bottom=2*zf}
    bottom<-y+zf
    if(bottom>=img.dim){top=img.dim-(2*zf);bottom=img.dim}
    
    par(xpd=TRUE)
    ymax<-par("usr")[4]
    xmax<-par("usr")[2]
    if(!is.null(dat$img1)){
        img1<-dat$img1
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax
        ybottom<-ymax-yinch(.8)
        rasterImage(img1[top:bottom,left:right,],xleft,ybottom,xright,ytop)
    }
    
    if(!is.null(dat$img2)){
        img2<-dat$img2
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax
        ybottom<-ymax-yinch(.8)
        rasterImage(img2[top:bottom,left:right,],xleft,ybottom,xright,ytop)
    }

    if(!is.null(dat$img3)){
        img3<-dat$img3
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax-yinch(.8)
        ybottom<-ymax-yinch(1.6)
        rasterImage(img3[top:bottom,left:right,],xleft,ybottom,xright,ytop)
    }
    
    if(!is.null(dat$img4)){
        img4<-dat$img4
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax-yinch(.8)
        ybottom<-ymax-yinch(1.6)
        rasterImage(img4[top:bottom,left:right,],xleft,ybottom,xright,ytop)
    }
}

# Y axis self adjusting	works with trac.click.3
#select trace added to select trace to plot
#yvar: logical.  If true y axis will vary
#ylim.max how to set top y limits.  Single value only
#zf added 170127
PeakFunc6 <- function(dat,n.names,t.type="t.dat",Plotit.trace=T,Plotit.both=F, info=T,lmain=NULL, bcex=.7, yvar=F, ylim.max=NULL, zf=40){
    if(class(t.type)=="character"){
        dat.select<-t.type
        dat.t<-dat[[dat.select]]
    }
    else{
        dat.select<-menu(names(dat))
        dat.t<-dat[[dat.select]]
    }

    if(yvar){
        ymax<-max(dat.t[,n.names])*1.05
        ymin<-min(dat.t[,n.names])*.95
        yrange<-ymax-ymin
    }else{		
        if(is.null(ylim.max)){ylim.max<-1.4}else{ylim.max<-ylim.max}
        if(Plotit.trace){ylim <- c(-.1,ylim.max)}
        if(Plotit.both){ylim <- c(-.5,ylim.max)}
        ymin<-min(ylim)
        ymax<-max(ylim)
        yrange<-ymax-ymin
    }

    if(Plotit.trace){ylim <- c(ymin,ymax)}
    if(Plotit.both){ymin<- -.5;ylim <- c(ymin,ymax)}
    par(xpd=FALSE)
    xlim <- range(dat.t[,1]) # use same xlim on all plots for better comparison
    
    #   ylim <- range(intensity(s1))
    levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    par(mar=c(6,4.5,3.5,11))
    plot(dat.t[,n.names]~dat.t[,1], main=paste(lmain,n.names),xlim=xlim,ylim=ylim,xlab="", ylab="",pch="", cex=.5)
    #axis(1, at=seq(0, length(dat.t[,1]), 5),tick=TRUE )  
    
    # Tool for labeling window regions
    wr<-dat$w.dat[,"wr1"]
    levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs]
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs]
    y1s <- rep(par("usr")[4],length(x1s))
    y2s <- rep(par("usr")[3],length(x1s))
    rect(x1s,y1s,x2s,y2s,col="grey95")
    #text(dat.t[match(levs,wr),"Time"],rep(-.1,length(levs)),levs,pos=4,offset=0,cex=.5)
    
    
    # Tool for labeling cellular aspects, gfp.1, gfp.2, tritc, area
    legend(x=max(xlim)*.95, y=ymax+(.45*yrange), xpd=TRUE, inset=c(0,-.14), legend=c(
        if(!is.null(dat$c.dat[n.names, "CGRP"])){paste("CGRP","",round(dat$c.dat[n.names,"CGRP"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp"])){paste("GFP","",round(dat$c.dat[n.names,"mean.gfp"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.1"])){paste("GFP.1","",round(dat$c.dat[n.names,"mean.gfp.1"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.2"])){paste("GFP.2","",round(dat$c.dat[n.names,"mean.gfp.2"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.dapi"])){paste("DAPI","",round(dat$c.dat[n.names,"mean.dapi"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "IB4"])){paste("IB4","",round(dat$c.dat[n.names,"IB4"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.tritc"])){paste("IB4","",round(dat$c.dat[n.names, "mean.tritc"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "area"])){paste("area","", round(dat$c.dat[n.names, "area"], digits=4))},
        if(!is.null(dat$c.dat[n.names, "ROI.Area"])){paste("area","", round(dat$c.dat[n.names, "ROI.Area"], digits=4))},
        #if(!is.null(dat$c.dat[n.names, "perimeter"])){paste("perimeter","", round(dat$c.dat[n.names, "perimeter"], digits=0))},
        if(!is.null(dat$c.dat[n.names, "circularity"])){paste("circularity","", round(dat$c.dat[n.names, "circularity"], digits=4))}
        )
    ,bty="n", cex=.7)
    
    #Adding binary scoring for labeling to plot
    par(xpd=TRUE)
    if(!is.null(dat$bin[n.names, "gfp.bin"])){text(y=ymax+(.25*yrange), x=max(dat.t[,1])*1.09, paste("mean.gfp :",dat$bin[n.names,"gfp.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "tritc.bin"])){text(y=ymax+(.25*yrange), x=max(dat.t[,1])*1.19, paste("IB4 :",dat$bin[n.names,"tritc.bin"]), cex=.7)}


    # Tool for lableing window region information
    if(info){
        x.name<-n.names
        levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])), "")
        levs.loc<-tapply(dat$w.dat[,"Time"],as.factor(wr),mean)[levs]
        mtext(c("max","tot"), side=1, at=-max(dat.t[,1])*.05, line=c(1.4, 2.1), cex=.6)
        for(i in levs){
            max.name<-paste(i,".max", sep="")
            max.val<-round(dat$scp[x.name, max.name], digits=3)
            mtext(max.val, side=1, at=levs.loc[i], line=1.4, cex=.6)
            
            tot.name<-paste(i,".tot", sep="")
            tot.val<-round(dat$scp[x.name, tot.name], digits=3)
            mtext(tot.val, side=1, at=levs.loc[i], line=2.1, cex=.6)
        }
        
    # Tool for labeling the binary score
        levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
        z<-t(dat$bin[n.names,levs])
        zz<-z==1
        zi<-attributes(zz)
        zzz<-which(zz, arr.ind=T)
        #levs<-zi$dimnames[[2]][zzz[,2]]
        levs<-unique(as.character(row.names(zzz)))
        x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs]
        x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs]
        y1s <- rep(par("usr")[4],length(x1s))
        y2s <- rep(par("usr")[3],length(x1s))
        rect(x1s,y1s,x2s,y2s,col="grey69")
        levs <- setdiff(unique(wr),"")
    }
    text(dat.t[match(levs,wr),"Time"],c(ymin, ymin+(yrange*.2)),levs,pos=4,offset=0,cex=bcex)	
    
    if(Plotit.both){
        if(!is.null(dat$der)){lines(dat$der[,n.names]~dat.t[-1,1], lwd=.01, col="paleturquoise4")}
        par(xpd=T)
        abline(h=0)
        lines(dat.t[,n.names]~dat.t[,1])
        points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
        par(xpd=F)
    }
    
    if(Plotit.trace){
        par(xpd=T)
        lines(dat.t[,n.names]~dat.t[,1])
        points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
        par(xpd=F)
    }
    
    ## Tool for adding rasterImages to plot
    
    ###Finding the picture loaction of the cells
    if(is.null(zf)){zf<-20
    }else{zf<-zf}

    img.dim<-dim(dat$img1)[1]
    x<-dat$c.dat[n.names,"center.x"]
    left<-x-zf
    if(left<=0){left=0; right=2*zf}
    right<-x+zf
    if(right>=img.dim){left=img.dim-(2*zf);right=img.dim}
    
    y<-dat$c.dat[n.names,"center.y"]
    top<-y-zf
    if(top<=0){top=0; bottom=2*zf}
    bottom<-y+zf
    if(bottom>=img.dim){top=img.dim-(2*zf);bottom=img.dim}
    
    par(xpd=TRUE)
    
    ### Where to plot pictures
    #ymax<-max(dat.t[,n.names])*1.05
    #ymin<-min(dat.t[,n.names])*.95
    #yrange<-ymax-ymin
    

    
    ymax<-par("usr")[4]
    xmax<-par("usr")[2]
    if(!is.null(dat$img1)){
        img1<-dat$img1
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax+yinch(.8)
        ybottom<-ymax
        rasterImage(img1[top:bottom,left:right,],xleft,ybottom,xright,ytop)
    }
    
    if(!is.null(dat$img2)){
        img2<-dat$img2
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax+yinch(.8)
        ybottom<-ymax
        rasterImage(img2[top:bottom,left:right,],xleft,ybottom,xright,ytop)
    }

    if(!is.null(dat$img3)){
        img3<-dat$img3
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax
        ybottom<-ymax-yinch(.8)
        rasterImage(img3[top:bottom,left:right,],xleft,ybottom,xright,ytop)
    }
    
    if(!is.null(dat$img4)){
        img4<-dat$img4
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax
        ybottom<-ymax-yinch(.8)
        tryCatch(rasterImage(img4[top:bottom,left:right,],xleft,ybottom,xright,ytop),error=function(e) rasterImage(img4[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }

    if(!is.null(dat$img5)){
        img5<-dat$img5
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax-yinch(.8)
        ybottom<-ymax-yinch(1.6)
        rasterImage(img5[top:bottom,left:right,],xleft,ybottom,xright,ytop)
    }
    
    if(!is.null(dat$img6)){
        img6<-dat$img6
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax-yinch(.8)
        ybottom<-ymax-yinch(1.6)
        rasterImage(img6[top:bottom,left:right,],xleft,ybottom,xright,ytop)
    }
    
    if(!is.null(dat$img7)){
        img7<-dat$img7
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax-yinch(1.6)
        ybottom<-ymax-yinch(2.4)
        rasterImage(img7[top:bottom,left:right,],xleft,ybottom,xright,ytop)
    }
    
    if(!is.null(dat$img8)){
        img8<-dat$img8
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax-yinch(1.6)
        ybottom<-ymax-yinch(2.4)
        rasterImage(img8[top:bottom,left:right,],xleft,ybottom,xright,ytop)
    }
    

    
    
}	

#This peak func allows for multiple t.types to be plotted
#170515: added pts and lns: (logical)
#added dat.n for insertation of the name for the rd file
PeakFunc7 <- function(dat,n.names,t.type="t.dat",Plotit.trace=T,Plotit.both=F, info=T,lmain=NULL, bcex=.7, yvar=T, ylim.max=NULL, zf=40, pts=T, lns=T, levs=NULL, underline=T, dat.n=""){
    dat.name<-deparse(substitute(dat))
    if(dat.name=="dat"){dat.name<-dat.n
    }else{dat.name<-dat.name}	
    
    if(is.null(lmain)){
        lmain=n.names
    }else{lmain=lmain}
    if(class(t.type)=="character")
    {
        dat.select<-t.type
        dat.t<-dat[[dat.select]]
    }else{
        dat.select<-select.list(names(dat), multiple=T)
        dat.t<-dat[[dat.select]]
    }

    if(yvar){
        ymax<-max(dat.t[,n.names])*1.05
        ymin<-min(dat.t[,n.names])*.95
        yrange<-ymax-ymin
    }else{		
        if(is.null(ylim.max)){ylim.max<-1.4}else{ylim.max<-ylim.max}
        if(Plotit.trace){ylim <- c(-.1,ylim.max)}
        if(Plotit.both){ylim <- c(-.5,ylim.max)}
        ymin<-min(ylim)
        ymax<-max(ylim)
        yrange<-ymax-ymin
    }

    if(Plotit.trace){ylim <- c(ymin,ymax)}
    if(Plotit.both){ymin<- -.5;ylim <- c(ymin,ymax)}
    par(xpd=FALSE)
    xlim <- range(dat.t[,1]) # use same xlim on all plots for better comparison
    
    #   ylim <- range(intensity(s1))
    if(is.null(levs)){levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    }else{levs<-levs}
    par(mar=c(9,6.2,3.5,13), bty="n")
    plot(dat.t[,n.names]~dat.t[,1], main=lmain,xlim=xlim,ylim=ylim,xlab="", ylab="",pch="", cex=.5)

    #axis(3,tick=TRUE, outer=F )
    axis(1, at= seq(0, max(dat.t[,1]),10), tick=TRUE)
    
    # Tool for labeling window regions
    wr<-dat$w.dat[,"wr1"]
    #levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs]
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs]
    y1s <- rep(par("usr")[4],length(x1s))
    y2s <- rep(par("usr")[3],length(x1s))
    rect(x1s,y1s,x2s,y2s,col="grey95")
    
    
    # Tool for labeling cellular aspects, gfp.1, gfp.2, tritc, area 
    legend(x=par("usr")[1]-xinch(1.45), y=par("usr")[3]-yinch(.25), xpd=TRUE, inset=c(0,-.14),bty="n", cex=.7, legend=c(
        if(!is.null(dat$c.dat[n.names, "CGRP"])){paste("CGRP","",round(dat$c.dat[n.names,"CGRP"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp"])){paste("GFP","",round(dat$c.dat[n.names,"mean.gfp"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.1"])){paste("GFP.1","",round(dat$c.dat[n.names,"mean.gfp.1"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.2"])){paste("GFP.2","",round(dat$c.dat[n.names,"mean.gfp.2"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.start"])){paste("mean.gfp.start","",round(dat$c.dat[n.names,"mean.gfp.start"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.end"])){paste("mean.gfp.end","",round(dat$c.dat[n.names,"mean.gfp.end"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.immuno"])){paste("CGRP immunostain","",round(dat$c.dat[n.names,"mean.gfp.immuno"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.dapi"])){paste("DAPI","",round(dat$c.dat[n.names,"mean.dapi"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "IB4"])){paste("IB4","",round(dat$c.dat[n.names,"IB4"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.tritc"])){paste("IB4","",round(dat$c.dat[n.names, "mean.tritc"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.tritc.start"])){paste("IB4.start","",round(dat$c.dat[n.names, "mean.tritc.start"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.tritc.end"])){paste("IB4.end","",round(dat$c.dat[n.names, "mean.tritc.end"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.tritc.immuno"])){paste("NF200 immunostain","",round(dat$c.dat[n.names, "mean.tritc.immuno"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.cy5.start"])){paste("IB4.start","",round(dat$c.dat[n.names, "mean.cy5.start"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.cy5.end"])){paste("IB4.end","",round(dat$c.dat[n.names, "mean.cy5.end"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "area"])){paste("area","", round(dat$c.dat[n.names, "area"], digits=4))},
        if(!is.null(dat$c.dat[n.names, "ROI.Area"])){paste("area","", round(dat$c.dat[n.names, "ROI.Area"], digits=4))},
        #if(!is.null(dat$c.dat[n.names, "perimeter"])){paste("perimeter","", round(dat$c.dat[n.names, "perimeter"], digits=0))},
        if(!is.null(dat$c.dat[n.names, "circularity"])){paste("circularity","", round(dat$c.dat[n.names, "circularity"], digits=4))}
        )
    )
    
    legend(x=par("usr")[2]+xinch(.8), y=par("usr")[3]-yinch(.9), xpd=TRUE, inset=c(0,-.14), bty="n", cex=.7, legend=dat.name)

    
    #Adding binary scoring for labeling to plot
    par(xpd=TRUE)
    if(!is.null(dat$bin[n.names, "gfp.bin"])){text(y=par("usr")[4]+yinch(.5), x=par("usr")[2]+xinch(1.8), paste("GFP:",dat$bin[n.names,"gfp.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "tritc.bin"])){text(y=par("usr")[4]+yinch(.25), x=par("usr")[2]+xinch(1.8), paste("IB4 :",dat$bin[n.names,"tritc.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "cy5.bin"])){text(y=par("usr")[4]+yinch(.25), x=par("usr")[2]+xinch(1.8), paste("IB4 :",dat$bin[n.names,"cy5.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "drop"])){text(y=par("usr")[4]+yinch(0), x=par("usr")[2]+xinch(1.8), paste("Drop :",dat$bin[n.names,"drop"]), cex=.7)}


    # Tool for lableing window region information
    levs.loc<-tapply(dat$w.dat[,"Time"],as.factor(wr),mean)[levs]
    if(info){
        x.name<-n.names
        #levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])), "")
        mtext(c("max","snr"), side=3, at=-max(dat.t[,1])*.05, line=c(0, .7), cex=.6)
        for(i in 1:length(levs)){
            max.name<-paste(levs[i],".max", sep="")
            max.val<-round(dat$scp[x.name, max.name], digits=3)
            mtext(max.val, side=3, at=levs.loc[ levs[i] ], line=0, cex=.6)
            
            tot.name<-paste(levs[i],".snr", sep="")
            tot.val<-round(dat$scp[x.name, tot.name], digits=3)
            mtext(tot.val, side=3, at=levs.loc[ levs[i] ], line=.7, cex=.6)
        }
        
    # Tool for labeling the binary score
        #levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
        z<-t(dat$bin[n.names,levs])
        zz<-z==1
        zi<-attributes(zz)
        zzz<-which(zz, arr.ind=T)
        #levs<-zi$dimnames[[2]][zzz[,2]]
        levs1<-unique(as.character(row.names(zzz)))
        x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs1]
        x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs1]
        y1s <- rep(par("usr")[4],length(x1s))
        y2s <- rep(par("usr")[3],length(x1s))
        rect(x1s,y1s,x2s,y2s,col="grey80")
        #levs <- setdiff(unique(wr),"")
    }
    
    #text(dat.t[match(levs,wr),"Time"],c(ymin, ymin+(yrange*.2)),levs,pos=4,offset=0,cex=bcex)	
    #text(dat.t[match(levs,wr),"Time"],par("usr")[3],levs,pos=3,offset=-4.2,cex=bcex, srt=90)    
    levs_cex <- nchar(levs)
	levs_cex[ levs_cex <= 12*1.3  ] <- 1
	levs_cex[ levs_cex > 12*1.3  ] <- 12/levs_cex[ levs_cex>12*1.3  ]*1.3

    text(levs.loc,par("usr")[3],levs,pos=3,offset=-4.3,cex=levs_cex, srt=90)	

    if(Plotit.both){
        if(!is.null(dat$der)){lines(dat$der[,n.names]~dat.t[-1,1], lwd=.01, col="paleturquoise4")}
        par(xpd=T)
        abline(h=0)
        if(lns){lines(dat.t[,n.names]~dat.t[,1])
        }else{}
        if(pts){points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
        }else{}
        par(xpd=F)
    }
    
    if(Plotit.trace){
        par(xpd=T)
        if(lns){lines(dat.t[,n.names]~dat.t[,1])
        }else{}
        
        if(pts){points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
        }else{}
        
        par(xpd=F)
    }
    
    ##Tool for adding underline to plot
    if(underline){
        par(xpd=F)
        abline(h=min(dat.t[,n.names]), col="black")
        par(xpd=T)
    }else{}

    
    ## Tool for adding rasterImages to plot
    
    ###Finding the picture loaction of the cells
    if(!is.null(dat$img1)){
        if(is.null(zf)){zf<-20
        }else{zf<-zf}

        img.dim<-dim(dat$img1)[1]
        x<-dat$c.dat[n.names,"center.x"]
        left<-x-zf
        if(left<=0){left=0; right=2*zf}
        right<-x+zf
        if(right>=img.dim){left=img.dim-(2*zf);right=img.dim}
        
        y<-dat$c.dat[n.names,"center.y"]
        top<-y-zf
        if(top<=0){top=0; bottom=2*zf}
        bottom<-y+zf
        if(bottom>=img.dim){top=img.dim-(2*zf);bottom=img.dim}
        
        par(xpd=TRUE)
    }
    ### Where to plot pictures
    #ymax<-max(dat.t[,n.names])*1.05
    #ymin<-min(dat.t[,n.names])*.95
    #yrange<-ymax-ymin
    

    
    ymax<-par("usr")[4]
    xmax<-par("usr")[2]
    if(!is.null(dat$img1)){
        img1<-dat$img1
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax+yinch(.8)
        ybottom<-ymax
        tryCatch(
            rasterImage(img1[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img1[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img2)){
        img2<-dat$img2
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax+yinch(.8)
        ybottom<-ymax
        tryCatch(
            rasterImage(img2[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img2[top:bottom,left:right],xleft,ybottom,xright,ytop))

    }

    if(!is.null(dat$img3)){
        img3<-dat$img3
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax
        ybottom<-ymax-yinch(.8)
        tryCatch(
            rasterImage(img3[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img3[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img4)){
        img4<-dat$img4
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax
        ybottom<-ymax-yinch(.8)
        tryCatch(
            rasterImage(img4[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img4[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }

    if(!is.null(dat$img5)){
        img5<-dat$img5
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax-yinch(.8)
        ybottom<-ymax-yinch(1.6)
        tryCatch(
            rasterImage(img5[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img5[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img6)){
        img6<-dat$img6
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax-yinch(.8)
        ybottom<-ymax-yinch(1.6)
        tryCatch(
            rasterImage(img6[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img6[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img7)){
        img7<-dat$img7
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax-yinch(1.6)
        ybottom<-ymax-yinch(2.4)
        tryCatch(
            rasterImage(img7[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img7[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img8)){
        img8<-dat$img8
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax-yinch(1.6)
        ybottom<-ymax-yinch(2.4)
        tryCatch(
            rasterImage(img8[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img8[top:bottom,left:right],xleft,ybottom,xright,ytop))

    }
    

    
    
}	

PeakFunc8 <- function(dat, n.names, t.type="t.dat", xlim=NULL, Plotit.trace=T,Plotit.both=F, info=T, lmain=NULL, bcex=.7, yvar=T, ylim=NULL, zf=40, pts=T, lns=T, levs=NULL, underline=T, dat.n="", lwd=1, img_plot=T){
    #How to Add the name of the experiment to plot
    dat.name<-deparse(substitute(dat))
    if(dat.name=="dat"){dat.name<-dat.n
    }else{dat.name<-dat.name}	
    #How to add a name to the plot automatically
    if(is.null(lmain)){
        lmain=n.names
    }else{lmain=lmain}
    #Choose the trace to display on the plot
    if(class(t.type)=="character")
    {
        dat.select<-t.type
        dat.t<-dat[[dat.select]]
    }else{
        dat.select<-select.list(names(dat), multiple=T)
        dat.t<-dat[[dat.select]]
    }
    #y limit plots 
    if(yvar){
        ymax<-max(dat.t[,n.names])*1.05
        ymin<-min(dat.t[,n.names])*.95
        yrange<-ymax-ymin
    }else{		
        ylim<-ylim
    }
    
    #if(Plotit.trace){ylim <- c(ymin,ymax)}
    #if(Plotit.both){ymin<- -.5;ylim <- c(ymin,ymax)}
    par(xpd=FALSE)
   
    #Tool to chagne the display of the xlimits
    if(is.null(xlim)){
        xlim <- range(dat.t[,1])# use same xlim on all plots for better comparison
    }else{xlim<-xlim}
    
    #   ylim <- range(intensity(s1))
    if(is.null(levs)){levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    }else{levs<-levs}
    if(img_plot){
        par(mar=c(9,6.2,3.5,13), bty="n")
    }else{
        par(mar=c(9,6.2,3.5,5), bty="n")
    }
    plot(dat.t[,n.names]~dat.t[,1], main=lmain,xlim=xlim,ylim=ylim,xlab="", ylab=expression(paste(Delta,"F/F")),pch="", cex=.5)

    #axis(3,tick=TRUE, outer=F )
    axis(1, at= seq(0, max(dat.t[,1]),10), tick=TRUE)
    
    # Tool for labeling window regions
    wr<-dat$w.dat[,"wr1"]
    #levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs]
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs]
    y1s <- rep(par("usr")[4],length(x1s))
    y2s <- rep(par("usr")[3],length(x1s))
    rect(x1s,y1s,x2s,y2s,col="grey95")
    
    # Tool for labeling cellular aspects, gfp.1, gfp.2, tritc, area 
    legend(x=par("usr")[1]-xinch(1.45), y=par("usr")[3]-yinch(.25), xpd=TRUE, inset=c(0,-.14),bty="n", cex=.7, legend=c(
        if(!is.null(dat$c.dat[n.names, "CGRP"])){paste("CGRP","",round(dat$c.dat[n.names,"CGRP"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp"])){paste("GFP","",round(dat$c.dat[n.names,"mean.gfp"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.1"])){paste("GFP.1","",round(dat$c.dat[n.names,"mean.gfp.1"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.2"])){paste("GFP.2","",round(dat$c.dat[n.names,"mean.gfp.2"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.start"])){paste("mean.gfp.start","",round(dat$c.dat[n.names,"mean.gfp.start"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.end"])){paste("mean.gfp.end","",round(dat$c.dat[n.names,"mean.gfp.end"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.immuno"])){paste("CGRP immunostain","",round(dat$c.dat[n.names,"mean.gfp.immuno"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.dapi"])){paste("DAPI","",round(dat$c.dat[n.names,"mean.dapi"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "IB4"])){paste("IB4","",round(dat$c.dat[n.names,"IB4"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.tritc"])){paste("IB4","",round(dat$c.dat[n.names, "mean.tritc"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.tritc.start"])){paste("IB4.start","",round(dat$c.dat[n.names, "mean.tritc.start"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.tritc.end"])){paste("IB4.end","",round(dat$c.dat[n.names, "mean.tritc.end"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.tritc.immuno"])){paste("NF200 immunostain","",round(dat$c.dat[n.names, "mean.tritc.immuno"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.cy5.start"])){paste("IB4.start","",round(dat$c.dat[n.names, "mean.cy5.start"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.cy5.end"])){paste("IB4.end","",round(dat$c.dat[n.names, "mean.cy5.end"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "area"])){paste("area","", round(dat$c.dat[n.names, "area"], digits=4))},
        if(!is.null(dat$c.dat[n.names, "ROI.Area"])){paste("area","", round(dat$c.dat[n.names, "ROI.Area"], digits=4))},
        #if(!is.null(dat$c.dat[n.names, "perimeter"])){paste("perimeter","", round(dat$c.dat[n.names, "perimeter"], digits=0))},
        if(!is.null(dat$c.dat[n.names, "circularity"])){paste("circularity","", round(dat$c.dat[n.names, "circularity"], digits=4))}
        )
    )
    
    legend(x=par("usr")[2]+xinch(.8), y=par("usr")[3]-yinch(.9), xpd=TRUE, inset=c(0,-.14), bty="n", cex=.7, legend=dat.name)

    #Adding binary scoring for labeling to plot
    par(xpd=TRUE)
    if(!is.null(dat$bin[n.names, "gfp.bin"])){text(y=par("usr")[4]+yinch(.5), x=par("usr")[2]+xinch(1.8), paste("GFP:",dat$bin[n.names,"gfp.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "tritc.bin"])){text(y=par("usr")[4]+yinch(.25), x=par("usr")[2]+xinch(1.8), paste("IB4 :",dat$bin[n.names,"tritc.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "cy5.bin"])){text(y=par("usr")[4]+yinch(.25), x=par("usr")[2]+xinch(1.8), paste("IB4 :",dat$bin[n.names,"cy5.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "drop"])){text(y=par("usr")[4]+yinch(0), x=par("usr")[2]+xinch(1.8), paste("Drop :",dat$bin[n.names,"drop"]), cex=.7)}

    # Tool for lableing window region information
    levs.loc<-tapply(dat$w.dat[,"Time"],as.factor(wr),mean)[levs]
    if(info){
        x.name<-n.names
        #levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])), "")
        mtext(c("max","snr"), side=3, at=-max(dat.t[,1])*.05, line=c(0, .7), cex=.6)
        for(i in 1:length(levs)){
            max.name<-paste(levs[i],".max", sep="")
            max.val<-round(dat$scp[x.name, max.name], digits=3)
            mtext(max.val, side=3, at=levs.loc[ levs[i] ], line=0, cex=.6)
            
            tot.name<-paste(levs[i],".snr", sep="")
            tot.val<-round(dat$scp[x.name, tot.name], digits=3)
            mtext(tot.val, side=3, at=levs.loc[ levs[i] ], line=.7, cex=.6)
        }
        
    # Tool for labeling the binary score
        #levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
        z<-t(dat$bin[n.names,levs])
        zz<-z==1
        zi<-attributes(zz)
        zzz<-which(zz, arr.ind=T)
        #levs<-zi$dimnames[[2]][zzz[,2]]
        levs1<-unique(as.character(row.names(zzz)))
        x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs1]
        x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs1]
        y1s <- rep(par("usr")[4],length(x1s))
        y2s <- rep(par("usr")[3],length(x1s))
        rect(x1s,y1s,x2s,y2s,col="grey80")
        #levs <- setdiff(unique(wr),"")
    }
    
    #text(dat.t[match(levs,wr),"Time"],c(ymin, ymin+(yrange*.2)),levs,pos=4,offset=0,cex=bcex)	
    #text(dat.t[match(levs,wr),"Time"],par("usr")[3],levs,pos=3,offset=-4.2,cex=bcex, srt=90)
	levs_cex <- nchar(levs)
	levs_cex[ levs_cex<=12 ] <- 1
	levs_cex[ levs_cex > 12 ] <- 12/levs_cex[ levs_cex>12 ]*1.3
    text(levs.loc,par("usr")[3],levs,pos=3,offset=-4.2,cex=levs_cex, srt=90)	

    if(Plotit.both){
        if(!is.null(dat$der)){lines(dat$der[,n.names]~dat.t[-1,1], lwd=.01, col="paleturquoise4")}
        par(xpd=T)
        abline(h=0)
        if(lns){lines(dat.t[,n.names]~dat.t[,1], lwd=lwd)
        }else{}
        if(pts){points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
        }else{}
        par(xpd=F)
    }
    
    if(Plotit.trace){
        par(xpd=F)
        if(lns){lines(dat.t[,n.names]~dat.t[,1], lwd=lwd)
        }else{}
        
        if(pts){points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
        }else{}
        
        #par(xpd=F)
    }
    
    ##Tool for adding underline to plot
    if(underline){
        par(xpd=F)
        abline(h=min(dat.t[,n.names]), col="black")
        par(xpd=T)
    }else{}

    ## Tool for adding rasterImages to plot
    ######################################
    #Adding the image plots to the traces
    ######################################
    if(img_plot){
        ###Finding the picture loaction of the cells
        if(!is.null(dat$img1)){
            if(is.null(zf)){zf<-20
            }else{zf<-zf}

            img.dim<-dim(dat$img1)[1]
            x<-dat$c.dat[n.names,"center.x"]
            left<-x-zf
            if(left<=0){left=0; right=2*zf}
            right<-x+zf
            if(right>=img.dim){left=img.dim-(2*zf);right=img.dim}
            
            y<-dat$c.dat[n.names,"center.y"]
            top<-y-zf
            if(top<=0){top=0; bottom=2*zf}
            bottom<-y+zf
            if(bottom>=img.dim){top=img.dim-(2*zf);bottom=img.dim}
            
            par(xpd=TRUE)
        }
        ### Where to plot pictures
        #ymax<-max(dat.t[,n.names])*1.05
        #ymin<-min(dat.t[,n.names])*.95
        #yrange<-ymax-ymin

        ymax<-par("usr")[4]
        xmax<-par("usr")[2]
        if(!is.null(dat$img1)){
            img1<-dat$img1
            xleft<-xmax
            xright<-xmax+xinch(.8)
            ytop<-ymax+yinch(.8)
            ybottom<-ymax
            tryCatch(
                rasterImage(img1[top:bottom,left:right,],xleft,ybottom,xright,ytop),
                error=function(e) rasterImage(img1[top:bottom,left:right],xleft,ybottom,xright,ytop))
        }
        
        if(!is.null(dat$img2)){
            img2<-dat$img2
            xleft<-xmax+xinch(.8)
            xright<-xmax+xinch(1.6)
            ytop<-ymax+yinch(.8)
            ybottom<-ymax
            tryCatch(
                rasterImage(img2[top:bottom,left:right,],xleft,ybottom,xright,ytop),
                error=function(e) rasterImage(img2[top:bottom,left:right],xleft,ybottom,xright,ytop))

        }

        if(!is.null(dat$img3)){
            img3<-dat$img3
            xleft<-xmax
            xright<-xmax+xinch(.8)
            ytop<-ymax
            ybottom<-ymax-yinch(.8)
            tryCatch(
                rasterImage(img3[top:bottom,left:right,],xleft,ybottom,xright,ytop),
                error=function(e) rasterImage(img3[top:bottom,left:right],xleft,ybottom,xright,ytop))
        }
        
        if(!is.null(dat$img4)){
            img4<-dat$img4
            xleft<-xmax+xinch(.8)
            xright<-xmax+xinch(1.6)
            ytop<-ymax
            ybottom<-ymax-yinch(.8)
            tryCatch(
                rasterImage(img4[top:bottom,left:right,],xleft,ybottom,xright,ytop),
                error=function(e) rasterImage(img4[top:bottom,left:right],xleft,ybottom,xright,ytop))
        }

        if(!is.null(dat$img5)){
            img5<-dat$img5
            xleft<-xmax
            xright<-xmax+xinch(.8)
            ytop<-ymax-yinch(.8)
            ybottom<-ymax-yinch(1.6)
            tryCatch(
                rasterImage(img5[top:bottom,left:right,],xleft,ybottom,xright,ytop),
                error=function(e) rasterImage(img5[top:bottom,left:right],xleft,ybottom,xright,ytop))
        }
        
        if(!is.null(dat$img6)){
            img6<-dat$img6
            xleft<-xmax+xinch(.8)
            xright<-xmax+xinch(1.6)
            ytop<-ymax-yinch(.8)
            ybottom<-ymax-yinch(1.6)
            tryCatch(
                rasterImage(img6[top:bottom,left:right,],xleft,ybottom,xright,ytop),
                error=function(e) rasterImage(img6[top:bottom,left:right],xleft,ybottom,xright,ytop))
        }
        
        if(!is.null(dat$img7)){
            img7<-dat$img7
            xleft<-xmax
            xright<-xmax+xinch(.8)
            ytop<-ymax-yinch(1.6)
            ybottom<-ymax-yinch(2.4)
            tryCatch(
                rasterImage(img7[top:bottom,left:right,],xleft,ybottom,xright,ytop),
                error=function(e) rasterImage(img7[top:bottom,left:right],xleft,ybottom,xright,ytop))
        }
        
        if(!is.null(dat$img8)){
            img8<-dat$img8
            xleft<-xmax+xinch(.8)
            xright<-xmax+xinch(1.6)
            ytop<-ymax-yinch(1.6)
            ybottom<-ymax-yinch(2.4)
            tryCatch(
                rasterImage(img8[top:bottom,left:right,],xleft,ybottom,xright,ytop),
                error=function(e) rasterImage(img8[top:bottom,left:right],xleft,ybottom,xright,ytop))

        }
    }
}	

# USe to sort based on features from bin and c.dat
c.sort<-function(dat,char=NULL){
    tmp<-cbind(dat$c.dat, dat$bin)
    bob<-row.names(tmp[order(tmp[,char], decreasing=T),])
    return(bob)
    }

# Click through a set of selected cell and create a stack plot
# Could use labeling improvements
Trace.Click.1<-function(dat, cells=NULL){
    graphics.off()
    dev.new(width=14,height=4)    
    dev.new(width=12,height=8)  
    if(is.null(cells)){cnames <- names(dat$t.dat[,-1])}
    else{cnames<-cells}
    lines.flag <- 0
    cell.i <- 1
    g.names<-NULL
    click.i <- 1
    #group.names<-NULL
    linefunc <- function(dat,m.names,snr=NULL,lmain="",cols=NULL,m.order=NULL,rtag=NULL,rtag2=NULL,rtag3=NULL, sf=.25,lw=3,bcex=1,p.ht=7,p.wd=10)
    {
    t.dat<-dat$t.dat
    wr<-dat$w.dat[,2]
    levs<-unique(as.character(dat$w.dat[,2]))[-1]
    m.names <- intersect(m.names,names(t.dat))
    xseq <- t.dat[,1]
    
    library(RColorBrewer)
    if(length(m.names) > 0)
    {        
        if(!is.null(m.order)){	
        dat<-dat$c.dat[m.names,]
        n.order<-dat[order(dat[,m.order]),]
        m.names <- row.names(n.order)
        }
        #else{
            #m.pca <- prcomp(t(t.dat[,m.names]),scale=F,center=T)
            #morder <- m.pca$x[,1] * c(1,-1)[(sum(m.pca$rot[,1]) < 0)+1]
            #m.names <- m.names[order(m.pca$x[,1],decreasing=sum(m.pca$rot[,1]) < 0)]
            #um.names <- m.names[order(morder)]
        #}
        
        
        if(is.null(cols)){
        #cols <- rainbow(length(m.names),start=.55)
        cols <-brewer.pal(8,"Dark2")
        cols <- rep(cols,ceiling(length(m.names)/length(cols)))
        cols <- cols[1:length(m.names)]
        } 
        else { cols<-cols
         cols <- rep(cols,ceiling(length(m.names)/length(cols)))
         cols <- cols[1:length(m.names)]
        }
        
        hbc <- length(m.names)*sf+max(t.dat[,m.names])
        hb <- ceiling(hbc)
        par(mar=c(4,1,4,1))
        plot(xseq,t.dat[,m.names[1]],ylim=c(0,hbc),xlab="Time (min)",main=lmain,type="n", xaxt="n",yaxt="n",xlim=c(min(xseq)-1.5,max(xseq)+1.5))#-sf
        axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
        if(length(wr) > 0)
        {
            if(!is.null(levs))
            {
            #levs <- setdiff(unique(wr),"")
            x1s <- tapply(xseq,as.factor(wr),min)[levs]
            x2s <- tapply(xseq,as.factor(wr),max)[levs]
            y1s <- rep(-.3,length(x1s))
            y2s <- rep(hbc+.2,length(x1s))
            rect(x1s,y1s,x2s,y2s,col=NA,border="darkgrey")
            cpx <- xseq[match(levs,wr)+round(table(wr)[levs]/2,0)]
            offs <- nchar(levs)*.5
            text(cpx,rep(c(sf/2,sf),length=length(levs)),levs,pos=1,cex=bcex)#,offset=-offs
            }
        }
        for(i in 1:length(m.names))
        {
            lines(xseq,t.dat[,m.names[i]]+i*sf,col=cols[i],lwd=lw)
            if(!is.null(snr))
            {
            pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
            pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
                                        #                pp3 <- dat$crr[,m.names[i]] > 0
            points(xseq[pp1],t.dat[pp1,m.names[i]]+i/10,pch=1,col=cols[i])
            points(xseq[pp2],t.dat[pp2,m.names[i]]+i/10,pch=0,col=cols[i])
                                        #                points(xseq[pp3],t.dat[pp3,m.names[i]]+i/10,pch=2,col=cols[i],cex=.5)
                                        }    
        }
        text(rep(0,length(m.names)),seq(1,length(m.names))*sf+t.dat[1,m.names],m.names,cex=.8*bcex,col=cols,pos=2)
        
        if(is.null(rtag)){
        if(!is.null(m.order)){
            rtag <- dat$c.dat[m.names,m.order]
            text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag,cex=.8*bcex,col=cols,pos=4)
        }}
        else{
            rtag <- dat$c.dat[m.names,rtag]
            text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag2,cex=.8*bcex,col=cols,pos=4)
         }

        if(!is.null(rtag2)){
            rtag2 <- dat$c.dat[m.names,rtag2]
            text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag2,cex=.8*bcex,col="green4",pos=3)
            text(rep(max(xseq),length(n.names)),seq(1,length(n.names))*sf+t.dat[nrow(t.dat),n.names],rtag2,cex=.8*bcex,col="green4",pos=3)

        }
        if(!is.null(rtag3)){
            rtag3 <- dat$c.dat[m.names,rtag3]
            text(rep(max(xseq),length(m.names)),seq(1,length(m.names))*sf+t.dat[nrow(t.dat),m.names],rtag3,cex=.8*bcex,col="Red",pos=1)
        }

        } 
    }


    while(click.i!=4)
    {
        cell.pick <- cnames[cell.i]
        dev.set(dev.list()[1])
        p1 <- PeakFunc2(dat$mp,cell.pick,shws=2,phws=20,Plotit=T,wr=dat$w.dat$wr1,SNR.lim=2,bl.meth="SNIP")
        p1.par<-par()
        if(lines.flag==1){dev.set(dev.list()[2]);linefunc(dat, g.names);lines.flag <- 0}
        if(lines.flag==0){dev.set(dev.list()[1])}
        #title(sub=paste("Group ",group.i," n=",g.num," Cell ",cell.i,sep=""))
        xs <- rep(dat$t.dat[50,"Time"],4)
        points(x=xs,y=c(1.2,1.1,1.0,.9),pch=16)
        text(x=xs,y=c(1.2,1.1,1.0,.9),labels=c("Cell +","Cell -","Stack", "off"),pos=2,cex=.5)
        click.i <- identify(x=xs,y=c(1.2,1.1,1.0,.9),n=1,plot=F)
        
        if(click.i==1)
        {cell.i <- cell.i + 1;if(cell.i>length(cnames)){cell.i<-1}}
        if(click.i==2)
        {cell.i <- cell.i - 1;if(cell.i<1){cell.i<-length(cnames)}}
        if(click.i==3)
        {g.names<-union(g.names,cnames[cell.i]);lines.flag<-1}
        if(click.i==4){graphics.off()}
    }
    print(g.names)
}

# Click Throug cells, and zoom on cell of interest
Trace.Click.2<-function(dat, cells=NULL,img=NULL, plotit=T){
    graphics.off()
    dev.new(width=14,height=4)    
    dev.new(width=10,height=6)  
    dev.new(width=8, height=8)
    if(is.null(cells)){cnames <- names(dat$t.dat[,-1])}
    else{cnames<-cells}
    lines.flag <- 0
    cell.i <- 1
    g.names<-NULL
    click.i <- 1
    #group.names<-NULL

    while(click.i!=5)
    {
        cell.pick <- cnames[cell.i]
        dev.set(dev.list()[1])
        p1 <- PeakFunc5(dat,cell.pick,ylim.max=1.6)
        p1.par<-par()
        if(lines.flag==2){dev.set(dev.list()[3]);cell.veiw.2048(dat, img=img, cell=cell.pick, cells=cells,cols="red",plot.new=F,cell.name=T);lines.flag <- 0}
        if(lines.flag==1){dev.set(dev.list()[2]);LinesEvery.2(dat,g.names,plot.new=FALSE);lines.flag <- 0}
        if(lines.flag==0){dev.set(dev.list()[1])}
        #title(sub=paste("Group ",group.i," n=",g.num," Cell ",cell.i,sep=""))
        #xs <- -(rep(dat$t.dat[50,"Time"],5)*1.08)
        xs<- rep(par("usr")[1]-yinch(.2), 5)
        ys<-seq(par("usr")[4],by=-yinch(.5), length.out=5)
        points(x=xs,y=ys,pch=16)
        text(x=xs,y=ys,labels=c("Cell +","Cell -","Veiw","Stack","off"),pos=2,cex=.5)
        
        ## How many cells are you looking at
        maxy<-par("usr")[4]
        text(par("usr")[1], par("usr")[4]+yinch(.3),paste(cell.i, ":",length(cnames)))
        click.i <- identify(x=xs,y=ys,n=1,plot=F)
        
        if(click.i==1)
        {cell.i <- cell.i + 1;if(cell.i>length(cnames)){cell.i<-1};lines.flag<-0}
        if(click.i==2)
        {cell.i <- cell.i - 1;if(cell.i<1){cell.i<-length(cnames)};lines.flag<-0}
        if(click.i==3)
        {lines.flag<-2}
        if(click.i==4)
        {g.names<-union(g.names,cnames[cell.i]);lines.flag<-1}
        if(click.i==5){graphics.off()}
}
print(g.names)}

Trace.Click<-function(dat, cells=NULL,img=dat$img1, yvar=FALSE, t.type="t.dat", plot.new=F, info=T, pts=T, lns=T, bcex=1){
    if(plot.new){graphics.off()}
    dev.new(width=14,height=4)
    click.window<-dev.cur()
    
    dev.new(width=10,height=6) 
    lines.window<-dev.cur()
    
    dev.new(width=8, height=8)
    view.window<-dev.cur()
    
    if(is.null(cells)){cnames <- names(dat$t.dat[,-1])}
    else{cnames<-cells}

    lines.flag <- 0
    cell.i <- 1
    g.names<-NULL
    click.i <- 1
    #group.names<-NULL

    while(click.i!=9)
    {
        cell.pick <- cnames[cell.i]
        #dev.set(dev.list()[1])
        dev.set(which=click.window)
        p1 <- PeakFunc7(dat,cell.pick, t.type=t.type,yvar=yvar, info=info, bcex=bcex, pts=pts, lns=lns)
        p1.par<-par()
        if(lines.flag==1){
            #dev.set(dev.list()[2])
            dev.set(which=lines.window)
            LinesEvery.5(dat,g.names,plot.new=F, img=img, t.type=t.type, col="black")
            lines.flag <- 0
        }

        if(lines.flag==2){
            #dev.set(dev.list()[3])
            dev.set(which=view.window)
            cell.view(dat,cell=cell.pick, img=img,cols="red",plot.new=F,cell.name=T, zoom=FALSE)
            lines.flag <- 0
        }
        
        if(lines.flag==0){
            #dev.set(dev.list()[1]) 
            dev.set(which=click.window)
        }		
        
        #title(sub=paste("Group ",group.i," n=",g.num," Cell ",cell.i,sep=""))
        xs<- rep(par("usr")[1]-xinch(.5), 9)
        ys<-seq(par("usr")[4],by=-yinch(.2), length.out=9)
        points(x=xs,y=ys,pch=16)
        text(x=xs,y=ys,labels=c("Cell +","Cell -","Veiw","Stack","yvar","Select Trace","Points","Lines","off"),pos=2,cex=.5)
        
        ## How many cells are you looking at
        text(par("usr")[1], par("usr")[4]+yinch(.3),paste(cell.i, ":",length(cnames)))
        click.i <- identify(x=xs,y=ys,n=1,plot=F)
        
        if(click.i==1)
        {cell.i <- cell.i + 1;if(cell.i>length(cnames)){cell.i<-1};lines.flag<-0}
        
        if(click.i==2)
        {cell.i <- cell.i - 1;if(cell.i<1){cell.i<-length(cnames)};lines.flag<-0}
        
        if(click.i==3)
        {lines.flag<-2}
        
        if(click.i==4)
        {g.names<-union(g.names,cnames[cell.i]);lines.flag<-1}
        
        if(click.i==5){
            if(yvar){yvar<-FALSE}else{yvar<-TRUE}
        }
        
        if(click.i==6){
            t.type<-select.list(names(dat))
        }
        
        if(click.i==7){
            if(pts){pts<-FALSE}else{pts<-TRUE}
        }
        
        if(click.i==8){
            if(lns){lns<-FALSE}else{lns<-TRUE}
        }
        
        if(click.i==9){
            #graphics.off()
            dev.off(which=click.window)
            dev.off(which=lines.window)
            dev.off(which=view.window)
        }
    }
    print(g.names)
}

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

#170606 Added: 
#up arrow: move through list specified in entry
#down arrow: Move down through list spcified in entry
#c: add cells to g.names
#r: reset g.names
#1-0 : add cells to g.names1 through g.name10
#shift+# removes cell from that group
#s: stack g.names
#y: Zoom yaxis automatically
#t: brings up list of RD file. Select Trace
#o: order cells in a new way
#p: Toggles points on graph
##d: changes drop collumn to 1 automatically.  Remeber to save RD file at end of experiment
##k: changes drop collumn to 0 automatically.  Remeber to save RD file at end of experiment
#l: choose window region to display on stack trace plot
#i: Select image to display on multi view options
#I: image for stacked traces
#u: Underlines the Trace
#f: New trace fitting for pottassium pulses
#F: New smoothing factor for imputer
#z: zoom factor to apply to the view of the cell on the right side of the trace, and to the view window
tcd<-function(dat, cells=NULL,img=dat$img1, l.img=c("img1"), yvar=FALSE, t.type="t.dat", plot.new=F, info=T, pts=T, lns=T, bcex=1, levs=NULL, klevs=NULL, sft=NULL, underline=T, zf=20, lw=2, sf=1, dat.name=NULL, view_func_description=F, save_question = T){
    graphics.off()
    print(environment())
    if(is.null(dat.name)){
        dat.name<-deparse(substitute(dat))
    }else{dat.name<-dat.name}
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
    #y: Zoom yaxis automatically
    #z: image zoom
    ")
    }else{}

    dat.tmp<-dat
    if(plot.new){graphics.off()}
    if(is.null(sft)){sft<-7}
    
    windows(width=14,height=4,xpos=0, ypos=50)
    click.window<-dev.cur()
    
    windows(width=10,height=6,xpos=0, ypos=450) 
    lines.window<-dev.cur()
    
    dimx<-dim(img)[2]
    dimy<-dim(img)[1]
    haight<-10*dimy/dimx
    windows(width=haight*dimx/dimy, height=haight,xpos=1130, ypos=200)
    view.window<-dev.cur()
    
    windows(width=8, height=8,xpos=1130, ypos=0)
    multipic.window<-dev.cur()
    
    windows(width=12, height=2,xpos=0, ypos=550)
    traceimpute.window<-dev.cur()
    
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
    if(is.null(levs)){levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    }else{levs<-levs}
    
    if(is.null(klevs)){klevs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    }else{klevs<-levs}
    
    while(keyPressed!="q")
    {
        cell.pick <- cnames[cell.i]
        
        dev.set(which=click.window)
        p1 <- PeakFunc7(dat,cell.pick, t.type=t.type, yvar=yvar, info=info, bcex=bcex, pts=pts, lns=lns, levs=levs, underline=underline, dat.n=dat.name, zf=zf)
        p1.par<-par()
        
        ##LinesEvery
        if(lines.flag==1){
            #if(length(p.names)<100){
                if(length(p.names)>11){
                    dev.off(which=lines.window)
                    windows(width=10,height=12,xpos=0, ypos=100) 
                    lines.window<-dev.cur()
                }else{
                    dev.off(which=lines.window)
                    windows(width=10,height=7,xpos=0, ypos=250) 
                    lines.window<-dev.cur()
                }
                dev.set(which=lines.window)
                tryCatch(LinesEvery.5(dat,p.names,plot.new=F, img=l.img,lmain=paste(gsub("[$]","",p.namez), 'n=',length(p.names)), t.type=t.type, lw=lw, col=lines.color, lns=lns, levs=levs, bcex=1, underline=underline, dat.n=dat.name, zf=zf, sf=sf, values=values),error=function(e) print("You haven't stacked traces yet, yo."))
                lines.flag <- 0
            #}
        }
        
        if(lines.flag==2){
            sample.to.display<-as.numeric(select.list(as.character(c(5,10,20,50,70,100))),title='Sample Number?')
            tryCatch(dev.off(which=lines.window.2), error=function(e) print("this windows hasn't been opened yet"))
            windows(width=10,height=12,xpos=0, ypos=250) 
            lines.window.2<-dev.cur()
            dev.set(which=lines.window.2)
            
            tryCatch(
                LinesEvery.5(
                    dat,
                    sample(p.names)[1:sample.to.display],
                    plot.new=F,
                    lmain=paste("Sample",sample.to.display,"out of",length(p.names)), 
                    img=l.img, lw=lw, t.type=t.type, col="black", lns=lns, levs=levs, bcex=1, underline=underline, dat.n=dat.name, zf=zf,sf=sf, values=values)
            ,error=function(e) print("You haven't stacked traces yet, yo."))
            lines.flag <- 0
        }

        ##Pulse Imputer
        if(lines.flag==3){

            #dev.off(which=traceimpute.window)
            #windows(width=2*length(klevs),height=2,xpos=0, ypos=550) 
            #traceimpute.window<-dev.cur()
            
            dev.set(which=traceimpute.window)
            tryCatch(PulseImputer(dat,cell.pick,levs,sf=sf),error=function(e) print("You haven't stacked traces yet, yo."))
            lines.flag <- 0
        }
        
        ##Pic zoom
        if(window.flag==1){
            dev.set(which=view.window)
            tryCatch(cell.view(dat,cell=p.names, img=img,cols="Yellow",plot.new=F,cell.name=T, lmain=paste(gsub("[$]","",p.namez)), zoom=FALSE),error=function(e) print("You haven't collected cells to view"))
            
            dev.set(which=multipic.window)
            tryCatch(multi.pic.zoom(dat,p.names,img, plot.new=F, zf=zf, labs=F) ,error=function(e) print("You haven't collected cells to view"))
            window.flag <- 0
        }
        
        if(lines.flag==0){
            #dev.set(dev.list()[1]) 
            dev.set(which=click.window)
        }		
        
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
            cellTypeId <- grep('^cell',names(dat), value=T)
            if(length(cellTypeId) > 0){
                #get only cell types that we want to reassign
                cellTypeNames <- names(dat[[cellTypeId]])
                cellTypesToNotClean <- c('neurons', 'glia')
                cellTypesToClean <- setdiff(cellTypeNames, cellTypesToNotClean)
                
                #remove it from all groups
                dat[[cellTypeId]][cellTypesToClean] <- lapply(dat[[cellTypeId]][cellTypesToClean], function(X) setdiff(X,cellToReassign))
                
                ###now that we have indicated that we would like to place this cell into a new group
                #First lets find all groups we can assign to\
                bringToTop(-1)
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
            bringToTop(-1)
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
            # levsMiddle <- tapply(            levsMiddle <- tapply(dat$w.dat[,1], as.factor(dat$w.dat$wr1),mean)[levs]
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
            bringToTop(-1)
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
            tryCatch(p.names<-c.sort.2(dat,p.names),error=function(e) print("You have not stacked traces yet."))
            lines.flag<-1
            window.flag<-1
        }	
    #p: Toggles points on graph
        if(keyPressed=="p"){
            if(pts){pts<-FALSE}else{pts<-TRUE}
            lines.flag<-1
        }
    #P: Pick a group/cells to click through
        if(keyPressed=="P"){
            bringToTop(-1)
            print("Pick a Group of cells or a single cell to observe \nIf you Click cancel, all cells will be returned")
            selection<-select.list(c("group","cells"))
            if(selection=="group"){
                gt.to.click<-select.list(names(gt.names), multiple=F)
                if( is.null(gt.names[[gt.to.click]]) | is.logical( gt.names[[gt.to.click]]) ){
                    bringToTop(-1)
                    print("Nothing is in this Group")
                }else{
                    cell.i<-1
                    print(gt.to.click)
                    cnames<-gt.names[[gt.to.click]]
                    tryCatch(
                        cnames<-c.sort.2(dat[[datfram]],cnames,collumn=collumn),
                    error=function(e) print("Something went wrong try again") )
                    print(cnames)
                }
            }
            if(selection=="cells"){
                cell.i<-1
                cnames<-select.list(as.character(dat$c.dat$id), multiple=T)
                tryCatch(cnames<-c.sort.2(dat[[datfram]],cnames,collumn=collumn),error=function(e) print("Something went wrong try again"))

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
            bringToTop(-1)
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
            print(p.namez)
            p.names<-gt.names[[p.namez]]
            #p.names<-get(ls(pattern=p.namez))
            print(p.names)
            lines.flag<-1
        }
        
    #S: Sample selected groups
        if(keyPressed=="S"){
            p.namez<-paste(select.list(names(gt.names)),sep="")
            print(p.namez)
            p.names<-gt.names[[p.namez]]
            #p.names<-get(ls(pattern=p.namez))
            print(p.names)
            lines.flag<-2
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
            bringToTop(-1)
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
                drops <- row.names(dat$bin[dat$bin$drop,])
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
            bringToTop(-1)
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
            #first open a new window
            #after undergoing a logical test to see if it exists
            if(length(ls(pattern='bp.selector.window'))==0){
                dev.new(width=14, height=8)
                #give this window a name
                bp.selector.window<-dev.cur()
            }else{}
            #give the focus to the new window
            dev.set(bp.selector.window)
            #empty gt.names[[12]]
            gt.names[[12]]<-NA
            #remove the NA, which will be repalced with a logical(0)
            gt.names[[12]]<-lapply(gt.names[[12]], function(x) x[!is.na(x)])
            #do the function bp.selector to gather data
            bringToTop(-1)
            cat("This function allows you to create statistics based on the statistic you select. 
            \n This Function finds a represention of peak amplification and or block 
            \n This function will take in what ever you are currently scrolling through
            \n You have the option to localize your boxplot. This means, select cells
            \n specifically based on where you click on the boxplot. Two clicks means you need
            \n to specigy the lower range followed by the upper range.
            \n One click will take everything greater than your click
            \n The Other option that will arise is, would you like the save the stat.
            \n If you do, the console will prompt you to enter a name. Ensure no spaces in the name
            \n The next option will be whether you would like to make another statistic.")
            cat(" \n Would you like to localize your boxplot? \n")
            print("T=yes, F=no")
            localize_log<-scan(n=1,what="character")
            print(localize_log != "T")
            if(localize_log != "T"){localize_log<-"F"}
            print(cnames[cell.i])
            dev.set(bp.selector.window)
            gt.names[[12]]<-bp.selector(dat,cnames[cell.i],cnames,plot.new=F,dat.name=NULL,env=environment(),localize=localize_log)
            #Now fill TCD with the cells just selected.
            cnames<-gt.names[[12]]	
            cell.i<-1
            lines.flag<-1
            windows.flag<-1
        }
        
        #F2: Advanced Statistic maker This function uses the function (After-Before)/(After+Before)
        #this function allows you to save the stat.  This will be added to the scp dataframe at the bottom.
        #if you have created statistics, be sure to save your RD file before you close
        if(keyPressed=="F2"){
        
            #first open a new window
            #after undergoing a logical test to see if it exists
            if(length(ls(pattern='bp.selector.window'))==0){
                dev.new(width=14, height=8)
                #give this window a name
                bp.selector.window<-dev.cur()
            }else{}
            #give the focus to the new window
            dev.set(bp.selector.window)
            #empty gt.names[[12]]
            gt.names[[12]]<-NA
            #remove the NA, which will be repalced with a logical(0)
            gt.names[[12]]<-lapply(gt.names[[12]], function(x) x[!is.na(x)])
            #do the function bp.selector to gather data
            bringToTop(-1)
            cat("This function allows you to create statistics based on the statistic you select. 
            \n This Function finds a represention of peak amplification and or block 
            \n This function will take in what ever you are currently scrolling through
            \n You have the option to localize your boxplot. This means, select cells
            \n specifically based on where you click on the boxplot. Two clicks means you need
            \n to specigy the lower range followed by the upper range.
            \n One click will take everything greater than your click
            \n The Other option that will arise is, would you like the save the stat.
            \n If you do, the console will prompt you to enter a name. Ensure no spaces in the name
            \n The next option will be whether you would like to make another statistic.")
            cat(" \n Would you like to localize your boxplot? \n")
            print("T=yes, F=no")
            localize_log<-scan(n=1,what="character")
            print(localize_log != "T")
            if( length(localize_log) == 0 ){
                localize_log<-"F"
            }else{ 
                if(localize_log != "T"){
                    localize_log<-"F"
                }
            }
            print(cnames[cell.i])
            dev.set(bp.selector.window)
            gt.names[[12]]<-bp.selector.advanced(dat,cnames[cell.i],cnames,plot.new=F,dat.name=NULL,env=environment(),localize=localize_log)
            #Now fill TCD with the cells just selected.
            cnames<-gt.names[[12]]	
            cell.i<-1
            lines.flag<-1
            windows.flag<-1
        }
        

        #F3: Plotting the Density plots.  There are many options for this plot
        if(keyPressed=="F3"){
            if(length(ls(pattern="density_win"))==0){
                dev.new(width=10,height=10)
                density_win<-dev.cur()
            }else{}
            bringToTop(-1)
            cat("What dataframe wil contain your stat? \n")
            dense_df_q<-select.list(names(dat))
            cat("What attribute would you like to see the distribution? \n")
            dense_df_att<-menu(names(dat[[dense_df_q]]))
            statz<-dat[[dense_df_q]][dense_df_att]

            #define the top xlim value
            cat("Define Top xlim value \n")
            cat("Enter 0 to allow default Max value \n")
            xlim_top<-scan(n=1)
            if(xlim_top==0){
                xlim_top<-max(dat[[dense_df_q]][dense_df_att])
            }
            
            cat("Define bottom xlim value \n")
            xlim_bottom<-scan(n=1)
            if(xlim_bottom==0){
                xlim_bottom<-min(dat[[dense_df_q]][dense_df_att])
            }
            
            dev.set(density_win)
            density_ct_plotter(dat,g.names,cell_types=NULL, stat=statz,overlay=T, dense_sep=TRUE,plot_new=F,xlim_top=xlim_top,xlim_bottom=xlim_bottom,dat.name=dat.name)
            lines.flag<-1
        }
        
        #F4: Utilizing Topview
        if(keyPressed=="F4"){
            p.namez<-paste(select.list(names(gt.names)),sep="")
            p.names<-gt.names[[p.namez]]

            aux_var<-c('area')

            #What i need to do is selectively import gfp and tritc variables into the 
            #topview function
            #this means search in the bin data frame for ib4 and gfp
            add_vars <- grep('mcherry|cy5|gfp|drop', names(dat$bin),value=T)
            aux_var<-c(aux_var, add_vars)
            TopView(dat, p.names, 12, 6, dat_name=dat.name, aux.var=aux_var)
        }

        #F5: Censusus Viewer
        if(keyPressed=="F5"){
			cnames_orig <- cnames
			cells_to_view <- census_viewer(dat)
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
                p.namez <- cells_to_view$name
				p.names <- gt.names[[12]]
				lines.flag<-1
			}
		}    
        
        #F6: Censusus Viewer
        if(keyPressed=="F6"){
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
                bringToTop(-1)             
                print("\nI have filled in your cell_types to choose by pressing \'P\' ENJOY!\n")
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
            dev.off(which=click.window)
            dev.off(which=lines.window)
            tryCatch(dev.off(which=lines.window.2), error=function(e) print("this windows hasn't been opened yet"))			
            dev.off(which=view.window)
            dev.off(which=multipic.window)
            dev.off(which=traceimpute.window)

        }
    }
    #rd.name <- as.character(substitute(dat))
    #print(rd.name)
    #assign(rd.name, dat, envir=.GlobalEnv)
    #gt.names<-list(g.names1=g.names1, g.names2=g.names2, g.names3=g.names3, g.names4=g.names4, g.names5=g.names5, g.names6=g.names6, g.names7=g.names7, g.names8=g.names8,g.names9=g.names9, g.names10=g.names10, g.names11=g.names11, g.names12=g.names12, g.names=g.names)
    BACKUP<<-gt.names 
    
    assign(dat.name,dat, envir=.GlobalEnv)
    bringToTop(-1)
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

#create a trace.click that allows for scoring while clicking
Trace.Click.repair<-function(dat, cells=NULL,img=dat$img1, yvar=FALSE, t.type="t.dat", plot.new=F, info=T, bcex=1, save.bp=F,view.cells=F){
    if(plot.new){graphics.off()}
    dev.new(width=14,height=4)
    click.window<-dev.cur()
    
    dev.new(width=10,height=6) 
    lines.window<-dev.cur()
    
    dev.new(width=8, height=8)
    view.window<-dev.cur()
    
    if(is.null(cells)){cnames <- names(dat$t.dat[,-1])}
    else{cnames<-cells}

    lines.flag <- 0
    cell.i <- 1
    g.names<-NULL
    click.i <- 1
    #group.names<-NULL

    while(click.i!=7)
    {
        cell.pick <- cnames[cell.i]
        #dev.set(dev.list()[1])
        dev.set(which=click.window)
        p1 <- PeakFunc6(dat,cell.pick, t.type=t.type,yvar=yvar, info=info, bcex=bcex)
        p1.par<-par()
        if(lines.flag==1){
            #dev.set(dev.list()[2])
            dev.set(which=lines.window)
            LinesEvery.5(dat,g.names,plot.new=F, img=img, t.type=t.type, col="black")
            lines.flag <- 0
        }

        if(lines.flag==2){
            #dev.set(dev.list()[3])
            dev.set(which=view.window)
            cell.view(dat,cell=cell.pick, img=img,cols="red",plot.new=F,cell.name=T, zoom=FALSE)
            lines.flag <- 0
        }
        
        if(lines.flag==0){
            #dev.set(dev.list()[1]) 
            dev.set(which=click.window)
        }		
        
        #title(sub=paste("Group ",group.i," n=",g.num," Cell ",cell.i,sep=""))
        xs<- rep(par("usr")[1]-xinch(.5), 7)
        ys<-seq(par("usr")[4],by=-yinch(.2), length.out=7)
        points(x=xs,y=ys,pch=16)
        text(x=xs,y=ys,labels=c("Cell +","Cell -","Veiw","Stack","yvar","Select Trace","off"),pos=2,cex=.5)
        
        ## How many cells are you looking at
        text(par("usr")[1], par("usr")[4]+yinch(.3),paste(cell.i, ":",length(cnames)))
        click.i <- identify(x=xs,y=ys,n=1,plot=F)
        
        if(click.i==1)
        {cell.i <- cell.i + 1;if(cell.i>length(cnames)){cell.i<-1};lines.flag<-0}
        if(click.i==2)
        {cell.i <- cell.i - 1;if(cell.i<1){cell.i<-length(cnames)};lines.flag<-0}
        if(click.i==3)
        {lines.flag<-2}
        if(click.i==4)
        {g.names<-union(g.names,cnames[cell.i]);lines.flag<-1}
        if(click.i==5){
            if(yvar){yvar<-FALSE}else{yvar<-TRUE}
        }
        if(click.i==6){
            t.type<-select.list(names(dat))
        }
        if(click.i==7){
        #graphics.off()
        dev.off(which=click.window)
        dev.off(which=lines.window)
        dev.off(which=view.window)}
    }

    print(g.names)
}

bp.selector<-function(dat,cell=NULL,cells=NULL,dat.name=NULL,plot.new=T,save.bp=F,view.cells=F, env=NULL, localize=T){
    print(environment())
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
    
    
    ###################################################################
    #This region needs significant work to improve to all data aspects
    ###################################################################
    
    ## Selcet eith Area or Peak Height
    type<-select.list(c("Peak Height", "Area"), multiple=F, title="Parameter?")
    if(type=="Peak Height"){type<-".max"
    }else{type<-".tot"}


    
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
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    par(bg="gray90")
    #define the open window
    peakfunc.window<-dev.cur()
    #plot the trace specified at the beigning
    PeakFunc7(dat,cell, lmain="  ",bcex=1.5, info=F)
    title(expression("RED"* phantom("/BLUE")), col.main="red")
    title(expression(phantom("RED/")*"BLUE"),col.main="blue")
    title(expression(phantom("RED")*"/"*phantom("BLUE")),col.main="black")

    # add point to the plot to define buttons
    ys<-rep(par("usr")[3],length(levs))

    points(levs.mean, ys, pch=16, cex=2)
    #label each point with levs text
    #text(levs.mean,ys,labels=names(levs.mean),pos=c(1,3),cex=1, srt=90)
    
    ###Selecting Control Windows
    bringToTop(-1)
    cat("Choose one or more window regions for the denominator in the equations,
    
    Amplification-or-block = active.window / control.window
    
    CLICK LARGE BLACK DOTS to select
    Click stop in the top left.
    
    "
    )
    #Select windows to define numerator
    controlwindows <- identify(x=levs.mean,y=ys,labels="X",plot=T, col="red", cex=1.5)
    #collect the names of what you have selected
    controlwindows<- levs[controlwindows]
    
    ###Selecting Active Windows
    bringToTop(-1)
    cat("Choose one or more window regions for the numerator in the equations,
    
    Amplification-or-block = active.window / control.window
    
    Click stop in the top left, and then STOP LOCATOR from the drop down
    
    "
    )
    #change focus back to the peakwindow for active window selection
    dev.set(peakfunc.window)
    activewindows <- identify(x=levs.mean,y=ys,labels="X",plot=T, col="blue",cex=1.5)
    activewindows<-levs[activewindows]
    
    #now if there are multiple control windows selected, 
    if(length(controlwindows)>1){
        #create the name for scp collumn lookup
        controlmax<-paste(controlwindows, type, sep="")
        #add that name to scp, and do a row mean
        controlmaxmean<-data.frame(rowMeans(dat$scp[controlmax]))
    }else{
        controlmax<-paste(controlwindows, type, sep="")
        controlmaxmean<-dat$scp[controlmax]
    }
    #same as above!
    if(length(activewindows)>1){
        activemax<-paste(activewindows, type, sep="")
        activemaxmean<-data.frame(rowMeans(dat$scp[activemax]))
    }else{
        activemax<-paste(activewindows, type, sep="")
        activemaxmean<-dat$scp[activemax]
    }

    max_amp_mean<-activemaxmean/controlmaxmean
    max_amp_mean[,2]<-seq(from=1,to=dim(max_amp_mean)[1],by=1)
    max_amp_mean_cells<-data.frame(activemaxmean[cells,])/data.frame(controlmaxmean[cells,])
    max_amp_mean_cells[,2]<-seq(from=1,to=dim(max_amp_mean_cells)[1],by=1)
    
    row.names(max_amp_mean_cells)<-cells
    
    # Calculate percent change and select for cells
    print("Would you like to save this statistic to scp?")
    save_stat_op<-select.list(c("yes","no"), title="Save Stat?")
    if(save_stat_op=="yes"){
        print("Enter the name of the statistic to be added to scp")
        stat.name<-scan(n=1, what='character')
        dat$scp[stat.name]<-max_amp_mean
        assign(dat.name,dat, envir=env)
    }
    
    density_ct_plotter(dat, cells, NULL, max_amp_mean[1],xlim_top=3,xlim_bottom=0, overlay=T,dense_sep=F,plot_new=F)
    #dev.new(width=15, height=5)
    #par(mfrow=c(1,2), bty="l")
    
    #hist(max_amp_mean_cells[,1], breaks=length(max_amp_mean_cells[,1])/2, xlim=c(0,2))
 
    boxplot(max_amp_mean_cells[,1], outline=F, ylim=c(0,2),width=10, lty=1, lwd=2, main=paste(activewindows,"Amplification Cutoff"), ylab="Active.Max/Control.Max", horizontal=T)
    text(
        jitter(
            rep(
                1, 
                length(max_amp_mean_cells[,1])
            ),10
        )~max_amp_mean_cells[,1], 
        labels=row.names(max_amp_mean_cells), 
        cex=.5, 
        col=rgb(1,1,1,4, maxColorValue=10)
    )#,ylim=c(0,2.5), add=T, vertical=T, method="jitter", jitter=.2)
    
    #170131 adding 2 point localization
    if(localize=="T"){
        selector<-select.list(c("one", "two"), title="Left side first!")
        
        if(selector=="one"){loc<-locator(n=1, type="p", pch=15, col="red")}
        if(selector=="two"){loc<-locator(n=2, type="p", pch=15, col="red")}

        abline(v=loc$x,col="red")
        
        if(length(loc$x)==1){
            x.names<-row.names(which(max_amp_mean[1]>loc$x, arr.ind=T, useNames=T))
            x.names<-row.names(max_amp_mean[order(max_amp_mean[,1],decreasing=T),])
        }

        if(length(loc$x)==2){
            x.names<-which(max_amp_mean[1]>loc$x[1] & max_amp_mean[1]<loc$x[2], arr.ind=T,useNames=T)
            x.names<-row.names(max_amp_mean[order(max_amp_mean[,1],decreasing=T),])
        }
    }else{
        x.names<-row.names(max_amp_mean_cells[order(max_amp_mean_cells[,1],decreasing=T),])
        print(x.names)
    }
    if(view.cells){
        continue<-select.list(c("Yes", "No"), multiple=F, title="View Selected Cells?")
    }else{continue<-"No"}
    
    if(continue=="Yes"){
        print(length(x.names))
        #graphics.off()
        real.cells<-tcd(dat, x.names,dat.name=dat.name)
        return(real.cells)
    }else{
        return(x.names)
    }
}

bp.selector.advanced<-function(dat,cell=NULL,cells=NULL,dat.name=NULL,plot.new=T,save.bp=F,view.cells=F, env=NULL, localize=T){
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
    
    
    ###################################################################
    #This region needs significant work to improve to all data aspects
    ###################################################################
    
        
    ## Selcet eith Area or Peak Height
    type<-select.list(c("Peak Height", "Area"), multiple=F, title="Parameter?")
    if(type=="Peak Height"){type<-".max"
    }else{type<-".tot"}

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
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    par(bg="gray90")
    #define the open window
    peakfunc.window<-dev.cur()
    #plot the trace specified at the beigning
    PeakFunc7(dat,cell, lmain="  ",bcex=1.5, info=F)
    
    title("(After-Before)/(After+Before)")
    
    # add point to the plot to define buttons
    ys<-rep(par("usr")[3],length(levs))

    points(levs.mean, ys, pch=16, cex=2)
    #label each point with levs text
    #text(levs.mean,ys,labels=names(levs.mean),pos=c(1,3),cex=1, srt=90)
    continue<-"yes"
    
    while(continue=="yes"){

        ###Selecting Control Windows
        bringToTop(-1)
        cat("
        Choose the Pulse Following the compound of interest. 
        
        This is the AFTER pulse
        
        CLICK LARGE BLACK DOTS to select
        You Only Get one shot.
        
        CLICK ANY KEY TO CONTINUE
        "
        )
        
        scan(n=1)
        #Select windows to define numerator
        afterwindows <- identify(x=levs.mean,y=ys,labels="X",plot=T, col="red", cex=1.5,n=1)
        #collect the names of what you have selected
        afterwindows<- levs[afterwindows]
        
        ###Selecting Active Windows
        bringToTop(-1)
        cat("
        ###############################################
        Choose the Pulse Before the compound of interest.
        
        This is the BEFORE pulse
        You only get one click.
        
        PRESS ANY KEY TO CONTINUE
        "
        )
        scan(n=1)
        #change focus back to the peakwindow for active window selection
        dev.set(peakfunc.window)
        beforewindows <- identify(x=levs.mean,y=ys,labels="X",plot=T, col="blue",cex=1.5,n=1)
        beforewindows<-levs[beforewindows]
        
        #Find the scp collumn to provide the best stat
        aftermax<-paste(afterwindows, type, sep="")
        aftermaxmean<-dat$scp[aftermax]
        
        beforemax<-paste(beforewindows, type, sep="")
        beforemaxmean<-dat$scp[beforemax]
        
        max_amp_mean<-(aftermaxmean-beforemaxmean)/(aftermaxmean+beforemaxmean)
        max_amp_mean[,2]<-seq(from=1,to=dim(max_amp_mean)[1],by=1)
        
        max_amp_mean_cells<-(
        ( data.frame(aftermaxmean[cells,])-data.frame(beforemaxmean[cells,]) )/
        ( data.frame(aftermaxmean[cells,])+data.frame(beforemaxmean[cells,]) ) )
        
        max_amp_mean_cells[,2]<-seq(from=1,to=dim(max_amp_mean_cells)[1],by=1)
        
        row.names(max_amp_mean_cells)<-cells
        
        # Calculate percent change and select for cells
        cat("Would you like to save this statistic to scp? \n")
        save_stat_op<-select.list(c("yes","no"), title="Save Stat?")
        if(save_stat_op=="yes"){
            cat("Enter the name of the statistic to be added to scp \n")
            stat.name<-scan(n=1, what='character')
            dat$scp[stat.name]<-max_amp_mean
            assign(dat.name,dat, envir=env)
        }
        cat("
        Make another stat?")
        continue<-select.list(c("yes","no"))
    }

    density_ct_plotter(dat, cells, NULL, max_amp_mean[1],xlim_bottom=-1,xlim_top=1, overlay=T,dense_sep=F,plot_new=F)
    #dev.new(width=15, height=5)
    #par(mfrow=c(1,2), bty="l")
    
    #hist(max_amp_mean_cells[,1], breaks=length(max_amp_mean_cells[,1])/2, xlim=c(0,2))
 
    boxplot(max_amp_mean_cells[,1], outline=F, ylim=c(-1,1),width=10, lty=1, lwd=2, main="Amplification Cutoff", ylab="Active.Max/Control.Max", horizontal=T)
    text(
        jitter(
            rep(
                1, 
                length(max_amp_mean_cells[,1])
            ),10
        )~max_amp_mean_cells[,1], 
        labels=row.names(max_amp_mean_cells), 
        cex=.5, 
        col=rgb(1,1,1,4, maxColorValue=10)
    )#,ylim=c(0,2.5), add=T, vertical=T, method="jitter", jitter=.2)
    
    #170131 adding 2 point localization
    if(localize){
        selector<-select.list(c("one", "two"), title="Left side first!")
        
        if(selector=="one"){loc<-locator(n=1, type="p", pch=15, col="red")}
        if(selector=="two"){loc<-locator(n=2, type="p", pch=15, col="red")}

        abline(v=loc$x,col="red")
        
        if(length(loc$x)==1){
            #now we need to 
            #1.select cells based on the first click on the boxplot graphic
            
            x.names<-row.names(which(max_amp_mean_cells[1]>loc$x, arr.ind=T, useNames=T))
            #now that we have fouynd the cells which respond in these ways we will 
            #sort the dataframe based on these stats
            new_max_amp_mean_cells<-max_amp_mean_cells[x.names,]
            x.names<-row.names(new_max_amp_mean_cells[order(new_max_amp_mean_cells[1], decreasing=T),])
        }

        if(length(loc$x)==2){
            x.names<-row.names(which(max_amp_mean_cells[1]>loc$x[1] & max_amp_mean_cells[1]<loc$x[2], arr.ind=T,useNames=T))
            new_max_amp_mean_cells<-max_amp_mean_cells[x.names,]
            x.names<-row.names(new_max_amp_mean_cells[order(new_max_amp_mean_cells[1], decreasing=T),])
        }
    }else{
        x.names<-row.names(max_amp_mean_cells[order(max_amp_mean_cells[,1],decreasing=T),])
        print(x.names)
    }
    if(view.cells){
        continue<-select.list(c("Yes", "No"), multiple=F, title="View Selected Cells?")
    }else{continue<-"No"}
    
    if(continue=="Yes"){
        print(length(x.names))
        #graphics.off()
        real.cells<-tcd(dat, x.names,dat.name=dat.name)
        return(real.cells)
    }else{
        return(x.names)
    }
}

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

#Repairs score from levs only
# Uses peakfunc5
bin.repair<-function(dat, n.names=NULL){
    if(is.null(n.names)){n.names<-names(dat$t.dat[,-1])}
    cell.i<-1
    cell<-n.names[cell.i]
    levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    levs.mean<-sort(tapply(dat$t.dat[,"Time"], as.factor(dat$w.dat$wr1), mean))
    levs<-setdiff(names(levs.mean),"")
    levs.mean<-levs.mean[levs]
    xs <- c(levs.mean,rep(dat$t.dat[50,"Time"],4))
    ys<-c(rep(1.4, length(levs.mean)),1.2, 1.1, 1.0, 0.9)
    dev.new(width=14, height=5)
    dev.set(dev.list()[1])
    PeakFunc6(dat,cell, Plotit.both=F)
    linesflag<-0
    click.i<-0
    
    while(click.i!=length(levs.mean)+4){
        points(x=xs,y=ys,pch=16)
        text(x=xs,y=c(rep(1.4, length(levs.mean)),1.2,1.1,1.0,0.9),labels=c(names(levs.mean),"Cell +","Cell -","drop","off"),pos=2,cex=.5)
        click.i <- identify(x=xs,y=ys,n=1,plot=T)
        cell<-n.names[cell.i]
        if(click.i<=length(levs.mean)){
            if(dat$bin[cell, levs[click.i]]==1){dat$bin[cell, levs[click.i]]=0;dat$bin[cell,"drop"]=0;linesflag<-0}
            else{dat$bin[cell, levs[click.i]]=1;dat$bin[cell,"drop"]=0;linesflag<-0}
            dev.set(dev.list()[1]);PeakFunc6(dat, cell, Plotit.both=F)
        }
        
        if(click.i==length(levs.mean)+1){cell.i <- cell.i + 1;if(cell.i>length(n.names)){cell.i<-1};linesflag<-1}
        if(click.i==length(levs.mean)+2){cell.i <- cell.i - 1;if(cell.i<1){cell.i<-length(n.names)};linesflag<-1}
        if(click.i==length(levs.mean)+3){dat$bin[cell, "drop"]=1;dev.set(dev.list()[1]);PeakFunc6(dat, cell, Plotit.both=F)} #dat$bin[cell,levs]=0;
        if(linesflag==1){PeakFunc6(dat, n.names[cell.i], Plotit.both=F)}

    }
            graphics.off()
            return(dat$bin)
}
    
### Repairs GFP and TRITC score from label bin
# uses peakfunc5	
bin.repair.2<-function(dat, n.names=NULL){
    if(is.null(n.names)){n.names<-names(dat$t.dat[,-1])}
    cell.i<-1
    cell<-n.names[cell.i]
    levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    levs.mean<-sort(tapply(dat$t.dat[,"Time"], as.factor(dat$w.dat$wr1), mean))
    levs<-setdiff(names(levs.mean),"")
    levs.mean<-levs.mean[levs]
    
    rep(max(dat$t.dat[,"Time"])-(max(dat$t.dat[,"Time"])*1.095),4)
    xs <- c(levs.mean,c(max(dat$t.dat[,1])*1.09, max(dat$t.dat[,1])*1.19),rep(max(dat$t.dat[,"Time"])-(max(dat$t.dat[,"Time"])*1.095),4))
    ys<-c(rep(1.8, length(levs.mean)+2),1.2, 1.0, 0.8, 0.6)
    dev.new(width=14, height=5)
    dev.set(dev.list()[1])
    PeakFunc5(dat,cell, Plotit.both=T)
    linesflag<-0
    click.i<-0
    
    while(click.i!=length(levs.mean)+2+4){
        points(x=xs,y=ys,pch=16)
        text(x=xs,y=ys,labels=c(names(levs.mean),"mean.gfp", "tritc","Cell +","Cell -","drop","off"),pos=3,cex=.5)
        click.i <- identify(x=xs,y=ys,n=1,plot=T)
        cell<-n.names[cell.i]
        
        if(click.i<=length(levs.mean)){
            if(dat$bin[cell, levs[click.i]]==1){dat$bin[cell, levs[click.i]]=0;dat$bin[cell,"drop"]=0;linesflag<-0}
            else{dat$bin[cell, levs[click.i]]=1;dat$bin[cell,"drop"]=0;linesflag<-0}
            dev.set(dev.list()[1]);PeakFunc5(dat, cell, Plotit.both=T)
        }
        
        if(click.i==length(levs.mean)+1){
            if(dat$bin[cell, "gfp.bin"]==1){dat$bin[cell, "gfp.bin"]=0;dat$bin[cell,"drop"]=0;linesflag<-0}
            else{dat$bin[cell, "gfp.bin"]=1;dat$bin[cell,"drop"]=0;linesflag<-0}
            dev.set(dev.list()[1]);PeakFunc5(dat, cell, Plotit.both=T)
        }

        if(click.i==length(levs.mean)+2){
            if(dat$bin[cell, "tritc.bin"]==1){dat$bin[cell, "tritc.bin"]=0;dat$bin[cell,"drop"]=0;linesflag<-0}
            else{dat$bin[cell, "tritc.bin"]=1;dat$bin[cell,"drop"]=0;linesflag<-0}
            dev.set(dev.list()[1]);PeakFunc5(dat, cell, Plotit.both=T)
        }
        
        if(click.i==length(levs.mean)+3){cell.i <- cell.i + 1;if(cell.i>length(n.names)){cell.i<-1};linesflag<-1}
        if(click.i==length(levs.mean)+4){cell.i <- cell.i - 1;if(cell.i<1){cell.i<-length(n.names)};linesflag<-1}
        if(click.i==length(levs.mean)+5){dat$bin[cell, "drop"]=1;dev.set(dev.list()[1]);PeakFunc5(dat, cell, Plotit.both=T)} #dat$bin[cell,levs]=0;
        if(linesflag==1){PeakFunc5(dat, n.names[cell.i], Plotit.both=T)}
        }

        graphics.off()
        neuron.response<-select.list(levs, title="What defines Neurons?", multiple=T)
        neurons<-cellz(dat$bin,neuron.response, 1)
        drop<-cellz(dat$bin, "drop", 1)
        neurons<-setdiff(neurons,drop)
        pf<-apply(dat$bin[,c("gfp.bin", "tritc.bin")],1,paste, collapse="")
        dat$bin["lab.pf"]<-as.factor(pf)
        lab.groups<-unique(dat$bin$lab.pf)
        
        cells<-list()
        for(i in lab.groups){
            x.names<-cellz(dat$bin[neurons,], "lab.pf", i)
            cells[[i]]<-x.names
        }
        
        glia.response<-select.list(c(levs, "none"), title="What defines glia?", multiple=T)
        if(glia.response!="none"){
            drop<-cellz(dat$bin, "drop", 1)
            glia<-cellz(dat$bin,glia.response, 1)
            glia<-setdiff(glia,drop)
            cells[["000"]]<-setdiff(glia, neurons)
        } 
        else {cells[["000"]]<-setdiff(row.names(dat$c.dat), neurons)}
        dat$cells<-cells
        return(dat)

}

bin.rep.cells<-function(dat){
    
    cells<-dat$cells
    
    for(i in 1:length(cells)){
        dat<-bin.repair.2(dat, cells[[i]])
        }
    return(dat)
}

# Creates Binary socring for labeling
# Input RD list, and # of cells to observe for sampling
# Outuput bin dataframe with added intensity scoring
label.bin<-function(dat, cells=10){
    rand.names<-attributes(sample(dat$c.dat$id))$levels
    n.names<-rand.names[1:cells]

    cell.i<-1
    dev.new(width=15, height=3)
    yes.green<-vector()
    no.green<-vector()
    yes.red<-vector()
    no.red<-vector()

    for(i in 1:length(n.names)){
        
        par(mfrow=c(1,5))
        multi.pic.zoom(dat, n.names[i], dat$img1, plot.new=F)
        multi.pic.zoom(dat, n.names[i], dat$img2, plot.new=F)
        multi.pic.zoom(dat, n.names[i], dat$img3, plot.new=F)
        multi.pic.zoom(dat, n.names[i], dat$img4, plot.new=F)


        par(mar=c(0,0,0,0))
        xloc<-c(2,2,2,2)
        yloc<-c(3.5,2.5,1.5,0.5)
        loc<-cbind(xloc, yloc)
        plot(loc,xlim=c(0,4), pch=15, ylim=c(0,4), xaxt="n", yaxt="n", cex=1.5)
        text(loc, c("+GFP","+TRITC", "+GFP & +TRITC","No Label") ,pos=4, cex=1.5)
        click.i<-identify(loc, n=1, plot=T)	
        
        if(click.i==1){yes.green[i]<-dat$c.dat[n.names[i],"mean.gfp"];no.red[i]<-dat$c.dat[n.names[i],"mean.tritc"]}
        if(click.i==2){yes.red[i]<-dat$c.dat[n.names[i],"mean.tritc"];no.green[i]<-dat$c.dat[n.names[i],"mean.gfp"]}
        if(click.i==3){yes.red[i]<-dat$c.dat[n.names[i],"mean.tritc"];yes.green[i]<-dat$c.dat[n.names[i],"mean.gfp"]}
        if(click.i==4){no.red[i]<-dat$c.dat[n.names[i],"mean.tritc"];no.green[i]<-dat$c.dat[n.names[i],"mean.gfp"]}
    }
    graphics.off()
    
    if(length(yes.green)>=1){yes.green<-setdiff(yes.green,c("NA",NA))}
    if(length(no.green)>=1){no.green<-setdiff(no.green,c("NA",NA))}
    if(length(yes.red)>=1){yes.red<-setdiff(yes.red,c("NA",NA))}
    if(length(no.red)>=1){no.red<-setdiff(no.red,c("NA",NA))}

    dat$bin["gfp.bin"]<-0
    dat$bin["tritc.bin"]<-0

    if(length(yes.green)>=1){green.names<-row.names(dat$c.dat)[dat$c.dat$mean.gfp>min(yes.green)]}
    if(length(yes.red)>=1){red.names<-row.names(dat$c.dat)[dat$c.dat$mean.tritc>min(yes.red)]}

    if(length(yes.green)>=1){dat$bin[green.names,"gfp.bin"]<-1}
    if(length(yes.red)>=1){dat$bin[red.names,"tritc.bin"]<-1}
    
    print(paste("Green Cells : ",min(yes.green)))
    print(paste("Red Cells : ",min(yes.red)))
    print(paste("No label Green : ",max(no.green),"No label Red", max(no.red)))
    
    pf<-apply(dat$bin[,c("gfp.bin", "tritc.bin")],1,paste, collapse="")
    dat$bin["lab.pf"]<-as.factor(pf)

    return(dat$bin)
}	

##############################################################################################
# Cell Group Review
##############################################################################################
#Group summarry
#generate pdfs with line graphs
#table of means and frequencies for all c.dat
#THIS MUST BE CLEANED UP 040314
GroupSummary <- function(dat,snr,c.dat,wr,levs,groups,pref="Group"){
    g.levs <- unique(groups)
    for(i in g.levs)
    {
        cnames <- names(groups[groups==i])
        pdf.name <- paste(pref,i,".pdf",sep="")
        lmain <- paste(pref,i,sep="")
        LinesEvery(dat,snr,cnames,wr,levs,lmain,pdf.name)
        dev.off()
    }
    res.tab <- data.frame(mean=apply(c.dat[names(groups),],2,mean))
    res.tab["sd"] <- apply(c.dat[names(groups),],2,sd)
    for(i in g.levs)
    {
        cnames <- names(groups[groups==i])
        res.tab[paste(pref,i,".mean",sep="")] <- apply(c.dat[cnames,],2,mean)
        res.tab[paste(pref,i,".sd",sep="")] <- apply(c.dat[cnames,],2,sd)
    }
    tab.name <- paste(pref,".table.csv",sep="")
    write.csv(res.tab,file=tab.name)
    #lines figure similar to boxplot
    ## tmp <- scale(c.dat[names(groups),],center=T,scale=T)
    ## tmp.mn <- data.frame(t(apply(tmp,2,function(x){tapply(x,as.factor(groups),mean)})))
    ## tmp.sd <- data.frame(t(apply(tmp,2,function(x){tapply(x,as.factor(groups),sd)})))
    ## tmp.se <- t(t(tmp.sd)/sqrt(summary(as.factor(groups))))
    ## ylim <- c(min(tmp.mn)-2,max(tmp.mn))
    ## miny <- min(ylim)+1
    ## dev.new()
    ## par(xaxt="n",mar=c(2,4,4,2))
    ## plot(seq(1,nrow(tmp.mn)),tmp.mn[,1],ylim=ylim,xlim=c(0,(nrow(tmp.mn)+1)),type="n",ylab="Normalized Mean +- SE",xaxt="n")
    ## cols <- rainbow(ncol(tmp.mn),start=.3)
    ## names(cols) <- names(tmp.mn)
    ## nudge <- 0
    ## ## for(i in names(tmp.mn))
    ## {
    ##     xseq <- seq(1,nrow(tmp.mn))
    ##     rect(nudge+seq(1,nrow(tmp.mn))-.05,tmp.mn[,i]-tmp.se[,i],nudge+seq(1,nrow(tmp.mn))+.05,tmp.mn[,i]+tmp.se[,i],col=cols[i],border=NA)
    ##     points(nudge+seq(1,nrow(tmp.mn)),tmp.mn[,i],pch=16,col=cols[i],lwd=2,type="b")        
    ##     nudge <- nudge+.1
    ## }
    ## text(rep(nrow(tmp.mn),ncol(tmp.mn)),tmp.mn[nrow(tmp.mn),],paste(pref,names(tmp.mn),sep=""),cex=.8,col=cols,pos=4)
    ## text(seq(1,nrow(tmp.mn))+.25,miny,names(c.dat),srt=90,pos=3)
    c.mn <- data.frame(t(apply(c.dat,2,function(x){tapply(x,as.factor(groups),mean)})))
    c.sd <- data.frame(t(apply(c.dat,2,function(x){tapply(x,as.factor(groups),sd)})))
    c.se <- t(t(c.sd)/sqrt(summary(as.factor(groups))))
    return(list(mean=c.mn,sd=c.sd,se=c.se))   
}

# Fucntion plotting cell locations, barplots of labeled intensities, stacked traces, and
# single traces of all scored groups.
# Needs work on click funcitons, and recognition of NULL intensities from experiemnts
GroupReview.2 <- function(dat,bp.plot=T,shws=2,phws=20,wr.i=2,bl.meth="TopHat"){
    library(cluster)
    graphics.off()
    #peakfunc window= dev.list()[1]
    windows(width=8,height=4, xpos=0, ypos=0)    
    #linefunc window= dev.list()[2]
    windows(width=8,height=5, xpos=0, ypos=360) 
    #bpfunc window= dev.list()[3]
    windows(width=5,height=4, xpos=800, ypos=420) 
    #cell.locate window= dev.list[4]
    windows(width=12,height=12, xpos=820, ypos=0) 
    #gui window= dev.list[5]
    windows(width=2,height=2, xpos=1400, ypos=620) 
    # Plotting all traces ontop of each other	
    # Could attempt something like a LinesEvery function 
    # Should replace linesfunce with linesevery.2.  If there are more than 15 cells
    # then i need to plot traces like tracechase. Needs window plotting.
    # shade windows according to scoring


    #Cell locate still needs to be able to move through images.  \
    # New data set will have 4-5 images
    # Also, this function needs have all click features available, including click 
    # cells for peakfunc selections

    # Create a table with binary groups as rows
    # collumn 1=total cells in group
    # collumn 2=group number
    total.cell<-sort(summary(dat$c.dat[,"pf"]))
    group.sum<-cbind(total.cell, seq(1,length(total.cell), by=1))
    as.table(group.sum)
    colnames(group.sum)<-c("c.tot", "g.num")
    #make clust (which is the definition of clusters) be equal to the group numbers
    #in group.sum
    #clust<-group.sum[,"g.num"]
    levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    pf<-apply(dat$bin[,levs],1,paste,collapse="")
    pf.sum<-summary(as.factor(pf),maxsum=500)
    pf.sum<-pf.sum[order(pf.sum,decreasing=T)]
    pf.ord<-pf.sum
    pf.ord[]<-seq(1,length(pf.sum))
    dat$c.dat["pf"]<-as.factor(pf)
    dat$c.dat["pf.sum"]<-pf.sum[pf]
    dat$c.dat["pf.ord"]<-pf.ord[pf]
    clust<-dat$c.dat[,"pf.ord"]
    clust.name <- unique(clust)
    
    levs<-setdiff(unique(as.character(dat$w.dat[,2])),"")

    dev.set(dev.list()[5])
    par(mar=c(0,0,0,0))
    plot(2,2, pch=NA)
    points(x=c(rep(1.75,5),rep(2.5,6)),y=c(2.5,2.25,2.0,1.75,1.5,2.5,2.25,2,1.75,1.5,1.25),pch=16)
    text(x=c(rep(1.75,5),rep(2.5,6)),y=c(2.5,2.25,2.0,1.75,1.5,2.5,2.25,2,1.75,1.5,1.25),
    labels=c("Group +","Group -","Cell +","Cell -","Done", "Image 1", "Image 2", "Image 3", "Zoom", "+ Pulse", "- Pulse"),pos=2,cex=.8)
    
    
    img<-dat$img1
    #an intiator of the linesfunc if lines.flag=1
    lines.flag <- 1
    #this is a list of all cell names
    g.names <- names(dat$t.dat[,-1])
    #highest group #
    pam.k <- max(clust)
    #initial group and cell to start analysis
    group.i <- 1
    cell.i <- 1
    peak.i<-1
    
    # define first click
    click.i <- 1
    while(click.i)
    {
    #initiate the single peak plot, but only if the group exists
        g.num <- sum(clust==group.i)
        if(g.num > 0)
        {
        #first group defined above, but can be further defined below
        group.names <- g.names[clust==group.i]
        #first cell is defined above, but can be further defined below
        cell.pick <- group.names[cell.i]
                # Intial setting for image changer
        
        #move to next plot and start peakfunc2
        #p1 <- PeakFunc2(dat,cell.pick,shws=shws,phws=phws,Plotit=T,wr=dat$w.dat[,wr.i],SNR.lim=2,bl.meth=bl.meth)
        }
        #start boxplot of color intensities
        if(lines.flag==1){
            dev.set(dev.list()[1]);PeakFunc5(dat, cell.pick)
            dev.set(dev.list()[2]);if(length(group.names)>10){LinesStack(dat, group.names, plot.new=F)}else{LinesEvery.2(dat, group.names,plot.new=F)}
            dev.set(dev.list()[3]);bpfunc(dat,group.names)
            dev.set(dev.list()[4]);cell.zoom.2048(dat,img, group.names, plot.new=F);lines.flag <- 0
        }
        
        dev.set(dev.list()[5])
    
        click.i <- identify(x=c(rep(1.75,5),rep(2.5,6)),y=c(2.5,2.25,2.0,1.75,1.5,2.5,2.25,2,1.75,1.5,1.25),n=1,plot=F)
        # syntax for first click on peakfumc2. if click group+ group.i+1
        
        if(click.i==1)
        {group.i <- group.i + 1;if(group.i > pam.k){group.i <- 1};cell.i<-1;lines.flag <- 1}
        
        if(click.i==2)
        {group.i <- group.i - 1;if(group.i < 1){group.i <- pam.k};cell.i<-1;lines.flag <- 1}
        
        if(click.i==3){
            cell.i <- cell.i + 1
            if(cell.i > g.num){cell.i <- 1}
            dev.set(dev.list()[1]);PeakFunc5(dat, cell.pick)
        }
        
        if(click.i==4){
            cell.i <- cell.i - 1
            if(cell.i < 1){cell.i <- g.num}
            dev.set(dev.list()[1]);PeakFunc5(dat, cell.pick)
        }
        
        if(click.i==5)
        {graphics.off();stop()}
        
        if(click.i==6){if(!is.null(dat$img1)){img<-dat$img1};lines.flag<-1}
        
        if(click.i==7){if(!is.null(dat$img2)){img<-dat$img2};lines.flag<-1}
        
        if(click.i==8){if(!is.null(dat$img3)){img<-dat$img3};dev.set(dev.list()[4]);lines.flag<-1}
        
        if(click.i==9){}#cell.pick<-group.names[cell.i];cell.locate(cell.pick, zoom=5)}
        
        if(click.i==10){peak.i<-peak.i+1;group.names<-row.names(dat$bin)[dat$bin[,levs[peak.i]]==1];lines.flag <- 1}
        
        if(click.i==11){group.names<-p.names[[peak.i-1]]; cell.i<-1 ;lines.flag<-1}
    }
    dev.off()
}
##############################################################################################
# Trace Searching
##############################################################################################
#topdown parsing of all traces
TraceChase <- function(dat,blc=NULL,levs=NULL,x.names=NULL,scale=T){
    library(cluster)
    if(is.null(blc)){
    if(is.element("blc",names(dat))){blc <- dat$blc}
    else
    {tmp.pcp <- ProcConstPharm(dat);blc <- tmp.pcp$blc}}
    if(is.null(levs))
    {
        levs <- unique(dat$w.dat[,"wr1"])
        levs <- select.list(levs,multiple=T,title="Select Regions for clustering")
    }
    dmat <- t(scale(blc[is.element(dat$w.dat[,"wr1"],levs),-1],scale=scale,center=scale))
    a.names <- names(blc)[-1]
    if(!is.null(x.names)){a.names <- intersect(x.names,names(blc))}
    done=FALSE
    while(!done)
    {
        if(length(a.names) < 21)
        {
            x.names <- TraceSelect(dat,a.names,dat$w.dat[,"wr1"],levs, "Final Select")
            done=TRUE
        }
        else
        {
            #pam20 <- pam(dmat[a.names,],k=20)
            clmb20 <- ClimbTree(dmat[a.names,],k=20)
            lmain <- paste("Select Traces (all or none to end) n=",length(a.names))
            #x.names <- SmashSelect(blc[c("Time",a.names)],pam20$clustering,row.names(pam20$medoids),dat$w.dat[,"wr1"],levs,lmain=lmain)				
            x.names <- SmashSelect(blc[c("Time",a.names)],clmb20,names(clmb20)[match(1:length(unique(clmb20)),clmb20)],dat$w.dat[,"wr1"],levs,lmain=lmain)							
            if(length(a.names)==length(x.names)){done = TRUE}
            if(length(x.names)==0){done= TRUE}
            a.names <- x.names
        }
    
    }
    return(x.names)	
}

#given a set of traces (or trace seqments)
#calculate the distances and group into K groups
#by height of tree cutting. One of the K groups will
#be a catch-all for all small groups
ClimbTree <- function(x,k=20){
    tabstat <- function(x){return(list(mean=mean(x),length=length(x),median=median(x),sd=sd(x),gt5c=sum(x>5)))}
    library(cluster)
    d1 <- dist(x)
    h1 <- hclust(d1)
    q1 <- quantile(h1$height,probs=1:10/10)
    clust <- cutree(h1,h=q1[5])
    clust.tab <- table(clust)
    clust.tab <- clust.tab[order(clust.tab,decreasing=T)]	
    new.num <- clust.tab
    new.num[] <- seq(1,length(new.num))
    clust[] <- new.num[as.character(clust)]
    clust.tab <- table(clust)
    if(length(clust.tab) > k)
    {
        in.grp <- names(clust.tab[1:(k-1)])
        out.grp <- setdiff(names(clust.tab),in.grp)
        clust[is.element(clust,out.grp)] <- k	
    }
    return(clust)
    # clust.stat <- data.frame(tabstat(clust.tab))
    # for(i in 2:length(q1))
    # {
        # clust <- cutree(h1,h=q1[i])
        # clust.tab <- table(clust)
        # clust.stat[i,] <- tabstat(clust.tab)
    # }
    # return(clust.stat)
    
}

#smash select plot the smashes and return the selected.
#all data in t.dat is ploted (1st col must be time)
#m.names are taken to be the medoids of the clusters
SmashSelect <- function(t.dat,clust,m.names,wr,levs=NULL,lmain=""){	
    rtag <- table(clust)
    names(rtag) <- m.names[order(clust[m.names])]
    sf <- 1
    gcol <- rgb(10,10,10,alpha=120,max=255)
    #gcol <- "grey"
    x <- t.dat[,-1]
    xm <- apply(x,2,max)
    xn <- scale(x,center=F,scale=xm)
    for(i in 1:nrow(xn)){xn[i,] <- xn[i,]+clust}
    
    library(RColorBrewer)
    lwds <- 2
    
    xseq <- t.dat[,1]
    cols <-brewer.pal(8,"Dark2")
    cols <- rep(cols,ceiling(length(m.names)/length(cols)))
    cols <- cols[1:length(m.names)]
    dev.new(width=14,height=9)
    op <- par(yaxt="n",bty="n",mar=c(4,0,2,1),cex=1)
    plot(xseq,xn[,m.names[1]],ylim=c((min(xn)-2),max(xn)),xlab="Time (min)",ylab="Ratio with shift",main=lmain,type="n", xaxt="n")
    axis(1, at=seq(0, length(t.dat[,1]), 5))
    apply(xn,2,lines,x=xseq,col=gcol,lwd=2)
    hbc <- 1
    if(length(wr) > 0)
    {
        if(is.null(levs)){levs <- setdiff(unique(wr),"")}
        x1s <- tapply(xseq,as.factor(wr),min)[levs]
        x2s <- tapply(xseq,as.factor(wr),max)[levs]
        y1s <- rep(-.3,length(x1s))
        y2s <- rep(hbc+.2,length(x1s))
        #rect(x1s,y1s,x2s,y2s,col="lightgrey")
        text(xseq[match(levs,wr)],rep(c(.2,-.2),length.out=length(levs)),levs,pos=4,offset=0,cex=1)
    }
    x.sel <- NULL
    xs <-c(rep(0,length(m.names)),c(.1,.1,.1))
    ys <- xn[1,m.names]
    ys <- as.vector(c(ys,c(sf*.9,0,-sf*.9)))
    #    xs[(length(xs)-2):length(xs)] <- c(0,5,10)
    p.names <- c(rep(" ",length(m.names)),"ALL","NONE","FINISH")
    done.n <- length(p.names)
    none.i <- done.n-1
    all.i <- none.i-1
    p.cols <- c(cols,c("black","black","black"))
    for(i in 1:length(m.names))
    {
          #lines(xseq,xn[,m.names[i]],col="black",lwd=lwds*.5)
        lines(xseq,xn[,m.names[i]],col=cols[i],lwd=lwds)
    }
    text(x=rep(max(xseq),length(m.names)),y=xn[nrow(xn),m.names],cex=.9,rtag,pos=4,col=p.cols)
    text(x=xs,y=ys,labels=p.names,pos=2,cex=.7,col=p.cols)
    points(x=xs,y=ys,pch=16,col=p.cols,cex=1.5)
    click.i <- 1    
    while(click.i != done.n)
    {
        click.i <- identify(xs,ys,n=1,plot=F)
        if(click.i < (length(m.names)+1) & click.i > 0)
        {
            i <- click.i
            if(is.element(i,x.sel))
            {
                lines(xseq,xn[,m.names[i]],col=cols[i],lwd=lwds)
                x.sel <- setdiff(x.sel,i)
            }
                else
                {
                lines(xseq,xn[,m.names[i]],col="black",lwd=lwds)
                x.sel <- union(x.sel,i)
            }
        }
        if(click.i == none.i)
        {
            x.sel <- NULL
            for(i in 1:length(m.names))
            {
                lines(xseq,xn[,m.names[i]],col=cols[i],lwd=lwds)
            }
        }
        if(click.i == all.i)	
        {
            x.sel <- seq(1,length(m.names))
            for(i in 1:length(m.names))
            {
                lines(xseq,xn[,m.names[i]],col="black",lwd=lwds)
            }
            
        }
    }
    c.sel <- clust[m.names[x.sel]]
    x.ret <- names(clust[is.element(clust,c.sel)])
    dev.off()
    return(x.ret)
}

#this simply finds the traces in t.dat that are similar to targs
#note this is "complete" similarity other options may be
#"average" and "best"
GetCloser <- function(t.dat,targs,k=20){
    x.names <- setdiff(names(t.dat),targs)
    ct <- cor(t.dat[,x.names],t.dat[,targs])
    x.max <- apply(ct,1,min)
    y.names <- x.names[order(x.max,decreasing=T)[1:k]]
    return(y.names)	
}

#this is a bit raw still
#Given a set of traces (t.dat) and a list of targets (targs)
#identify the 20 most similar traces using wr and the select levs.
#allow the user to select from those to add to the master list.
SimilarSelect <- function(t.dat,targs,wr,levs=NULL){
    plot(t.dat[,1],t.dat[,targs[1]],type="n",ylim=c(min(t.dat[-1]),(length(targs)+50)*.2))
    sf <- 0
    for(i in targs){lines(t.dat[,1],t.dat[,i]+sf);sf<-sf+.2}
    a.names <- setdiff(names(t.dat)[-1],targs)
    rjct <- rep(0,length(a.names))
    names(rjct) <- a.names
    done=FALSE
    tps <- seq(1:nrow(t.dat))
    if(!is.null(levs)){tps <- tps[is.element(wr,levs)]}
    while(!done)
    {	
        if(sum(rjct==0) < 21)
        {done=TRUE}
        else
        {
            x.names <- GetCloser(t.dat[tps,c(a.names[rjct==0],targs)],targs)
            rjct[x.names] <- 1
            y.names <- TraceSelect(t.dat,,x.names,wr)
    
            if(length(y.names)==0){done=TRUE}
            else
            {
                targs <- c(targs,y.names)
                for(i in y.names){lines(t.dat[,1],t.dat[,i]+sf);sf<-sf+.2}
            }
        }
    }
    return(targs)
    #plot targs and allow user to
    #paint region of interest if you can do this it makes a very good window adjust function.
    #find matches within t.dat
    #show matches in trace select allow user to choose.
    #merge all selected and return that list.
    
    
}
##############################################################################################
# Interactive Image analysis
##############################################################################################
#Function to Automatically Rename images from the time i used a different naming scheme	
ImageRenamer<-function(dat){	

    image.names<-grep("img",names(dat))	
    for(i in 1:length(image.names)){
        names(dat)[image.names[i]]<-paste("img",i,sep="")
    }
    
return(dat)
}

#How to create a function to select and add images to the specified experiment files 
#from the current working directory.
ImageFiller<-function(dat){
    require(png)
    potential.images<-list.files(pattern='png')
    print(potential.images)
    bringToTop(-1)
    print("######################")
    print("These are the images you have the option of selecting")
    print("Now select the images to fill in for image 1 to 8")
    
    for(i in 1:8){
        image.to.add<-select.list(list.files(pattern='png'),title=paste('img',i,sep=''))
        
        if(image.to.add==""){dat[[paste('img',i,sep='')]]<-NULL
        }else{
            dat[[paste('img',i,sep='')]]<-readPNG(image.to.add)
        }
    }
    return(dat)
}

ImageFillerv2 <- function(dat, img_name_vec){

    if( is.null(img_name_vec) ){
        img_name_vec<-c(
            "bf.gfp.tritc.start.png",
            "gfp.tritc.start.ci.ltl.rs.png",
            "tritc.start.ci.ltl.rs.png",
            "gfp.start.ci.ltl.rs.png",
            "bf.start.lab.png",
            "fura2.png",
            "fura2.divide.start.png",
            "roi.img.png")
        image_question <- T
    }else{
        image_question <- F
    }
    # Add images
    img_list<-list()
    for( j in 1:length(img_name_vec) ){
        dat[[ paste0("img",j) ]] <- tryCatch(readPNG(img_name_vec[j]), error=function(e)NULL)
    }
    if(image_question == T){
        cat('\nThese are the images I have attempted to load for you\nIf any are NULL, and want to add different images say yes to the \nnext question. You will be asked to select a png image for each loaction.\n\n')
        cat(img_name_vec, sep='\n')
        cat(str(img_list))

        cat('\nDO YOU WANT DIFFERENT IMAGES[y,n]?\n')
        img_reselect <- scan(n=1,what='character')
        if( img_reselect=='y' ){
            cat("\nAlright buddy I am going to give you options if you don't\nwant any image there just go ahead and press 0\n\n")
            png_imgs <- list.files(pattern='png')
            for( j in 1:8 ){
                cat("\nWhat do you want for image ", j, '\n')
                selection <- menu(png_imgs)
                if(selection==0){
                    dat[[paste0("img",j)]] <- NULL
                }else{
                    dat[[paste0("img",j)]] <- readPNG(png_imgs[selection])
                }
                cat('\nI have added ', png_imgs[selection],' to position ',j,'\n')
            }
            }
    }
    return(dat)    
}




# Fucntion locates single cell or groups of cells on plot.  
# Needs more optional assignments
cell.veiw.2048<-function(dat, img=NULL, cell=NULL, cells=NULL, cols=NULL,lmain="", bcex=.5, plot.new=T, cell.name=T){
    if(plot.new){dev.new()}
    require(png)
    require(zoom)
    par(mar=c(0,0,1,0))
    cells.x<-dat$c.dat[cells,"center.x"]
    cells.y<-dat$c.dat[cells,"center.y"]
    cell.x<-dat$c.dat[cell,"center.x"]
    cell.y<-dat$c.dat[cell,"center.y"]

    if(is.null(img)){img<-dat$img1}
    else{img<-img}
    if(is.null(cols)){cols="white"}
    else{cols=cols}

    plot(0, 0, xlim=c(0,2048),ylim=c(2048,0), main=lmain,xaxs="i", yaxs="i", xlab="Pixels", ylab="Pixels")
    rasterImage(img, 0, 2048, 2048, 0)

    points(cell.x, cell.y, col=cols, pch=4, cex=1)
    text(cell.x, cell.y, labels=cell, col=cols, pos=2, cex=1)
    
    points(cells.x, cells.y, col="white", pch=4, cex=bcex)
    text(cells.x, cells.y, labels=dat$c.dat[cells,1], col="white", pch=4, pos=2, cex=bcex)
}

cell.view<-function(dat, cell=NULL,img=NULL,  zoom=TRUE, cols=NULL,lmain="", bcex=.8, labs=T, plot.new=T, cell.name=T){
    if(plot.new){dev.new()}
    require(png)
    par(mar=c(0,0,1,0))
    x<-dat$c.dat[,"center.x"]
    y<-dat$c.dat[,"center.y"]
    cell.x<-dat$c.dat[cell,"center.x"]
    cell.y<-dat$c.dat[cell,"center.y"]

    if(is.null(img)){img<-dat$img1}
    else{img<-img}
    if(is.null(cols)){cols="white"}
    else{cols=cols}
    img.dimy<-dim(img)[1]
    img.dimx<-dim(img)[2]

    plot(0, 0, xlim=c(0,img.dimx),ylim=c(img.dimy,0), main=lmain,xaxs="i", yaxs="i", xlab="Pixels", ylab="Pixels")
    rasterImage(img, 0, img.dimy, img.dimx, 0)

    if(labs){
        if(!is.null(cell)){
            points(cell.x, cell.y, col=cols, pch=0, cex=2)
            text(cell.x, cell.y, labels=cell, col=cols, pos=2, cex=bcex)
        }
        else{
            points(x, y, col=cols, pch=4, cex=2)
            text(x, y, labels=dat$c.dat[,1], col=cols, pch=0, pos=2, cex=bcex)
        }
    }

    if(zoom==TRUE & length(cell)>1){
        cell.1<-row.names(dat$c.dat[order(dat$c.dat$center.x),])
        cell<-intersect(cell,cell.1)
        multi.pic.zoom(dat,cell,img) 
    }
    }

    cell.zoom.640.480<-function(dat, img=NULL, cell=NULL, zoom=NULL, cols=NULL, labs=T, plot.new=T, cell.name=T)
    {
    if(plot.new){dev.new()}
    require(png)
    require(zoom)
    par(mar=c(0,0,0,0))
    x<-dat$c.dat[,"center.x"]
    y<-dat$c.dat[,"center.y"]
    cell.x<-dat$c.dat[cell,"center.x"]
    cell.y<-dat$c.dat[cell,"center.y"]

    if(is.null(img)){img<-dat$img1}
    else{img<-img}
    if(is.null(cols)){cols="white"}
    else{cols=cols}

    plot(0, 0, xlim=c(0,640),ylim=c(480,0), xaxs="i", yaxs="i", xlab="Pixels", ylab="Pixels")
    rasterImage(img, 0, 480, 640, 0)

    if(labs){
    if(!is.null(cell)){
        points(cell.x, cell.y, col=cols )
        text(cell.x, cell.y, labels=cell, col=cols, pos=2, cex=.8)
        }
    else{
        points(x, y, col=cols)
        text(x, y, labels=dat$c.dat[,1], col=cols, pos=2, cex=.5)
    }}

    if(!is.null(zoom)){
    zoomplot.zoom(x=cell.x, y=cell.y, fact=zoom)
    }
    else{zm()}
}

XYtrace.640.480 <- function(dat, img=NULL, cols=NULL, labs=T){
    x.coor<-grep("\\.x",names(dat$c.dat), value=T, ignore.case=T)
    y.coor<-grep("\\.y",names(dat$c.dat), value=T, ignore.case=T)
    area<-grep("area",names(dat$c.dat), value=T, ignore.case=T)
    
    lab1<-grep("cgrp",names(dat$c.dat), value=T, ignore.case=T)
    if(length(lab1)==0){lab1<-grep("gfp.1",names(dat$c.dat), value=T, ignore.case=T)}
    
    lab1.1<-grep("cgrp",names(dat$c.dat), value=T, ignore.case=T)
    if(length(lab1.1)==0){lab1.1<-grep("gfp.2",names(dat$c.dat), value=T, ignore.case=T)}
    
    lab2<-grep("ib4",names(dat$c.dat), value=T, ignore.case=T)
    if(length(lab2)==0){lab2<-grep("tritc",names(dat$c.dat), value=T, ignore.case=T)}
    
    cell.coor<-dat$c.dat[,c(x.coor, y.coor)]
        
    # select the names of the collumns containing coordinates
    levs <- unique(dat$w.dat[,"wr1"])
    levs<-setdiff(levs, "")
    if(labs==TRUE){
    if(is.null(cols)){cols="grey5"} else{cols=cols}}
    pch=16
    
    dev.new(height=4,width=12)
    dev.new(width=10, height=8)
    dev.new(height=8,width=12)
    lmain<-"XY ROI"

    dev.set(dev.list()[2])
    par(mar=c(0,0,0,0))
    plot(0, 0, xlim=c(0,640),ylim=c(480,0),xaxs="i", yaxs="i",col=cols,pch=".")
    
    if(is.null(img)){img<-dat$img1}
    if(!is.null(img)){rasterImage(img, 0, 480, 640, 0);points(cell.coor[,1],cell.coor[,2],col=cols,pch=0,cex=2.4)}
    else{
    points(cell.coor[,1],cell.coor[,2], col=cols, cex=dat$c.dat[,area]/200)
    points(cell.coor[,1],cell.coor[,2],col=cols, pch=4)}
    
    

    i <- identify(cell.coor[,1],cell.coor[,2],n=1,plot=F, col=NA, tolerance=0.05)
    i.names<-row.names(dat$c.dat)[i]
    while(length(i) > 0)
    {	#selected name of cell
        s.names <- row.names(dat$c.dat)[i]
        dev.set(dev.list()[1])
        PeakFunc2(dat,s.names,3,30,TRUE,,lmain=lmain)
        dev.set(dev.list()[2])
        # If a cell is selected, that has already been selected, 
        # then remove that cell from the list
        if(length(intersect(i.names,s.names))==1){
        i.names<-setdiff(i.names,s.names)
        points(cell.coor[s.names,1],cell.coor[s.names,2],col="grey90",pch=0,cex=2.4)
        points(cell.coor[i.names,1],cell.coor[i.names,2],col="red",pch=0,cex=2.4)}
        # If it han't been selected, then add it to the list
        else{i.names<-union(i.names,s.names)
        points(cell.coor[i.names,1],cell.coor[i.names,2],col="red",pch=0,cex=2.4)}
        
        if(length(i.names)>=2){dev.set(dev.list()[3]);LinesEvery.2(dat,m.names=i.names, plot.new=F)}		
        
        dev.set(dev.list()[2])
        i <- identify(cell.coor[,1],cell.coor[,2],labels=dat$c.dat[,1],n=1,plot=T, pch=0,col="grey90", tolerance=0.05)
    }
    dev.off()
    graphics.off()
    return(dat$c.dat[i.names,1])
       
}

# Function allows for selection and deselection of cells to build stacked traces
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

XYtrace.2<-function(dat, cells=NULL, img=NULL, cols=NULL, zoom=T, labs=T, yvar=F, zf=40, t.type=NULL, sf=1,plot.labs=T){
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

# View Individual cell picture
multi.pic.zoom<-function(dat, m.names, img, labs=T,plot.new=T, zf=20){
    col.row<-ceiling(sqrt(length(m.names)))
    
    if(plot.new)
    {
        dev.new()
        par(mfrow=c(col.row, col.row))
        par(mar=c(0,0,0,0))
    }
    else{
        par(mfrow=c(col.row, col.row))
        par(mar=c(0,0,0,0))
    }	
    
    m.names<-rev(m.names)
    for(i in 1:length(m.names)){		

        img.dim<-dim(img)
        x<-dat$c.dat[m.names[i],"center.x"]
        left<-x-zf
        if(left<=20){left=0; right=zf}
        right<-x+zf
        if(right>=img.dim-zf){left=img.dim-zf;right=img.dim}
                    
        y<-dat$c.dat[m.names[i],"center.y"]
        top<-y-zf
        if(top<=20){top=0; bottom=zf}
        bottom<-y+zf
        if(bottom>=img.dim-zf){top=img.dim-zf;bottom=img.dim}

        par(xpd=TRUE)
        xleft<-0
        xright<-20
        ytop<-0
        ybottom<-20
        plot(c(xright, xleft), c(ytop, ybottom), ylim=c(20,0) ,xaxs="i", yaxs="i", axes=F)
        rasterImage(img[top:bottom,left:right,],xleft,ybottom,xright,ytop)
        text(4,1.5, m.names[i], col="white", cex=.8)
        box(lty = 1, col = "white",lwd=2)
        text(16.5, 2, labels=dat$c.dat[m.names[i], "area"], col="white")

        if(labs){
            points(x=10,y=10, type="p", pch=3, cex=2,col="white")
            text(16.5, 2, labels=dat$c.dat[m.names[i], "ROI.Area"], col="white")
            text(16.5, 3.5, labels=dat$c.dat[m.names[i], "mean.gfp.1"], col="green")
            text(16.5, 3.5, labels=dat$c.dat[m.names[i], "mean.gfp"], col="green")
            text(16.5, 3.5, labels=dat$c.dat[m.names[i], "CGRP"], col="green")
            text(16.5, 5, labels=dat$c.dat[m.names[i], "mean.tritc"], col="red")
            text(16.5, 5, labels=dat$c.dat[m.names[i], "IB4"], col="red")
            text(16.5, 6.5, labels=dat$c.dat[m.names[i], "mean.dapi"], col="blue")
        }
    }
        
}

# View Individual cell picture creates a png image
# must assgin multi.pic.zoom to a variable name
# For use in linesEvery.4
multi.pic.zoom.2<-function(dat, m.names, img, labs=F, zf=NULL, cols=NULL){
    if(is.null(cols)){cols<-rep("white", length(m.names))}else{cols<-cols}
    col.row<-ceiling(sqrt(length(m.names)))
    #png("tmp.png",width=6,height=6,units="in",res=72,bg="transparent", type="cairo")
    #dev.new()
    png('tmp.png', res=70)	
            par(mfrow=c(col.row, col.row))
            par(mar=c(0,0,0,0))
        #else{par(mar=c(0,0,0,0))}	
        m.names<-rev(m.names)
        img.dim<-as.numeric(dim(img)[1])
        
    for(i in 1:length(m.names)){		

            if(is.null(zf)){zf<-20}else{zf<-zf}
            #zf<-20
            x<-dat$c.dat[m.names[i],"center.x"]
            left<-x-zf
            right<-x+zf
            if(left<=zf){left=0; right=zf}
            
            if(right>=img.dim){left=img.dim-zf;right=img.dim
            }else{right=right}

            y<-dat$c.dat[m.names[i],"center.y"]
            top<-y-zf
            if(top<=zf){top=0; bottom=zf}
            bottom<-y+zf
            if(bottom>=img.dim-zf*2){top=img.dim-zf;bottom=img.dim}

            par(xpd=TRUE)
            xleft<-0
            xright<-20
            ytop<-0
            ybottom<-20
            plot(c(xright, xleft), c(ytop, ybottom), ylim=c(20,0) ,xaxs="i", yaxs="i", axes=F)
            
            if(length(dim(img))>2){rasterImage(img[top:bottom,left:right,],xleft,ytop,xright,ybottom)
            }else{rasterImage(img[top:bottom,left:right],xleft,ytop,xright,ybottom)}
            points(x=10,y=10, type="p", pch=3, cex=2,col="white")
            text(4,1.5, labels=m.names[i], col="white", cex=1.3)
            #text(4,1.5, labels=m.names[i], col=cols[i], cex=1.2)
            box(lty = 1, col = "white",lwd=2)
        if(labs){
            #label.names<-c("ROI.Area", "mean.gfp.1", "CGRP", "IB4")
            label.names<-c("area","mean.gfp","mean.tritc", "mean.dapi")
            label.y.location<-c(2,3.5,5,6.5)
            label.cols<-c("white", "green", "red", "blue")
            for(j in 1:length(label.names)){
                text(16.5, label.y.location[j], labels=tryCatch(round(dat$c.dat[m.names[i],label.names[j]],digits=5),error=function(e) NULL), col=label.cols[j])
            }
        }
    }	
        dev.off()
        tmp.png <- readPNG("tmp.png")
        unlink("tmp.png")
        return(tmp.png)			
}

#multipiczoom
multi.pic.zoom.3<-function(dat, m.names, img, labs=T,plot.new=T, zf=20){
    col.row<-ceiling(sqrt(length(m.names)))
    
    if(plot.new){
        dev.new()
        par(mfrow=c(col.row, col.row))
        par(mar=c(0,0,0,0))
    }
    else{par(mar=c(0,0,0,0))}	
        m.names<-rev(m.names)
    for(i in 1:length(m.names)){		
        x<-dat$c.dat[m.names[i],"center.x"]
        left<-x-zf
        if(left<=20){left=0; right=zf}
        right<-x+zf
        if(right>=1004){left=2048-zf;right=2048}
                    
        y<-dat$c.dat[m.names[i],"center.y"]
        top<-y-zf
        if(top<=20){top=0; bottom=zf}
        bottom<-y+zf
        if(bottom>=1004){top=2048-zf;bottom=2048}

        par(xpd=TRUE)
        xleft<-0
        xright<-20
        ytop<-0
        ybottom<-20
        plot(c(xright, xleft), c(ytop, ybottom), ylim=c(20,0) ,xaxs="i", yaxs="i", axes=F)
        rasterImage(img[top:bottom,left:right,],xleft,ytop,xright,ybottom)
        points(x=10,y=10, type="p", pch=3, cex=2,col="white")
        box(lty = 1, col = "white",lwd=2)
        if(labs){
            text(4,1.5, labels=m.names[i], col="white", cex=1.2)
            text(16.5, 2, labels=dat$c.dat[m.names[i], "area"], col="white")
            text(16.5, 2, labels=dat$c.dat[m.names[i], "ROI.Area"], col="white")
            text(16.5, 3.5, labels=dat$c.dat[m.names[i], "mean.gfp.1"], col="green")
            text(16.5, 3.5, labels=dat$c.dat[m.names[i], "mean.gfp"], col="green")
            text(16.5, 3.5, labels=dat$c.dat[m.names[i], "CGRP"], col="green")
            text(16.5, 5, labels=dat$c.dat[m.names[i], "mean.tritc"], col="red")
            text(16.5, 5, labels=dat$c.dat[m.names[i], "IB4"], col="red")
            text(16.5, 6.5, labels=dat$c.dat[m.names[i], "mean.dapi"], col="blue")
        }
    }
        
}

image.selector<-function(tmp.rd, multi=T){
    img.names<-grep(names(tmp.rd),pattern="img", value=T)
    
    null.images<-vector()
    for(i in 1:length(img.names)){null.images[i]<-!is.null(tmp.rd[[img.names[i]]])}
    img.logical<-cbind(img.names,null.images)
    real.imgs<-which(img.logical[,2]=="TRUE")
    
    img.names<-img.logical[real.imgs, 1]
    
    
    dev.new(width=ceiling(sqrt(length(img.names)))*4, height=ceiling(sqrt(length(img.names)))*4)
    img.sel<-dev.cur()
    par(mfrow=c(ceiling(sqrt(length(img.names))),ceiling(sqrt(length(img.names)))))
    
    for(i in 1:length(img.names)){
        par(mar=c(0,0,0,0))
        img<-tmp.rd[[img.names[[i]]]]
        img.dim.y<-dim(img)[1]
        img.dim.x<-dim(img)[2]	
        
        top<-img.dim.y*.25
        bottom<-img.dim.y*.75
        left<-img.dim.x*.25
        right<-img.dim.x*.75
        
        plot(0, 0, 
            xlim=c(img.dim.x*.4, img.dim.x*.6),
            ylim=c(img.dim.y*.4,img.dim.y*.6),xaxt="n", yaxt="n",pch="."
            )
        rasterImage(img[top:bottom,left:right,], 0, img.dim.y, img.dim.x, 0)
        text(img.dim.x*.45,img.dim.y*.45,labels=paste(i), cex=2, col="white")
    }
    img<-select.list(img.names, title="Images", multiple=multi)
    dev.off(img.sel)
    return(img)
}

PointTrace <- function(lookat,png=F,col=rep("black",nrow(lookat)),pch=16,cex=1,lmain="PointTrace",ylim=c(-2,2),x.trt=NULL,y.trt=NULL,wr="wr1",t.names=NULL){
    if(!is.null(x.trt)){lookat["x"] <- lookat[,x.trt]}
    if(!is.null(y.trt)){lookat["y"] <- lookat[,y.trt]}	
    dev.new(height=4,width=14)
    rr.dev <- dev.cur()
    dev.new(height=4,width=4)
    
    plot(lookat[,"x"],lookat[,"y"],col=col,pch=pch,cex=cex,main=lmain,xlab=x.trt,ylab=y.trt, ylim=ylim)
    ret.list <- NULL
    i <- identify(lookat[,"x"],lookat[,"y"],n=1,plot=F)
    my.dev <- dev.cur()
    while(length(i) > 0)
    {
        x.names <- lookat[i,"trace.id"]
        #points(lookat[i,"x"],lookat[i,"y"],pch=8,cex=.5)
        rn.i <- row.names(lookat)[i]
        tmp <- get(lookat[i,"rd.name"])
        levs <- unique(tmp$w.dat[,"wr1"])
        lmain <- paste(i,lookat[i,"rd.name"])
        #LinesEvery(tmp$t.dat,,x.names,tmp$w.dat[,"wr1"],levs,lmain=lmain)
        dev.set(which=rr.dev)
        PeakFunc5(tmp,x.names,lmain=lookat[i,"rd.name"])
        if(!is.null(t.names)){mtext(paste(t.names,tmp$c.dat[x.names,t.names],collapse=":"))}
        if(png==TRUE)
        {
            f.name <- paste(lookat[i,"rd.name"],lookat[i,"trace.id"],"png",sep=".")
            png(f.name,heigh=600,width=1200)
            PeakFunc2(tmp$t.dat,x.names,3,30,TRUE,tmp$w.dat[,wr],lmain=lookat[i,"rd.name"])
            dev.off()
        }

        dev.set(which=my.dev)
        if(is.element(rn.i,ret.list))
            {points(lookat[i,"x"],lookat[i,"y"],col=col[i],pch=pch,cex=cex);ret.list <- setdiff(ret.list,rn.i)}		
        else
            {points(lookat[i,"x"],lookat[i,"y"],col="red",pch=pch,cex=cex);ret.list <- union(rn.i,ret.list)}		
        i <- identify(lookat[,"x"],lookat[,"y"],n=1,plot=F)
    }
    return(ret.list)

}

PointTrace.2 <- function(lookat,png=F,col=rep("black",nrow(lookat)),pch=16,cex=1,lmain="PointTrace",x.trt=NULL,y.trt=NULL,wr="wr1",t.names=NULL){
    graphics.off()
    if(!is.null(x.trt)){lookat["x"] <- lookat[,x.trt]}else{lookat["x"] <- lookat[,select.list(names(lookat))]}
    if(!is.null(y.trt)){lookat["y"] <- lookat[,y.trt]}else{lookat["y"] <- lookat[,select.list(names(lookat))]}
    dev.new(height=4,width=14)
    rr.dev <- dev.cur()
    dev.new(height=4,width=4)
    
    plot(lookat[,"x"],lookat[,"y"],pch=pch,cex=cex,main=lmain,xlab=x.trt,ylab=y.trt, col="white")
    text(lookat[,"x"],lookat[,"y"],labels=lookat$trace.id)
    ret.list <- NULL
    i <- identify(lookat[,"x"],lookat[,"y"],n=1,plot=F)
    my.dev <- dev.cur()
    while(length(i) > 0)
    {
        x.names <- lookat[i,"trace.id"]
        #points(lookat[i,"x"],lookat[i,"y"],pch=8,cex=.5)
        rn.i <- row.names(lookat)[i]
        tmp <- get(lookat[i,"rd.name"])
        levs <- unique(tmp$w.dat[,"wr1"])
        lmain <- paste(i,lookat[i,"rd.name"])
        #LinesEvery(tmp$t.dat,,x.names,tmp$w.dat[,"wr1"],levs,lmain=lmain)
        dev.set(which=rr.dev)
        #PeakFunc5(tmp,x.names,lmain=lookat[i,"rd.name"])
        rtpcr.multi.plotter(tmp,x.names,pdf=F,bcex=1, melt.plot=T, plot.new=F)
        if(!is.null(t.names)){mtext(paste(t.names,tmp$c.dat[x.names,t.names],collapse=":"))}
        if(png==TRUE)
        {
            f.name <- paste(lookat[i,"rd.name"],lookat[i,"trace.id"],"png",sep=".")
            png(f.name,heigh=600,width=1200)
            PeakFunc2(tmp$t.dat,x.names,3,30,TRUE,tmp$w.dat[,wr],lmain=lookat[i,"rd.name"])
            dev.off()
        }

        dev.set(which=my.dev)
        if(is.element(rn.i,ret.list))
            {points(lookat[i,"x"],lookat[i,"y"],col=col[i],pch=pch,cex=cex);ret.list <- setdiff(ret.list,rn.i)}		
        else
            {points(lookat[i,"x"],lookat[i,"y"],col="red",pch=pch,cex=cex);ret.list <- union(rn.i,ret.list)}		
        i <- identify(lookat[,"x"],lookat[,"y"],n=1,plot=F)
    }
    return(ret.list)

}

##############################################################################################
# Multi Experiment Analysis
##############################################################################################
#calculate means and sems for all cnames of dat
#divided by the levels of fac.name
#make a bargraph of these
MeanSemGraph <- function(dat,cnames,fac.name,t.cols=NULL,ylab=NULL,main.lab=NULL,x.labs=NULL,bt=.1,lgc="topleft",ylim=NULL){
    semfunc <- function(x)
    {
        n <- sum(!is.na(x))
        if(n < 3){return(NA)}
        return(sd(x,na.rm=T)/sqrt(n))
    }
    x <- as.factor(dat[,fac.name])
    x.levs <- levels(x)
    if(1/length(x.levs) < bt){bt <- 1/length(x.levs)}
    sem.levs <- paste(x.levs,"sem",sep=".")
    x.res <- data.frame(apply(dat[x==x.levs[1],cnames,drop=F],2,mean,na.rm=T))
    for(i in x.levs)
    {
        x.res[i] <- apply(dat[x==i,cnames,drop=F],2,mean,na.rm=T)
        x.res[paste(i,"sem",sep=".")] <- apply(dat[x==i,cnames,drop=F],2,semfunc)
    }
    xlim <- c(1,length(cnames)+length(x.levs)*bt)
    if(is.null(ylim)){ylim <- c(-.02,max(x.res[,x.levs]+x.res[,sem.levs]*2)*1.2)}
    
    if(is.null(t.cols)){t.cols <- rainbow(length(x.levs));names(t.cols) <- x.levs}
    plot(x.res[,x.levs[1]],xlim=xlim,ylim=ylim,type="n",xaxt="n",xlab="",ylab=ylab,main=main.lab)
    for(i in 1:length(x.levs))
    {
        x1 <- seq(1,length(cnames))+(i-1)*bt
        y1 <- x.res[,x.levs[i]]
        rect(x1,rep(0,length(x1)),x1+bt,y1,col=t.cols[x.levs[i]])
    }
    for(i in 1:length(x.levs))
    {
        x1 <- seq(1,length(cnames))+(i-1)*bt+(bt)/2
        y1 <- x.res[,x.levs[i]] + x.res[,sem.levs[i]]*2
        y2 <- x.res[,x.levs[i]] - x.res[,sem.levs[i]]*2		
        arrows(x1,y2,x1,y1,angle=90,col="black",length=bt*.25,code=3)
    }
    
    if(is.null(x.labs)){x.labs <- row.names(x.res)}
    text(seq(1,length(cnames)),rep(-.02,length(cnames)),x.labs,pos=4,cex=.8,offset=0)
    legend(lgc,col=t.cols,names(t.cols),pch=15)
    return(x.res[,-1])
}

cells.plotter<-function(dat, tmp.names, subset.n=5,multi=TRUE, pic=TRUE){
    rd.names<-unique(dat$rd.name)
    rd.list<-list()
    
    for(i in 1:length(rd.names)){
        x.names<-row.names(dat[tmp.names,])[dat[tmp.names,"rd.name"]==rd.names[i]]
        x.names<-dat[x.names,"id"]
        x.names<-setdiff(x.names, "NA")
        x.name<-setdiff(x.names, NA)
        rd.list[[i]]<-x.names
        names(rd.list)[i]<-rd.names[i]
    }
    
    if(multi){
        for(i in 1:length(rd.list)){
            tmp<-load(paste(names(rd.list)[i],".rdata",sep=""))
            LinesStack.2(get(tmp), rd.list[[i]], names(rd.list[i]), subset.n=subset.n)
            rm(tmp)		
            }
    }
    if(pic){
        for(i in 1:length(rd.list)){
            LinesEvery.3(get(names(rd.list)[i]), rd.list[[i]],img=get(names(rd.list)[i])$img1, lmain=names(rd.list[i]))
        }
    }
    return(rd.list)
}

bg.plotter<-function(gid.bin, dat, subset.n=5,multi=TRUE, pic=TRUE){
    tmp.names<-row.names(dat)[dat$gid.bin==gid.bin]
    rd.names<-unique(dat$rd.name)
    rd.list<-list()
    
    for(i in 1:length(rd.names)){
        x.names<-row.names(dat[tmp.names,])[dat[tmp.names,"rd.name"]==rd.names[i]]
        x.names<-dat[x.names,"id"]
        x.names<-setdiff(x.names, "NA")
        x.name<-setdiff(x.names, NA)
        rd.list[[i]]<-x.names
        names(rd.list)[i]<-rd.names[i]
    }
    
    if(multi){
        for(i in 1:length(rd.list)){
            LinesStack.2(get(names(rd.list)[i]), rd.list[[i]], names(rd.list[i]), subset.n=subset.n)
        }
    }
    if(pic){
        for(i in 1:length(rd.list)){
            LinesEvery.3(get(names(rd.list)[i]), rd.list[[i]],img=get(names(rd.list)[i])$img1, lmain=names(rd.list[i]))
        }
    }
    return(rd.list)
}

pf.plotter<-function(dat,pf, subset.n=5,multi=TRUE, pic=TRUE){
    tmp.names<-row.names(dat)[dat$pf==pf]
    rd.names<-unique(dat$rd.name)
    rd.list<-list()
    
    for(i in 1:length(rd.names)){
        x.names<-row.names(dat[tmp.names,])[dat[tmp.names,"rd.name"]==rd.names[i]]
        x.names<-dat[x.names,"id"]
        #x.names<-na.exclude(x.names)
        rd.list[[i]]<-x.names
        names(rd.list)[i]<-rd.names[i]
    }
    
    if(multi){
        for(i in 1:length(rd.list)){
            LinesStack.2(get(names(rd.list)[i]), rd.list[[i]], names(rd.list[i]), subset.n=subset.n)
        }
    }
    if(pic){
        for(i in 1:length(rd.list)){
            LinesEvery.3(get(names(rd.list)[i]), rd.list[[i]],img=get(names(rd.list)[i])$img3, lmain=names(rd.list[i]))
        }
    }
    return(rd.list)
}

#Updated with Linesevery3
levs.plotter<-function(dat,levs,levs.no, subset.n=5,multi=F, pic=T, click=F){
    tmp.names<-row.names(dat)[dat[,levs]==1]
    
    rd.names<-unique(dat$rd.name)
    rd.list<-list()
    
    for(i in 1:length(rd.names)){
        x.names<-row.names(dat[tmp.names,])[dat[tmp.names,"rd.name"]==rd.names[i]]
        x.names<-dat[x.names,"id"]
        #x.names<-na.exclude(x.names)
        rd.list[[i]]<-x.names
        names(rd.list)[i]<-rd.names[i]
    }
    
    if(multi){
        for(i in 1:length(rd.list)){
            tmp<-load(paste(names(rd.list)[i],".rdata",sep=""))
            LinesStack(get(tmp), rd.list[[i]], names(rd.list[i]), subset.n=subset.n)
        }
    }
    if(pic){
        for(i in 1:length(rd.list)){
            tmp<-load(paste(names(rd.list)[i],".rdata",sep=""))
            LinesEvery.3(get(tmp), rd.list[[i]],img=get(names(rd.list)[i])$img3, lmain=names(rd.list[i]), pic.plot=F, XY.plot=T)
        }
    }
    if(click){
        for(i in 1:length(rd.list)){
            tmp<-load(paste(names(rd.list)[i],".rdata",sep=""))
            Trace.Click.3(get(tmp), rd.list[[i]],img=get(names(rd.list)[i])$img3, lmain=names(rd.list[i]), pic.plot=F, XY.plot=T)
        }
    }

    rm(list=ls(rd.list))
    return(rd.list)
}

all.plotter<-function(dat, subset.n=5,multi=F, pic=F, click=T){
    tmp.names<-row.names(dat)
    
    rd.names<-unique(dat$rd.name)
    rd.list<-list()
    
    for(i in 1:length(rd.names)){
        x.names<-row.names(dat[tmp.names,])[dat[tmp.names,"rd.name"]==rd.names[i]]
        x.names<-dat[x.names,"id"]
        #x.names<-na.exclude(x.names)
        rd.list[[i]]<-x.names
        names(rd.list)[i]<-rd.names[i]
    }
    
    if(multi){
        for(i in 1:length(rd.list)){
            tmp<-load(paste(names(rd.list)[i],".rdata",sep=""))
            LinesStack(get(tmp), rd.list[[i]], names(rd.list[i]), subset.n=subset.n)
        }
    }
    if(pic){
        for(i in 1:length(rd.list)){
            tmp<-load(paste(names(rd.list)[i],".rdata",sep=""))
            LinesEvery.3(get(tmp), rd.list[[i]],img=get(names(rd.list)[i])$img3, lmain=names(rd.list[i]), pic.plot=F, XY.plot=T)
        }
    }
    selected.cells<-list()
    if(click){
        for(i in 1:length(rd.list)){
            tmp<-load(paste(names(rd.list)[i],".rdata",sep=""))
            selected.cells[[i]]<-Trace.Click.3(get(tmp), rd.list[[i]])
            names(selected.cells)[i]<-rd.names[i]
        }
    }

    rm(list=ls(rd.list))
    
    if(multi==T | pic==T){return(rd.list)}
    if(click==T){return(selected.cells)}
    
}

noci.plotter<-function(dat,type, subset.n=5,multi=F, pic=T){
    tmp.names<-row.names(dat)[dat$noci.type==type]
    tmp.names<-setdiff(tmp.names, "NA")
    rd.names<-unique(dat$rd.name)
    rd.list<-list()
    
    for(i in 1:length(rd.names)){
        x.names<-row.names(dat[tmp.names,])[dat[tmp.names,"rd.name"]==rd.names[i]]
        x.names<-dat[x.names,"id"]
        x.names<-setdiff(x.names, c("NA",NA))
        #x.names<-na.exclude(x.names)
        rd.list[[i]]<-x.names
        names(rd.list)[i]<-rd.names[i]
    }
    
    if(multi){
        for(i in 1:length(rd.list)){
        tmp<-load(paste(names(rd.list)[i],".rdata",sep=""))
        LinesStack.2(get(tmp), rd.list[[i]], names(rd.list[i]), subset.n=subset.n)
        rm(tmp)
        }
    }
    if(pic){
        for(i in 1:length(rd.list)){
            tmp<-load(paste(names(rd.list)[i],".rdata",sep=""))
            LinesEvery.3(get(tmp), rd.list[[i]],img=get(names(rd.list)[i])$img3, lmain=names(rd.list[i]))
            rm(tmp)
        }
    }
    rm(list=ls(rd.list))
    return(rd.list)
}

### Function to select rows based on collumn parameters
# dat can be either a raw RD object or an RD dataframe
# ex dat -or- dat$bin

cellzand<-function(dat,collumn=NULL, parameter=1,cells=NULL){
    
    bob<-list()
     if(is.null(cells)){cells<-dat$c.dat$id}else{cells<-cells}
    if(class(dat)=="list"){
        dat.select<-select.list(names(dat), title="Select DataFrame")
        dat<-dat[[dat.select]]
        if(is.null(cells)){cells<-row.names(dat)}else{cells<-cells}

    }else{
        dat<-dat
        if(is.null(cells)){cells<-row.names(dat)}else{cells<-cells}
    }
    
    if(is.null(collumn)){
        collumn<-select.list(names(dat), multiple=T, title="Select Collumn")
    }else(collumn<-collumn)
    
    if(is.null(parameter)){
        parameter<-1
    }else(parameter<-parameter)
    
    for(i in collumn){
        bob[[i]]<-row.names(dat)[dat[,i]>=parameter]
    }
    
    bob<-Reduce(union, bob)
    #bob<-intersect(bob, cells)
    
    bob<-intersect(bob,cells)
    return(bob)
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


cellzor<-function(dat,collumn=NULL, parameter=1,cells=NULL){
    
    bob<-list()
     if(is.null(cells)){cells<-dat$c.dat$id}else{cells<-cells}
    if(class(dat)=="list"){
        dat.select<-select.list(names(dat), title="Select DataFrame")
        dat<-dat[[dat.select]]
        if(is.null(cells)){cells<-row.names(dat)}else{cells<-cells}

    }else{
        dat<-dat
        if(is.null(cells)){cells<-row.names(dat)}else{cells<-cells}
    }
    
    if(is.null(collumn)){
        collumn<-select.list(names(dat), multiple=T, title="Select Collumn")
    }else(collumn<-collumn)
    
    if(is.null(parameter)){
        parameter<-1
    }else(parameter<-parameter)
    
    for(i in collumn){
        bob[[i]]<-row.names(dat)[dat[,i]>=parameter]
    }
    
    bob<-Reduce(intersect, bob)
    #bob<-intersect(bob, cells)
    
    bob<-intersect(bob,cells)
    return(bob)
}
    
cellz<-function(dat,collumn=NULL, parameter){
    bob<-list()
    if(class(dat)=="list"){
        dat.select<-select.list(names(dat))
        dat<-dat[[dat.select]]
    }else{dat<-dat}
    
    if(is.null(collumn)){
        collumn<-select.list(names(dat), multiple=T)}
    else(collumn<-collumn)
    if(is.null(parameter)){
        parameter<-1}
    else(parameter<-parameter)

    for(i in collumn){
        bob[[i]]<-row.names(dat)[dat[,i]==parameter]
    }
    bob<-Reduce(union, bob)
    

    
    #bob<-intersect(bob,cells)
    return(bob)
}


# function to obtained sorted cell names based off 
# collumn names from c.dat and bin
c.sort<-function(dat,char=NULL){
    char<-select.list(names(dat))
    sort.dir<-select.list(c("TRUE", "FALSE"), title="Decreasing?")
    bob<-row.names(dat[order(dat[,char], decreasing=sort.dir),])
    return(bob)
}

c.sort.2<-function(dat,cells=NULL,collumn=NULL){
    if(class(dat)=="list"){
        dat.selector<-select.list(intersect(names(dat), c("c.dat","bin", "scp")), title="Select DataFrame")
        dat<-dat[[dat.selector]]
    }else{dat<-dat}
    
    if(is.null(collumn)){
        collumn<-select.list(names(dat), title="Select Variable to Sort")
    }else{collumn=collumn}
    
    sort.dir<-select.list(c("TRUE", "FALSE"), title="Decreasing?")
    bob<-row.names(dat[order(dat[,collumn], decreasing=sort.dir),])
    if(!is.null(cells)){bob<-intersect(bob,cells)}
    return(bob)
}

#create a list that uses the names input for the names in the list
named.list<-function(...){
    bob<-list(...)
    names(bob)<-as.character(substitute((...)))[-1]
    return(bob)
    }


    cell.ti<-function(dat, x.names, img=NULL){
    graphics.off()
    dev.new(width=15, height=5)
    PeakFunc5(dat, x.names)
    if(is.null(img)){img<-dat$img1}else{img<-img}
    cell.view(dat,x.names,img)
    multi.pic.zoom(dat, x.names, img, zf=80)
}

#given a list of file names collect and merge all bin scp and c.dat data
CollectMulti <- function(f.names,rd.names=NULL){
    if(is.null(rd.names))
    {
        rd.names <- sub("\\.rdata$","",sub(".*\\/","",f.names),ignore.case=T)
        for(i in f.names){load(i)}
    }
    
    b.names <- NULL
    s.names <- NULL
    cnames <- NULL
    for(i in rd.names)
    {
        tmp <- get(i)
        names(tmp$bin) <- make.names(names(tmp$bin))
        names(tmp$scp) <- make.names(names(tmp$scp))		
        names(tmp$c.dat) <- make.names(names(tmp$c.dat))		
        b.names <- union(b.names,names(tmp$bin))
        s.names <- union(s.names,names(tmp$scp))
        cnames <- union(cnames,names(tmp$c.dat))
    }
    cnames <- setdiff(cnames,b.names)
    s.names <- setdiff(s.names,b.names)
    cnames <- setdiff(cnames,s.names)
    
    tot.names <- c(b.names,s.names,cnames,"rd.name","trace.id")
    ret.dat <- data.frame(matrix(rep(1,length(tot.names)),ncol=length(tot.names)))
    names(ret.dat) <- tot.names
    for(i in rd.names)
    {
        tmp <- get(i)
        names(tmp$bin) <- make.names(names(tmp$bin))
        names(tmp$scp) <- make.names(names(tmp$scp))		
        names(tmp$c.dat) <- make.names(names(tmp$c.dat))			
        ret.tmp <- data.frame(cbind(tmp$bin,tmp$scp,tmp$c.dat))
        ret.tmp["rd.name"] <- i
        ret.tmp["trace.id"] <- row.names(tmp$bin)
    #		ret.dat <- merge(ret.dat,ret.tmp)
        i.names <- setdiff(tot.names,names(ret.tmp))

        for(j in i.names)
        {
            ret.tmp[j] <- NA
        }

        ret.add <- ret.tmp[,tot.names]
        ret.dat <- rbind(ret.dat,ret.add)
    }
    ret.dat <- ret.dat[-1,]
    return(ret.dat)	
}

census.brewer<-function(dat){
    cell.types<-dat$cell.types


    dev.new(width=10, height=5)
    stacked.traces<-dev.cur()
    LinesEvery.5.1(dat, sample(row.names(dat$c.dat)[1:5]), plot.new=F, lmain="WAZZZUPPPP", t.type="t.dat", img=dat$img1)

    print("How Many groups to census?")
    group.number<-scan(n=1, what='numeric')

    print("enter the names of your census groups seperated by '.' (6)")
    census.names<-scan(n=as.numeric(group.number), what='character')
    dev.off(stacked.traces)


    selected.cell.groups<-select.list(names(cell.types), title="Select groups to census", multiple=T)
    print("These are the cells you have chosen")
    print(selected.cell.groups)
    census<-list()

    for(i in 1:length(selected.cell.groups))
    {
        print(selected.cell.groups[i])
        
        if(length(cell.types[[selected.cell.groups[i]]])>1){
            census[[i]]<-tcd(dat, Reduce(union,cell.types[[selected.cell.groups[i]]]))
            names(census[[i]])<-census.names
        }else{
            census[[i]]<-NA
        }	
    }
    print(names(census))
    print(selected.cell.groups)

    names(census)<-selected.cell.groups

    dat$census<-census
    return(dat)
}

census.brewer.2<-function(dat){
    cell.types<-dat$cell.types

    if(is.null(dat$census)){

        dev.new(width=12, height=5)
        stacked.traces<-dev.cur()
        LinesEvery.5.1(dat, sample(row.names(dat$c.dat)[1:5]), plot.new=F, lmain="Reference Plot", t.type="t.dat", img=dat$img1)

        print("How Many groups to census?")
        group.number<-scan(n=1, what='numeric')

        print("enter the names of your census groups seperated by '.' (6)")
        census.names<-scan(n=as.numeric(group.number), what='character')
        dev.off(stacked.traces)
    }else{
    census.names<-names(dat$census[[1]])[!is.na(names(dat$census[[1]]))]
    }


    selected.cell.groups<-select.list(names(cell.types), title="Select groups to census", multiple=T)

    if(is.null(dat$census)){dat$census<-list()}

    for(i in selected.cell.groups)
    {
        print(i)
        
        if(length(cell.types[[i]])>1){
            dat$census[[i]]<-Trace.Click.dev(dat, Reduce(union,cell.types[[i]]))
            names(dat$census[[i]])<-census.names
        }else{
            dat$census[[i]]<-NULL
        }	
    }
    #names(dat$census)<-selected.cell.groups

    #dat$census<-census

    census.df<-dat$bin
    census.df.cn<-names(dat$census[[1]])[!is.na(names(dat$census[[1]]))]#census data frame column names

    for(a in 1:length(dat$census)){
        for(b in 1:length(census.df.cn)){
            census.df[dat$census[[a]][[b]],census.df.cn[b]]<-1
        }
    }

    census.df[is.na(census.df)]<-0#convert all NA to 0
    dat$bin<-census.df
        
    return(dat)
}

census_to_table<-function(dat){
    census.df<-dat$bin
    i<-1
    while( is.na(dat$census[[i]]) ){
        i<-i+1
    }

    (census.df.cn<-names(dat$census[[i]])[!is.na(names(dat$census[[i]]))])#census data frame column names)

    for(i in 1:length(census.df.cn)){
        census.df[census.df.cn[i]]<-0
    }
    
    for(a in 1:length(dat$census)){
        if(!is.na(dat$census[[a]])){
            for(b in 1:length(census.df.cn)){
                census.df[ dat$census[[a]][[b]],census.df.cn[b]]<-1
            }
        }
    }

    #census.df[is.na(census.df)]<-0#convert all NA to 0
    dat$bin<-census.df
        
    return(dat)
}

multi.plotter<-function(dat,cells,levs=NULL, values=NULL){
    dat<-dat
    tmp.rd<-dat
    rd.name<-sub(".Rdata|.rdata", "", list.files(pattern="RD.1"))

    if(is.null(levs)){levs<-select.list(names(dat$bin), multiple=TRUE)
    }else{levs<-levs}

    #if(is.null(values)){values<-select.list(names(dat$c.dat), multiple=T)
    if(is.null(values)){values<-"area"
    }else{values<-values}

    #(img<-image.selector(dat))
    img<-"img1"
    channel<-list(c(1:3))
    cell.length<-length(cells)
    print(cell.length)
    cseries1<-seq(1,2000,10)
    cseries2<-seq(10,2000,10)
    cell.groups<-max(which(cseries1-round(cell.length)<=0, arr.ind=T))

    for(k in 1:cell.groups){
    LinesEvery.5(tmp.rd,cells[cseries1[k]:cseries2[k]], plot.new=T, levs=levs, t.type="mp2", values=values, img=img, channel=channel, bcex=1.2, lmain=paste(rd.name,"LU.ide"))
    }
}

TraceImpute.2 <- function(x,ts=seq(0,(length(x)-1))*.03,xspan=5/length(x),time.step=1/120,plotit=F, lmain=NULL){
    if(is.null(lmain)){lmain="traceImpute.2"}else{lmain=lmain}
    targ <- data.frame(ts=seq(min(ts),max(ts),by=time.step))
    xloe <- loess(x ~ ts,span=xspan)
    xp <- predict(xloe,newdata=targ)
    cols<-rainbow(n=length(ts))
    if(plotit)
    {
        plot(ts,x,pch=16,cex=1.2,xlab="time",ylab="Response", col="black", main=lmain)
        points(targ[,1],xp,type="p", col="red", pch=16, cex=.8)
    }
    return(xp)
}

PulseImputer<-function(tmp,cell,pulse.names=NULL,plot.new=F,sf=8){

    if(is.null(pulse.names)){pulse.names<-intersect(grep("^K",names(tmp$bin),ignore.case=T,value=T),tmp$w.dat[,"wr1"])
    }else{pulse.names<-pulse.names}

    if(plot.new){dev.new(width=1.5*length(pulse.names), height=2)}
    par(mfrow=c(1,length(pulse.names)), mar=c(1,1,1,1))
    for(i in 1:length(pulse.names)){
        cell.pulse<-tmp$mp[tmp$w.dat[,"wr1"]==pulse.names[i],cell]
        cell.time<-tmp$mp[tmp$w.dat[,"wr1"]==pulse.names[i],1]
        alpha<-sf/length(cell.pulse)
        TraceImpute.2(cell.pulse,cell.time,plotit=T,lmain=pulse.names[i], xspan=alpha)
    }
}

#function to build a table with defined cell types, and selected collumns
TableBrewer<-function(dat, ct.names=NULL){
    require(xlsx)
    dat.name<-deparse(substitute(dat))
    pulse<-select.list(names(dat$bin), multiple=T, title="select variables for table")
    ct.sum<-data.frame()
    
    
    
    if(is.null(ct.names)){
        cell.type.names<-names(dat$cell.types)
        cell.types<-dat$cell.types
    }else{
        cell.type.names<-names(ct.names)
        cell.types<-ct.names
    }
    
    for(z in 1:length(pulse)){
            for(x in 1:length(cell.type.names)){ 
                #first count the number of cells in the cell type group
                ct.sum[as.character(dat.name),cell.type.names[x]]<-length(cell.types[[ cell.type.names[x] ]])
                #sum the collumn with only the cell.types defined rows based on the current selected collumn
                ct.sum[pulse[z],cell.type.names[x]]<-sum(dat$bin[cell.types[[ cell.type.names[x] ]],pulse[z]])
            }
    }
    print('Endter you file name without sapces')
    save.names<-scan(n=1, what='character')
    print(paste(save.names,'xlsx',sep=''))
    write.xlsx(ct.sum, file=paste(save.names,'.xlsx',sep=''))
    return(ct.sum)
}

#########################################
##############################################################
#Function with 3 options. Edit_ct, classify UL , classify thermos
#This follows Marios scheme for classifying our cell types
Cell_Typer<-function(tmp.rd, edit_ct=T, UL_classify=T, thermos_classify=T){

    dropped<-cellzand(tmp.rd$bin,"drop",1)
    #selected bin and dropped
    print("Select The response that coorespond to Neurons,
    ex.
    K+.40mM, and capsaicin.300nM")
    neurons<-cellzand(tmp.rd$bin, , 1)
    neurons<-setdiff(neurons, dropped)
    greeR.Cells<-cellzand(tmp.rd$bin,"gfp.bin" ,1) #selected bin then gfp.bin
    red.cells<-cellzand(tmp.rd$bin,"cy5.bin" ,1)
    caG.Cells<-cellzand(tmp.rd$bin,  grep("caps",names(tmp.rd$bin),ignore.case=T, value=T), 1)
    aitc.cells<-cellzand(tmp.rd$bin,  grep("aitc",names(tmp.rd$bin),ignore.case=T, value=T), 1)
    menth.cells<-cellzand(tmp.rd$bin, grep("menth",names(tmp.rd$bin),ignore.case=T, value=T), 1)
    menth.only<-setdiff(menth.cells, aitc.cells)
    large.cells.330<-cellzand(tmp.rd$c.dat,"area" ,330)#selected c.dat then area

    glia<-setdiff(tmp.rd$c.dat$id, neurons)
    glia<-setdiff(glia, dropped)

    peptidergic<-greeR.Cells
    not.cgrp<-setdiff(neurons, greeR.Cells)

    #Sort green cells first by capsaicin then aitc
    G.C<-intersect(greeR.Cells, caG.Cells)
    G.0<-setdiff(greeR.Cells, G.C)

    G.C.A<-intersect(G.C, aitc.cells)
    G.C<-setdiff(G.C, G.C.A)

    G.A<-intersect(G.0, aitc.cells)
    G.A<-setdiff(G.A, G.C.A)
    G.0<-setdiff(G.0, G.A)
    G.M<-intersect(G.0, menth.cells)
    G.0<-setdiff(G.0, G.M)

    #This gives us G.C, G.C.A, G.A, G.M, and G.0 (zero) under the green cells
    #next we seperate red from unlabeled

    nonpep<-intersect(not.cgrp, red.cells)
    nonpep<-setdiff(nonpep, menth.only)
    unlabeled<-setdiff(not.cgrp, nonpep)

    #Chase down the red classes, the two that are pretty unambiguous are R.A, R.C and R.other

    R.A<-intersect(nonpep, aitc.cells)
    R.A<-setdiff(R.A, caG.Cells)

    R.other<-setdiff(nonpep, R.A)
    R.C<-intersect(R.other, caG.Cells)
    R.C<-setdiff(R.C, aitc.cells)

    R.other<-setdiff(R.other, R.C)

    #This gives us our red groups: R.A, R.C and R.other
    #Finally we chase down our unlabeled groups (unlabeled)

    thermos<-menth.only
    thermos<-intersect(thermos, neurons)
    unlabeled<-setdiff(unlabeled, thermos)
    UL<-intersect(large.cells.330,unlabeled)
    UL<-intersect(UL,neurons)
    UL<-setdiff(UL, caG.Cells)
    UL<-setdiff(UL, aitc.cells)

    US<-setdiff(unlabeled, UL)
    US<-intersect(US,neurons)

    US.A<-intersect(US, aitc.cells)
    US.A<-setdiff(US.A, caG.Cells)
    US.C<-intersect(US, caG.Cells)
    US.C<-setdiff(US.C, US.A)
    US.0<-setdiff(US, US.A)
    US.0<-setdiff(US.0, US.C)

    #review the autosorted cell classes. Remove the cells that are not part of each class and put into discard pile, press "1" to move cells to the discard pile
    if(edit_ct){
        cat(
        "Review the autosorted cell classes. Remove the cells that are not part of each class 
        and put into discard pile, press '1' to move cells to the discard pile (button 1)
        ")
        cat("Editing UL
        Remove anything that responds to capsaicin, menthol or aitc. 
        Also remove cells labeled with IB4 or CGRP-GFP
        ")
        UL.edit<-tcd(tmp.rd, c(UL))
            discard<-UL.edit[[1]]
        UL<-setdiff(UL, discard)

        cat("Editing G.M
        Remove cells that are not green or do not respond to menthol. 
        Also remove red cells and cells responding to AITC	
        ")
        G.M.edit<-tcd(tmp.rd, c(G.M))
            discard1<-G.M.edit[[1]]
        G.M<-setdiff(G.M, discard1)
        discard<-union(discard, discard1)

        cat("Editing G.0
        Remove cells that respond to menthol, capsaicin or aitc, 
        cells that are not green and cells that are red.
        ")
        G.0.edit<-tcd(tmp.rd, c(G.0))
            discard1<-G.0.edit[[1]]
        G.0<-setdiff(G.0, discard1)
        discard<-union(discard, discard1)

        
        cat("Editing G.A
        Remove cells that respond to Capsaicin or have 'noisy' variable baseline, 
        or cells that AITC response does not return to baseline"
        )
        G.A.edit<-tcd(tmp.rd, c(G.A))
            discard1<-G.A.edit[[1]]
        G.A<-setdiff(G.A, discard1)
        discard<-union(discard, discard1)
        
        cat("Editing G.C")
        G.C.edit<-tcd(tmp.rd, c(G.C))
            discard1<-G.C.edit[[1]]
        G.C<-setdiff(G.C, discard1)
        discard<-union(discard, discard1)

        cat("Editing G.C.A")
        G.C.A.edit<-tcd(tmp.rd, c(G.C.A))
            discard1<-G.C.A.edit[[1]]
        G.C.A<-setdiff(G.C.A, discard1)
        discard<-union(discard, discard1)
    
        cat("Editing R.A")
        R.A.edit<-tcd(tmp.rd, c(R.A))
            discard1<-R.A.edit[[1]]
        R.A<-setdiff(R.A, discard1)
        discard<-union(discard, discard1)

        cat("Editing R.C
        Remove cells that respond to aitc or are green
        ")
        R.C.edit<-tcd(tmp.rd, c(R.C))
            discard1<-R.C.edit[[1]]
        R.C<-setdiff(R.C, discard1)
        discard<-union(discard, discard1)
        
        cat("Editing R.Other")
        R.other.edit<-tcd(tmp.rd, c(R.other))
            discard1<-R.other.edit[[1]]
        R.other<-setdiff(R.other, discard1)
        discard<-union(discard, discard1)
        
        cat("Editing Thermosensors
        Remove cells that either do not respond to menthol or 
        have an aitc response larger than the menthol response
        ")
        
        thermos.edit<-tcd(tmp.rd, c(thermos))
            discard1<-thermos.edit[[1]]
        thermos<-setdiff(thermos, discard1)
        discard<-union(discard, discard1)
        
        cat("Editing US.A")
        US.A.edit<-tcd(tmp.rd, c(US.A))
            discard1<-US.A.edit[[1]]
        US.A<-setdiff(US.A, discard1)
        discard<-union(discard, discard1)
        
        cat("Editing US.C")
        US.C.edit<-tcd(tmp.rd, c(US.C))
            discard1<-US.C.edit[[1]]
        US.C<-setdiff(US.C, discard1)
        discard<-union(discard, discard1)
        
        cat("Editing US.0")
        US.0.edit<-tcd(tmp.rd, c(US.0))
            discard1<-US.0.edit[[1]]
        US.0<-setdiff(US.0, discard1)
        discard<-union(discard, discard1)
        cat("
        Sort the discard pile
        1 UL
        2 G.M
        3 G.0
        4 G.A
        5 G.C
        6 G.A.C
        7 R.A
        8 R.C
        9 R.other
        10 thermos
        11 US.A
        12 US.C
        ")
        if(length(discard)>0){
            hand.sort<-tcd(tmp.rd, c(discard))

            #Union between hand sorted and edited autosort. Sort UL into 4 groups based on R3J response (1=prop, 2=jagged, 3=ide, 4=NE), Create a few thermos groups
            UL<-union(UL, hand.sort[[1]])
            G.M<-union(G.M, hand.sort[[2]])
            G.0<-union(G.0, hand.sort[[3]])
            G.A<-union(G.A, hand.sort[[4]])
            G.C<-union(G.C, hand.sort[[5]])
            G.C.A<-union(G.C.A, hand.sort[[6]])
            R.A<-union(R.A, hand.sort[[7]])
            R.C<-union(R.C, hand.sort[[8]])
            R.other<-union(R.other, hand.sort[[9]])
            thermos<-union(thermos, hand.sort[[10]])
            US.A<-union(US.A, hand.sort[[11]])
            US.C<-union(US.C, hand.sort[[12]])
        }else{}
    }else{}
    
    if(UL_classify){
        cat(" Sort the Unlabled Large into
            1:Propriocepters
            2:Jagged
            3:IDE only
            4:No Effect
        ")
        UL.groups<-tcd(tmp.rd, c(UL))
        UL.1<-UL.groups[[1]] #proprioceptor	
        UL.2<-UL.groups[[2]] #jagged
        UL.3<-UL.groups[[3]] #IDE only
        UL.4<-UL.groups[[4]] #no effect
    }
        
    if(thermos_classify){
        thermos.groups<-tcd(tmp.rd, c(thermos))
        thermos.high<-thermos.groups[[1]]
        thermos.low<-thermos.groups[[2]]
        thermos.C<-intersect(thermos, caG.Cells)
    }

    cell.types<-named.list(
        neurons, 
        glia,
        UL,
        G.M,
        G.0,
        G.A,
        G.C,
        G.C.A,
        R.A,
        R.C,
        R.other,
        thermos,
        US.A, 
        US.C,
        US.0
    )
    if(UL_classify){
        UL_ct<-named.list(
            UL.1,
            UL.2,
            UL.3,
            UL.4)
        cell.types<-append(cell.types,UL_ct)
    }else{}
    
    if(thermos_classify){
        thermos_ct<-named.list(
            thermos.high,
            thermos.low,
            thermos.C)
    cell.types<-append(cell.types,thermos_ct)
    }else{}

    tmp.rd$cell.types<-cell.types
    return(tmp.rd)
}

##############################################################
#Function with 3 options. Edit_ct, classify UL , classify thermos
#This follows Marios scheme for classifying our cell types
#########################################
##############################################################
#Function with 3 options_ Edit_ct, classify UL , classify thermos
#This follows Marios scheme for classifying our cell types

#edit_ct=Logical, if true each cell class will be double checked
#UL_classify= If TRUE then classify large diameter cells
#GFP=logical, if TRUE then classify green cells
#cell_types=list input.  This is mainly used if the large cell types have already been classified.
#if so then then the large cell types are passed straight to the cell_types
#181016 If the large diameter cells have been classified, then do not score again.
Cell_Typer_2<-function(tmp_rd, edit_ct=F, UL_classify=T, GFP=T, cell_types=NULL){
    
    if(is.null(cell_types)){
        large_cell_types_names <- NULL
    }else{
        #If your cell_types is not null do large celltyping
        UL_classify <- T
        #perform a test on your cell_types to see if there are large ones
        #based on the names within cell_types
        cell_types_names<-names(cell_types)
        #find ones that have an L
        large_cell_types_names<-grep("^L",cell_types_names,value=T)
        #If you have any that have and L
        if(length(large_cell_types_names)>1){
            #Do not cell_type the large cells
            UL_classify <- F
            UL_ct<-cell_types[large_cell_types_names]
			UL_classes_logic <- T
        }
	}


    dropped<-cellzand(tmp_rd$bin,"drop",1)
    #selected bin and dropped
    cat("Select The response that coorespond to Neurons,
        ex_
        K+_40mM, and capsaicin_300nM
    ")
    #identfy Neurons
    neurons<-cellzand(tmp_rd$bin, , 1)
    #Remove dropped cells from he neuron class
    neurons<-setdiff(neurons, dropped)
    #Idenfy green cells_ Corrected with ROIReview
    if(GFP){
        green_cells<-cellzand(tmp_rd$bin,"gfp.bin" ,1) #selected bin then gfp_bin
    }
    #identify red cells
    ib4_label <- grep("cy5|tritc", names(tmp_rd$bin), value=T)
    red_cells<-cellzand(tmp_rd$bin,ib4_label ,1)
    #define Unlabeled cells as not green or red labeling
    if(GFP){
        unlabeled<-setdiff(neurons, green_cells)
    }else{unlabeled<-neurons}
    unlabeled<-setdiff(unlabeled, red_cells)

    #cells that respond to capsaicin_ These cells wer
    cap_cells<-cellzand(tmp_rd$bin,  
        grep("cap",names(tmp_rd$bin),ignore.case=T, value=T)[length(grep("cap",names(tmp_rd$bin),ignore.case=T, value=T))],
        1)
    #identify AITC responses
    aitc_cells<-cellzand(tmp_rd$bin,  
        grep("aitc",names(tmp_rd$bin),ignore.case=T, value=T)[length(grep("aitc",names(tmp_rd$bin),ignore.case=T, value=T))],
        1)
    #Indentify Menthol Responses
    menth_cells<-cellzand(tmp_rd$bin, 
        grep("men",names(tmp_rd$bin),ignore.case=T, value=T)[length(grep("men",names(tmp_rd$bin),ignore.case=T, value=T))], 
        1)
    #Remove aitc responders to find trpm8 only cells
    menth_only<-setdiff(menth_cells, aitc_cells)
    #Find AITC and capsaicin
    aitc_and_caps<-intersect(aitc_cells,cap_cells)
    #define large cells as larger that 330uM^2
    large_cells_330<-cellzand(tmp_rd$c.dat,"area" ,330)
    
    #define glia is a very weak way_ Antyhing that isnt a
    #neuron is considered glia
    glia<-setdiff(tmp_rd$c.dat$id, neurons)
    glia<-setdiff(glia, dropped)

    cell_types<-named.list(neurons, glia)
    discard<-c()



    ####################
    #GREEN Group
    #Sort green cells first by capsaicin then aitc
    if(GFP){
        #G8 gpf+, menthol negative, capsaicin only 
        G8<-intersect(green_cells, cap_cells)
        G8<-setdiff(G8, aitc_cells)
        G8<-setdiff(G8, menth_cells)
        #now clean the green group?

        #G9 gfp+, menthol negative, aitc and capsaicin
        #first discover cells that are positive for caps and aitc
        #now intersect the green cells with a+ c+
        G9<-intersect(green_cells, aitc_and_caps)
        
        #G10 gfp+, AITC positive only
        G10<-intersect(green_cells, aitc_cells)
        #remove capsaicin from this group
        G10<-setdiff(G10, cap_cells)
        
        #G7 gfp+, Menthol + only
        G7<-intersect(green_cells, menth_cells)
        #remove aitc responders
        G7<-setdiff(G7, aitc_cells)
        #Create G7_capsaicin cells
        G7_c<-intersect(G7,cap_cells)
        
        #now create a group of green responding cell that are
        #not classified by th previous green groups
        #This groups contains miscored Menthol responses and the large cell groups
        G_0<-setdiff(green_cells, c(G7,G8,G9,G10))
		print(G_0)
    }	
    ########################################
    #RED ONLY GROUP
    ########################################
    #remove any green from red cells
    if(GFP){
        red_only<-setdiff(red_cells,green_cells)
    }else{red_only<-red_cells}
    #Chase down the red classes, the two that are pretty unambiguous are R_A, R_C and R_other
    
    #R13 IB4 only,AITC only
    R13<-intersect(red_only, aitc_cells)
    #remove capsaisin responses from this group
    R13<-setdiff(R13, cap_cells)
    
    #R11 IB4 only, Capsaicin only
    R11<-intersect(red_only, cap_cells)
    #remove AITC from this group
    R11<-setdiff(R11, aitc_cells)
    
    #R12 IB4 only, Capsaicin and AITC
    R12<-intersect(red_only, aitc_and_caps)
    
    #R_0 Where the unclassified Red only cells are stored
    R_0<-setdiff(red_only, c(R11,R12,R13))
    #This gives us our red groups: R_A, R_C and R_other
    #Finally we chase down our unlabeled groups (unlabeled)

    #######################################
    #Unlabeled Cell Types
    #######################################
    #N15 no-label, Menthol sensitive
    
    #How do we find menthol responses larger or equal to the
    #aitc response.
    #1 find the cells that respond to menthol 
    #2 find cells taht respond to aitc
    #3 Compare peak heights of these two responses.
    #4 if the AITC response is >= 90% of the menthol response
    #4a add to a new group
    
    aitc_stat<-intersect(
        grep(".max", names(tmp_rd$scp), ignore.case=T, value=T),
        grep("aitc",names(tmp_rd$scp),ignore.case=T, value=T)
    )
    
    menth_stat<-intersect(
        grep(".max", names(tmp_rd$scp), ignore.case=T, value=T),
        grep("men",names(tmp_rd$scp),ignore.case=T, value=T)
    )
    
    menth_stat<-menth_stat[length(menth_stat)]
    
    #find trpm8 and trpa1 containing neurons
    menth_and_aitc_cells<-intersect(menth_cells, aitc_cells)
    if(length(menth_and_aitc_cells) > 0 ){
		menth_and_aitc_neurons<-intersect(neurons, menth_and_aitc_cells)
		trpm8_trpa1<-c()
		for(i in 1:length(menth_and_aitc_neurons)){
			if(tmp_rd$scp[menth_and_aitc_neurons[i],menth_stat] >= 
				((tmp_rd$scp[menth_and_aitc_neurons[i],aitc_stat])*1.1)
			){
				trpm8_trpa1<-c(trpm8_trpa1,menth_and_aitc_neurons[i])
			}
		}
		
		N15_a<-trpm8_trpa1
	}else{
		N15_a <- NULL
	}
    ################################
    #Unlabeled
    ################################

    #Unlabeled smaller neurons responding to menthol and not AITC
    N15<-menth_only
    if(GFP){
        N15<-setdiff(N15, G7)
    }
    #ensure these are neurons
    N15<-intersect(N15, neurons)
    #remove these cells from the unlabeled group
    unlabeled<-setdiff(unlabeled, N15)	
    
    #N15_c Menthol capsaicin
    N15_c<-intersect(N15, cap_cells)
    
    #Now create an unlabeled large group of cells
    UL<-intersect(large_cells_330,unlabeled)
    #ensure they are neurons
    UL<-intersect(UL,neurons)
    #remoce any capsaicin or aitc responders
    UL<-setdiff(UL, c(cap_cells,aitc_cells))
    
    #Create N13, and N16, US is a super category
    US<-setdiff(unlabeled, UL)
    US<-intersect(US,neurons)

    #N14 unlabeled  capsaicin negative
    N14<-intersect(US, aitc_cells)
    N14<-setdiff(N14, cap_cells)
    #N16 unlabeled, capsaicin positive
    N16<-intersect(US, cap_cells)
    N16<-setdiff(N16, aitc_cells)
    
    #create a US_0 class where these additional values are stored
    US_0<-setdiff(US, c(N14,N16))
    
    #N14 is a miscellaneous class that stores addiional unclassified cells
    N14<-union(N14, R_0)
    UC<- union(R_0, US_0)
    #N14<-union(N14, US_0)

    #######################################
    #UL
    #######################################
    if(UL_classify){
        cat(" Sort the Unlabled Large into
            1:Propriocepters
            2:Jagged
            3:IDE only
            4:No Effect
            5:Discard
            PRESS ANY KEY TO CONTINUE
        ")
        scan(n=1)
        UL_groups<-tcd(tmp_rd, c(UL), save_question=F)
        L1<-UL_groups[[1]] #proprioceptor	
        L2<-UL_groups[[2]] #jagged
        L3<-UL_groups[[3]] #IDE only
        L4<-UL_groups[[4]] #no effect
        
        if(edit_ct){discard<-union(discard, UL_groups[[5]])}

        if(GFP){
            cat(" Sort the Unlabled Large into
            1:R3J IDE
            2:no Effect
            3:Discard
            
            PRESS ANY KEY TO CONTINUE
            ")
            scan(n=1)
            
            G_0_sort<-tcd(tmp_rd, c(G_0), save_question=F)
            L5<-G_0_sort[[1]]
            L6<-G_0_sort[[2]]
            if(edit_ct){discard<-union(discard,G_0_sort[[3]])}
        }
    }else{
        UL<-large_cells_330
        UL<-setdiff(UL, c(cap_cells, aitc_cells, menth_cells) )
        #print(UL)
    }

    
    #review the autosorted cell classes_ Remove the cells that are not part of each class and put into discard pile, press "1" to move cells to the discard pile
    if(edit_ct){
        cat(
        "Review the autosorted cell classes_ Remove the cells that are not part of each class 
        and put into discard pile, press '1' to move cells to the discard pile (button 1)
        ")
        if(GFP){
            cat("G7: GFP+, Menthol Only")
            G7.edit<-tcd(tmp.rd, c(G7))
                discard1<-G7.edit[[1]]
            G7<-setdiff(G7, discard1)
            discard<-union(discard, discard1)
            
            cat("G8 GFP+, Capsaicin Only")
            G8.edit<-tcd(tmp.rd, c(G8))
                discard1<-G8.edit[[1]]
            G8<-setdiff(G8, discard1)
            discard<-union(discard, discard1)

            cat("G9 GFP+,AITC AND Capsaicin+")
            G9.edit<-tcd(tmp.rd, c(G9))
                discard1<-G9.edit[[1]]
            G9<-setdiff(G9, discard1)
            discard<-union(discard, discard1)

            cat("G10 GFP+, AITC+ only")
            G10.edit<-tcd(tmp.rd, c(G10))
                discard1<-G10.edit[[1]]
            G10<-setdiff(G10, discard1)
            discard<-union(discard, discard1)
        }
        
        cat("R11 IB4+ Only, Capsaicin Only")
        R11.edit<-tcd(tmp.rd, c(R11))
            discard1<-R11.edit[[1]]
        R11<-setdiff(R11, discard1)
        discard<-union(discard, discard1)

        cat("R12 IB4 only, Capsaicin and AITC only")
        R12.edit<-tcd(tmp.rd, c(R12))
            discard1<-R12.edit[[1]]
        R12<-setdiff(R12, discard1)
        discard<-union(discard, discard1)
        
        cat("R13 IB4 only, AITC only")
        R13.edit<-tcd(tmp.rd, c(R13))
            discard1<-R13.edit[[1]]
        R13<-setdiff(R13, discard1)
        discard<-union(discard, discard1)
    
        cat("N14 No-label, capsaicin only")
        N14.edit<-tcd(tmp.rd, c(N14))
            discard1<-N14.edit[[1]]
        N14<-setdiff(N14, discard1)
        discard<-union(discard, discard1)
        
        cat("N15 No-label, Menthol")
        N15.edit<-tcd(tmp.rd, c(N15))
            discard1<-N15.edit[[1]]
        N15<-setdiff(N15, discard1)
        discard<-union(discard, discard1)
        
        #Remove cells that either do not respond to menthol or have an aitc response larger than the menthol response
        cat("N16 No-label, Capsaicin only")
        N16.edit<-tcd(tmp.rd, c(N16))
            discard1<-N16.edit[[1]]
        N16<-setdiff(N16, discard1)
        discard<-union(discard, discard1)
        
        cat("N17 No-label, Menthol and AITC+ CONVINCING TRPM8 AND TRPA1")
        N17.edit<-tcd(tmp.rd, N17)
            discard1<-N17.edit[[1]]
        N17<-setdiff(N17,discard1)
        discard<-union(discard, discard1)
        
        cat("
            Sort the discard pile
            #1 Large
            #2 Large.green
            #3 G7 (G.M)
            #4 G8 (G.C)
            #5 G9 (G.A.C)
            #6 G10 (G.A)
            #7 R11 (R.C)
            #8 R12 (R.A.C)
            #9 R13 (R.A)
            #10 N14 (US)
            #11 N15 (thermos)
            #12 N16 (US.C)
            #c N17 (trpm8 trpa1)
        ")
        
        if(length(discard)>0){
            hand_sort<-tcd(tmp_rd, c(discard))
            #Union between hand sorted and edited autosort_ Sort UL into 4 groups based on R3J response (1=prop, 2=jagged, 3=ide, 4=NE), Create a few thermos groups
            
            Large.sort<-tcd(tmp.rd, hand.sort[[1]])
            if(GFP){
                Large.green.sort<-tcd(tmp.rd, hand.sort[[2]])
                G7<-union(G7, hand.sort[[3]])
                G8<-union(G8, hand.sort[[4]])
                G9<-union(G9, hand.sort[[5]])
                G10<-union(G10, hand.sort[[6]])
            }
            R11<-union(R11, hand.sort[[7]])
            R12<-union(R12, hand.sort[[8]])
            R13<-union(R13, hand.sort[[9]])
            N14<-union(N14, hand.sort[[10]])
            N15<-union(N15, hand.sort[[11]])
            N16<-union(N16, hand.sort[[12]])
        }else{}#discard option
    }else{}#edit_ct

    if(UL_classify){
        UL_ct<-named.list(
            L1,
            L2,
            L3,
            L4
        )
        cell_types<-append(cell_types,UL_ct)
        if(GFP){
            UL_gfp_ct<-named.list(
                L5,
                L6
            )
			cell_types<-append(cell_types,UL_gfp_ct)
        }else{
        }
    }else{
        if(UL_classes_logic){
            cell_types<-append(cell_types,UL_ct)
        }else{
            cell_types<-append(cell_types,named.list(UL))
        }
    }
    if(GFP){
        gfp_ct<-named.list(
            G7,
            G7_c,
            G8,
            G9,
            G10
        )
        cell_types<-append(cell_types,gfp_ct)
    }else{}
    
    red_ul_ct<-named.list(
        R11,
        R12,
        R13,
        N14,
        N15,
        N15_c,
        N15_a,
        N16,
        UC
    )
    
    cell_types<-append(cell_types,red_ul_ct)
    
    tmp_rd$cell_types<-cell_types
    
    for(i in 1:length(cell_types)){
        print(
            paste( 
                names(tmp_rd$cell_types)[i],
                "=",
                length( tmp_rd$cell_types[[i]] ) 
            )
        )
    }
    return(tmp_rd)
}

# I have a series of pdf files
gif_maker<-function(dense=200, fps=2, file.name=NULL, type='png'){
    require(magick)
    
    #select the reader for 
    if(type=='pdf'){
        reader <- get( paste0('image_read_', "pdf") )
    }
    if(type=='png'){
        reader <- get('image_read')
    }

    
    #MAKE FILE NAME
    if(is.null(file.name)){
        cat("\nThis function will create a gif for either png's or pdfs.\nPlease Enter the name of the file you want to create.\nex. pdfs_in_gif.png\n")
        file.name<-scan(n=1,what="character")
    }
    
    #ASKING AND ANSWERING QUESTIONS
    cat("\nLets create a gif with this data, below are all",type,"s in your experiment \n")
    cat(list.files(pattern=type),sep="\n")
    pdf_imgs<-list.files(pattern=type)
    cat("How many images would you like in your gif? \n")
    imgs_for_gif<-scan(n=1)
    
    #SELECT EACH PDF FOR 
    pdfs_for_gif<-c()
    for(i in 1:imgs_for_gif){
        img_selection<-menu(list.files(pattern=type),title=paste("Select image ",i))
        pdfs_for_gif[i]<-pdf_imgs[as.numeric(img_selection)]
        cat("These are the selected images \n")
        cat(pdfs_for_gif,sep="\n")
    }

    #BEGIN MAKING PDFs, FIRST HAS RED BORDER
    gif<-reader(pdfs_for_gif[1],density=dense)
    gif<-image_border(gif,"red","10x10")
    for(i in 2:length(pdfs_for_gif)){
        gifz<-reader(pdfs_for_gif[i],density=dense)
        gifz<-image_border(gifz,"black","10x10")
        gif<-c(gif,gifz)
    }   

    animation<-image_animate(gif,fps=fps)
    image_write(animation,paste0(file.name,'.gif'))
}
# I have a series of pdf files
gif_png_maker<-function(dense=200,fps=2,file.name=NULL){
    require(magick)
    if(is.null(file.name)){
        cat("Write the name fo the file you would like the ned image to be \n")
        file.name<-scan(n=1,what="character")
        file.name <- paste0(file.name,'.gif')
    }
    
    cat("Lets create a gif with this data, below are all pngs in your experiment \n")
    cat(list.files(pattern="[pP][nN][gG]"),sep="\n")
    imgs<-list.files(pattern="[pP][nN][gG]")

    cat("How many images would you like in your gif? \n")
    imgs_to_add <-scan(n=1)
    
    imgs_for_gif<-c()
    for(i in 1:imgs_to_add){
        img_selection <- menu(list.files(pattern="[pP][nN][gG]"),title=paste("Select image ",i))
        imgs_for_gif[i] <- imgs[as.numeric(img_selection)]
        cat("These are the selected images \n")
        cat(imgs_for_gif,sep="\n")
    }

    #dense<-200
    gif<-image_read(imgs_for_gif[1], density=dense)
    gif<-image_border(gif, "red", "10x10")
    for(i in 2:length(imgs_for_gif)){
        gifz<-image_read(imgs_for_gif[i], density=dense)
        gifz<-image_border(gifz,"black","10x10")
        gif<-c(gif,gifz)
    }   

    #fps=2
    animation<-image_animate(gif,fps=fps)

    image_write(animation,file.name)
}

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
        png(png.name[i], 40,40, res=20, bg="transparent")
        par(bty="n",mai=c(0,0,0,0))
        plot(dat[[t.type]][ levs_min:levs_max, x.names[i] ],type='l',lwd=2,xaxt='n',yaxt='n',col="white", ylim=c(-0.2,peak_max))
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

dice <- function(x, n,min.n=10){
    x.lst <- split(x, as.integer((seq_along(x) - 1) / n))
    x.i <- length(x.lst)
    if(length(x.lst[[x.i]]) < min.n & x.i > 1)
    {
        x.lst[[x.i-1]] <- c(x.lst[[x.i-1]],x.lst[[x.i]])
        x.lst <- x.lst[1:(x.i-1)]
    }
    return(x.lst)
}

#################################
#Welcome to a new method to score cells
#################################
RDView_2<-function(dat, cells=NULL, levs=NULL){
    dat.name<-deparse(substitute(dat))
    cat(
    "HOWDY partner, we R bout to score some rowdy responses \n
    from your cells. Please selact what we should score \n
    and how we should initially sort this data. \n")
        if(is.null(levs)){
            levs<-setdiff(unique(dat$w.dat$wr1),"")
        }else{levs<-levs}
        if(is.null(cells)){
            cells<-dat$c.dat$id
        }else{
            cells<-cells
        }
        cat(
    "\nWitch window region would you like to score????\n \n What do you say?\n")
    
    lev<-levs[menu(levs)]
    #lev<-levs[26]
    #how would you like to sor this variable?
    cat("#############\nAnd how shall we sort? \n ############### \n")
    sorted.cells<-c.sort.2(dat$scp[grep(lev, names(dat$scp),value=T)],cells)
    #sorted.cells
    subset.list<-dice(sorted.cells, 300, 300/4)
    #subset.list


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

osmo_correct<-function(vol_des_mL=50, osmo_original=281, osmo_desired=300){
    osmo_no_glucose <- 270
    gluc_original_M <- osmo_original - osmo_no_glucose
    gluc_desired_M <- osmo_desired - osmo_no_glucose
    
    glucose_fw <- 180
    
    gluc_orignal_grams<-(gluc_original_M)/1000 * glucose_fw * (vol_des_mL)/1000 
    gluc_desired_grams<-(gluc_desired_M)/1000 * glucose_fw * (vol_des_mL)/1000 
    
    glucose_to_add <- gluc_desired_grams-gluc_orignal_grams
    cat(paste("\nGood Day Sir, to correct your orignal osmolarity from",osmo_original,"to", osmo_desired, "please add, \n"))
    cat(paste(glucose_to_add*1000, "mg D-glucose"))
    cat("\nto your orignal solution.\n")
}

#######################################################
#TOPVIEW DEVELOPED BY KEVIN CHASE
########################################################
LinesEvery.TV <- function(dat,m.names, img=dat$img1, pic.plot=TRUE, multi.pic=T, zf=NULL, t.type="mp", snr=NULL, lmain="", cols=NULL, levs=NULL, levs.cols="grey90",  m.order=NULL,rtag=NULL, rtag2=NULL, rtag3=NULL, plot.new=F, sf=1, lw=2, bcex=.6, p.ht=7, p.wd=10, lns=T, pts=F){
    require(png)
    if(class(t.type)=="character"){t.dat<-dat[[t.type]]}# if trace type is empty select the data, you would like your trace to be
    else{t.type<-menu(names(dat));t.dat<-dat[[t.type]]}
    wr<-dat$w.dat[,"wr1"]
    if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")}
    else{levs<-levs}
    m.names <- intersect(m.names,names(t.dat))
    xseq <- t.dat[,1]
    hbc <- length(m.names)*sf+max(t.dat[,m.names])
    hb <- ceiling(hbc)
    library(RColorBrewer)
    
    if(length(m.names) > 0)
    {
        if(!is.null(m.order))
        {	
            tmp<-dat$c.dat[m.names,]
            n.order<-tmp[order(tmp[,m.order]),]
            m.names <- row.names(n.order)
        }
        ### Picture Plotting!		
        #if(XY.plot==T){cell.zoom.2048(dat, cell=m.names,img=img, cols="white",zoom=F, plot.new=T)}
        ## Tool for color labeleing
        if(is.null(cols)){
            cols <-brewer.pal(8,"Dark2")
            cols <- rep(cols,ceiling(length(m.names)/length(cols)))
            cols <- cols[1:length(m.names)]
        } 
        ## Tool for single color labeling
        else {cols<-cols
            cols <- rep(cols,ceiling(length(m.names)/length(cols)))
            cols <- cols[1:length(m.names)]
        }
        
        if(multi.pic){
            if(plot.new){
                if(length(m.names)>10){dev.new(width=16,height=10);layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(10,6), heights=c(6,6))}
                else(dev.new(width=12,height=8))
            }
            else{
                if(length(m.names)>10){layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(10,6), heights=c(6,6))}
                }
        }else{dev.new(width=12,height=8)}
        par(xpd=TRUE,mar=c(4,3,2,2), bty="l")
        plot(xseq,t.dat[,m.names[1]],ylim=c(0,hbc),xlab="Time (min)",main=lmain,type="n", xaxt="n",yaxt="n",xlim=c(min(xseq)-1.5,max(xseq)))#-sf
        bob<-dev.cur()
        axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
        #axis(2, 1.4, )
        text(rep(0,length(m.names)),seq(1,length(m.names))*sf+t.dat[1,m.names],m.names, cex=.5,col=cols,pos=2)

        ## Tool for adding window region labeling
        if(length(wr) > 0){
            #levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
            x1s <- tapply(dat$w.dat[,1],as.factor(wr),min)[levs]
            x2s <- tapply(dat$w.dat[,1],as.factor(wr),max)[levs]
            y1s <- rep(-.3,length(x1s))
            y2s <- rep(hbc+.2,length(x1s))
            rect(x1s,y1s,x2s,y2s,col=levs.cols,border="black")
            cpx <- xseq[match(levs,wr)+round(table(wr)[levs]/2,0)]
            offs <- nchar(levs)*.5
            #			par(xpd=TRUE)
            text(dat$t.dat[match(levs,wr),"Time"],rep(c((sf*.7)*.5,(sf*.7),(sf*.7)/5),length=length(levs)),levs,pos=4,offset=0,cex=bcex*.8)#,offset=-offs}
            #par(xpd=FALSE)
        }
    
        ## Tool for adding line, point and picture to the plot
        for(i in 1:length(m.names)){
            ypos<-t.dat[,m.names[i]]+i*sf
            if(lns){lines(xseq,ypos, lty=1,col=cols[i],lwd=lw)}
            if(pts){points(xseq,ypos,pch=16,col=cols[i],cex=.3)}
            if(!is.null(snr)){
                pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
                pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
                points(xseq[pp1],t.dat[pp1,m.names[i]]+i/10,pch=1,col=cols[i])
                points(xseq[pp2],t.dat[pp2,m.names[i]]+i/10,pch=0,col=cols[i])
            }
        }

    if(is.null(img)){
    img.p<-dat[[select.list(grep("img",names(dat), value=T))]]
    if(is.null(img.p)){img.p<-dat$img1}
    }else{img.p<-img}
    if(is.null(zf)){zf<-20}else{zf<-zf}
    
    #if(pic.plot==TRUE & length(m.names)<=10){
    if(pic.plot==TRUE){
    if(length(m.names)<=100){

        pic.pos<-list()
        for(i in 1:length(m.names)){
            ypos<-t.dat[1,m.names[i]]+i*sf
            pic.pos[[i]]<-ypos}

            for(i in 1:length(m.names)){		
                #if(dat$bin[m.names[1],"mean.gfp.bin"]!=1 & dat$bin[m.names[1],"mean.tritc.bin"]!=1){img.p<-dat$img.gtd #if the cell is neither red or green, then make the img to plot img.gtd
                #}else{img.p<-img} 
                #img.p<-img
                
                img.dim<-dim(dat$img1)[1]

                x<-dat$c.dat[m.names[i],"center.x"]
                left <- x-zf
                if(left<=0){left=0; right=2*zf}
                right<- x+zf
                if(right>=img.dim){left=img.dim-(2*zf);right=img.dim}
                
                y<-dat$c.dat[m.names[i],"center.y"]
                top<-y-zf
                if(top<=0){top=0; bottom=2*zf}
                bottom<-y+zf
                if(bottom>=img.dim){top=img.dim-(2*zf);bottom=img.dim}
                
                #par(xpd=TRUE)
                xleft<-min(dat$t.dat[,1])-xinch(1)
                xright<-min(dat$t.dat[,1])-xinch(.5)
                ytop<-pic.pos[[i]]+yinch(.25)
                ybottom<-pic.pos[[i]]-yinch(.25)
                
                tryCatch(rasterImage(img.p[top:bottom,left:right,],xleft,ybottom,xright,ytop),error=function(e) rasterImage(img.p[top:bottom,left:right],xleft,ybottom,xright,ytop))
            }
        }
    else{
        par(mar=c(0,0,0,0))
        plot(0,0,xlim=c(0,6), ylim=c(0,6), xaxs="i",yaxs="i", xaxt='n', yaxt='n')
        tmp.img<-multi.pic.zoom.2(dat, m.names,img=img.p, labs=T, zf=zf, cols=cols)
        dev.set(bob) # FUCK THIS!
        rasterImage(tmp.img, 0,0,6,6)
    }
    }
    }
        #if(!is.null(pdf.name))
        #{dev.off()}

    #return(pic.pos)
}

#matrix image
PlotHeatMat <- function(mat,wt=14,ht=10,new.dev=T,title="TOPVIEW", dat_name = NULL){
    if(is.null(dat_name)){
        dat_name<-""
    }
    mat <- mat-min(mat)
    mat <- mat/max(mat)
    if(new.dev)
    {
        dev.new(width=wt,height=ht,family="mono",canvas="black",title=title)
        par(fg="darkgrey",col.axis="white",col.lab="grey",col.main="grey",mar=c(1,3,6,1))
        plot(c(0,1),c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",type="n")
    }
    
    rasterImage(mat,0,0,1,1,interpolate=F)
    gx <- !grepl("gap",dimnames(mat)[[2]])
    ux <- unique(dimnames(mat[,gx])[[2]])
    gi <- match(ux,dimnames(mat)[[2]])/ncol(mat)
    xi <- seq(0,(ncol(mat)-1))/(ncol(mat)-1)
    xi.mat <- data.frame(min=tapply(xi,dimnames(mat)[[2]],min))
    xi.mat[,"max"] <- tapply(xi,dimnames(mat)[[2]],max)
    xi.mat[,"med"] <- (xi.mat[,"max"]+xi.mat[,"min"])/2
    xi.lab <- row.names(xi.mat)
    xi.lab[grep("gap",xi.lab)] <- ""
    #xi.mat <- xi.mat[!grepl("gap",row.names(xi.mat)),]
    axis(side=3,at=xi.mat[xi.lab != "","med"],labels=xi.lab[xi.lab != ""],las=2,col.axis="darkgrey")
    axis(side=3,at=xi.mat[xi.lab=="","min"],labels=NA,tck=1,lwd=.2)
    axis(side=3,at=xi.mat[xi.lab=="","max"],labels=NA,tck=1,lwd=.2)
    yi <- seq(0,(nrow(mat)-1))/(nrow(mat)-1)
    y.names <- dimnames(mat)[[1]]
    y.names <- sub("^w[0987654321]*\\.","",y.names)
    yi.mat <- data.frame(min=tapply(yi,y.names,min))
    yi.mat[,"max"] <- tapply(yi,y.names,max)
    yi.mat[,"med"] <- (yi.mat[,"min"]+yi.mat[,"max"])/2
    yi.mat <- yi.mat[!grepl("blank",row.names(yi.mat)),]
    axis(side=2,at=1-yi.mat[,"med"],labels=row.names(yi.mat),las=3,cex.axis=.5,col.axis="darkgrey")

    par(xpd=T)
    points(par('usr')[1], par('usr')[3], col='white', cex=6, pch = 16)
    text(par('usr')[1], par('usr')[3], "END", col='black', cex=1, pch = 16)
    text(par('usr')[2]-xinch(.5), par('usr')[3]-yinch(.1), dat_name, col='white', cex=.7 )
    

    #dev.off()
}

#tmp is an RD object
#x.names defines the cells to display
#wt = device window width
#ht = device window height
#scale.var is the variable used to scale each cell.  If not in the RD object defaults to a log scale transform then each row scale 0-1 min to max.
#aux.var is a list of auxillary variables to be displayed to the right of the traces.  If there are missing values, no variation or variables not in the RD object they are not shown.
#img is the img name sent to LinesEvery.TV
#t.type is the trace type data sent to LinesEvery.TV
#title is the device window title.
#190508, Added a stat val. This sorts the traces based on the trace statistic you want it to

TopView <- function(tmp, x.names=NULL, wt=7, ht=4, scale.var="mean.sm.sd", aux.var=c("diameter","IB45.bin","gfp5.bin"), img="img1", t.type="blc", title="TOPVIEW", stat_val='max', dat_name=NULL){
    if( is.null(dat_name) ){
        dat_name <- deparse(substitute(tmp))
    }
    
    #vet the vars
    #m.tot <- CollectMulti(rd.names=c(deparse(substitute(tmp))))
    m.tot <- CollectMulti(rd.names=dat_name)
    
    scale.var <- intersect(scale.var,names(m.tot))
    aux.var <- intersect(aux.var,names(m.tot))
    
    if(is.null(x.names)){
        x.names <- row.names(tmp$bin)
    }

    blc.s <- as.matrix( t(tmp$blc[,x.names]) )
    name2 <- make.names( tmp$w.dat[,"wr1"], unique=F)
    name2[is.element(name2,c("X","epad"))] <- "gap"	
    dimnames(blc.s)[[2]] <- name2

    if(length(scale.var)==1)
    {
        blc.s <- sweep(blc.s,1,m.tot[x.names,scale.var],'/')
    }
    else
    {
        blc.s <- log10(blc.s+1)
        med <- apply(blc.s,1,min)
        blc.s <- sweep(blc.s,1,med,'-')
        blc.s[blc.s<0] <- 0
    }
    blc.s <- blc.s-min(blc.s)
    blc.s <- blc.s/max(blc.s)	

    if(length(aux.var) > 0)
    {
        aux.mat <- NULL
        n <- ceiling(nrow(tmp$w.dat)/50)
        for(i in aux.var)
        {
            mat1 <- matrix(rep(m.tot[x.names,i],n),ncol=n)
            dimnames(mat1)[[1]] <- x.names
            dimnames(mat1)[[2]] <- rep(i,n)
            aux.mat <- cbind(aux.mat,mat1)
        }
        aux.mat[is.na(aux.mat)]<-0
        aux.min <- apply(aux.mat,2,min)
        aux.mat <- sweep(aux.mat,2,aux.min,'-')
        aux.max <- apply(aux.mat,2,max)
        aux.mat <- sweep(aux.mat,2,aux.max,'/')
        aux.mean <- apply(aux.mat,2,mean)
        #aux.mat <- aux.mat[,is.na(aux.mean)]

        blc.s <- cbind(blc.s,aux.mat)
    }

    name2 <- dimnames(blc.s)[[2]]	
    name2[!is.element(name2,names(m.tot))] <- "gap"
    seqi <- seq(0,(ncol(blc.s)-1))/(ncol(blc.s)-1)
    click.id <- data.frame(x=tapply(seqi,as.factor(dimnames(blc.s)[[2]]),median))
    click.id[,"y"] <- 1
    click.id <- click.id[intersect(row.names(click.id),names(m.tot)),]
    click.id <- click.id[order(click.id$x),]
    click.id.rn.torn<-row.names(click.id)[nrow(click.id)-length(aux.var)]
    vals_to_click <- paste0(click.id.rn.torn,'.', stat_val)
    click.vals <- m.tot[ x.names,  row.names(click.id)]
    
    PlotHeatMat(blc.s,wt=wt,ht=ht,title=title, dat_name=dat_name)
    xy.click <- list(x=1,y=1)
    while(xy.click$x > 0 | xy.click$y > 0)
    {
        xy.click <- locator(n=1,type="n")
        if(xy.click$y > 1)
        {
            sort.trt <- row.names(click.id)[which.min(abs(xy.click$x-click.id[,"x"]))]
            sval <- click.vals[x.names,sort.trt]
            x.names <- x.names[order(sval)]
            blc.s <- blc.s[x.names,]
            PlotHeatMat(blc.s,new.dev=F,wt=wt,ht=ht, dat_name=dat_name)
        }
        if(xy.click$y < 1)
        {
            len1 <- length(x.names)-1
            y.i <- abs(seq(0,len1)/len1 - (1-xy.click$y))
            names(y.i) <- x.names			
            sort.i <- names(sort(y.i)[1:10])
            di <- dev.cur()
            if( length( ls(pattern="^lines_tv$") )<1 ){
                dev.new(width=10, height=12)
                lines_tv<-dev.cur()
            }else{ 
                dev.set(lines_tv) 
            }
            LinesEvery.TV(tmp,sort.i,lw=3,levs.cols=grey(.95),img=tmp[[img]],t.type=t.type,m.order <- seq(1,length(sort.i)),rtag="diameter",rtag2="gfp5.bin",rtag3="IB45.bin",zf=15,cols="black")
            dev.set(di)
        }
        
    }
    dev.off(di)
    dev.off(lines_tv)
}

census_viewer  <- function(dat){
cat(

"This function will essecially return cells from a specified cell\nin the census table
\n1. Select all of cells from a specific cell class. 
\n1a. If you click cancel all cells will be returned.
\n2. bin >> collumn >> cell class cells scored as one.
\n3. returns a vector of cell names ex c(X.3, X.30)
"
)

	(cell_list_name <- grep("^cell", names(dat), value=T))
	(cell_types <- names( dat[[ cell_list_name ]] ))
	(cell_type_name <- select.list( cell_types, title="Select the cell_type" ))
	
	#Tool to return all cells if cancel is selected.
	if(cell_type_name == ''){
		cell_type <- dat$c.dat$id
	}else{
		(cell_type <-dat[[ cell_list_name ]] [[ cell_type_name ]])
	}
	
	(bin_col <- select.list(names(dat$bin), title="Select bin collumn"))
	(cells <- cell_type[ dat$bin[cell_type , bin_col] == 1 ])
	
	if( length(cells) == 0 ){
		return(NA)
	}else{
		cells_to_view <- list()
		cells_to_view[[ 'name' ]]<- paste0(cell_type_name,"__", bin_col)
		cells_to_view[[ 'cells' ]] <- cells
		return(cells_to_view)
	}
}	
	
#Funciton to save the work along with create a unique savehistory

saveRD <- function(dat){
    cat("\nDO NOT CLOSE UNTIL I SAY YOU CAN!\nWait for the sound...")
    flush.console()
    bringToTop(-1)
    Sys.sleep(1)

    #History Saver
    experimentorsName <- strsplit(getwd(),'/')[[1]][2]
    historyName <- paste(experimentorsName, Sys.time(), 'History.r')
    historyName <- gsub(":", '_',historyName)
    savehistory(historyName)

    #Exp Saver
    expName <- deparse(substitute(dat))
    #expToSave <- get(expName, envir = .GlobalEnv)
    assign(expName, dat)
    save(list=expName, file=paste0(expName,".Rdata") )
    alarm()
    cat('\nYou can now close. Please consider cleaning up the file,\n',historyName,'\n')
}

