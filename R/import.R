
#' Main robust import function
#' @param area_conversion conversion factor from the object used from the experiment
#' @param img_name_vec vector of image names ex c('bf.gfp.cy5', 'bf.gfp')
#' @param image_question will ask you to fill in images if you don't have an img_name_vec
#' @export 
pharming_harvest <- function(main_dir=NULL, area_conversion=1.625, img_name_vec = NULL, image_question=T){
    cat('#########################################\nPHARM HARVEST\n#########################################\n')
    #tryCatch(cat(readLines("Y:/Box Sync/procpharm/farm.txt"),sep='\n'), error=function(e) cat('\nPHARM ON PHARM ANIMAL\n'))
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
        cell_data_name <- list.files(pattern= '^cell.*txt$')
        #cell_data_name<-"celldatacells_filtered.txt"

        ################################################
        #IMAGE IMPORT
        ################################################
        #Names of the files you want loaded into the RD file
        if( is.null(img_name_vec) ){
            img_name_vec<-c(
                "bf.gfp.tritc.start.png",
                "gfp.tritc.end.ci.ltl.rs.png",
                "tritc.end.ci.ltl.rs.png",
                "gfp.end.ci.ltl.rs.png",
                "bf.start.lab.png",
                "fura2.png",
                "fura2.divide.start.png",
                "dapi.end.lab.png")
        }
        # Add images
        img_list<-list()
        for( j in 1:length(img_name_vec) ){
            img_list[[ paste0("img",j) ]] <- tryCatch(png::readPNG(img_name_vec[j]), error=function(e)NULL)
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
                        img_list[[paste0("img",j)]] <- png::readPNG(png_imgs[selection])
                    }
                    cat('\nI have added ', png_imgs[selection],' to position ',j,'\n')
                }
             }
        }
            
            
            
        ########################################################
        #VIDEO PROCESSING
        ########################################################
        c_dat_make <- file.info(cell_data_name)$mtime
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
			}else{
                video <- list.files(pattern="^[Vv]ideo.*nd2$")
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
		#video <- list.files(pattern="^[Vv]ideo.*nd2$")
		
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
		# If it is a wr1.docx, continue
        if( length(wrdef_logic) == 1 ){
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
            #Wr<-length(wr[,1])#complete and revise this section
            wr['at'] <- wr['at'] - (10/60)

            w.dat <- MakeWr(t.dat,wr)
		}
        ########################################################
        #CELL DATA PROCESSING
        ########################################################\
        c.dat<-read.delim(cell_data_name, header=T, sep="\t")
        #Now use the neewly created roi_checker i built in python
        roiToRemove <- Reduce(c,py_pharm$roi_checker()) + 1
        if(length(roiToRemove)>0){
            cat("\nremoving ROI:\n")
            print(roiToRemove)
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
        tmp.rd <- list(t.dat=t.dat, t.340=t.340, t.380=t.380, w.dat=w.dat, c.dat=c.dat)

        levs<-setdiff(unique(as.character(w.dat[,2])),"")
        snr.lim=5; hab.lim=.05; sm=2; ws=3; blc="SNIP"
        pcp <- ProcConstPharm(tmp.rd, sm, ws, blc)
        tmp.rd$scp <- ScoreConstPharm(tmp.rd, pcp$blc, pcp$snr, pcp$der, snr.lim, hab.lim, sm)
        tmp.rd$blc <- pcp$blc
        #bin <- bScore(pcp$blc, pcp$snr, snr.lim, hab.lim, levs, tmp.rd$w.dat[,"wr1"])
        #bin <- bin[,levs]
        tmp.rd <- bscore2(tmp.rd)
        tmp.rd$bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp 

        tmp.rd <- TraceBrewer(tmp.rd) 
        
        tmp.rd <- c(tmp.rd, img_list)

        # Add the models
        tmp.rd <- traceProbMaker(tmp.rd)
        tmp.rd <- imageProbMaker(tmp.rd, verbose = T)

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

#' Old school way to import data
#' 
#' First step of this function is to change the working directory \code{setwd("Z:/Lee Leavitt/experimentFolder")}
#' Next use \code{ReadDataDump.lee.2("RD.date.age.gender.microscope.protocol", "bf.gfp.cy5")}
#' To make this work, you already need to have a \code{'Data (full).txt', 'ROI Data.txt', 'wr1.docx'}
#' @export
ReadDataDump.lee.2 <- function(rd.name=NULL,img1=NULL,img2=NULL,img3=NULL,img4=NULL,img5=NULL, img6=NULL, img7=NULL, img8=NULL, fancy=F,fname="Data (full).txt",wrdef="wr1.docx", Wr=NULL, c.dat="ROI Data.txt" ,sep="\t"){
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
        c.dat.name <- c.dat
        c.dat<-read.delim(file=c.dat.name,fileEncoding="UCS-2LE", sep=sep)
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
        
        c.dat <- c.dat[cnames]#create c.dat with specified collumns from cnames
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

        }else{
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
    c.dat <- cbind(c.dat, read.delim(file=c.dat.name,fileEncoding="UCS-2LE", sep=sep))
    
    #####################################################
    # Window Region Definition
    #####################################################
    wrdef <- list.files(pattern = '^wr1')
    wrdef_logic <- grep(".docx", wrdef, ignore.case=T, value=T)
    # If it is a wr1.docx, continue
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
        #Wr<-length(wr[,1])#complete and revise this section
        wr['at'] <- wr['at'] - (10/60)

        w.dat <- MakeWr(t.dat,wr)
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
    # Initial and simple Data processing
    tmp.rd <- list(t.dat=t.dat,w.dat=w.dat,c.dat=c.dat)
    levs<-setdiff(unique(as.character(w.dat[,2])),"")
    snr.lim=5; hab.lim=.05; sm=2; ws=3; blc="SNIP"
    pcp <- ProcConstPharm(tmp.rd,sm,ws,blc)
    scp <- ScoreConstPharm(tmp.rd,pcp$blc,pcp$snr,pcp$der,snr.lim,hab.lim,sm)
    tmp.rd <- bscore2(tmp.rd)
    bin <- tmp.rd$bin

    #bin <- bScore(pcp$blc, pcp$snr, snr.lim, hab.lim, levs, tmp.rd$w.dat[,"wr1"])
    bin <- bin[,levs,drop=F]
    bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp 

    tmp.rd <- list(
        t.dat=t.dat,
        w.dat=w.dat,
        c.dat=c.dat, 
        bin=bin, 
        scp=scp, 
        #snr=pcp$snr, 
        blc=pcp$blc)
        #der=pcp$der) 
    
    tmp.rd <- TraceBrewer(tmp.rd) 
    tmp.rd <- list(
        t.dat=t.dat,
        w.dat=w.dat,
        c.dat=c.dat, 
        bin=bin, 
        scp=tmp.rd$scp, 
        #snr=pcp$snr, 
        blc=tmp.rd$blc)
        #der=pcp$der) 


        # tmp.rd <- list(t.dat=t.dat,w.dat=w.dat,c.dat=c.dat)
        # # #####################################################
        # # #Create Despiked data
        # # #####################################################
        # wts <- tmp.rd$t.dat
        # for(i in 1:5) #run the despike 5 times.
        # {
        #     wt.mn3 <- Mean3(wts)
        #     wts <- SpikeTrim2(wts,1,-1)
        #     print(sum(is.na(wts))) #this prints out the number of points removed should be close to 0 after 5 loops.
        #     wts[is.na(wts)] <- wt.mn3[is.na(wts)]
        # }
        # tmp.rd$mp <- wts

        # # Initial Data processing
        # levs<-setdiff(unique(as.character(w.dat[,2])),"")
        # snr.lim=4;hab.lim=.05;sm=2;ws=20;blc="SNIP"
        
        # pcp <- ProcConstPharm(tmp.rd$mp,sm,ws,blc)
        # scp <- ScoreConstPharm(tmp.rd,pcp$blc,pcp$snr,pcp$der,snr.lim,hab.lim,sm)
        # bin <- bScore(pcp$blc,pcp$snr,snr.lim,hab.lim,levs,tmp.rd$w.dat[,"wr1"])
        # bin <- bin[,levs]
        # bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
        # bin<-pf.function(bin,levs)

        # tmp.rd <- TraceBrewer(tmp.rd)

        
        # tmp.rd$t.dat<-t.dat
        # tmp.rd$w.dat<-w.dat
        # tmp.rd$c.dat<-c.dat
        # tmp.rd$bin<-bin
        # tmp.rd$scp<-tmp.rd$scp
        # tmp.rd$blc<-tmp.rd$blc
        

        # Add images
        if(!is.null(img1)){tmp.rd$img1<-png::readPNG(img1)}
        if(!is.null(img2)){tmp.rd$img2<-png::readPNG(img2)}
        if(!is.null(img3)){tmp.rd$img3<-png::readPNG(img3)}
        if(!is.null(img4)){tmp.rd$img4<-png::readPNG(img4)}
        if(!is.null(img5)){tmp.rd$img5<-png::readPNG(img5)}
        if(!is.null(img6)){tmp.rd$img6<-png::readPNG(img6)}
        if(!is.null(img7)){tmp.rd$img7<-png::readPNG(img7)}
        if(!is.null(img8)){tmp.rd$img8<-png::readPNG(img8)}


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
    assign(rd.name, tmp.rd, envir=.GlobalEnv)
	assign(rd.name,tmp.rd)
    save(list=rd.name,file=f.name)
    return(paste(nrow(tmp.rd$c.dat),"traces read saved to ",f.name))
    #save as RD file
}

#########################################################################
# Window Import

#IMPORT DOCX tables
docx.wr1.importer<-function(file.name='wr1.docx'){
    #if( !library(docxtractr, logical.return=T) ){install.packages('docxtractr')}else{}
    #read in docx
    wr1<- invisible(docxtractr::read_docx(file.name))
    #Extract each table
    wr1 <- invisible(docxtractr::docx_extract_all_tbls(wr1, guess_header=F))

    #out table is the third one
    wr1<-Reduce(c,wr1[[length(wr1)]])
    #split up each vaue based on a single space
    wr1<-strsplit(wr1, ' ')
    
    #Now perform a test and provide a wait if there is an error where the window 
    #region has to little information
    error<-0
    for(i in 1:length(wr1)){
        if(length(wr1[[i]])>1 & length(wr1[[i]])<3){
            cat("\nThere is an >1 <2 info error at\n",wr1[[i]])
            error<-error+1
        }
    }
    
    for(i in 1:length(wr1)){
        if(length(wr1[[i]])>3){
            cat("\nThere is an >3 info error at\n",wr1[[i]])
            error<-error+1
        }
    }

    if(error>0){
        cat("\nYou have a total of ",error," errors\n")
        alarm()
        cat("\nFix these Errors\n")
        cat("Once fixed PRESS ENTER\n")
        scan(n=1)
        cat("\nThese are your window region definitions
        If you would like to make anymore changes do so now
        ")
        
        wr1<-docxtractr::read_docx(file.name)
        #Extract each table
        wr1<-docxtractr::docx_extract_all_tbls(wr1, guess_header=F)
        #out table is the third one
        wr1<-Reduce(c,wr1[[length(wr1)]])
        #split up each vaue based on a single space
        wr1<-strsplit(wr1, ' ')
    }

    wr1.logic <- unlist(lapply(wr1, function(x) length(x)>1 ))
    wr1.locations <- which(wr1.logic, arr.ind=T)
    wr1 <- wr1[wr1.locations]
    #wr1<-Reduce(rbind,wr1)
    wr1 <- do.call(cbind, lapply(wr1, data.frame, stringsAsFactors=F))
    row.names(wr1) <- c('at','treatment','duration')
    wr1['at',] <- as.numeric(wr1['at',])-(10/60)
    return(wr1)
}

# Creates Window region from the docx file
MakeWr.docx <- function(t.dat,wr1,padL=0,padR=0){
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

# Function that creates the w.dat
MakeWr <- function(t.dat,wr1,padL=0,padR=0){
    w.dat <- t.dat[,1:2]
    names(w.dat)[2] <- "wr1"
    w.dat["wr1"] <- ""
    wr1["treatment"] <- make.names(wr1[,"treatment"],unique=T)
    for(i in 1:nrow(wr1))
    {
        x1 <- which.min(abs(wr1[i,"at"]-t.dat[,"Time"]))
        x2 <- which.min(abs((wr1[i,"at"] + wr1[i,"duration"]) - t.dat[,"Time"]))
        w.dat[max((x1-padL),1):min((x2+padR),nrow(t.dat)),"wr1"] <- as.character(wr1[i,"treatment"])
    }
    return(w.dat)
}

#' Made a mistake in you window regions?
#' If it is only the time sequence, but all other info is corrected
#' then make complete F.  This will allow you to select the windows that need reapri
#' if the naming is off, then make complete=F.  You will need to do a complete reapri
#' you will lose all information from RDView
#' @export 
WindowRepair <- function(dat, rescore = F){
    # Get rid of that dam epad!
    if('epad' %in% dat$w.dat$wr1){
        epadRM <- dat$w.dat$wr1 != "epad"
        
        dat$blc <- dat$blc[epadRM,]
        dat$w.dat <- dat$w.dat[epadRM,]
        dat$t.dat <- dat$t.dat[epadRM,]
    }
    
    datOrig <- dat
    #Now do all of this to tmp.rd
    wrdef <- list.files(pattern = '^wr1')
    if(length(wrdef) > 1){
        cat("\nWhich file contains the updated information?\n")
        wrdef <- wrdef[menu(wrdef, title=)]
    }

    # Is the selection a docx?
    wrdef_logic <- grep("[.]docx$", wrdef, ignore.case=T, value=T)
    # If it is a wr1.docx, continue
    if( length(wrdef_logic) == 1 ){
        wr <- docx.wr1.importer(wrdef)
        dat$w.dat <- MakeWr.docx(dat$t.dat, wr)		

        ## Check for duplicated rows	
        if(length(which(duplicated(row.names(dat$t.dat))))>=1){
            dup<-which(duplicated(row.names(dat$t.dat)))
            dat$t.dat <- dat$t.dat[-dup,]
            dat$w.dat <- dat$w.dat[-dup,]
        }
    }else{
        tryCatch(wr <- ReadResponseWindowFile(wrdef)
            ,error=function(e) print("YOU NEED A wr1.docx or wr1.csv")
        )
        # Since this is the csv we need to move the window back ten seconds
        wr['at'] <- wr['at'] - (10/60)

        dat$w.dat <- MakeWr(dat$t.dat, wr)
    }

    tmp <- TraceBrewer(dat, F, F)
    # magic i need todo to keep the additional scp stats if there are any
    tmpscp <- cbind(tmp$scp, dat$scp[ (dim(tmp$scp)[2] + 1) : dim(dat$scp)[2] ])
    tmp$scp <- tmpscp
    dat <- tmp

    # Now we need to understand a few aspects of this series of data
    # 1: Are there more or less window regions?
    ## if more window regions we need to do a 
    newPulses <- setdiff(unique(as.character(dat$w.dat$wr1)),"")
    origPulses <- setdiff(unique(as.character(datOrig$w.dat$wr1)), "")

    # If the number ofPulses have not changed, then determine which window regions have changed
    if(length(newPulses) == length(origPulses)){
        # Here we look at if the starting points have changes
        origStarts <- tapply(datOrig$w.dat$Time, as.factor(datOrig$w.dat$wr1),min)[origPulses]
        newStarts <- tapply(dat$w.dat$Time, as.factor(dat$w.dat$wr1),min)[newPulses]
        startsLog <- origStarts %in% newStarts
        
        # Here we look if the ending points have changed
        origEnds <- tapply(datOrig$w.dat$Time, as.factor(datOrig$w.dat$wr1),max)[origPulses]
        newEnds <- tapply(dat$w.dat$Time, as.factor(dat$w.dat$wr1),max)[newPulses]
        endsLog <- origEnds %in% newEnds

        # Here we look if the window regions have the same names

        # This shows me which windows have changed
        changedWindowsLogic <- apply(cbind(!startsLog, !endsLog), 1, any)
        newPulsesToUpdate <- newPulses[changedWindowsLogic]

        # Now we need to see if the windows can enter specific neural networks
        toCatch <- c("^[aA][iI][Tt][Cc]","^[mM][eE][nN][tT][hH]", "^[cC][aA][pP][sS]","[kK].*40")
        pulsesToAddInNN <- unlist(sapply(toCatch, function(x) grep(x, newPulsesToUpdate, value = T)))
        if(length(pulsesToAddInNN) > 0){
            cat("These pulses will enter neuralNets (smart score):\n", pulsesToAddInNN,"\n")
            dat <- traceProbMaker(dat, T, names(pulsesToAddInNN))
        }

        pulsesToBscore <- setdiff(newPulsesToUpdate, pulsesToAddInNN)
        if(length(pulsesToBscore) > 0){
            cat("These pulses will enter bScore (dumb score):\n", pulsesToBscore,"\n")
            dat <- bscore2(dat, pulsesToBscore)
        }

        # Now see if there has been any name changes
        nameChange <- newPulses %in% origPulses
        if(any(!nameChange)){
            changedNameLogic <- !nameChange
            
            # This shows me which windows have changed
            changedWindowsLogic <- apply(cbind(changedNameLogic), 1, any)
            newPulsesToUpdate <- newPulses[changedWindowsLogic]
            originalPulsesToUpdate <- origPulses[changedWindowsLogic]

            for(i in 1:length(newPulsesToUpdate)){
                # the w.dat and the scp have the correct names
                # bin
                names(dat$bin)[names(dat$bin) == originalPulsesToUpdate[i]] <- newPulsesToUpdate[i]
                
                # uncMat
                tryCatch({
                    names(dat$uncMat)[names(dat$uncMat) == originalPulsesToUpdate[i]] <- newPulsesToUpdate[i]
                }, error = function(e) NULL)
            }

            if(rescore){
                # Now we need to see if the windows can enter specific neural networks
                toCatch <- c("^[aA][iI][Tt][Cc]","^[mM][eE][nN][tT][hH]", "^[cC][aA][pP][sS]","[kK].*40")
                pulsesToAddInNN <- unlist(sapply(toCatch, function(x) grep(x, newPulsesToUpdate, value = T)))
                if(length(pulsesToAddInNN) > 0){
                    cat("These pulses will enter neuralNets (smart score):\n", pulsesToAddInNN,"\n")
                    dat <- traceProbMaker(dat, T, names(pulsesToAddInNN))
                }

                pulsesToBscore <- setdiff(newPulsesToUpdate, pulsesToAddInNN)
                if(length(pulsesToBscore) > 0){
                    cat("These pulses will enter bScore (dumb score):\n", pulsesToBscore,"\n")
                    dat <- bscore2(dat, pulsesToBscore)
                }
            }
        }
        
    } else if(length(newPulses) != length(origPulses)){
        cat("There seem to be too many changes for us to deal with, all scores will be\n wiped. Do you want to continue?\n")
        
        alarm()
        
        sel <- c('yes', 'no')
        sel <- sel[menu(sel)]

        if(sel == 'yes'){
            dat <- bscore2(dat)
            dat <- traceProbMaker(dat)
        }else if(sel == 'no'){
            stop("No changes were made")
        }
    }

    dat$SETTINGS <- NULL
    return(dat)
}

#' Something happened within you experiment and you need to delete a region of row fro the traces
#' to recover the data and clean up the display.
#' dat: This is the RD object
#' cell: Cell the display on the plot, if left empty, X.1 will be used
#' complete: logical, If you say true, all window regions will be reassessed. If False window region can be selected.
#' @export
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

    tmp.rd<-TraceBrewer(tmp.rd)
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
    
    tryCatch(bringToTop(-1), error=function(e)NULL)
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

ReadResponseWindowFile <- function(fname){
    dat <- read.csv(fname)
    return(dat)
}

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

#each point is replaced with the mean of the two neighboring points
Mean3 <- function(wt){
	wt.mn <- (wt[-c(1,2),]+wt[-c(nrow(wt),(nrow(wt)-1)),])/2
	wt[2:(nrow(wt)-1),] <- wt.mn
	return(wt)
}

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
