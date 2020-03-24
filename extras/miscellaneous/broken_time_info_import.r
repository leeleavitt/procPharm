main_dir <- "Y:/Sabrina Kozel"
area_conversion <- 1.625
img_name_vec <- c(
       "bf.start.png",
       "fura2.png",
       "bf.start.lab.png",
       "roi.img.png"     
)
image_question=T

cat('#########################################\nPHARM HAVEST\n#########################################\n')
cat(readLines( "Z:/farm.txt"),sep='\n')
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
exp_dir_1 <- sub("./","",exp_dir)
exp_dir_2 <- lapply(strsplit(exp_dir_1,' '), function(x) x[1])
exp_dir_3 <- Reduce(c, exp_dir_2)
exp_dir_4 <- lapply(strsplit(exp_dir_3, '/'), function(x) x[length(x)] )
exp_dir_5 <- Reduce(c, exp_dir_4)
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
    if(length(list.files(pattern="video_data.txt")) < 1){
        py_pharm <- import('python_pharmer')
        video <- list.files(pattern="^video.*nd2$")
        if( length(video) > 1 ){
            cat("\nSelect your video to process\n")
            video_num <- menu(video)
            video <- video[video_num]
        }
        # Now read in video Data
        py_pharm$video_roi_extractor_faster(video)
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
    time_info <- read.delim("time.info.txt", fileEncoding="UCS-2LE")    
    time_min<-seq(from=0, to=57.7, length.out=923)
    #time_min <- round(time_info[2]/1000/60, digits=3)
    # Create row.names
    cell.names <- paste("X.", colnames(t.340)[-1], sep="")
    traces <- ls(pattern = "^[t.]{2}.*[0-9a-z]{3}")
    for( j in 1:length(traces) ){
        traze <- get( traces[j] )
        traze[,1] <- time_min
        traze <- as.data.frame(traze)
        colnames(traze) <- c("Time",cell.names)
        row.names(traze)<-time_min
        assign(traces[j], traze)
    }

    ########################################################
    # wr1 import
    ########################################################
    wrdef<-"wr1.docx"
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

    ########################################################
    #CELL DATA PROCESSING
    ########################################################\
    c.dat<-read.delim(cell_data_name, header=T, sep="\t")
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
    bin["drop"] <- 0 #maybe try to generate some drop criteria from the scp file.
    bin<-pf.function(bin,levs)

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
