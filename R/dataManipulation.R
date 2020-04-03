#' Now lets create the function to put premade groups into a binary collumns
#' @export
col_binner<-function(dat,cells){
    cell.names<-select.list(names(cells),multiple=T)
    cells<-cells[cell.names]
    for(i in 1:length(cell.names)){
        dat$bin[cell.names[i]]<-0
        dat$bin[ cells[[i]], cell.names[i] ]<-1
    }
    return(dat)
}

#' Create Binary Classes of cells
#' dat IS THE rd. INPUT
#' LEVS, IS AN OPTIONAL ARGUEMENT, IF LEFT BLACK THE FUNCTION WILL LOOK IN THE BIN DATA.FRAME COLLUMN ANMES
#' TO ALLOW TO MANUALLY SECT THE COLLUMNS YOU WANT TO COMBINE
#' @export
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

#' Create Binary Classes of cells
#' @export
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

#' @export
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

#' This takes a pf and allows you to create a binarry table based on the barcode
#' Created in pf.function
#' @export
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

#' Function to select rows based on collumn parameters
#' dat can be either a raw RD object or an RD dataframe
#' ex dat -or- dat$bin
#' @export
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

#' @export
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

#' given a list of file names collect and merge all bin scp and c.dat data
#' @export
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

#' function to obtained sorted cell names based off 
#' collumn names from c.dat and bin
#' @export
c.sort<-function(dat,char=NULL){
    char<-select.list(names(dat))
    sort.dir<-select.list(c("TRUE", "FALSE"), title="Decreasing?")
    bob<-row.names(dat[order(dat[,char], decreasing=sort.dir),])
    return(bob)
}

#' @export
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

#' function to build a table with defined cell types, and selected collumns
#' @export
TableBrewer<-function(dat, ct.names=NULL, save=T, xlsx=T){
    dat.name<-deparse(substitute(dat))
    pulse<-select.list(names(dat$bin), multiple=T, title="select variables for table")
    ct.sum<-data.frame()

    if(is.null(ct.names)){
        #F7: Load cell Types into the groups to pick with 'P'
        cellTypeId <- grep('^cell', names(dat), value=T)
        if(length(cellTypeId)>0){
			if(length(cellTypeId)>1){
				bringToTop(-1)
				cat('\n Select the cell type to load in \n')
				cellTypeId <- select.list(cellTypeId, title="Select Cell Type")
            }
        }
        cell.type.names <- names(dat[[cellTypeId]])
        cell.types <- dat[[cellTypeId]]
    }else{
        cell.type.names <- names(ct.names)
        cell.types <- ct.names
    }
    
    for(z in 1:length(pulse)){
            for(x in 1:length(cell.type.names)){ 
                #first count the number of cells in the cell type group
                ct.sum[as.character(dat.name),cell.type.names[x]]<-length(cell.types[[ cell.type.names[x] ]])
                #sum the collumn with only the cell.types defined rows based on the current selected collumn
                ct.sum[pulse[z],cell.type.names[x]]<-sum(dat$bin[cell.types[[ cell.type.names[x] ]],pulse[z]])
            }
    }

    if(save){
        print('Enter you file name without spaces')
        save.names <- scan(n=1, what='character')
        print(paste(save.names,'xlsx',sep=''))
        if(xlsx){
            require(xlsx)
            tryCatch(
                write.xlsx(ct.sum, file=paste(save.names,'.xlsx',sep='')), 
                error=function(e) print("You Forgot to input cells.")
            )
        }else{
            write.csv(ct.sum, file=paste(save.names,'.csv',sep=''))
        }
    }
    return(ct.sum)
}

#' I have a series of pdf files
#' @export
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
#' I have a series of pdf files
#' @export
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

#' Funciton to save the work along with create a unique savehistory
#' @export
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

#' @export
census_viewer <- function(dat){
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

#' Function to rename experiment to the name of the folder it
#' resides
#' @export
renamer <- function(){
    expName <- rev(strsplit(getwd(), "/")[[1]])[1]
    expName <- strsplit(expName, " ")[[1]][1]
    expName <- paste0("RD.", expName)

    print(ls(pattern = "^RD[.]"))
    expToRename <- get(ls(pattern = "^RD[.]", envir=.GlobalEnv))

    assign(expName, expToRename)

    save(list=expName, file=paste0(expName,".Rdata") )
}

