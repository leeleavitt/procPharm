#' create a list that uses the names input for the names in the list
#' @export
named.list<-function(...){
    bob<-list(...)
    names(bob)<-as.character(substitute((...)))[-1]
    return(bob)
}


cell.ti<-function(dat, x.names, img=NA){
    graphics.off()
    dev.new(width=15, height=5)
    PeakFunc5(dat, x.names)
    if(is.null(img)){img<-dat$img1}else{img<-img}
    cell.view(dat,x.names,img)
    multi.pic.zoom(dat, x.names, img, zf=80)
}

#' This is the old way to create a census! Recently updated though
#' @export
census.brewer<-function(dat){

    cellTypesOptions <- grep('^cell([_]|[.])types$', names(dat), value=TRUE)
    if(length(cellTypesOptions) > 1){
        correctCellTypesName <- select.list(cellTypesOptions, title='Select the correct Cell Types')
    }else{
        correctCellTypesName <- cellTypesOptions
    }
    cell.types <- dat[[correctCellTypesName]]


    dev.new(width=10, height=8)
    stacked.traces<-dev.cur()
    LinesEvery.5(dat, sample(row.names(dat$c.dat))[1:10], plot.new=F, lmain="WAZZZUPPPP", t.type="t.dat", img='img1')

    cat("HOWDY PARTNER How Many groups to census?\n")
    tryCatch(bringToTop(-1), error=function(e)NA)
    group.number<-scan(n=1, what='numeric')

    cat("\nEnter the names of your census groups seperated by '.'\n")
    census.names<-scan(n=as.numeric(group.number), what='character')
    dev.off(stacked.traces)


    selected.cell.groups<-select.list(names(cell.types), title="Select groups to census", multiple=T)
    cat("\nThese are the cells you have chosen\n")
    print(selected.cell.groups)
    census<-list()

    for(i in 1:length(selected.cell.groups)){
        print(selected.cell.groups[i])
        
        if(length(cell.types[[selected.cell.groups[i]]])>1){
            census[[i]]<-tcd(dat, Reduce(union,cell.types[[selected.cell.groups[i]]]), save_question=F, track = F)
            names(census[[i]])<-census.names
        }else{
            census[[i]]<-NA
        }	
    }
    print(names(census))
    print(selected.cell.groups)

    names(census)<-selected.cell.groups

    dat$census<-census

    dat <- census_to_table(dat)
    return(dat)
}

##############################################################
#Function with 3 options. Edit_ct, classify UL , classify thermos
#This follows Marios scheme for classifying our cell types
#########################################
##############################################################
#' Function with 3 options_ Edit_ct, classify UL , classify thermos
#' This follows Marios scheme for classifying our cell types

#' @param edit_ct Logical, if true each cell class will be double checked
#' @param UL_classify default is T If TRUE then classify large diameter cells
#' @param GFP logical, default is T if TRUE then classify green cells
#' @param cell_types list input.  This is mainly used if the large cell types have already been classified. If so then then the large cell types are passed straight to the cell_types
#' @export
Cell_Typer_2<-function(tmp_rd, edit_ct=F, UL_classify=T, GFP=T, cell_types=NA){
    
    if(is.na(cell_types)){
        large_cell_types_names <- NA
    }else{
        #If your cell_types is not NA do large celltyping
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

    #selected bin and dropped
    dropped<-cellzand(tmp_rd$bin,"drop",1)

    cat("Select The response that coorespond to Neurons,
        ex_
        K+_40mM, and capsaicin_300nM
    ")
    #identfy Neurons
    neurons <- cellzand(tmp_rd$bin, , 1)
    #Remove dropped cells from he neuron class
    neurons <- setdiff(neurons, dropped)
    
    #Identify green cells_ Corrected with ROIReview
    if(GFP){
        green_cells<-cellzand(tmp_rd$bin, "gfp.bin" ,1, neurons) #selected bin then gfp_bin
    }
    
    #identify red cells
    ib4_label <- grep("cy5|tritc", names(tmp_rd$bin), value=T)
    red_cells <- cellzand(tmp_rd$bin,ib4_label ,1, neurons)
    
    #define Unlabeled cells as not green or red labeling
    if(GFP){
        unlabeled <- setdiff(neurons, green_cells)
    }else{
        unlabeled <- neurons
    }
    unlabeled<-setdiff(unlabeled, red_cells)

    # Identify capsaicin Response
    cap_cells<-cellzand(tmp_rd$bin,  
        grep("cap",names(tmp_rd$bin),ignore.case=T, value=T)[length(grep("cap",names(tmp_rd$bin),ignore.case=T, value=T))],
        1,
        neurons 
    )

    # Identify AITC Response
    aitc_cells<-cellzand(tmp_rd$bin,  
        grep("aitc",names(tmp_rd$bin),ignore.case=T, value=T)[length(grep("aitc",names(tmp_rd$bin),ignore.case=T, value=T))],
        1,
        neurons
    )

    # Identify Menthol Responses
    menth_cells<-cellzand(tmp_rd$bin, 
        grep("men",names(tmp_rd$bin),ignore.case=T, value=T)[length(grep("men",names(tmp_rd$bin),ignore.case=T, value=T))], 
        1,
        neurons
    )

    #Remove aitc responders to find trpm8 only cells
    menth_only<-setdiff(menth_cells, aitc_cells)
    
    #Find AITC and capsaicin
    aitc_and_caps <- intersect(aitc_cells, cap_cells)
    
    #define large cells as larger that 330uM^2
    large_cells_330<-cellzand(tmp_rd$c.dat,"area" ,330, neurons)
    
    #define glia is a very weak way_ Antyhing that isnt a
    #neuron is considered glia
    glia<-setdiff(tmp_rd$c.dat$id, neurons)
    glia<-setdiff(glia, dropped)

    cell_types <- named.list(neurons, glia)
    discard<-c()

    ####################
    #GREEN Group
    #Sort green cells first by capsaicin then aitc
    if(GFP){
        # G7: gfp+, Menthol + only
        G7<-intersect(green_cells, menth_cells)
        #remove aitc responders
        G7<-setdiff(G7, aitc_cells)
        #Create G7_capsaicin cells
        G7_c<-intersect(G7,cap_cells)
		G7_m <- setdiff(G7, G7_c)

        # G8: gpf+, menthol negative, capsaicin only 
        G8<-intersect(green_cells, cap_cells)
        G8<-setdiff(G8, aitc_cells)
        G8<-setdiff(G8, menth_cells)

        # G9: gfp+, menthol negative, aitc and capsaicin
        #first discover cells that are positive for caps and aitc
        #now intersect the green cells with a+ c+
        G9<-intersect(green_cells, aitc_and_caps)
        
        # G10: gfp+, AITC positive only
        G10<-intersect(green_cells, aitc_cells)
        #remove capsaicin from this group
        G10<-setdiff(G10, cap_cells)
        
        #now create a group of green responding cell that are
        #not classified by th previous green groups
        #This groups contains miscored Menthol responses and the large cell groups
        G_0<-setdiff(green_cells, c(G7,G8,G9,G10))
    }

    ########################################
    #RED ONLY GROUP
    ########################################
    # remove any green from red cells
    if(GFP){
        red_only <- setdiff(red_cells,green_cells)
    }else{
        red_only <- red_cells
    }
    
    # R11: IB4 only, Capsaicin only
    R11 <- intersect(red_only, cap_cells)
    #remove AITC from this group
    R11 <- setdiff(R11, aitc_cells)

    # R12: IB4 only, Capsaicin and AITC
    R12 <- intersect(red_only, aitc_and_caps)

    # R13: IB4 only,AITC only
    R13 <- intersect(red_only, aitc_cells)
    # remove capsaisin responses from this group
    R13<-setdiff(R13, cap_cells)
    
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
    
    # find aitc max stat
    aitc_stat<-intersect(
        grep(".max", names(tmp_rd$scp), ignore.case=T, value=T),
        grep("aitc",names(tmp_rd$scp),ignore.case=T, value=T)
    )
    
    # Find menthol max stat
    menth_stat<-intersect(
        grep(".max", names(tmp_rd$scp), ignore.case=T, value=T),
        grep("men",names(tmp_rd$scp),ignore.case=T, value=T)
    )
    
    # N15_a: N15 with aitc. Menthol response must be larger than AITC
    # Sloppy way of obtaining the menthol stat if something else has men in it
    menth_stat <- menth_stat[length(menth_stat)]
    
    #find trpm8 and trpa1 containing neurons
    menth_and_aitc_cells <- intersect(menth_cells, aitc_cells)

    # If these exist
    if(length(menth_and_aitc_cells) > 0 ){
        # extract them
		menth_and_aitc_neurons <- intersect(neurons, menth_and_aitc_cells)
		trpm8_trpa1<-c()
		for(i in 1:length(menth_and_aitc_neurons)){
            # Is the menthol Response bigger than the AITC response
            menthGreaterAitcLogic <- tmp_rd$scp[menth_and_aitc_neurons[i],menth_stat] >= ((tmp_rd$scp[menth_and_aitc_neurons[i],aitc_stat])*1.1)
			if(menthGreaterAitcLogic){
				trpm8_trpa1 <- c(trpm8_trpa1, menth_and_aitc_neurons[i])
			}
		}
		N15_a <- trpm8_trpa1
	}else{
		N15_a <- NA
	}
    ################################
    #Unlabeled
    ################################
    # Unlabeled smaller neurons responding to menthol and not AITC
    N15 <- menth_cells
    if(GFP){
        N15 <- setdiff(N15, G7)
    }
    #ensure these are neurons
    N15 <- intersect(N15, neurons)
    #remove these cells from the unlabeled group
    unlabeled <- setdiff(unlabeled, N15)	
    
    #N15_c Menthol capsaicin
    N15_c <- intersect(N15, cap_cells)
    R11 <- setdiff(R11, N15_c)
    
    #Now create an unlabeled large group of cells
    UL<- intersect(large_cells_330, unlabeled)
    #ensure they are neurons
    UL <- intersect(UL, neurons)
    #remoce any capsaicin or aitc responders
    UL <- setdiff(UL, c(cap_cells,aitc_cells))
    
    #Create N13, and N16, US is a super category
    US <- setdiff(unlabeled, UL)
    US <- intersect(US, neurons)

    #N14 unlabeled  capsaicin negative
    US_noRep <- setdiff(US, c(aitc_cells, menth_cells, cap_cells))
    US_aitc_cells <- intersect(US, aitc_cells)
    US_cap_cells <- intersect(US, cap_cells)
    US_menth_cells <- intersect(US, menth_cells)

    # Unlabeled with no responses and aitc responses
    N14 <- c(US_noRep, US_aitc_cells)
    
    # N15: menthol cells
    N15 <- union(N15, US_menth_cells)
    N14 <- setdiff(N14, menth_cells)
	N15_m <- setdiff(N15, c(N15_c, N15_a))

    #N16 unlabeled, capsaicin positive
    N16 <- US_cap_cells
    N16 <- setdiff(N16, aitc_cells)
    N16 <- setdiff(N16, menth_cells)
    
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

        if(length(UL) > 0){
            UL_groups <- tcd(tmp_rd, c(UL), save_question=F, track = F)
            L1<-UL_groups[[1]] #proprioceptor	
            L2<-UL_groups[[2]] #jagged
            L3<-UL_groups[[3]] #IDE only
            L4<-UL_groups[[4]] #no effect
            if(edit_ct){discard<-union(discard, UL_groups[[5]])}
         }else{
            L1 <- NA
            L2 <- NA
            L3 <- NA
            L4 <- NA
         }

        if(GFP){
            cat(" Sort the Unlabled Large into
            1:R3J IDE
            2:no Effect
            3:Discard
            
            PRESS ANY KEY TO CONTINUE
            ")
            scan(n=1)
            if(length(G_0)>0){
                G_0_sort<-tcd(tmp_rd, c(G_0), save_question=F, track = F)
                L5<-G_0_sort[[1]]
                L6<-G_0_sort[[2]]
                if(edit_ct){discard<-union(discard,G_0_sort[[3]])}
            }else{
                L5 <- NA
                L6 <- NA
            }
        }
    }else{
        UL<-large_cells_330
        UL<-setdiff(UL, c(cap_cells, aitc_cells, menth_cells) )
        #print(UL)
    }

    
    #review the autosorted cell classes_ Remove the cells that are not part of each class and put into discard pile, press "1" to move cells to the discard pile
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
			G7_m,
            G7_c,
            G8,
            G9,
            G10
        )
        cell_types <- append(cell_types,gfp_ct)
    }else{}
    
    red_ul_ct<-named.list(
        R11,
        R12,
        R13,
        N14,
        N15,
		N15_a, 
		N15_m,
        N15_c,
        N16
    )
    
    cell_types<-append(cell_types,red_ul_ct)

    #UC: Cells in neurons that couldn't find a home
    allCells <- Reduce(union, cell_types)
    UC <- setdiff(neurons, allCells)
    UC <- named.list(UC)
    cell_types <- append(cell_types, UC)

    tmp_rd$cell_types <- cell_types
    
    for(i in 1:length(cell_types)){
        print(
            paste( 
                names(tmp_rd$cell_types)[i],
                "=",
                length( tmp_rd$cell_types[[i]] ) 
            )
        )
    }

    # # This adds the cells to c.dat
    # if(!any(names(tmp_rd) == "cell_type_model")){
    #     tmp_rd <- cellTypeAdder(tmp_rd)
    #     # This models the cells
    #     tmp_rd <- cell_type_modeler(tmp_rd)

    #     # Now lets use the neuralnetworks to correctly place R13, N14, N15
    #     toReClassify <- Reduce(union, list(tmp_rd$cell_types$R13, tmp_rd$cell_types$N14))

    #     x <- tmp_rd$cell_type_model[toReClassify]
    #     for(i in 1:length(x)){
    #         x[[i]]['gfp',] <- x[[i]]['gfp',] * 0.8
    #         x[[i]]['k40',] <- x[[i]]['k40',] * 2
    #         x[[i]]['menth',] <- x[[i]]['menth',] * .8
    #         x[[i]]['caps',] <- x[[i]]['caps',] * 0.6
    #         x[[i]]['ib4',] <- x[[i]]['ib4',] * 1.2

    #     }
    #     tmp_rd$cell_type_model[toReClassify] <- x
    # }
    # # Here i noticed we need to understand the uncertain vs certain
    # # cells. To do this we will use the kurtosis function
    # # This is now added to the scp
    # kurtosisinfo <- lapply(tmp_rd$cell_type_model, function(x){kurtosis(apply(x,2,sum))})
    # tmp_rd$scp['cell_type_kurtosis'] <- Reduce(c,kurtosisinfo)

    # # Correct the cell type
    # selectCT <- c("G8", "R13", "N14")
    # correctedCT <- sapply(tmp_rd$cell_type_model[toReClassify], function(x){
    #     rowSums <- apply(x[selectCT], 2, sum)
    #     rowMaxLogic <- rowSums == max(rowSums)
    #     names(rowSums)[rowMaxLogic]
    # })

    # # Now remove all of the names from the selectCT
    # tmp_rd$cell_types[selectCT] <- lapply(tmp_rd$cell_types[selectCT], function(x) setdiff(x,names(correctedCT)))

    # # Now add the cells to the correct class
    # toCorrectTo <- unique(correctedCT)
    # for(i in 1:length(toCorrectTo)){
    #     tmp_rd$cell_types[[ toCorrectTo[i] ]]<- union(tmp_rd$cell_types[[ toCorrectTo[i] ]], names(which(correctedCT == toCorrectTo[i] , arr.ind = T)))
    # }

    return(tmp_rd)
}

#' Updated cell_typer which incorporates neural networks to distinguish n14 vs R13
#' Now only the proprioceptors and jagged need to be classified. Everything else is automated.
#' @param dat RD.experiment
#' @param UL_classify boolean decision to classify large diameter cells
#' @export
Cell_Typer_3 <- function(dat){
    # This logic assesess whether to reclassify
    largeClassifiedAlready <- any(c("L1","L2", "L3", "L4", "L5", "L6") %in% names(dat$cell_types))

    # Is gfp present?
    gfpLogic <- grep("gfp", names(dat$c.dat))
    if(length(gfpLogic) > 0){
        gfpLogic <- T
    }else{
        gfpLogic <- F
    }

    # Is there R3J?
    r3jLogic <- grep('^[rR](3|[I]{3})[jJ].*[mM]$', names(dat$bin), value = T)
    if(length(r3jLogic) > 0){
        UL_classify <- T
    }else{
        UL_classify <- F
    }
    
    # Begin Cell typing
    cell_types <- list()
    
    allCells <- dat$c.dat$id
    # find the windowNames
    menthName <- grep("^[mM][eE][nN][tT]", names(dat$bin), value = T)[1]
    aitcName <- grep("^[aA][iI][tT][cC]", names(dat$bin), value = T)[1]
    capsName <- grep("^[cC][aA][pP][sS]", names(dat$bin), value = T)[1]
    k40Name <- grep("[kK].*40", names(dat$bin), value = T)[1]

    # Define drops
    dropLogic <- dat$bin$drop == 1
    drops <- allCells[dropLogic]

    # First define neurons 
    resps <- c(aitcName, capsName, k40Name)
    #neuronResps <- sapply(resps, function(x) grep(x,names(dat$bin), value  = T))
    neuronLogic <- apply(dat$bin[,resps] == 1, 1, any)
    neurons <- allCells[neuronLogic]
    neurons <- setdiff(neurons, drops)

    if(gfpLogic){
        # G7
        g7Logic <-  dat$bin[neurons, menthName] == 1 &
            dat$bin[neurons, "gfp.bin"] == 1
        cell_types$G7 <- neurons[g7Logic]

        g7.cLogic <- dat$bin[cell_types$G7, capsName] == 1
        cell_types$G7_c <- cell_types$G7[g7.cLogic]
        
        g7.mLogic <- dat$bin[cell_types$G7, capsName] != 1
        cell_types$G7_m <- cell_types$G7[g7.mLogic]

        # G8
        g8Logic <- dat$bin[neurons, menthName] != 1 &
            dat$bin[neurons, aitcName] != 1 &
            dat$bin[neurons, capsName] == 1 &
            dat$bin[neurons, "gfp.bin"] == 1
        cell_types$G8 <- neurons[g8Logic]

        # G9 
        g9Logic <- dat$bin[neurons, menthName] != 1 &
            dat$bin[neurons, aitcName] == 1 &
            dat$bin[neurons, capsName] == 1 &
            dat$bin[neurons, "gfp.bin"] == 1
        cell_types$G9 <- neurons[g9Logic]

        # G10
        g10Logic <- dat$bin[neurons, menthName] != 1 &
            dat$bin[neurons, aitcName] == 1 &
            dat$bin[neurons, capsName] != 1 &
            dat$bin[neurons, "gfp.bin"] == 1
        cell_types$G10 <- neurons[g10Logic]
    }

    # R11
    r11Logic <- dat$bin[neurons, menthName] != 1 &
        dat$bin[neurons, aitcName] != 1 &
        dat$bin[neurons, capsName] == 1 &
        dat$bin[neurons, "gfp.bin"] != 1 &
        dat$bin[neurons, "cy5.bin"] == 1
    cell_types$R11 <- neurons[r11Logic]

    # R12 
    r12Logic <- dat$bin[neurons, menthName] != 1 &
        dat$bin[neurons, aitcName] == 1 &
        dat$bin[neurons, capsName] == 1 &
        dat$bin[neurons, "gfp.bin"] != 1 &
        dat$bin[neurons, "cy5.bin"] == 1
    cell_types$R12 <- neurons[r12Logic]

    # R13 
    r13Logic <- dat$bin[neurons, aitcName] == 1 &
        dat$bin[neurons, capsName] != 1 &
        dat$bin[neurons, "gfp.bin"] != 1 &
        dat$bin[neurons, "cy5.bin"] == 1 &
        (dat$scp[neurons, paste0(aitcName,".max")]  >= (dat$scp[neurons, paste0(menthName,".max")] * .8))
    cell_types$R13 <- neurons[r13Logic]

    # N14 
    n14Logic <- dat$bin[neurons, menthName] != 1 &
        dat$bin[neurons, k40Name] == 1 & 
        dat$bin[neurons, capsName] != 1 & 
        dat$bin[neurons, "gfp.bin"] != 1 &
        dat$c.dat[neurons, 'area'] <= 400
        #dat$bin[neurons, "cy5.bin"] == 1
    cell_types$N14 <- neurons[n14Logic]

    # N15 
    n15Logic <- dat$bin[neurons, menthName] == 1 &
        dat$bin[neurons, "gfp.bin"] != 1 &
        (dat$scp[neurons, paste0(menthName,".max")] > (dat$scp[neurons, paste0(aitcName,".max")] * .8))
    cell_types$N15 <- neurons[n15Logic]

    n15.cLogic <-   dat$bin[cell_types$N15, capsName] == 1 &
                    dat$bin[cell_types$N15, aitcName] != 1
    cell_types$N15_c <- cell_types$N15[n15.cLogic]

    n15.aLogic <-   dat$bin[cell_types$N15, capsName] != 1 &
                    dat$bin[cell_types$N15, aitcName] == 1
    cell_types$N15_a <- cell_types$N15[n15.aLogic]

    n15.acLogic <-   dat$bin[cell_types$N15, capsName] == 1 &
                    dat$bin[cell_types$N15, aitcName] == 1
    cell_types$N15_ac <- cell_types$N15[n15.acLogic]

    n15.mLogic <-   dat$bin[cell_types$N15, capsName] != 1 &
                    dat$bin[cell_types$N15, aitcName] != 1
    cell_types$N15_m <- cell_types$N15[n15.mLogic]

    # N16
    n16Logic <- dat$bin[neurons, "cy5.bin"] != 1 &
                dat$bin[neurons, "gfp.bin"] != 1 &
                dat$bin[neurons, aitcName] != 1 &
                dat$bin[neurons, menthName] != 1 &
                dat$bin[neurons, capsName] == 1
    cell_types$N16 <- neurons[n16Logic]

    if(!largeClassifiedAlready){
        if(UL_classify){
            cat(" Sort the Unlabled Large into
                1:Propriocepters
                2:Jagged
                
                Use F1/F2 for God's sake
                press X to discard
                Press ENTER to continue
            ")
            scan(n=1)
            uLLogic <- dat$c.dat[neurons, "area"] > 400 &
                        dat$bin[neurons, "gfp.bin"] != 1 &
                        dat$bin[neurons, "cy5.bin"] != 1 &
                        dat$bin[neurons, aitcName] != 1 &
                        dat$bin[neurons, menthName] != 1 &
                        dat$bin[neurons, capsName] != 1
            UL <- neurons[uLLogic]
            levs <- unique(dat$w.dat$wr1)
            r3jLocation <- grep('[rR](3|[I]{3})[jJ].*[mM]$', levs)[1]
            beforeR3j <- r3jLocation - 1
            afterR3j <- r3jLocation + 1


            if(length(UL) > 1){
                UL <- UL[!is.na(UL)]
                largeCells <- tcd(dat, UL, save_question = F)[1:2]
                names(largeCells) <- c("L1", "L2")

                tot <- Reduce(c, largeCells)
                L3L4 <- setdiff(UL, tot)
                
                l3Logic <-  (dat$scp[L3L4, paste0(levs[afterR3j], ".max")] * 0.7) > dat$scp[L3L4, paste0(levs[beforeR3j], ".max")] 
                L3 <- L3L4[l3Logic]

                l4Logic <-  (dat$scp[L3L4, paste0(levs[afterR3j], ".max")] * 0.7) <= dat$scp[L3L4, paste0(levs[beforeR3j], ".max")]
                L4 <- L3L4[l4Logic]

                largeCells <- c(largeCells, list(L3 = L3, L4 = L4))
            }else{
                cat("There are no unlabeled Large Diameter cells to cell type.\n")
                largeCells <- list(L1 = NA, L2 = NA, L3 = NA, L4 = NA, L5 = NA, L6 = NA)
            }

            if(gfpLogic){
                # L5   
                l5Logic <-  (dat$scp[neurons, paste0(levs[afterR3j], ".max")] * 0.7) > dat$scp[neurons, paste0(levs[beforeR3j], ".max")] &
                            dat$bin[neurons, "gfp.bin"] == 1 &
                            dat$bin[neurons, aitcName] != 1 &
                            dat$bin[neurons, menthName] != 1 &
                            dat$bin[neurons, capsName] != 1
                L5 <- neurons[l5Logic]

                # L6
                l6Logic <-  (dat$scp[neurons, paste0(levs[afterR3j], ".max")] * 0.7) <= dat$scp[neurons, paste0(levs[beforeR3j], ".max")] &
                            dat$bin[neurons, "gfp.bin"] == 1 &
                            dat$bin[neurons, aitcName] != 1 &
                            dat$bin[neurons, menthName] != 1 &
                            dat$bin[neurons, capsName] != 1
                L6 <- neurons[l6Logic]
                
                largeCells <- c(largeCells, list(L5 = L5, L6 = L6, L5.split = NA, L6.split = NA))
            }
        }else{
            uLLogic <- dat$c.dat[neurons, "area"] > 300 &
                    dat$bin[neurons, "gfp.bin"] != 1 &
                    dat$bin[neurons, "cy5.bin"] != 1 &
                    #dat$c.dat[neurons, 'area'] > 400 &
                    dat$bin[neurons, aitcName] != 1 &
                    dat$bin[neurons, menthName] != 1 &
                    dat$bin[neurons, capsName] != 1
            UL <- neurons[uLLogic]
            largeCells <- list(UL = UL)
            if(gfpLogic){
                uLgLogic <- dat$c.dat[neurons, "area"] > 300 &
                    dat$bin[neurons, "gfp.bin"] == 1 &
                    dat$bin[neurons, "cy5.bin"] != 1 &
                    #dat$c.dat[neurons, 'area'] > 400 &
                    dat$bin[neurons, aitcName] != 1 &
                    dat$bin[neurons, menthName] != 1 &
                    dat$bin[neurons, capsName] != 1

                UL.g <- neurons[uLgLogic]
                largeCells <- c(largeCells, list(UL.g = UL.g))
            }
        }
    }else{
        largeNames <- grep("^L", names(dat$cell_types), value = T)
        largeCells <- dat$cell_types[largeNames]
    }

    cell_types <- c(largeCells, cell_types)

    # UC
    assignedCellTypes <- Reduce(c, cell_types)
    UC <- setdiff(neurons, assignedCellTypes)

    cell_types <- c(cell_types, list(UC = UC))
    cell_types <- lapply(cell_types, function(x){x = x[!is.na(x)]; return(x)})
    dat$cell_types <- cell_types

    # Send the R13 and N14 through the neural networks to get the correct class
    selectCT <- c("R13", "N14")
    toReClassify <- Reduce(union, cell_types[selectCT])
    if( !("cellTypeModel" %in% names(dat)) ){
        # Correct my models to improve things a bit
        dat <- cell_type_modeler(dat)
        
        x <- dat$cellTypeModel[toReClassify]
        for(i in 1:length(x)){
            x[[i]]['gfp',] <- x[[i]]['gfp',] * 0.8
            x[[i]]['k40',] <- x[[i]]['k40',] * 2
            x[[i]]['menth',] <- x[[i]]['menth',] * .8
            x[[i]]['caps',] <- x[[i]]['caps',] * 0.6
            x[[i]]['ib4',] <- x[[i]]['ib4',] * 1.2
        }
        dat$cellTypeModel[toReClassify] <- x
    }

    # Correct the cell type
    correctedCT <- sapply(dat$cellTypeModel[toReClassify], function(x){
        rowSums <- apply(x[selectCT], 2, sum)
        rowMaxLogic <- rowSums == max(rowSums)
        names(rowSums)[rowMaxLogic]
    })
    # Now remove all of the names from the selectCT
    cell_types[selectCT] <- lapply(cell_types[selectCT], function(x) setdiff(x,names(correctedCT)))

    # Now add the cells to the correct class
    toCorrectTo <- unique(correctedCT)
    for(i in 1:length(toCorrectTo)){
        cell_types[[ toCorrectTo[i] ]] <- union(cell_types[[ toCorrectTo[i] ]], names(which(correctedCT == toCorrectTo[i] , arr.ind = T)))
    }

    cell_types <- c(list(neurons = neurons, glia = setdiff(allCells, neurons)), cell_types)
    dat$cell_types <- cell_types
    dat <- cellTypeAdder(dat)
    
    # Print the cell types at the end
    for(i in 1:length(cell_types)){
        print(
            paste( 
                names(dat$cell_types)[i],
                "=",
                length( dat$cell_types[[i]] ) 
            )
        )
    }

    dat$cell_types <- lapply(dat$cell_types, function(x) {setdiff(x,drops)})

    return(dat)
}

