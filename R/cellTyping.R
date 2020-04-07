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

#' This is the old way to create a census! Recently updated though
#' @export
census.brewer<-function(dat){
    cell.types<-dat$cell.types


    dev.new(width=10, height=5)
    stacked.traces<-dev.cur()
    LinesEvery.5.1(dat, sample(row.names(dat$c.dat)[1:5]), plot.new=F, lmain="WAZZZUPPPP", t.type="t.dat", img=dat$img1)

    cat("HOWDY PARTNER How Many groups to census?\n")
    tryCatch(bringToTop(-1), error=function(e)NULL)
    group.number<-scan(n=1, what='numeric')

    cat("\nEnter the names of your census groups seperated by '.'\n")
    census.names<-scan(n=as.numeric(group.number), what='character')
    dev.off(stacked.traces)


    selected.cell.groups<-select.list(names(cell.types), title="Select groups to census", multiple=T)
    cat("\nThese are the cells you have chosen\n")
    print(selected.cell.groups)
    census<-list()

    for(i in 1:length(selected.cell.groups))
    {
        print(selected.cell.groups[i])
        
        if(length(cell.types[[selected.cell.groups[i]]])>1){
            census[[i]]<-tcd(dat, Reduce(union,cell.types[[selected.cell.groups[i]]]), save_question=F)
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

        if(length(UL) > 0){
            UL_groups <- tcd(tmp_rd, c(UL), save_question=F)
            L1<-UL_groups[[1]] #proprioceptor	
            L2<-UL_groups[[2]] #jagged
            L3<-UL_groups[[3]] #IDE only
            L4<-UL_groups[[4]] #no effect
            if(edit_ct){discard<-union(discard, UL_groups[[5]])}
         }else{
            L1 <- NULL
            L2 <- NULL
            L3 <- NULL
            L4 <- NULL
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
                G_0_sort<-tcd(tmp_rd, c(G_0), save_question=F)
                L5<-G_0_sort[[1]]
                L6<-G_0_sort[[2]]
                if(edit_ct){discard<-union(discard,G_0_sort[[3]])}
            }else{
                L5 <- NULL
                L6 <- NULL
            }
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
