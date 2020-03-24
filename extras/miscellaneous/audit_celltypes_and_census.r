#FIND COLLUMN OF INTEREST and CELLS OF INTEREST
cells_responded <- cellzand(tmp_rd$bin)
cell_types_names <- select.list(names(tmp_rd$cell_types),multiple=T)
for( j in 1:length(cells_responded) ){
    counter <- 0
    cell_types_combine <- c()
    for(i in 1:length(cell_types_names)){
        present_in_ct <- tmp_rd$cell_types[[ cell_types_names[i] ]][ cells_responded[j]==tmp_rd$cell_types[[ cell_types_names[i] ]] ]
        if( length(present_in_ct) > 0 ){
            counter <- counter + 1
            cell_types_combine <- c(cell_types_combine, cell_types_names[i])
        }
     }
     if(counter > 1){
        print(cells_responded[j])
        print(cell_types_combine)
     } 
}     
            