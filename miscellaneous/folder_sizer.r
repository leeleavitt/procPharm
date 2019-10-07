main_dir<-"E:/Data"
setwd(main_dir)
experimentors <- list.dirs(".", recursive=F)

folder_size<-c()
for(i in 1:length(experimentors)){
	print(i)
	setwd( experimentors[i] )
	experimentor_name <- sub("./", '', experimentors[i])
	all_files <- list.files(".", all.files=T, recursive=T, include.dirs=T)
	file_sizes <- file.info(all_files)$size / 1e9
	folder_size[i] <- sum(file_sizes)
	names(folder_size)[i] <- experimentor_name
	setwd(main_dir)
}

folder_size <- round(folder_size, digits=0)
cbind(names(folder_size),folder_size)

as.data.frame(folder_size)