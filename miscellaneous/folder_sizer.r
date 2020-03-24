main_dir<-"Y:/"
setwd(main_dir)
experimentors <- list.dirs(".", recursive=F)

print(experimentors)

folder_size<-c()
file_infos <- list()
for(i in 1:length(experimentors)){
	#Show me where i am
	print(experimentors[i] )
	setwd( experimentors[i] )
	#make sure i can see what is going on
	flush.console()
	
	#Create a name
	experimentor_name <- sub("./", '', experimentors[i])
	#grab all the files from this arena
	all_files <- list.files(".", all.files=T, recursive=T, include.dirs=T)
	
	# Compute all File info
	file_infos[[experimentor_name]] <- file.info(all_files)
	#This is the size of each file in there
	file_sizes <- file_infos[[experimentor_name]]$size / 1e9
	
	#print the sum of the file sizes
	print(sum(file_sizes))
	#sum up the file sizes
	folder_size[experimentor_name] <- sum(file_sizes[!is.na(file_sizes)])
	setwd(main_dir)
}

# folder_size <- c()
# for(i in 1:length(file_infos)){
	# file_sizes <- file_infos[[i]]$size / 1e9
	# folder_size[names(file_infos[i])] <- sum(file_sizes[!is.na(file_sizes)])
# }


#folder_size <- round(folder_size, digits=0)
cbind(names(folder_size),folder_size)

as.data.frame(folder_size)