#copy a bunch of files for practice on the experiments

main.dir<-"Z:/Lee Leavitt/Pulicatin"
setwd(main.dir)

dirs_to_copy_from<-select.list(list.dirs(), multiple=T)

setwd(main.dir)
for(i in 1:length(dirs_to_copy_from)){
    setwd(dirs_to_copy_from[i])
    print(dirs_to_copy_from[i])
    print(list.files(pattern="roi.1024"))
    print(intersect(list.files(pattern="video"),list.files(pattern="nd2")))
    setwd(main.dir)
    }

#Now lets go to the harddrive i wanna copy to
second_dir<-"I:/hpc_project"
    setwd(second_dir)
for(i in 1:length(dirs_to_copy_from)){
    dir.create(dirs_to_copy_from[i])
    }

#now that i have created the directories in the new file locations
#i will begin copying the data to these new locations
setwd("Z:/Lee leavitt")
setwd(main.dir)
for(i in 1:length(dirs_to_copy_from)){
    start<-Sys.time()
    setwd(dirs_to_copy_from[i])
    file.copy(    
        paste(
            main.dir,
            sub("[.]","",dirs_to_copy_from[i])
            ,"/",
            intersect(list.files(pattern="video"),list.files(pattern="nd2"))
            ,sep="")
         ,   
        paste(
            second_dir,
            sub("[.]","",dirs_to_copy_from[i])
            ,"/",
            intersect(list.files(pattern="video"),list.files(pattern="nd2"))
            ,sep="")
        ,overwrite=T
    )
    end<-Sys.time()
print(paste("completed", end-start))
setwd(main.dir)   
}


 
#setwd(dir_exp)
roi.file<-list.files(pattern="roi.1024")



file.copy(