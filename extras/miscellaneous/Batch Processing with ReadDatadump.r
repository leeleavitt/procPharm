## Today i will attempt to creat an aplication that allows users to 

#1 Select experiment Folders to process
main_dir<-"Z:/Cell Picking"
setwd(main_dir)
cat("Select Experiments to Process with 
ReadDataDump.lee.2

ENTER ANY KEY TO CONTINUE")

exps_to_process<-select.list(list.dirs(), multiple=T)


#quick assesment of required docx files
for(i in 1:length(exps_to_process)){
    setwd(main_dir)
    setwd(exps_to_process[i])
    
    cat(
exps_to_process[i],
"
",
list.files(pattern=".docx"),
"
"  
    )
    setwd(main_dir)
}



for(i in 1:length(exps_to_process)){
    setwd(main_dir)
    setwd(exps_to_process[i])
    
    rd_name<-paste("RD.",strsplit(sub("./","",exps_to_process[i])," ")[[1]][1], sep="")
    
    #Select images to upload
    img_names<-list.files(pattern=".png")
    cat(
        "These are the images you can upload.

    ")
    print(img_names)
    cat("
    
        Would you like to continue? =1
        Would you like to remake your images? =2
    ")
    
    decision_to_continue<-scan(n=1)
    
    if(decision_to_continue==2){
        print("When Read press ANY key")
        scan(n=1)
        #Select images to upload
        img_names<-list.files(pattern=".png")
        cat(
        "These are the images you can upload.

        ")
        print(img_names)
        decision_to_continue=1
    }

    if(decision_to_continue==1){
        cat(
        "how many would you like to upload?
        MAXIMUN IS 8
        ")
        images_to_upload<-scan(n=1)
        
        img_list<-c()
        for(i in 1:images_to_upload){
            cat("These are the current images selected, in order
            ")
            print(img_list)
            print(paste("Select Img",i,sep=""))
            img_list[i]<-menu(img_names)
        }

        
    
    
    
