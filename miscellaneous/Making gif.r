# I have a series of pdf files
require(magick)

pdfs_for_gif<-select.list(list.files(pattern="pdf"),multiple=T)
dense<-200
gif<-image_read_pdf(pdfs_for_gif[1],density=dense)
gif<-image_border(gif,"red","10x10")
for(i in 2:length(pdfs_for_gif)){
    
    gifz<-image_read_pdf(pdfs_for_gif[i],density=dense)
    gifz<-image_border(gifz,"black","10x10")
    gif<-c(gif,gifz)
}   

animation<-image_animate(gif,fps=2)

image_write(animation,"yo.gif")