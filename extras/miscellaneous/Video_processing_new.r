#what are the names to import
img_name_vec<-c(
	"bf.gfp.tritc.start.png",
	"gfp.tritc.dapi.end.ci.ltl.rs.png",
	"gfp.start.ci.ltl.rs.png",
	"tritc.start.ci.ltl.rs.png",
	"bf.start.lab.png",
	"fura2.png",
	"fura2.divide.start.png",
	"roi.img.png")

pharming_harvest(
	main_dir="Y:/Cris Urcino/190505", #where the experiments are located
	img_name_vec = img_name_vec, #image names to upload
	image_question = F)
