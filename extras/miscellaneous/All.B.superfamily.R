fn1 <- list.files(path="~/Box Sync/Transcriptomes/transcriptome",pattern="anot.cono.pep.uniq.fa",recursive=T,full.names=T)
length(fn1)
tn1 <- sub(".*/","",fn1)
folder1 <- sub("^.*/","",sub("/[^/]*$","",fn1))
gstr <- "supFam:B\t"
gstr <- "conotoxin.B\t"
for(i in fn1)
{
	sysstr <- paste("cat \'",i,"\' | grep -A 1 \'",gstr,"\' >> B.pep.out.fa",sep="")
	system(sysstr)
}
#cat B.pep.out.fa | grep -v "^--" > B.pep2.out.fa

library(seqinr)
b.pep <- read.fasta("B.pep2.out.fa",seqtype="AA",strip=F)
b.str <- sapply(b.pep,c2s)
b.str[1:5]
length(unique(b.str))
length(b.str)
b.i <- match(unique(b.str),b.str)
length(b.i)
b.pep <- b.pep[b.i]

write.fasta(b.pep,names(b.pep),file.out="b.pep.unique.fa",open="w",nbchar=100000)

#clustalo -i b.pep.unique.fa -t Protein --infmt fasta --outfmt clustal --force -o tmp.align
tmp.align <- read.alignment("tmp.align","clustal",forceToLower=F)
tmp.dist <- dist.alignment(tmp.align)
tmp.dist[is.na(tmp.dist)] <- 1
tmp.hc <- hclust(tmp.dist)
plot(tmp.hc)
tmp.cut <- cutree(tmp.hc,h=0)
#this doesn't take into consideration the differences in peptide length.
#instead just order by hc and write it out.
b.pep <- b.pep[tmp.hc$order]
write.fasta(b.pep,names(b.pep),file.out="b.pep.order.unique.fa",open="w",nbchar=100000)

sort(table(tmp.cut))
sapply(b.pep[tmp.cut==76],c2s)
b.i <- match(unique(tmp.cut),tmp.cut)
length(b.i)
b.pep2 <- b.pep[b.i]
length(b.pep2)
length(tmp.cut)
names(tmp.cut)

write.fasta(b.pep2,names(b.pep2),file.out="b.pep2.unique.fa",open="w",nbchar=100000)
#clustalo -i b.pep2.unique.fa -t Protein --infmt fasta --outfmt clustal --force -o tmp2.align
tmp2.align <- read.alignment("tmp2.align","clustal",forceToLower=F)
tmp2.dist <- dist.alignment(tmp2.align)
tmp2.dist[is.na(tmp2.dist)] <- 1
tmp2.hc <- hclust(tmp2.dist)
plot(tmp2.hc)
tmp2.cut <- cutree(tmp2.hc,h=0)
sort(table(tmp2.cut))
sapply(b.pep2[tmp2.cut==13],c2s)
