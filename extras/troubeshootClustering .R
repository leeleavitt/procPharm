dat <- get(load("Y:/Jackson/THP experiments/200803.M.31.lavage.LPS.m3.p1/RD.200803.M.31.lavage.LPS.m3.p1.Rdata"))

dat$blc<- dat$blc[1:dim(dat$t.dat)[1],]

names(pam5)

pam5$clusinfo

objs <- c()
asw <- c()
for(i in 1:20){
    print(i)
    pam5 <- pam(t(t.dat[x.min:x.max,m.names]),k=i)
    objs[i] <- pam5$objective[1]

    asw[i] <- pam5$silinfo$avg.width
}

plot(objs)

k.best <- which.max(asw)
dev.new()
plot(asw, type= "b", main = "pam() clustering assessment",
     xlab= "k  (# clusters)", ylab = "average silhouette width")
axis(1, k.best, paste("best",k.best,sep="\n"), col = "red", col.axis = "red")

