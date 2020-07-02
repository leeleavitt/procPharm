#first find the middle region of each levs to correct the scoring
xLoc <- tapply(dat$w.dat[,1], as.factor(dat$w.dat$wr1),mean)[levs]
yLoc <- rep(par('usr')[3] + yinch(.1), length(levsMiddle))
par(xpd=T)

points(xLoc, yLoc, pch=19, cex = 1.5))
identify(levsMiddle, yLocation, n=1)

