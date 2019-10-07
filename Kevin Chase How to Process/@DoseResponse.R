#nls estimate of dose response
#using an initial fit of logit 
#right now this is a mess don't use it!!!!
DoseResponse <- function(dose,response,plotit=T,sf=0,cex=1,rtag=NULL)
{
	require(car)
	data <- data.frame(y=response,x=dose)
	asym <- max(data[,"y"])*1.2 #initial estimate of the asymptote.
	start.coef <- coef(lm(logit((y)/asym)~x,data=data)) #logit transform seems to work but maybe something better.	

nlt <-tryCatch(nls(y~phi1/(1+exp(-(phi2+phi3*x))),start=list(phi1=asym,phi2=start.coef[1],phi3=start.coef[2]),data=data,trace=F),error=function(e) NULL)#sometimes fails so catch the error and just report the fail.

	if(is.null(nlt)){return(c(NA,NA,NA))}
	if(plotit==T)
	{
		
		x1 <- seq(min(dose)*.9,max(dose)*1.1,by=(max(dose)-min(dose))/50)
		phi1<-coef(nlt)[1] #this is the asymptote
		phi2<-coef(nlt)[2]
		phi3<-coef(nlt)[3] #rate parameter
		y1 <-phi1/(1+exp(-(phi2+phi3*x1))) 
		points(dose,response+sf,cex=cex)		
		lines(x1,y1+sf)
		text(max(dose),max(response+sf),rtag)
	}
	return(coef(nlt))
}


# # x <- p.tab[,"V4"]
# p.tab[,"MP"] <- x
# x <- c(median(x[1:3]),x,median(x[(length(x)-2):length(x)])) #add median of first 3 to start and median of last 3 to end
# mx <- (x[-c(1,2)]+x[-c(length(x),(length(x)-1))])/2 #calculate the mean of left and right neighbors at each point
# dx <- x[-1] - x[-length(x)] #delta values
# #define spikes as a
# # change of 5 or more followed by a change in the opposite direction of 5 or more 
# #remove those points and replace with the mean of the left and right neighbors
# ud <- (dx[-1] < -5 & dx[-length(dx)] > 5)  
# du <- (dx[-1] > 5 & dx[-length(dx)] < -5)  
# if(sum(ud | du) > 0){x[2:(length(x)-1)][ud | du] <- mx[ud | du];p.tab[,"MP"] <- x[2:(length(x)-1)]}

# #set some colors and the base plot
# p.col <- tapply(p.tab[,"V2"],p.tab[,"V12"],min)
# p.col <- p.col[rank(p.col)] #order the flows by time.
# p.col[] <- rainbow(length(p.col)) #rainbow colors 
# plot(p.tab[,"V2"],p.tab[,"MP"],pch=1,col=p.col[as.character(p.tab[,"V12"])],main=file.name,ylab="response",xlab="Time")
# points(p.tab[,"V2"],p.tab[,"V4"],pch=3,cex=.5) #show the raw data for reference.
# legend("bottomright",legend=names(p.col),pch=16,col=p.col,cex=.5)

# #loop over all the flows and estimte nolinear fits 
# for(i in names(p.col)) #maybe leave off the last flow wash?
# {
	# data <- data.frame(x=p.tab[p.tab[,"V12"]==i,"V2"],y=p.tab[p.tab[,"V12"]==i,"MP"])
	# if(nrow(data) < 5){break} #not enough data points
	# min.val <- apply(data,2,min)
	# data <- sweep(data,2,min.val,"-") #baseline the data (e.g. values in the positive quandrant)
	# #linear starting estimate
	# asym <- max(data[,"y"])*1.2 #initial estimate of the asymptote.
	# start.coef <- coef(lm(logit((y+1)/asym)~x,data=data)) #logit transform seems to work but maybe something better.	
	# nlt <-tryCatch(nls(y~phi1/(1+exp(-(phi2+phi3*x))),start=list(phi1=asym,phi2=start.coef[1],phi3=start.coef[2]),data=data,trace=TRUE),error=function(e) NaN)#sometimes fails so catch the error and just report the fail.
	# if(length(nlt)==1)
	# {
		# text(min.val[1],min.val[2],"FAIL",col=p.col[i],cex=2,pos=4)
	# }
	# else
	# {
		
		# phi1<-coef(nlt)[1] #this is the asymptote
		# phi2<-coef(nlt)[2]
		# phi3<-coef(nlt)[3] #rate parameter
	
		# x1 <- c(min(data$x):(max(p.tab[,"V2"])))
		# y1 <-phi1/(1+exp(-(phi2+phi3*x))) 
		
		# predict<-data.frame(x+min.val[1],y+min.val[2])
		# lines(predict,col=p.col[i],lwd=2)
		# #save the results in the res.tab
		# nlt.sum <- summary(nlt)
		# ki <- match(k,res.tab[,"file.name"])
		# res.tab[ki,"phi1"] <- phi1
		# res.tab[ki,"phi2"] <- phi2
		# res.tab[ki,"phi3"] <- phi3				
		# res.tab[ki,"phi1.se"] <- coefficients(nlt.sum)[1,2]
		# res.tab[ki,"phi2.se"] <- coefficients(nlt.sum)[2,2]
		# res.tab[ki,"phi3.se"] <- coefficients(nlt.sum)[3,2]				
	# }
# }

# #close the pdf figure for output
# dev.off()
# }

