#performs a non-linear least squares calibration on a subset of data, then validation on a subset of witheld data.

#x is a matrix of dependent variables (for a single dependent variable, linearize and use cal.val)
#y is a vector of the independent variable
#prct.cal is the fraction of data to use in the calibration (between 0 and 1)
#n.it is the number of random assignments to calibration/validation datasets
#f is a function with the form of the model. For example, function(x1,x2,a,b1,b2) {a+exp(b1*x1)+b2*x2}
#s is a vector of starting values. The length of s equals ncol(x)

cal.val.nls<-function(x,y,prct.cal,n.it,f,s) {
	
	cal.valRMSE<-rep(0,n.it)
	cal.valCE<-rep(0,n.it)
	cal.valr2<-rep(0,n.it)
	cal.valBIC<-rep(0,n.it)
	cal.valint<-rep(0,n.it)
	cal.valcoeff<-matrix(,nrow=n.it,ncol=ncol(x))
	
	for (j in 1:n.it) {
	n.cal<-round(length(y)*prct.cal)															#number of records in calibration period
	n.val<-length(y)-n.cal
	rndm<-sample(seq(1:length(y)))																#randomly order data
	cal.i<-y[rndm[1:n.cal]]																		#select the first n.cal for calibration
	val.i<-y[rndm[(n.cal+1):length(y)]]
	cal.d<-x[rndm[1:n.cal],]
	val.d<-x[rndm[(n.cal+1):length(y)],]
	pred<-rep(0,length(val.i))																	#empty matrix to be filled with predicted values
	
	cal.dat<-data.frame(cal.d,cal.i)
	if (length(which(is.na(rowSums(cal.dat))=="TRUE"))>0) {
		cal.dat<-cal.dat[-which(is.na(rowSums(cal.dat))=="TRUE"),]
	}
	val.dat<-data.frame(val.d,val.i)
	if (length(which(is.na(rowSums(val.dat))=="TRUE"))>0) {
		val.dat<-val.dat[-which(is.na(rowSums(val.dat))=="TRUE"),]
	}
	
	if(ncol(x)==2) {
		fit<-nls(cal.i~f(X1,X2,a,b1,b2),data=cal.dat,start=c(a=s[1],b1=s[2],b2=s[3]))
		pred<-f(val.dat$X1,val.dat$X2,summary(fit)$coeff[1,1],summary(fit)$coeff[2,1],summary(fit)$coeff[3,1])
		} else {
			if(ncol(x)==3) {
				fit<-nls(cal.i~f(X1,X2,X3,a,b1,b2,b3),data=cal.dat,start=c(a=s[1],b1=s[2],b2=s[3],b3=s[4]))
				pred<-f(val.dat$X1,val.dat$X2,val.dat$X3,summary(fit)$coeff[1,1],summary(fit)$coeff[2,1],summary(fit)$coeff[3,1],summary(fit)$coeff[4,1])	
				} else {
					if(ncol(x)==4) {
						fit<-nls(cal.i~f(X1,X2,X3,X4,a,b1,b2,b3,b4),data=cal.dat,start=c(a=s[1],b1=s[2],b2=s[3],b3=s[4],b4=s[5]))
						pred<-f(val.dat$X1,val.dat$X2,val.dat$X3,val.dat$X4,summary(fit)$coeff[1,1],summary(fit)$coeff[2,1],summary(fit)$coeff[3,1],summary(fit)$coeff[4,1],summary(fit)$coeff[5,1])	
						} else {
							print("You picked more than 4 variables and need to edit the code")
						}
				}
		}
	
cal.valr2[j]<-cor(cal.dat$cal.i,predict(fit))
cal.valBIC[j]<-BIC(fit)
cal.valint[j]<-summary(fit)$coeff[1,1]
cal.valcoeff[j,]<-summary(fit)$coeff[2:(ncol(x)+1),1]
cal.valRMSE[j]<-sqrt(sum((pred-val.dat$val.i)^2)/nrow(val.dat))
cal.valCE[j]<-1-(sum((pred-val.dat$val.i)^2)/sum((mean(val.dat$val.i)-val.dat$val.i)^2))

}
cal.valsummary<-matrix(,nrow=9,ncol=2)
cal.valsummary[1,1]<-mean(cal.valr2)
cal.valsummary[1,2]<-sd(cal.valr2)
cal.valsummary[2,1]<-mean(cal.valBIC)
cal.valsummary[2,2]<-sd(cal.valBIC)
cal.valsummary[3,1]<-mean(cal.valint)
cal.valsummary[3,2]<-sd(cal.valint)
cal.valsummary[4:(ncol(x)+3),1]<-colMeans(cal.valcoeff)
cal.valsummary[4:(ncol(x)+3),2]<-apply(cal.valcoeff,2,sd)
cal.valsummary[8,1]<-mean(cal.valRMSE)
cal.valsummary[8,2]<-sd(cal.valRMSE)
cal.valsummary[9,1]<-mean(cal.valCE)
cal.valsummary[9,2]<-sd(cal.valCE)

cal.valsummary<<-cal.valsummary

}