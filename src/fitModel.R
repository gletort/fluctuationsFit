rm( list=ls() )
library(minpack.lm)

### Load the calculated fft by condition and fit the Heilrich model (see paper from Almonacid et al. 2019)

conds = c("MDA231", "MCF7", "MCF10A")
colors = c("red", "blue", "green")

datas = c()
fits = c()

png("fft_by_file.png")
plot( 0.1, 0.1, xlim=c(1,180), ylim=c(0.00001,1.1), xlab="Mode", ylab="|fft(r-R)|^2", type="n", log="xy" )
for (cond in conds )
{
    data = read.csv( paste( cond, "results_fft.csv", sep="/"), header=T, sep=";" )
    #data[,2:length(data[1,])] = data[,2:length(data[1,])]^2
    for (i in 2:length(data[1,]))
    {
        points( data[,1], data[,i], type="l", col=colors[which(conds==cond)] )
    }
    avg_data = cbind( data[,1], apply( data[,2:length(data[1,])], 1, mean ) )
    colnames(avg_data) = c("Mode", "Mean")
    avg_data = as.data.frame(avg_data)
    
    ## Fit the Heilrich model
    start = 2
    subavg_data = avg_data[start:180,]
    fit.c = nlsLM( Mean ~ b / (a * Mode^2 + 1* Mode^4 ) + c, data=subavg_data, start=list(a=0.01, b=0.001, c=0.01),, lower=c(0,0,0), upper=c(Inf,Inf,Inf), control=nls.lm.control(maxiter=1000, maxfev = 10000) )
    print(fit.c)
    fits = rbind( fits, coef(fit.c) )
    
    if (length(datas)==0)
    {
        datas = avg_data
    } else {
        datas = merge( datas, avg_data, by="Mode" )
    }
}
dev.off()


png("fft_avg.png")
plot( 0.1, 0.1, xlim=c(1,180), ylim=c(0.00001,0.5), xlab="Mode", ylab="|fft(r-R)|^2", type="n", log="xy"  )
for (i in 2:length(datas[1,]))
{
    fitdata = fits[i-1,2] / (fits[i-1,1] * datas[,1]^2 + 1* datas[,1]^4 ) + fits[i-1,3]
    points( datas[,1], datas[,i], type="l", col=colors[i-1] )
    points( datas[start:length(fitdata),1], fitdata[start:length(fitdata)], type="l", col=colors[i-1], lty=2 ) 
}
dev.off()

#coefs = conds
coefs = cbind(fits[,1]/fits[,2], 1/fits[,2], fits[,3], sqrt(1/fits[,1])) 

coefs = cbind( cond=conds, coefs )
colnames(coefs) = c("celltype", "sigma", "k", "Yinf", "sqrt(k/sigma)")
write.csv( coefs, "helfrich_coefs.csv", row.names=F, sep=";" )
