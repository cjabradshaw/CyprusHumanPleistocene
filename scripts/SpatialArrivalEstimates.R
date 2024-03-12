## Estimating arrival dates following the method of Solow (1993; http://www.jstor.org/stable/1940821);
## Principle = estimate a weighted likelihood, with weights depending on the distance between the focal point and the observations.
## The basic form is exponential, with the parameters regulated by calculations of bias and variance, as well as minimisation of the 
## integrated mean squared error (IMSE) using double kernel process
#
## full method is described in more mathematical detail in: Saltré et al. 2019 doi:10.1038/s41467-019-13277-0
#
## Biological hypothesis from Solow (1993): = age of the dated specimen in x is uniform between
## 0 and T(x) = arrival in x

## Explain main functions & their outputs
## solow.fct = function to estimate extinction/arrival dates at points (xx,yy)
##        outputs:  
##              - Test = estimated date
##              - bias = estimated bias
##		          - sd = estimated standard deviation
## sg.fct = function to calculate coordinates (xx,yy) of the points of interest
##              - nx = manages the departure grid step in longitude and latitude (grids nx x nx)
##              - dmmax = minimum distance between two points
##              - eliminates the points that fall on a grid cell with water and the points that are at least dmmax from another point
##              - nbcoeur = number of requested kernels

## required R libraries
library(maps)
library(doParallel)
library(foreach)
library(Imap)

# define functions
solow.fct <- function(xx=xx1,yy=yy1,Lon=Lon,Lat=Lat,Age=Age,SdAge=SdAge,nbcoeur=nbcoeur)
{
  cl <- makeCluster(nbcoeur,outfile="")
  registerDoParallel(cl)
  paramsig <- lm(log(SdAge)~Age)$coef ## pour simu plus loin
  distan <- matrix(0,length(xx),length(Lon))
  for( i in 1:length(xx)) {
  for(j in 1:length(Lon)) {
  distan[i,j] <- gdist(xx[i],yy[i],Lon[j],Lat[j],units="km")
   }}
  distan[is.na(distan)]<- max(distan,na.rm=T)
  iappartenance <- foreach (i = 1:ncol(distan),.combine='c') %dopar%
       { order(distan[,i])[1]}
  stopCluster(cl) 
  
  estim.fct <- function(yobs,sdobs,distan,pas)
   {
    lmv.fct <- function(tmax,yobs=yobs,sdobs=sdobs,pond=pond)
     {
      if(tmax<0){return(10**6)}
      yobs <- yobs/10000
      sdobs <- sdobs/10000
      lm0 <- pnorm(yobs-tmax,mean=0,sd=sdobs,lower.tail=F,log.p=F)
      lm1 <- pnorm(yobs,mean=0,sd=sdobs,lower.tail=F,log.p=F)
      lm2 <- lm0-lm1
      lm2[lm2<=0] <- 10**(-150)
      lm2[is.na(lm2)] <- 10**(-150)
      lm2[!is.finite(lm2)] <- 10**(-150)
      lm2 <- log(lm2)
      lmv <- sum(pond)*log(tmax)-sum(pond*lm2)
      lmv[is.na(lmv)] <- 10**6
      if(!is.finite(lmv)){lmv <- 10**6}
      return(lmv)
     }
   Test <- NULL
     for(i in 1:nrow(distan))
     {
      pond   <- pas*distan[i,]/max(distan)
      pond0 <- exp(-pond*pond)
      pond0[pond*pond>300] <-0
#      uu <- nlm(lmv.fct,max(yobs)/10000,yobs=yobs,sdobs=sdobs,pond=pond0,
#         stepmax=0.5,iterlim=1000)
#      #print(c(i,uu$est))
#      vv <- nlm(lmv.fct,uu$est,yobs=yobs,sdobs=sdobs,pond=pond0,iterlim=1000)
#      #print(c('done', i))
   vv <- optimize(lmv.fct,c(0,2*max(yobs)/10000),yobs=yobs,
          sdobs=sdobs,pond=pond0)
   vv$est <- vv$min	
      Test <- c(Test,vv$est)
     }
     Test <- Test*10**4
     return(Test)
    }

   simu.fct <- function(Test,iappartenance,paramsig)
    {
     ysim <- NULL
     for(i in 1:length(iappartenance))
       {
       Tm <- Test[iappartenance[i]]
       ysim <- c(ysim,Tm*runif(1))
       }
     sdsim <- exp(paramsig[1]+paramsig[2]*ysim)
     ysim <- ysim+rnorm(length(ysim),rep(0,length(ysim)),sdsim)
     return(list(ysim=ysim,sdsim=sdsim))
    }

  estpond.fct <- function(Test,iappartenance,paramsig,distan,pas,estim.fct)
   {
    Tsim <- NULL
    for(i in 1:100)
    {
     cdsim <- simu.fct(Test,iappartenance,paramsig)
     Tsim <- cbind(Tsim,estim.fct(cdsim$ysim,cdsim$sdsim,distan,pas))
    }
    msim <- apply(Tsim,1,mean)
    sdsim <- apply(Tsim,1,sd)
    biais <- msim-Test
    mse <- biais*biais+sdsim*sdsim
    imse <- sum(mse)
    return(list(Test=Test,biais=biais,sd=sdsim,imse=imse,pas=pas))
  }

  Test <- estim.fct(Age,SdAge,distan,5)
   cl <- makeCluster(nbcoeur,outfile="")
  registerDoParallel(cl)
  voir <- foreach(iteration=1:100) %dopar% {
    prop <- iteration
    #print(c(iteration,prop))
    estpond.fct(Test,iappartenance,paramsig,distan,prop,estim.fct)
   }
  stopCluster(cl)
  iv <- NULL
  for(i in 1:100){iv <- c(iv,voir[[i]]$imse)}
  pas <- voir[[order(iv)[1]]]$pas
  Testfinal <- estim.fct(Age,SdAge,distan,pas)
  result <- estpond.fct(Testfinal,iappartenance,paramsig,distan,pas,estim.fct)
  result$pas <- pas
  result$iv <- iv
  result$voir <- voir
  return(result)
}             ## fin solow.fct



sousgrille.fct <- function(nx,nbcoeur=nbcoeur)
 {
  xx1 <- rep(1:nx,nx)/(nx+1)
  yy1 <- sort(xx1)
  xx1 <- -180 + xx1*360
  yy1 <-  -59 +yy1*(74+59)
  
  #datxy<-data.frame(xx1,yy1)
  #idxy<-which((datxy$xx1 > 28) & (datxy$xx1 < 40) & (datxy$yy1 > 30) & (datxy$yy1 < 40)) 
  #xx1 <- xx1[idxy]
  #yy1 <- yy1[idxy]
  
  ## remove the point that fall into the sea, but it also remove the small islands
  ##jj <- map("world",fill=T,plot=F)
  ##selmap <- map.where(jj,xx1,yy1) 
  ##xx1 <- xx1[!is.na(selmap)]
  ##yy1 <- yy1[!is.na(selmap)]
  
  #cl <- makeCluster(nbcoeur,outfile="")
  #registerDoParallel(cl)

  ## chercher a eliminer les points tres proches
  #dd <- matrix(0,length(xx1),length(xx1))
  #dd <- foreach(i = 1:length(xx1),.combine='rbind',.packages='Imap') %dopar%  {
   #ifor <- NULL
   #for(j in 1:length(xx1))  {
    #ifor <- c(ifor, gdist(xx1[i],yy1[i],xx1[j],yy1[j],units="km"))
   #}
   #ifor
   #}
  #dd[is.na(dd)]<- max(dd,na.rm=T)
  #diag(dd) <- max(dd)
  #dm <- apply(dd,1,min)
  #dm[dm < dmmax] <- dmmax
  #dd[dd > dmmax] <- dmmax+10


  #ichang <- 1
  #while(ichang==1)
   #{
    #iu <- foreach (i = 1:length(unique(yy1)),.combine='c') %dopar% {
       #max(apply(dd[yy1 !=unique(yy1)[i],],2,min)/dm)
        #}
   #ichang <- 0
   #if(min(iu) <= 1.00000001) {
    #i1 <- order(iu)[1]
    #dd <- dd[yy1 != unique(yy1)[i1],]
    #xx1 <- xx1[yy1 != unique(yy1)[i1]]
    #yy1 <- yy1[yy1 != unique(yy1)[i1]]
    #ichang <- 1
    #}
    #}

  #ichang <- 1
  #while(ichang==1)
   #{
    #iu <- foreach (i = 1:length(xx1),.combine='c') %dopar% {
      #max(apply(dd[-i,],2,min)/dm)
      #} 
     #ichang <- 0
    #if(min(iu) <= 1.0000001) {
      #i1 <- order(iu)[1]
      #xx1 <- xx1[-i1]
      #yy1 <- yy1[-i1]
      #dd <- dd[-i1,]
      #ichang <- 1
    #}
  #}
#stopCluster(cl)
return(list(xx1=xx1,yy1=yy1))
}


Full.dat <- read.csv("Cyprus_data.csv",header=T,sep=",",dec=".",na.strings="na")

selobs <- !is.na(Full.dat$Lon) & !is.na(Full.dat$Lat)
Lon <- Full.dat$Lon[selobs]
Lat <- Full.dat$Lat[selobs]
Age <- Full.dat$Calibrated.age[selobs]
SdAge <- Full.dat$Calibrated.error[selobs]

## modification to readjust the fixed value of the minimum
minAge <- 1.01*min(Age)
Age <- Age-minAge

nx <-800
nbcoeur <- 10
sg <- sousgrille.fct(nx,nbcoeur=nbcoeur)
xx1 <- sg$xx1
yy1 <- sg$yy1

datxy<-data.frame(xx1,yy1)
idxy<-which((datxy$xx1 >= min(Lon)-2) & (datxy$xx1 <= max(Lon)+1) & (datxy$yy1 >= min(Lat)-1) & (datxy$yy1 <= max(Lat)+2)) 
xx1 <- xx1[idxy]
yy1 <- yy1[idxy]

est_solow <- solow.fct(xx=xx1,yy=yy1,Lon=Lon,Lat=Lat,Age=Age,SdAge=SdAge,nbcoeur=nbcoeur)

## save et finding the minimum
res <- as.data.frame(list(Lon=xx1,Lat=yy1,estim=est_solow$Test+minAge,biais=est_solow$biais,sd=est_solow$sd))

# write output
write.table(res,"Estimates_Cyprus(800pts)_v2.csv",sep=";")


