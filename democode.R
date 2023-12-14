library(dlnm)
library(tidyr)
library(lubridate)
library(readxl)
library(dplyr)
library(quantreg)
library(zoo)
library(mvmeta)
library(splines)
library(mice)
library(lmtest)
library(sandwich)

load("D:/straw demo/data.rdata")
citiesname   <- unique(demodata$city)
dlist  <- lapply(citiesname,function(x) demodata[demodata$city==x,])
names(dlist) <- citiesname
yall <- matrix(NA,length(dlist),3,dimnames=list(citiesname,paste("b",seq(3),sep="")))
Sall <- vector("list",length(dlist)) 
names(Sall) <- citiesname 

  pollutantname<-"PM2.5"
  
  for(i in seq(dlist)) {
    cat(i,"") 
    sub  <- dlist[[i]]
    time <- NULL
    for(j in 1:nrow(sub)){
      time<-cbind(time,t(difftime(sub$date[j],sub$date[1], units="days")+1))
    }
    sub$time<-as.numeric(time)
    cenrel    <- 0 
    
    suppressWarnings(cb <-      crossbasis(sub$ASOBcount,
                                           lag=7,
                                           argvar=list(fun="lin"), 
                                           arglag=list(fun="ns",knots=logknots(7,df=3,int=T))))
    
    suppressWarnings(cb.rain <- crossbasis(sub$cumu_rainfall,
                                           lag=7,
                                           argvar=list(fun="poly",degree=3),
                                           arglag=list(fun="ns",knots=logknots(7,df=3))))
    
    suppressWarnings(cb.wind <- crossbasis(sub$mean_windspeed,
                                           lag=7,
                                           argvar=list(fun="ns",knots=equalknots(sub$mean_windspeed,fun="ns",df=3)),
                                           arglag=list(fun="ns",knots=logknots(7,df=3))))
    
    suppressWarnings(cb.hum <-  crossbasis(sub$mean_humidity,
                                           lag=7,
                                           argvar=list(fun="ns",knots=equalknots(sub$mean_humidity,fun="ns",df=3)),
                                           arglag=list(fun="ns",knots=logknots(7,df=3))))
    
    suppressWarnings(cb.temp <- crossbasis(sub$mean_temperature,
                                           lag=7,
                                           argvar=list(fun="ns",knots=equalknots(sub$mean_temperature,fun="ns",df=3)),
                                           arglag=list(fun="ns",knots=logknots(7,df=3))))
    
    mod<-paste(pollutantname,"~","cb +cb.hum +cb.rain+cb.temp+cb.wind+ ns(time,42)+factor(dow)+factor(holiday)")
    model1 = glm(as.formula(mod),family = gaussian(), data=sub)
    
    yall[i,]  <- model1$coefficients[2:4]
    Sall[[i]] <- NeweyWest(model1, lag= 10, prewhite = F)[2:4, 2:4]
  }
  
  mvall  <- mvmeta(yall,Sall,method="reml")
  #summary(mvall)
  
  ## 所有城市整体火点序列中结点放置的位置
  ranges<-t(sapply(dlist,function(x) range(x$ASOBcount,na.rm=T)))

  
  cb <- crossbasis(1:10,
                   lag=7,
                   argvar=list(fun="lin"), 
                   arglag=list(fun="ns",knots=logknots(7,nk=2,int=T)))
  
  cpall <- crosspred(cb,
                     coef=coef(mvall),
                     vcov=vcov(mvall),
                     model.link="identity",
                     by=1,
                     from=0,
                     to=10,
                     cen=0,
                     cumul=T)
  
  
  cumu<-t(rbind(cpall$cumfit[2,],cpall$cumlow[2,],cpall$cumhigh[2,],cpall$cumfit[11,],cpall$cumlow[11,],cpall$cumhigh[11,]))

  per10point<-paste0(format(round(cumu[,4],2),nsmall=2),' (',
                  format(round(cumu[,5],2),nsmall=2),',',
                  format(round(cumu[,6],2),nsmall=2),')')
  #write.csv(ci,"RESULTS.csv"))

#plot##### 
  #tiff("3DPLOT.tiff"),width=12,height=12,unit="cm",res=300)
  plot(cpall, "3d", col="#0099FF", theta=400, phi=20 ,
       xlab="Firepoint counts", 
       ylab="Lag (days)",
       zlab="", 
       main="",
       cex.axis=1.5,
       cex.lab=1.5,
       cex.main=2.5)
  #dev.off()
  
  #tiff("SLICEPLOT.tiff"),width=15,height=12,unit="cm",res=300)
  plot(cpall,var=c(10),col="#E64B35E5",
       xlab="",ylab="",
       main="", lwd=3,cex.axis=1.5,cex.lab=1.5,cex.main=2.5)
  title(ylab="", line=2, cex.lab=1.5)
  title(xlab="Lag (days)", line=2, cex.lab=1.5)
  #dev.off()

