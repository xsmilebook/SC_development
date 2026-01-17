#Non-Linear LongCombat Function 
#Load Libraries as Dependencies (Also Depends on NeuroCombat.R script)
library(invgamma)
library(lme4)
library('pbkrtest')
library('mgcv')
library(gamm4)
library('matrixStats')
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)


nlongcombat=function(features,knot,ranef,X_nlin,X_lin,batch){
stand=c()
subid=ranef
maxi=matrix(0,nrow=nrow(features),ncol=ncol(features))
Fitted_site=matrix(0,nrow=nrow(features),ncol=ncol(features))
Residual=matrix(0,nrow=nrow(features),ncol=ncol(features))

for(i in 1:ncol(features)){
y=features[,i]
Xe=cbind(y,X_nlin,X_lin,batch)
Xe=as.data.frame(Xe)

e=""
for(j in 1:ncol(X_nlin)){
  e=paste(e,"+","s(",colnames(X_nlin)[j],",bs ='cr'",", k=",knot,")",sep="")
}
e=substring(e,2)

#Do the same for linear coefficients 
f=""
for(z in 1:ncol(X_lin)){
  f=paste(f,"+",colnames(X_lin)[z],sep="")
}
f=substring(f,2)

yur=gamm4(as.formula(paste('y','~',"as.factor(batch)",
                           "+",e,"+",f,"+1",sep='')),
          random=~(1|subid),data=data.frame(Xe))

maxi[,i]=as.matrix(yur$gam$fitted.values) 
Fitted_site[,i]=model.matrix(yur$gam)[,2]%*%as.matrix(yur$gam$coefficients[2])
stand[i]=as.data.frame(VarCorr(yur$mer))[3,5]
Residual[,i]=(y-(maxi[,i]-Fitted_site[,i]))/(as.data.frame(VarCorr(yur$mer))[3,5])
print(paste0('Fitting Feature',i,sep=' '))
}

data.harmonized=neuroCombat(t(Residual),as.factor(Xe$batch))
dat=t(data.harmonized$dat.combat)

for(i in 1:20){
  dat[,i]=dat[,i]*stand[i]
  #print(i)
}
harmonized=dat+(maxi-Fitted_site)

return(list(harmonized=harmonized))
}


#---------------------------------
#Example Input and Function
#---------------------------------
#knot = 45
#ranef=simdata$subid
#X_nlin=as.matrix(simdata$age)
#colnames(X_nlin)='age'
#X_lin=cbind(simdata$time)
#colnames(X_lin)=c('time')
#batch=simdata$batch
#features=simdata[,6:25]
#colnames(features)=colnames(simdata)[6:25]


#dr=nlongcombat(features,knot,ranef,X_nlin,X_lin,batch)
  
  
  
  

  
  
