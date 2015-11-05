  ##############################################################################
   sub.set=T
   M=1e+6
   setwd("F:\\Research\\MPO\\RPackage")
   source("BoostTree.R")
   begintime=strptime(date(), "%a %b %d %H:%M:%S %Y")
   sum.tree=c()
   nu.param=0.01
   nBoost.param=5
   n=2000
   x.vec=rnorm(n)
   p=10
   s.mat=matrix(runif(p*n),ncol=p,byrow=T)
   bet=cbind(2*(sin(2*pi*s.mat[,1]))^2+exp(2*s.mat[,2]-1),2*(cos(2*pi*s.mat[,1]))^2+8*s.mat[,2]*(1-s.mat[,2]))
   y.vec=apply(cbind(1,x.vec)*bet,1,sum)+rnorm(n,0,0.5)
   plot(x.vec,y.vec,type='p',pch='.',cex=2)
################################################################################
   names.part=paste("var",1:p,sep='')
   colnames(s.mat)=names.part
   index.train=(1:n) %in% sample(1:n,round(n*0.8))
   split.nd="1"
   Leaf=rep("1",length(y.vec))
   crit.root=sum((regression(y.vec[index.train], cbind(1,x.vec[index.train]),method="linear")$residuals)^2)
   res=Boost.Tree(s.mat[Leaf==split.nd &index.train,],as.matrix(cbind(1,x.vec)[Leaf==split.nd & index.train,]),y.vec[Leaf==split.nd & index.train],names.part,nBoost=nBoost.param,nu=nu.param,mini.size=20,fit.method="linear")
   #endtime=strptime(date(), "%a %b %d %H:%M:%S %Y")
   #difftime(endtime, begintime, units='secs')
   index.pred=!index.train
   res.pred=pred.Boost(s.mat[index.pred,],as.matrix(cbind(1,x.vec)[index.pred,]),y.vec[index.pred],res,3:4,nBoost.param,nu=nu.param)
   nBoost.tune=(1:nBoost.param)[res.pred$L2.error==min(res.pred$L2.error)]
   res.tune=Boost.Tree(s.mat[Leaf==split.nd &index.train,],as.matrix(cbind(1,x.vec)[Leaf==split.nd & index.train,]),y.vec[Leaf==split.nd & index.train],names.part,nBoost=nBoost.tune,nu=nu.param,mini.size=20,fit.method="linear",ncatelow=0,ncatehigh=1)
   save.image(paste("Sim1Boost",nBoost.param,".RData",sep=''))
   endtime=strptime(date(), "%a %b %d %H:%M:%S %Y")
   difftime(endtime, begintime, units='secs')
   ## Summarize the boosting results
   summary.Boost(res,names.part)
   ## Out-sample prediction
   index.pred=!index.train
   res.pred=pred.Boost(s.mat[index.pred,],as.matrix(cbind(1,x.vec)[index.pred,]),y.vec[index.pred],res,3:4,nBoost.param,nu=nu.param)
################################################################################

