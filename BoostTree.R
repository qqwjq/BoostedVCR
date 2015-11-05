library(MASS)
library(sfsmisc)
################################################################################
## Function for conducting linear regression and log-linear regression.
## method="linear"
## method="const.loglinear"
## given thresh.low and thresh.high
## Constant nest-level parameter.
regression=function(y,x,N.vec=0,n.vec=0,market=NA,method="linear",nEM=5000,tolern=1e-6,opt='exp',nests=NA,off.set=0,param.init=0)
{
  res=lm(y~x-1)
  if(method=="const.loglinear")
  {
   res.new=lm(log(y)~x-1)
   coefs=res.new$coef
   if(coefs[2]>thresh.high | coefs[2]<thresh.low)
   {
    #y.tmp=log(y)-x[,-1]*thresh
    #xmat.tmp=cbind(1,-x[,-1])
    #reg=constrp(xmat.tmp,y.tmp,t(c(0,1)))
    #coefs=reg$bhat
    #coefs=c(coefs[1],thresh-coefs[2])
    tmp.fn=function(bet,y.tmp,x.tmp)
    {
    return(sum((log(y.tmp)-x.tmp%*%bet)^2))
    }
    coefs=optim(c(0,0),tmp.fn,method="L-BFGS-B",lower=c(-Inf,thresh.low),upper=c(Inf,thresh.high),y.tmp=y,x.tmp=x)$par
   }
   resids=y-exp(x%*%coefs)
   res$residuals=resids
   res$coef=coefs
  }
  if(method=="const.linear")
  {
   res.new=lm(y~x-1)
   coefs=res.new$coef
   if(coefs[2]>0) coefs=c(mean(y),0)
   resids=y-x%*%coefs
   res$residuals=resids
   res$coef=coefs
  }
  if(method=="loglinear")
  {
   res.new=lm(log(y)~x-1)
   coefs=res.new$coef
   resids=y-exp(x%*%coefs)
   res$residuals=resids
   res$coef=coefs
  }
  ## choice model for a given market indicator.
  if(method=="choice")
  {
   gamm.cur=rep(0,ncol(x))
   crit.cur=NA
   ##seg,Nvec,nvec,zmat,gamm.param
   for(iEM in 1:nEM)
    {
    #cat(iEM,'\n')
    tmp.cur=exp(x%*%gamm.cur)
    p.cur=tapply(tmp.cur,market,sum)/(tapply(tmp.cur,market,sum)+0)
    xixi=mlogit.choice.IRLS(market,N.vec,n.vec,x,gamm.cur,oo=F,GAM=opt)
    cat(iEM,xixi$gamm,'\n')
    gamm.new=xixi$gamm
    gamm.new[1]=gamm.new[1]-mean(apply(x%*%gamm.new,1,sum))
    crit.cur=xixi$crit
    if(sum(abs(gamm.new-gamm.cur))<tolern) break;
    gamm.cur=gamm.new
   }
    res$residuals=NA
    res$coef=gamm.new
    res$crit=(-2)*crit.cur
    if(opt=='linear')
    {
    gamm.cur=c(mean(n.vec/N.vec),rep(0,ncol(x)-1))
    ## Likelihood criterion for choice model.
    crit.linear=function(gamm)
    {
     atr=Attraction(x%*%gamm,option=opt)$val
     prob.vec=atr/(group.fn(as.factor(market),atr,"sum")$results)
     crit=(-2)*sum(n.vec*log(prob.vec))
     return(crit)
    }
     res$coef=optim(gamm.cur,crit.linear)$par
     res$crit=optim(gamm.cur,crit.linear)$value
    }
  }
  ## Nested choice model for a given market indicator.
  if(method=="nest")
  {
   #gamm.cur=c(rep(0,ncol(x)-1))
   #tau.cur=rep(0,lengthunique(as.character(nests))-1)
   gamm.cur=as.vector(param.init[1:(ncol(x))])
   tau.cur=as.vector(param.init[-(1:(ncol(x)))])
   #if(homo) tau.cur=0
   #cat(length(c(gamm.cur,tau.cur)),'\n')
   ## Likelihood criterion for nested logit model.
   crit.nest=function(param)
   {
     cat(param,'\n')
     gamm=as.vector(param[1:(ncol(x))])
     tau=rep(as.vector(c(param[-(1:(ncol(x)))])),lengthunique(as.character(nests)))
     #tau[length(tau)]=0
     #tau=c(tau,0)
     tau.vec=rep(0,nrow(x))
     tau.vec=tau[nests]
     #cat(lengthunique(as.character(nests)),'\n')
     atr=Attraction((x%*%gamm+off.set)*exp(-tau.vec),option=opt)$val
     xixi=group.fn(as.factor(paste(market,nests,sep="-")),atr,"sum")
     #cat("okay1",'\n')
     prob.vec=atr/(xixi$results)
     I.vec=log(xixi$results)
     Qtemp.vec=exp(exp(tau.vec)*I.vec)
     Qtemp2.vec=Qtemp.vec/group.fn(as.factor(paste(market,nests,sep="-")),rep(1,length(atr)),"sum")$results
     Q.vec=Qtemp.vec/(group.fn(as.factor(market),Qtemp2.vec,"sum")$results)
     cat(summary(prob.vec),";",summary(Q.vec),'\n')
     crit=(-2)*sum(n.vec*log(prob.vec*Q.vec))
     #cat(crit,param,'\n')
     return(crit)
   }
   #nseq=20
   #param.mat=expand.grid(1,seq(-0.75,-0.25,length.out=nseq),seq(-0.1,.1,length.out=nseq),seq(-0.1,.1,length.out=nseq),0)
   #param.mat=cbind(param.mat,param.mat[,ncol(param.mat)])#,seq(-0.5,0.5,length.out=nseq))
   #cat(lengthunique(as.character(nests))-2,"\n")
   #value=c()
   #for(iseq in 1:nrow(param.mat))
   #{
   # cat(iseq/nrow(param.mat),'\n')
   # heihei=c(as.vector(unlist(param.mat[iseq,])))
    #cat(lengthunique(as.character(nests)),"\n")
    #cat(heihei,'\n')
   # value=c(value,crit.nest(heihei))
   #}
   #par(mfrow=c(2,3))
   #for(i in 1:ncol(param.mat))
   #{
   # plot(param.mat[,i],value,type='p',pch=19,cex=.5)
   #}
   cat(c(gamm.cur,tau.cur),'\n')
   res.cur=optim(c(gamm.cur,tau.cur),crit.nest,control=list(maxit=1000))
   res$coef=res.cur$par
   #res.new=res.cur$par+as.vector(rnorm(length(res.cur$par)))
   #cat("New param:",res.new,crit.nest(res.new),"\n")
   res$crit=res.cur$value
   res$conv=res.cur$convergence
   #res$val=value
  }
  return(res)
}
## examine the regression().
#x=cbind(1,rgamma(100,1,1))
#y=exp(x%*%c(2,4)+rnorm(100,0,0.1))
#thresh.high=3
#thresh.low=1
#regression(y,x,method="const.loglinear")
   #k=4
   #n=50
   #xs=sort(runif(n))
   #ys=-xs
   #ys[ys<(-0.31)]=-0.31
   #ys=ys+rnorm(n,0,.004)
   #const.fit(ys,xs,k=4)
   #plot(xs,ys,type='p')
   #lines(xs, const.fit(ys,xs,k=4)$fit,type='p',pch=19,col='red')
   #ys=c(1,2,10,11,13,17,18,19,10,5,1,1,2.48,2.17,5,3.79,3.79,3.48,2.48,1.55,2,1.24,1.24,1,0.93,0.93,0.93,0.93,3,8,3,6,10,11,5,7,7,3,2,4,1.44,1.26,1,1.62,1.62,1.44,1.44,1.44,0.9,1,0.72,0.72,1,0.54,0.54,0.54,0.54,2,4,4,8,17,19,13,14,4,5,3,1,4.72,2.38,3,3.06,2.72,2.72,2.72,1.7,2,1.36,1.36,1,1.02,1.02,1.02,1.02)
   #xs=c(1495,1455.5,1543.505,1501.183,1476.08,1449.88,1428.444,1426.934,1378.31,1383.402,1549.9,1337,1406.379,1406.378,1405.432,1444.248,1439.024,1450.466,1406.379,1406.381,1406.38,1265.742,1265.742,1406.38,1406.376,1406.376,1406.376,1265.742,1514,1458.255,1468,1513.172,1518.005,1450.191,1459.008,1408.579,1441.14,1402.34,1440.01,1392.5,1406.382,1406.381,1265.74,1406.383,1406.383,1406.382,1406.382,1406.382,1406.378,1406.38,1265.736,1265.736,1406.38,1406.389,1406.389,1406.389,1265.741,1575,1477.5,1523.25,1509.38,1483,1452.89,1456.003,1450.36,1476.5,1441.6,1272.01,1450,1467.235,1406.378,1265.74,1406.379,1406.379,1406.379,1406.379,1406.382,1406.38,1265.743,1265.743,1406.38,1406.382,1406.382,1406.382,1265.735)
   #const.fit(ys,xs,k=4)
################################################################################
### Given a training set and leaf indicator
### Create a table with each row representing a leaf.
summarize.tree=function(x,y,leaves)
{
  unik.leaves=unique(leaves)
  nleaves=length(unik.leaves)
  x.name=colnames(x)
  p=ncol(x)
  #cat(p,'\n')
  y.est=c()
  x.est=c()
  count=c()
  y.std=c()
  res.name=c()
  for(column in 1:p)
  {
   if(is.factor(x[,column]))
      res.name=c(res.name,x.name[column])
   if(is.numeric(x[,column]))
      res.name=c(res.name,paste(x.name[column],c("min","max"),sep='.'))
  }
  for(leaf in 1:nleaves)
  {
    #cat(leaf,'\n')
    sub.set=(leaves==unik.leaves[leaf])
    count=c(count,sum(sub.set))
    y.est=c(y.est,mean(y[sub.set]))
    y.std=c(y.std,sd(y[sub.set]))
    x.est.cur=c()
    for(column in 1:p)
    {
     if(is.factor(x[sub.set,column]))
        x.est.cur=c(x.est.cur,paste(as.character(unique(x[sub.set,column])),collapse=','))
     if(is.numeric(x[sub.set,column]))
        x.est.cur=c(x.est.cur,min(x[sub.set,column]),max(x[sub.set,column]))
    }
   x.est=rbind(x.est,t(x.est.cur))
  }
  res=cbind(x.est,y.est,y.std,count)
  #cat(c(x.name,"Fitted"),'\n')
  #cat(length(colnames(res)),'\n')
  colnames(res)=c(res.name,"Fitted","St.Dev","Size")
  return(res)
}
################################################################################
### Conduct a single split, given the vector and x represents the variable
### maxsplit represents the maximum number of splits, and minimum size of each node.
### returns a matrix.
### Error message: error=0; error=1: single level; error=2, all splits result in small nodes;
### error=3: no split will give enough reduction in SSE.
### crit.parent: the criterion for the parent node.
### crit.imp: the minimum improvement on the criterion in order to partition.
### where xmat and covariate must be matrices.
### Number of categories can't be too large for setparts().
#library(partitions)
################################################################################
### Error message: error=0; error=1: single level; error=2, all splits result in small nodes;
### error=3: no split will give enough reduction in SSE.
### crit.parent: the criterion for the parent node.
### crit.imp: the minimum improvement on the criterion in order to partition.
### where xmat and covariate must be matrices.
### Number of categories can't be too large for setparts().
split.node=function(xmat,node.cur,covariate,y,crit.parent,var.names,crit.min=0,maxsplit=8,minsize=20,reg.method="linear",ncate.low=5,ncate.high=40,F_test=F,ncate.cont=10)
{
  covariate=as.matrix(covariate)
  maxsplit.cur=maxsplit
  minsize.cur=minsize
  split.rule=c()
  split.est=c()
  split.var=c()
  ## number of columns for xmat.
  p=ncol(xmat)
  nx=nrow(xmat)
  ## The index of elements in the node in the current data set
  index.cur=rep(node.cur,nx)
  crit.imp=0
  crit.leftnode=c()
  crit.rightnode=c()
  coef.leftnode=c()
  coef.rightnode=c()
  split.crit=c()
  imp.v=rep(0,p)
  split.opt=var.names[1]
  #cat(p,nx,split.opt,'\n')
  if(F_test & ncol(xmat)>1)
  {
   #cat(dim(xmat),'\n')
   #cat(is.numeric(xmat[,1]),is.numeric(xmat[,2]),is.numeric(xmat[,3]),'\n')
   split.opt=F.stat(y,xmat,covariate,Lcont=ncate.cont)$opt
  }
  #cat(split.opt,":",var.names,'\n')
  #cat(node.cur,"size",nx,'\n')
  for(j in 1:p)
  {
  if(F_test) {if(split.opt!=var.names[j]) next}
  x=xmat[,j]
  #cat(var.names[j],is.numeric(x),is.factor(x),is.ordered(x),'\n')
  ## decide the variable type; splitting on a numerical variable.
  ## return the split points.
  if(is.numeric(x))
  {
   x.cut=sort(unique(x))
   if(length(x.cut)==1) next;
   x.cut=(x.cut[-1]+x.cut[-length(x.cut)])/2
   ## the vector of cutoff values.
   x.cut.crude=x.cut
   if(length(x.cut)>(3*maxsplit))
    {
    tmp=round(seq(1,length(x.cut),length.out=maxsplit))
    x.cut.crude=x.cut[tmp]
    }
    for(i in 1:length(x.cut.crude))
    {
      split.cur=(x <= x.cut.crude[i])
      if(sum(split.cur)< minsize| sum(split.cur)>(nx-minsize)) next;
      y.left=y[split.cur]
      covariate.left=as.matrix(covariate[split.cur,])
      ## use regression SSE
      reg.left=regression(y.left,covariate.left,method=reg.method)
      crit.left=sum((reg.left$residuals)^2)
      coef.left=reg.left$coef
      y.right=y[!split.cur]
      covariate.right=as.matrix(covariate[!split.cur,])
      ## use regression SSE
      reg.right=regression(y.right,covariate.right,method=reg.method)
      crit.right=sum((reg.right$residuals)^2)
      coef.right=reg.right$coef
      crit.imp.cur=crit.parent-crit.left-crit.right
      if(crit.imp.cur>crit.imp & crit.imp.cur>crit.min)
      {
        crit.imp=crit.imp.cur
        split.crit=paste(node.cur,"|",var.names[j],"|<=|{",x.cut.crude[i],"}",sep='')
        split.var=var.names[j]
        imp.v=rep(0,p)
        imp.v[j]=crit.imp.cur
        index.cur[split.cur]=paste(unique(node.cur),"1",sep='')
        index.cur[!split.cur]=paste(unique(node.cur),"0",sep='')
        crit.leftnode=crit.left
        crit.rightnode=crit.right
        coef.leftnode=coef.left
        coef.rightnode=coef.right
      }
    }
  }
  ## decide on the variable type; splitting on a categorical variable: no more than 15 categories;
  ## return the association matrix of category values.
  if((!is.numeric(x)) & (!is.ordered(x)))
  {
   ## if the number of levels of x is 1, then next variable.
   x=as.character(x)
   tmp=matrix(c(T,F),ncol=1,byrow=T)
   nleve=lengthunique(x)
   if(nleve==1) next;
   #cat(nleve, ncate.high, ncate.low,'\n')
   if(nleve > ncate.high | nleve <= ncate.low)
   {
    #cat("Wrong category!","\n")
   if(nleve > 2 & nleve <=  ncate.low)
   {
    hoho=(1:(2^(nleve-1)-1))-1
    tmp=(rbind(1,digitsBase(hoho,base=2,nleve-1))==1)
   }
   ## If the number of levels is larger than 40, treat it as ordinal.
   if(nleve > ncate.high)
   {
    #cat(nleve, ncate.high, ncate.low,"second time",'\n')
    tmp=c()
    #if(max(table(x))<= minsize) next;
    x.unik=names(sort(table(x),decreasing=T))
    covariate.mean=as.vector(apply(covariate,2,mean))
    crit.est=rep(0,length(x.unik))
    #print(table(x))
    for(kk in 1:length(x.unik))
    {
    crit.cur=(-Inf)
    if(sum(x==x.unik[kk])>= minsize)
     {
     crit.cur=sum(regression(y[x==x.unik[kk]],covariate[x==x.unik[kk],],method=reg.method)$coef*covariate.mean)
     }
    #if(sum(x==x.unik[kk])< minsize) crit.cur=0#sample(crit.est[crit.est!=0],1)
    crit.est[unique(x)==x.unik[kk]]=crit.cur
   }
   ## Take a random ordering if all the crit.est are -Inf.
   if(sum(crit.est != (-Inf)) <1) crit.est=rnorm(length(crit.est))
   crit.cut=(sort(crit.est)[-1]+sort(crit.est)[-length(crit.est)])/2
   #cat(minsize,":",crit.est,'\n')
   #cat(crit.cut,'\n')
   for(kk in 1:(length(x.unik)-1))
   {
     tmp=cbind(tmp,crit.est<crit.cut[kk])
   }
     #cat("Factor, ordered",'\n')
   }
   for(i in 1:ncol(tmp))
   {
      split.cur=(x %in% unique(x)[tmp[,i]])
      #cat(unique(x)[tmp[,i]],sum(split.cur),'\n')
      if(sum(split.cur)< minsize| sum(split.cur)>(nx-minsize)) next;
      y.left=y[split.cur]
      covariate.left=as.matrix(covariate[split.cur,])
      ## use regression SSE
      reg.left=regression(y.left,covariate.left,method=reg.method)
      crit.left=sum((reg.left$residuals)^2)
      coef.left=reg.left$coef
      y.right=y[!split.cur]
      covariate.right=as.matrix(covariate[!split.cur,])
      ## use regression SSE
      reg.right=regression(y.right,covariate.right,method=reg.method)
      crit.right=sum((reg.right$residuals)^2)
      coef.right=reg.right$coef
      crit.imp.cur=crit.parent-crit.left-crit.right
      #cat("xixi:",i,crit.imp.cur,'\n')
      if(crit.imp.cur>crit.imp & crit.imp.cur>crit.min)
      {
        crit.imp=crit.imp.cur
        split.crit=paste(node.cur,"|",var.names[j],"| in |{",paste(unique(x)[tmp[,i]],collapse=","),"}",sep='')
        split.var=var.names[j]
        imp.v=rep(0,p)
        imp.v[j]=crit.imp.cur
        index.cur[split.cur]=paste(unique(node.cur),"1",sep='')
        index.cur[!split.cur]=paste(unique(node.cur),"0",sep='')
        crit.leftnode=crit.left
        crit.rightnode=crit.right
        coef.leftnode=coef.left
        coef.rightnode=coef.right
      }
    }
   }
   if(nleve> ncate.low & nleve<= ncate.high)
   {
    #hoho=unique(floor(runif(15,0,1)*(2^(nleve-1)-1)))[1:5]
    tmp=(rbind(1,matrix(rbinom((nleve-1)*5,1,1/2),nrow=(nleve-1),byrow=T))==1)
    #cat(tmp,'\n')
    for(i in 1:ncol(tmp))
    {
      #cat(i,"th initial point","\n")
      split.cur=(x %in% unique(x)[tmp[,i]])
      if(sum(split.cur)< minsize| sum(split.cur)>(nx-minsize)) next;
      ## use regression SSE
      #cat("steep.descent",'\n')
      #cat(unique(x)[tmp[,i]],'\n')
      #cat(sum(split.cur),'\n')
      crit.left=sum((regression(y[split.cur],covariate[split.cur,],method=reg.method)$residuals)^2)
      y.right=y[!split.cur]
      covariate.right=as.matrix(covariate[!split.cur,])
      ## use regression SSE
      crit.right=sum((regression(y.right,covariate.right,method=reg.method)$residuals)^2)
      logic.v=steep.descent(tmp[,i],x,y,covariate,crit.parent-crit.left-crit.right,crit.parent,reg.method,minsize,crit.min)
      #cat("steep.descent:",logic.v,'\n')
      rm(split.cur)
      split.cur=(x %in% unique(x)[logic.v])
      ## use regression SSE
      reg.left=regression(y[split.cur],covariate[split.cur,],method=reg.method)
      crit.left=sum((reg.left$residuals)^2)
      coef.left=reg.left$coef
      ## use regression SSE
      reg.right=regression(y[!split.cur],covariate[!split.cur,],method=reg.method)
      crit.right=sum((reg.right$residuals)^2)
      coef.right=reg.right$coef
      crit.imp.cur=crit.parent-crit.left-crit.right
      #cat(xixi,'\n')
      if(crit.imp.cur>crit.imp & crit.imp.cur>crit.min)
      {
        crit.imp=crit.imp.cur
        split.crit=paste(node.cur,"|",var.names[j],"| in |{",paste(unique(x)[logic.v],collapse=","),"}",sep='')
        split.var=var.names[j]
        imp.v=rep(0,p)
        imp.v[j]=crit.imp.cur
        index.cur[split.cur]=paste(unique(node.cur),"1",sep='')
        index.cur[!split.cur]=paste(unique(node.cur),"0",sep='')
        crit.leftnode=crit.left
        crit.rightnode=crit.right
        coef.leftnode=coef.left
        coef.rightnode=coef.right
      }
    }
   }
  }
  if(is.ordered(x))
  {
   #cat("Ordinal",'\n')
   x.unik=(levels(x)[levels(x) %in% x])
   nleve=length(x.unik)
   if(nleve==1) next;
   tmp=matrix(rep(x,nleve-1),nrow=(nleve-1),byrow=T) < (x.unik[-1])
   for(i in 1:nrow(tmp))
    {
      #cat("Ordinal row:",i,'\n')
      split.cur=(tmp[i,])
      if(sum(split.cur)< minsize| sum(split.cur)>(nx-minsize)) next;
      y.left=y[split.cur]
      covariate.left=as.matrix(covariate[split.cur,])
      ## use regression SSE
      reg.left=regression(y.left,covariate.left,method=reg.method)
      crit.left=sum((reg.left$residuals)^2)
      coef.left=reg.left$coef
      y.right=y[!split.cur]
      covariate.right=as.matrix(covariate[!split.cur,])
      ## use regression SSE
      reg.right=regression(y.right,covariate.right,method=reg.method)
      crit.right=sum((reg.right$residuals)^2)
      coef.right=reg.right$coef
      crit.imp.cur=crit.parent-crit.left-crit.right
      #cat(crit.imp.cur,'\n')
      if(crit.imp.cur>crit.imp & crit.imp.cur>crit.min)
      {
        crit.imp=crit.imp.cur
        split.crit=paste(node.cur,"|",var.names[j],"| in |{",paste(x.unik[1:i],collapse=","),"}",sep='')
        split.var=var.names[j]
        imp.v=rep(0,p)
        imp.v[j]=crit.imp.cur
        index.cur[split.cur]=paste(unique(node.cur),"1",sep='')
        index.cur[!split.cur]=paste(unique(node.cur),"0",sep='')
        crit.leftnode=crit.left
        crit.rightnode=crit.right
        coef.leftnode=coef.left
        coef.rightnode=coef.right
      }
    }
  }
 }
 left.node=paste(unique(node.cur),"1",sep='')
 right.node=paste(unique(node.cur),"0",sep='')
 #cat(imp.v,'\n')
 modelest=rbind(t(c(as.numeric(left.node),sum(index.cur==left.node),coef.leftnode,crit.leftnode)),t(c(as.numeric(right.node),sum(index.cur==right.node),coef.rightnode,crit.rightnode)))
 return(list(imp.vec=imp.v,leaf.assoc=index.cur,split.krit=split.crit,Krit.imp=crit.imp,model.est=modelest,var.split=split.var))
}
################################################################################
## Fitting a simple tree with M nodes.
Simple.tree=function(x.mat,nodecur="1",covariate,y,crit.root,var.names,M=4,crit.min=0,maxsplit=8,minisize=20,regmethod="linear",ncate_low=5,ncate_high=40,F.test=F,L.cont=10)
{
 leaves.cur=nodecur
 split.rule=c()
 Leaf=rep(nodecur,nrow(x.mat))
 #cat(length(y),dim(x.mat),'\n')
 #cat(Leaf,'\n')
 node.summary=c()
 split.rule=c()
 spt.var=c()
 spt.imp=c()
 spt.impvec=rep(0,length(var.names))
 #cat(y,'\n')
 for(m in 1:(M-1))
 {
   crit.opt=0
   res.opt=c()
   spt.leaf=c()
   crit.opt=0
   for(ileaf in leaves.cur)
   {
    #cat(ileaf,'\n')
    crit.parent=node.summary[node.summary[,1]==as.numeric(ileaf),ncol(node.summary)]
    if(ileaf=="1") crit.parent=crit.root
    #cat("hehe",ileaf,'\n')
    #cat(ileaf,sum(Leaf==ileaf),'\n')
    #cat(ileaf,"before split",'\n')
    heihei=x.mat[Leaf==ileaf,]
    if(ncol(x.mat)<2)
     {
      heihei=as.matrix(heihei)
      #cat("Simple.tree",'\n')
      }
    #cat(ileaf,'\n')
    #cat(Leaf,'\n')
    #print(cbind(y,Leaf))
    res.tmp=split.node(heihei,ileaf,covariate[Leaf==ileaf,],y[Leaf==ileaf],crit.parent,var.names,minsize=minisize,reg.method=regmethod,ncate.low=ncate_low,ncate.high=ncate_high,F_test=F.test,ncate.cont=L.cont)
    #Krit.imp
    #cat(ileaf,"after split",'\n')
    #print(res.tmp$model.est)
    #cat(ileaf,'\n')
    #print(res.tmp$imp.vec)
    #cat(ileaf,",",res.tmp$Krit.imp,'\n')
    if(res.tmp$Krit.imp> crit.opt)
      {
       spt.leaf=ileaf
       crit.opt=res.tmp$Krit.imp
       res.opt=res.tmp
      }
   }
   #cat(m,is.null(res.opt$split.krit),res.opt$split.krit,'\n')
   if(is.null(res.opt$split.krit)) break;
   Leaf[Leaf==spt.leaf]=res.opt$leaf.assoc
   if(m>1) split.rule=paste(split.rule,res.opt$split.krit,sep=";")
   if(m==1) split.rule=res.opt$split.krit
   spt.var=c(spt.var,res.opt$var.split)
   spt.imp=c(spt.imp,res.opt$Krit.imp)
   if(is.null(res.opt$imp.vec)) res.opt$imp.vec=0
   #cat(m,"nodes,",,res.opt$imp.vec,'\n')
   spt.impvec=spt.impvec+res.opt$imp.vec
   #cat("here:")
   #print(res.opt$model.est)
   node.summary=rbind(node.summary,res.opt$model.est)
   leaves.cur=unique(Leaf)
 }
 mode.est=node.summary[node.summary[,1] %in% as.numeric(leaves.cur),]
 if(sum(node.summary[,1] %in% as.numeric(leaves.cur))==0)
 {
  #cat("no splits,",'\n')
  split.rule=NA
  mode.est=matrix(c(1,length(y),0,0),nrow=1,byrow=T)
 }
 #cat("Finished M nodes",'\n')
 return(list(split.impvec=spt.impvec,split.imp=spt.imp,split.var=spt.var,leaf.assoc=Leaf,split.krit=split.rule,model.est=mode.est))
}
################################################################################
## Fitting a simple tree with M nodes.
 #res=Simple.tree(xixi[,1:3],"1",as.matrix(cbind(1,xixi$z)),xixi$y,crit.root,name.split,M,crit.min=0,minsize=20,reg.method="linear")
 #res.pred=split.predict(xixi[,1:3],as.matrix(cbind(1,xixi$z)),res$split.krit,res$model.est,3:4)
Boost.Tree=function(varPart,varReg,varY,name.split,M=4,nBoost=100,nu=0.1,fit.method="linear",mini.size=20,ncatelow=5,ncatehigh=20)
{
 resid.cur=varY-mean(varY)
 coef.Boost=varReg-varReg
 coef.Boost[,1]=mean(varY)
 #cat(mean(varY))
 #print(coef.Boost)
 var.imp=rep(0,length(name.split))
 splt.crit=c()
 splt.modelest=c()
 R2.in=c()
 error.in=c()
 for(iBoost in 1:nBoost)
 {
  if(iBoost %% 9==1) {par(mfrow=c(3,3))}
  crit.root=sum((regression(resid.cur,varReg,method=fit.method)$residuals)^2)
  res=Simple.tree(varPart,"1",varReg,resid.cur,crit.root,name.split,M,crit.min=0,minisize=mini.size,regmethod=fit.method,ncate_low=ncatelow,ncate_high=ncatehigh)
  res.pred=split.estimate(res$leaf.assoc,varReg,res$model.est,3:(ncol(varReg)+2))
  coef.increment=nu*res.pred$coefs
  var.imp=var.imp+(nu^2)*res$split.impvec
  coef.Boost=coef.Boost+coef.increment
  resid.cur=varY-apply(varReg*coef.Boost,1,sum)
  splt.modelest[[iBoost]]=res$model.est
  splt.crit=c(splt.crit,res$split.krit)
  plot(varY,apply(varReg*coef.Boost,1,sum),type='p',pch='.',cex=2,main=iBoost,ylab='fits')
  abline(0,1,col='red')
  R2.in=rbind(R2.in,t(unlist(calc.R2(apply(varReg*coef.Boost,1,sum),varY))))
  error.in=rbind(error.in,t(c(mean((apply(varReg*coef.Boost,1,sum)-varY)^2),mean(abs(apply(varReg*coef.Boost,1,sum)-varY)))))
  cat(iBoost,":",sum(resid.cur^2),'\n')
 }
 #R2.in=unlist(R2.in)
 colnames(error.in)=c("L2","L1")
 return(list(variable.imp=var.imp,model.est=splt.modelest,split.krit=splt.crit,coefs=coef.Boost,R2in=R2.in,f0.hat=mean(varY),error.train=error.in))
}
################################################################################
## Summary of the boosting results.
summary.Boost=function(boost,var.names,len=12,nam="Variable Importance")
{
 var.imp=matrix(boost$variable.imp/max(boost$variable.imp),nrow=1,byrow=T)
 colnames(var.imp)=substr(var.names,1,len)[order(var.imp,decreasing=T)]
 var.imp[1,]=sort(var.imp,decreasing=T)
 par(mfrow=c(1,2))
 hist(boost$coefs[,1],main='Intercept',freq=F,prob=T,xlab="Intercept")
 hist(boost$coefs[,2],main='Slope',freq=F,prob=T,xlab="Slope")
 #x11()
 pdf("VarImp.pdf")
 barplot(var.imp[,var.imp>0.1]*100,horiz =F,main=nam,col="red",las=2,cex.axis=0.8)
 #barplot(var.freq[,var.imp>0.1]*100,horiz =F,main=nam,col="red",las=2,cex.axis=0.8)
 dev.off()
 return(list(imp.vec=var.imp))
}
################################################################################
## Out-sample prediction based on boosting.
pred.Boost=function(varPart,varReg,varY,obj.Boost,cols.coef,nBoost,nu)
{
  pred.mat=rep(obj.Boost$f0.hat,length(varY))
  pred=obj.Boost$f0.hat
  R2.out=c()
  errorL2.test=c()
  errorL1.test=c()
  coef.mat=varReg-varReg
  coef.mat[,1]=obj.Boost$f0.hat
  for(iBoost in 1:nBoost)
  {
   spt.pred=split.predict(varPart,varReg,obj.Boost$split.krit[iBoost],obj.Boost$model.est[[iBoost]],cols=cols.coef)
   pred=pred+spt.pred$pred*nu
   #pred.mat=cbind(pred.mat,pred)
   R2.out=c(R2.out,calc.R2(pred,varY)$R2)
   coef.mat=coef.mat+spt.pred$coefs*nu
   errorL2.test=c(errorL2.test,mean((pred-varY)^2))
   errorL1.test=c(errorL1.test,mean(abs(pred-varY)))
  }
  return(list(predictions=pred.mat,R2=R2.out,coefs=coef.mat,L2.error=errorL2.test,L1.error=errorL1.test))
  #return(list(R2=R2.out,coefs=coef.mat,L2.error=errorL2.test,L1.error=errorL1.test))
}
################################################################################
## Partial dependence plot.
 Partial.Boost=function(varPart,varReg,varY,name.var)
  {
   var.levels=sort(unique(varPart[,name.var]))
   #var.levels=sample(var.levels,min(length(var.levels),100),replace=F)
   #if(is.numeric(var.levels))
   #var.levels=as.factor(sort(var.levels))
   fits=varY-varY
   varPart.cur=varPart
   for(ilevels.no in 1:length(var.levels))
   {
    ilevels=var.levels[ilevels.no]
    varPart.cur[,name.var]=ilevels
    cat(ilevels,'\n')
    pred.cur=pred.Boost(varPart.cur,cbind(1,varReg),varY,res,3:4,nBoost.param,nu=nu.param)
    cat(apply(pred.cur$coefs,2,mean),'\n')
    fits[varPart[,name.var]==ilevels]=cbind(1,varReg)[varPart[,name.var]==ilevels,]%*%apply(pred.cur$coefs,2,mean)
  }
  return(fits)
  }
##############################################################################
## Given the vector of leaf associations, and model estimates, provide coefficients,
## as well as the predicted values.
split.estimate=function(leaves,covar,modelest,cols)
{
  #leaves=res$leaf.assoc
  #covar=as.matrix(cbind(1,xixi$z))
  #modelest=res$model.est
  #cols=3:4
  coef.mat=matrix(0,nrow=nrow(covar),ncol=length(cols),byrow=T)
  for(imod in 1:nrow(modelest))
  {
    imod.tmp=(as.numeric(leaves)==modelest[imod,1])
    coef.mat[imod.tmp,]=matrix(rep(modelest[imod,cols],sum(imod.tmp)),ncol=2,byrow=T)
  }
  return(list(pred=apply(coef.mat*covar,1,sum),coefs=coef.mat))
}
####################################################################################
## Compute R2 for given predicted and actual.
calc.R2=function(pred, act)
{
  return(list(R2=1-sum((act-pred)^2)/sum((act-mean(act))^2),R1=1-sum(abs(act-pred))/sum(abs(act-mean(act)))))
}
################################################################################
## Find the predicted responses given feature matrix xmat, predictor matrix covariate,
## the split rule split.rule, and summary of leaves (leaf ID and coefficients).
split.predict=function(xmat,covariate,split.rule,leaf.summary,cols=2:ncol(leaf.summary))
{
  #cat(split.rule,'\n')
  if(is.na(split.rule)) return(list(assoc=rep("1",nrow(covariate)),pred=0,coefs=0,spt.rules=NA))
  #leaf.summary=res$model.est[,c(1,3,4)]
  #xmat=xixi[,1:3]
  #covariate=cbind(1,xixi$z)
  #split.rule=res$split.krit
  ##############################################################################
  #split.rules=matrix(NA,ncol=4,nrow=length(split.rule),byrow=T)
  split.rules=matrix(unlist(lapply(unlist(strsplit(split.rule[!is.na(split.rule)],";")),strsplit,"\\|")),ncol=4,byrow=T)
  leaf.assoc=rep("1",nrow(xmat))
  leaf.pred=rep(0,nrow(xmat))
  coef.mat=matrix(0,nrow=nrow(xmat),ncol=ncol(covariate),byrow=T)
  for(i in 1:nrow(split.rules))
  {
    #cat(i,split.rules[i,2],'\n')
    logic.v=rep(F,nrow(xmat))
    if(is.numeric(xmat[,split.rules[i,2]])) logic.v[xmat[,split.rules[i,2]]<= as.numeric(substr(split.rules[i,4],2,nchar(split.rules[i,4])-1))]=T
    if(is.factor(xmat[,split.rules[i,2]]))
    {
     #cat(unlist(strsplit(substr(split.rules[i,4],2,nchar(split.rules[i,4])-1),'\\,')),'\n')
     logic.v[xmat[,split.rules[i,2]] %in% unlist(strsplit(substr(split.rules[i,4],2,nchar(split.rules[i,4])-1),'\\,'))]=T
    }
    leaf.assoc[substr(leaf.assoc,1,nchar(split.rules[i,1]))==split.rules[i,1] & logic.v]=paste(split.rules[i,1],"1",sep="")
    leaf.assoc[substr(leaf.assoc,1,nchar(split.rules[i,1]))==split.rules[i,1] & (!logic.v)]=paste(split.rules[i,1],"0",sep="")
   }
  #print(split.rules)
  #print(leaf.assoc[1:68])
  ## handle exceptions, assign to the leading category.
  leaf.assoc[(leaf.assoc %in% split.rules[,1])]=names(table(leaf.assoc)[table(leaf.assoc)==max(table(leaf.assoc))])
  for(i in 1:nrow(leaf.summary))
  {
   #cat(leaf.summary[i,1],'\n')
   logic.v=rep(F,nrow(xmat))
   logic.v[as.numeric(leaf.assoc)==leaf.summary[i,1]]=T
   if(sum(logic.v)==0) next;
   hehe=covariate[logic.v,]
   if(sum(logic.v)==1) hehe=as.matrix(t(hehe))
   if(length(leaf.pred[logic.v]) != nrow(as.matrix(hehe))) cat(length(leaf.pred[logic.v]),nrow(as.matrix(hehe)) ,"Wierd!",'\n')
   #cat(,dim(as.matrix(covariate[logic.v,])),'\n')
   leaf.pred[logic.v]=as.matrix(hehe) %*% as.vector(leaf.summary[i,cols])
   #cat("heihei",sum(logic.v),length(logic.v),'\n')
   #print(coef.mat[logic.v,])
   #cat("heihei",'\n')
   #print(matrix(rep(leaf.summary[i,cols],sum(logic.v)),nrow=sum(logic.v),byrow=T))
   coef.mat[logic.v,]=matrix(rep(leaf.summary[i,cols],sum(logic.v)),nrow=sum(logic.v),byrow=T)
  }
  return(list(assoc=leaf.assoc,pred=leaf.pred,coefs=coef.mat,spt.rules=split.rules))
}

