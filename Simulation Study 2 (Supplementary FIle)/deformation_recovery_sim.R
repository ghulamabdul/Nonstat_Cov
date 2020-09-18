###################################################
##### Simulation study : Deformation recovery #####
###################################################


#################################################
############ Loading required packages ##########
#################################################

library(fdasrvf)
library(geoR)
library(fields)
library(plot3D)
library(rgl)
library(scatterplot3d)
library(scoringRules)
library(mvtnorm)


##################################################
########## Setting working directory #############
##################################################
setwd("/Users/qadirga/Documents/Project 1 (Deformation)/Technometrics submission/Revision/Revised Final codes/Simulation 2 (supplementary)")


sim.index=7# original 7
  
##################################################
################ Defining grid ###################
##################################################
  
  N=70 #### Setting grid resolution ######
  x<-seq(0,2,length=N)
  y<-seq(0,2,length=N)
  
  coord.x<-rep(x,times=N)
  coord.y<-rep(y,each=N)
  
  grid<-expand.grid(x,y)
  colnames(grid)=c("x","y")
  par(mfrow=c(1,1),pty="s")
  plot(grid)
  
  ###################################################
  ####### Defining stationary Matérn function #######
  ###################################################
  
  my.matern<-function(h,a,sigma,nu,variog=FALSE)
  {
    h[h==0]<-1e-10
    num1<-(sigma^2)*(2^(1-nu))/gamma(nu)
    num2<-(h*a)^nu
    num3<-besselK(x=(h*a), nu=nu)
    temp1<-num1*num2*num3
    temp2<-sigma^2-temp1
    if(variog==FALSE)
      return(temp1)
    else{
      return(temp2)
    }
  }
  
  ##################################################
  ##### Seting deformation function    #############
  ##################################################
  #### Composition of radial basis functions #######
  ############## Elementary radial basis function ##############
  el.rbf<-function(s,a,b,c)
  {
    h<-s-c
    norm.h<-apply(h, 1, function(x) sqrt(sum(x^2)))
    norm.h<-1+(b*exp(-a*(norm.h^2)))
    tmpp<-cbind(norm.h,norm.h)
    return(c+h*norm.h)
  }
  
  
  myel1.1<-el.rbf(s=grid,a=2,b=1.0,c=cbind(rep(1.7,times=length(grid$x)),rep(0.5,times=length(grid$x))))
  plot(myel1.1)
  
  myel1.2<-el.rbf(s=myel1.1,a=2,b=1.0,c=cbind(rep(1.7,times=length(grid$x)),rep(1.5,times=length(grid$x))))
  plot(myel1.2)
  
  myel2.1<-el.rbf(s=myel1.2,a=2,b=-0.5,c=cbind(rep(0.5,times=length(grid$x)),rep(0.5,times=length(grid$x))))
  plot(myel2.1)
  
  myel2.2<-el.rbf(s=myel2.1,a=2,b=-0.5,c=cbind(rep(0.5,times=length(grid$x)),rep(1.5,times=length(grid$x))))
  plot(myel2.2)
  
  ####### Simulating nonstationary random field ######
  
  nonst.dist<-rdist(myel2.2)
  
  
  ########## Matern parameters #########
  a.nst<-1/(0.4^2)
  sigma.nst<-sqrt(5)
  nu.nst<-0.8
  my_Cns<-my.matern(h=nonst.dist,a=a.nst,sigma = sigma.nst,nu=nu.nst)
  ##################################################
  #### Assigning subregions to the locations #######
  ##################################################
  
  grid$subr<-ifelse(grid$x<1.2,"R1","R2")
  
  
  ######### Matrix reordering as per ``grid" locations #############
  R1locs<-grid[grid$subr=="R1",]
  R2locs<-grid[grid$subr=="R2",]
  
  Re_ord_grid<-rbind(R1locs,R2locs)
  ordering.index<-order(Re_ord_grid$y,Re_ord_grid$x)
  
  #### Setting seed ######
  set.seed(sim.index)
  
  ##### simulating nonstationary Gaussian random field  #######
  error<-rmvnorm(n=1,sigma = my_Cns,method = "chol")
  #########################################################################
  ################### Plotting Simulated Random Field #####################
  #########################################################################
  
  quilt.plot(grid[-3],error,nx=N,ny=N)
  
  
  ############# Drawing partitioning line ###########
  abline(v=1.2)
  
  
  
  ##################################################################
  ######## Splitting into training and testing datasets ############
  ##################################################################
  
  ### Choosing a random sample of data for estimating the deformation ###
  
  rmiss<-sample(1:(N^2),N^2-1200)
  
  
  
  
  ############################################################################
  ############### Creating dataframe for the complete simulated data #########
  ############################################################################
  
  
  
  data<-data.frame(x=grid$x,y=grid$y,z=c(error),subregion=grid$subr)
  data$split<-"Training"
  data$split[rmiss]<-"Test"
  
  
  ###########################################################################################################
  ############ Extracting data for subregion 1 and subregion and plotting them respectively #################
  ###########################################################################################################
  
  
  reg1<-data[data$subregion=="R1",]
  reg2<-data[data$subregion=="R2",]
  
  
  ##### extracting region wise training points ######
  
  reg1.tr<-reg1[reg1$split=="Training",]
  reg2.tr<-reg2[reg2$split=="Training",]
  
  ####################################################################
  ################ computing empirical variograms ####################
  ####################################################################
  
  vr1<-variog(coords = cbind(reg1.tr$x,reg1.tr$y),data=reg1.tr$z,option = "bin",uvec=30,max.dist=(sqrt(8))/2)
  plot(vr1,main="Empirical Variogram for Region 1")
  
  vr2<-variog(coords = cbind(reg2.tr$x,reg2.tr$y),data=reg2.tr$z,uvec=30,max.dist = (sqrt(8))/2,option = "bin")
  plot(vr2,main="Empirical Variogram for Region 2")
  
  ################################################################################
  ######## Fitting WLS Matérn variogram to be used as an MLE initial values ######
  ################################################################################
  
  #############################
  ###### For subregion 1 ######
  #############################
  
  
  myvario<-vr1
  
  my.var.loss<-function(p,vardist=myvario$u,gamma_hat=myvario$v,nbins=myvario$n)
  {
    a<-p[1]
    nu<-p[2]
    sigma<-p[3]
    if(sum(p<=0)!=0)
    {
      return(Inf)
    }
    else
    {
      gamma_theta<-(my.matern(h=0,a=a,nu=nu,sigma = sigma)-my.matern(h=vardist,a=a,nu=nu,sigma = sigma))
      nk<-nbins
      temp<-nk*(((gamma_hat-gamma_theta)/gamma_theta)^2)
      return(sum(temp))
    }
  }
  
  temp.r1<-optim(my.var.loss,par = c(runif(n=1,min=0.1,max=30),0.5,sd(reg1.tr$z)))
  temp.r1
  
  ####################################################################################
  ########## Now using MLE to estimate stationary Matérn covariance parameters #######
  ####################################################################################
  
  
  mle.comp.all<-function(locs,z,p)
  {
    #distmat<-rdist(locs)
    a<-p[1]
    nu<-p[2]
    sigma<-p[3]
    distmat<-rdist(locs)
    if(sum(c(p<=0))!=0)
    {
      return(list(mlv=Inf))
    }
    else
    {
      
      
      C<-my.matern(h=distmat,a=a,nu=nu,sigma = sigma)
      nlogl<--dmvnorm(x=z,mean=(rep(0,times=length(z))),sigma = C,log = T)
      
      return(list(mlv=nlogl))
    }
  }
  
  mle.comp_only_mlv.r1<-function(par)
  {
    return(mle.comp.all(locs=cbind(reg1.tr$x,reg1.tr$y),z=reg1.tr$z,p=par)$mlv
    )
  }
  optim_marg_comp.loglik.r1 <- function(par){
    optim(par=par,
          fn = mle.comp_only_mlv.r1,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=1000))
  }
  
  mle.comp_only_mlv.r2<-function(par)
  {
    return(mle.comp.all(locs=cbind(reg2.tr$x,reg2.tr$y),z=reg2.tr$z,p=par)$mlv
    )
  }
  optim_marg_comp.loglik.r2 <- function(par){
    optim(par=par,
          fn = mle.comp_only_mlv.r2,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=1000))
  }
  
  
  mle.r1<-optim_marg_comp.loglik.r1(par = temp.r1$par)
  mle.r1
  
  myvario<-vr2
  temp.r2<-optim(my.var.loss,par = c(runif(n=1,min=0.1,max=30),0.5,sd(reg2.tr$z)))
  temp.r2
  
  mle.r2<-optim_marg_comp.loglik.r2(par = temp.r2$par)
  mle.r2
  
  
  #############################################
  ######### Aligning Variograms ###############
  #############################################
  
  
  #### variable time represents the distances at which the regional variograms are evaluated for registration ###
  time<-seq(0,sqrt(8),length=1000)
  
  my_data<-data.frame(f1=my.matern(h=time,a=mle.r1$par[1],nu=mle.r1$par[2],sigma = mle.r1$par[3],variog = T),
                      f2=my.matern(h=time,a=mle.r2$par[1],nu=mle.r2$par[2],sigma = mle.r2$par[3],variog = T))
  my_out<-time_warping(f=as.matrix(my_data),time = time)
  
  ############################################
  ####### Plotting alignment results #########
  ############################################
  
  par(mfrow=c(3,2))
  plot(my_out)
  
  
  ##########################################################
  ####  Now we Inverting the distance warping function #####
  ##########################################################
  
  inv1<-invertGamma(my_out$gam[,1])
  inv2<-invertGamma(my_out$gam[,2])
  
  ###################################################
  #### Now we convert back to the original scale ####
  ###################################################
  
  scale_inv1<-max(time)*inv1
  scale_inv2<-max(time)*inv2
  
  
  ##############################################################################
  #### Plotting two fitted variograms and local distance warping functions #####
  ##############################################################################
  
  #par(mfrow=c(1,3))
  #plot(vr1$u,vr1$v,pch=19,col="black",ylim = c(0,max(vr1$v,vr2$v)),
  #    xlim = c(0,max(vr1$u,vr2$u)),
  #   main = "Empirical and Fitted variograms",xlab="||h||",ylab = bquote(gamma(h)))
  #points(vr2$u,vr2$v,pch=19,col="red")
  #lines(time,my.matern(h=time,a=mle.r1$par[1],nu=mle.r1$par[2],sigma = mle.r1$par[3],variog = T),
  #     col="black",lty=1)
  #lines(time,my.matern(h=time,a=mle.r2$par[1],nu=mle.r2$par[2],sigma = mle.r2$par[3],variog = T),
  #     col="red",lty=2)
  #legend("bottomright",lty=c(4,5),col=c("black","red"),c("Region1","Region2"))
  
  #plot(time,my_out$fn[,1],type="l",ylim=c(0,1.4),main="Aligned Variograms",ylab=expression(paste(gamma,"(h)")),xlab="h")
  #lines(time,my_out$fn[,2],type="l",col="red",lty=2)
  #legend("bottomright",lty=c(1,2),col=c("black","red"),c("Region1","Region2"))
  
  #plot(time,scale_inv1,type="l",main="Local Distance Warping Functions",xlab="h",ylab=expression(paste(phi['i'],"(h)")))
  #lines(time,scale_inv2,lty=2,col="red")
  #legend("bottomright",lty=c(1,2),col=c("black","red"),c("Region1","Region2"))
  
  
  #########################################################
  ######### Now we compute the warped distances ###########
  #########################################################
  
  
  ##############################################################################
  ########## b1 is the distance matrix for locations in subregion 1 ############
  ##############################################################################
  
  b1<-rdist(cbind(reg1$x,reg1$y))
  #max(b1)
  #####################################################################################################
  ############ Creating local distance warping function for subregion 1 with kernel smoothing #########
  ######### Gaussian kernel with bandwidth = 0.01 #####################################################
  #####################################################################################################
  
  dwarp1<-function(h)
  {
    
    h2<-matrix(NA,nrow=nrow(h),ncol=ncol(h))
    checkf<-data.frame(rows=c(row(h)), columns=c(col(h)),
                       values=c(h))
    checkf<-checkf[order(checkf$values),]
    h2[cbind(checkf$rows,checkf$columns)]<-ksmooth(x=time,y=scale_inv1,kernel = "normal",bandwidth = 0.01,x.points = checkf$values)$y
    return(h2)
    
  }
  
  
  ####################################################################################################
  ############# R11 is the warped distance matrix for subregion 1 ####################################
  ####################################################################################################
  #plot(time,scale_inv1)
  
  R11<-dwarp1(h=b1)
  diag(R11)<-0
  
  
  rm(b1) ####### Removing b1 matrix #######
  
  
  ########## b2 is the distance matrix for locations in subregion 2 ############
  
  b2<-rdist(cbind(reg2$x,reg2$y))
  
  #####################################################################################################
  ############ Creating local distance warping function for subregion 2 with kernel smoothing #########
  ######### Gaussian kernel with bandwidth = 0.01 #####################################################
  #####################################################################################################
  
  
  dwarp2<-function(h)
  {
    h2<-matrix(NA,nrow=nrow(h),ncol=ncol(h))
    checkf<-data.frame(rows=c(row(h)), columns=c(col(h)),
                       values=c(h))
    checkf<-checkf[order(checkf$values),]
    h2[cbind(checkf$rows,checkf$columns)]<-ksmooth(x=time,y=scale_inv2,kernel = "normal",bandwidth = 0.01,x.points = checkf$values)$y
    return(h2)
    
  }
  
  R22<-dwarp2(h=b2)
  diag(R22)<-0
  rm(b2)
  
  ####################################################################################################
  ############# R22 is the warped distance matrix for subregion 2 ####################################
  ####################################################################################################
  
  #####################################################################################################################################
  ########## Now computing the matrix R12, which is the warped distances for the two points in region 1 and region 2 respectively #####
  #####################################################################################################################################
  midpoint.f<-function(i,j,xmid=1.2)
  {
    x1<-reg1[i,]$x
    y1<-reg1[i,]$y
    x2<-reg2[j,]$x
    y2<-reg2[j,]$y
    ymid<-y1+((y2-y1)/(x2-x1))*(xmid-x1)
    return(ymid)
  }
  
  y.inter<-outer(1:length(reg1$x),1:length(reg2$x),midpoint.f)
  
  
  
  
  dwarp12<-function(region1=reg1,region2=reg2,y.intersect=y.inter)
  {
    lr1<-length(region1$x)
    lr2<-length(region2$x)
    h1<-t(apply(matrix(1:lr1,ncol=1),1,function(i) rdist(cbind(region1$x[i],region1$y[i]),cbind(1.2,y.intersect[i,]))))
    h2<-apply(matrix(1:lr2,ncol=1),1,function(i) rdist(cbind(1.2,y.intersect[,i]),cbind(region2$x[i],region2$y[i])))
    h12<-h1+h2
    p1<-dwarp1(h=h12)
    p2<-dwarp2(h=h12)
    R12<-(h1/(h12))*p1+(h2/(h12))*p2
    return(R12)
  }
  
  
  R12<-dwarp12()
  
  ################################################################################################
  ############# Constructing the global distance matrix wdist (Delta in the manuscript) ##########
  ################################################################################################
  
  row1<-cbind(R11,R12)
  row2<-cbind(t(R12),R22)
  wdist<-rbind(row1,row2)
  #rm(R11,R22,R12,row1,row2)
  sum(is.na(wdist))
  ################################################################################################################################
  ########### setting the diagonals of the distance matrix to be zero, as due to smoothing they might not exactly be zero ########
  ################################################################################################################################
  
  diag(wdist)<-0
  
  ############################################################################
  ####Constructing Deformed space coordinates in 30 dimensions ###############
  ############################################################################
  
  deform_coord_30d<-cmdscale(wdist,k=30)
  
  
  ######################################################
  ###### Choosing the optimal value of psi #############
  ######################################################
  
  ###################################################################
  ####### Defining the normalized mean squared error function #######
  ###################################################################
  
  nmse<-function(true,pred)
  {
    return(1-(sum((true-pred)^2))/(sum((true-mean(true))^2)))
  }
  
  nmse_wdist<-numeric(length = 30)
  nmse_wdist[1]<-NA
  for(i in 2:30)
  {
    tempdist<-rdist(deform_coord_30d[,1:i])
    nmse_wdist[i]<-nmse(true = wdist,pred = tempdist)
  }
  plot(1:30,nmse_wdist)
  
  
  ###### choosing the psi for which the NMSE is closest to 1 ######
  psi<-which.max(-abs(nmse_wdist-1))
  #nmse_wdist[psi]
  
  
  
  
  ##### first re_ordering deformed coordinates ######
  deform_ordered<-deform_coord_30d[ordering.index,]
  
  
  
  ##### Now, we plot the true deformation and the estimated deformation #####
  deform_ordered3<-deform_ordered[,1:2] #### Considering only the first two dimensions of maximum variation for visulaization
  par(mfrow=c(1,1),pty="s")
  quilt.plot(deform_ordered3[,1:2],error,nx=N,ny=N)
  par(pty="s")
  plot(deform_ordered3[,1:2])
  
  par(mfrow=c(1,2),pty="s")
  plot(myel2.2)
  ##### Rotating and translating the estimated deformed coordinates ####
  yscale<-min(myel2.2[,2])-min(deform_ordered3[,1])
  xscale<-min(myel2.2[,1])-min(deform_ordered3[,2])
  my_r.def<-cbind(deform_ordered3[,2],deform_ordered3[,1])
  my_r.def[,2]<-my_r.def[,2]+yscale
  my_r.def[,1]<-my_r.def[,1]+xscale
  par(mfrow=c(1,2))
  
  plot(myel2.2,pch=19,cex=0.2)
  plot(my_r.def,pch=19,cex=0.2)
  quilt.plot(myel2.2,error)
  quilt.plot(my_r.def,error)
  
  ####################################################################################
  ########### Now fitting a stationary Matérn in the deformed space ##################
  ########### on training data #######################################################
  ####################################################################################
  
  
  
  mle.comp_only_mlv.def<-function(par)
  {
    return(mle.comp.all(locs=deform_ordered[-rmiss,1:psi],z=error[-rmiss],p=par)$mlv
    )
  }
  optim_marg_comp.loglik.def <- function(par){
    optim(par=par,
          fn = mle.comp_only_mlv.def,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=1000))
  }
  
  mle.comp_only_mlv.geo<-function(par)
  {
    return(mle.comp.all(locs=grid[-rmiss,-3],z=error[-rmiss],p=par)$mlv
    )
  }
  optim_marg_comp.loglik.geo <- function(par){
    optim(par=par,
          fn = mle.comp_only_mlv.geo,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=1000))
  }
  
  
  f.ini<-rbind(temp.r1$par,temp.r2$par)
  f.ini<-colMeans(f.ini)
  mle.def<-optim_marg_comp.loglik.def(par = f.ini)
  #mle.def
  
  mle.geo<-optim_marg_comp.loglik.geo(par = f.ini)
  #mle.geo
  
  
  ####### computing fitted nonstationary/stationary covariance functions #######
  
  def_cov<-my.matern(h=rdist(deform_ordered[,1:psi]),a=mle.def$par[1],nu=mle.def$par[2],sigma = mle.def$par[3])
  geo_cov<-my.matern(h=rdist(cbind(grid$x,grid$y)),a=mle.geo$par[1],nu=mle.geo$par[2],sigma = mle.geo$par[3])
  
  par(mfrow=c(2,3),pty="s")
  # 4 is og
  set.seed(4)
  rlocs<-sample(1:N^2,3)
  quilt.plot(grid[,-3],def_cov[rlocs[1],]/def_cov[1,1],nx=N,ny=N)
  quilt.plot(grid[,-3],def_cov[rlocs[2],]/def_cov[1,1],nx=N,ny=N)
  quilt.plot(grid[,-3],def_cov[rlocs[3],]/def_cov[1,1],nx=N,ny=N)
  
  quilt.plot(grid[,-3],my_Cns[rlocs[1],],nx=N,ny=N)
  quilt.plot(grid[,-3],my_Cns[rlocs[2],],nx=N,ny=N)
  quilt.plot(grid[,-3],my_Cns[rlocs[3],],nx=N,ny=N)
  
  
  ##########################
  cov_vis<-data.frame(x=rep(grid$x,times=6),
                      y=rep(grid$y,times=6),corr=c(def_cov[rlocs[1],]/def_cov[1,1],def_cov[rlocs[2],]/def_cov[1,1],def_cov[rlocs[3],]/def_cov[1,1],
                             my_Cns[rlocs[1],]/my_Cns[1,1],my_Cns[rlocs[2],]/my_Cns[1,1],my_Cns[rlocs[3],]/my_Cns[1,1]),
                      locs=c(rep("Location 1",times=length(def_cov[rlocs[1],])),rep("Location 2",times=length(def_cov[rlocs[1],])),rep("Location 3",times=length(def_cov[rlocs[1],])),
                             rep("Location 1",times=length(def_cov[rlocs[1],])),rep("Location 2",times=length(def_cov[rlocs[1],])),rep("Location 3",times=length(def_cov[rlocs[1],]))),
                      model=rep(c("Estimated correlations","True correlations"),each=3*length(def_cov[rlocs[1],])))
  
  
  ggplot(cov_vis, aes(x=x, y=y))+geom_point(data = cov_vis, aes(x, y, colour = corr),pch=19, size=2) +
    facet_wrap(~model+locs)+scale_color_viridis("",option = "D")+labs(x="x",y="y")
  
  
  par(mfrow=c(1,1),pty="s")

  plot(grid[,-3],pch=19,cex=0.2)
  
  
  plot(myel2.2,pch=19,cex=0.2)  
  points(myel2.2,pch=19,cex=0.2,col=color.scale(error,col=viridis(n=5000)))
  
  plot(my_r.def,pch=19,cex=0.2,xlab="x",ylab="y") 
  points(my_r.def,pch=19,cex=0.2,col=color.scale(error,col=viridis(n=5000)))
  
  quilt.plot(grid[,-3],error,col = viridis(n=5000))
  abline(v=1)  
  
  
  ######################################################
  ############# Plots for the main manuscript ##########
  ######################################################
  
  par(mar=c(4,4,0.5,0.5),pty="s")

  plot(grid[,c(1,2)], pch=19,cex=0.2)  
  plot(myel2.2, pch=19,cex=0.2,xlab="x'",ylab="y'")  
  
  par(mar=c(4,4,0.5,0.5))
  quilt.plot(grid$x,grid$y,error,col = viridis(n=5000),horizontal = T,legend.width = 0.5,nx=N,ny=N,xlab="x",ylab="y")
  
  par(mar=c(7,4,0.1,0.1),pty="s")
  plot(myel2.2,pch=19,cex=0.2,xlab="x'",ylab="y'")  
  points(myel2.2,pch=19,cex=0.2,col=color.scale(error,col=viridis(n=5000)))
  quilt.plot(grid$x,grid$y,error,col = viridis(n=5000),horizontal = T,legend.width = 0.5,nx=N,ny=N,xlab="x",ylab="y",legend.only=T,legend.mar = 2)
  
  par(mar=c(7,4,0.1,0.1),pty="s")
  plot(my_r.def,pch=19,cex=0.2,xlab="x'",ylab="y'")  
  points(my_r.def,pch=19,cex=0.2,col=color.scale(error,col=viridis(n=5000)))
  quilt.plot(grid$x,grid$y,error,col = viridis(n=5000),horizontal = T,legend.width = 0.5,nx=N,ny=N,xlab="x",ylab="y",legend.only=T,legend.mar = 2)
  
  