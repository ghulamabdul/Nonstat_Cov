########################################################################
######## Empirical study for the robustness of deformation method ######
########################################################################

###########################################
####### Matern covariance function ########
###########################################

my.matern<-function(h,a,sigma,nu)
{
  h[h==0]<-1e-10
  num1<-(sigma^2)*(2^(1-nu))/gamma(nu)
  num2<-(2*sqrt(nu)*h*a)^nu
  num3<-besselK(x=(2*sqrt(nu)*h*a), nu=nu)
  return(num1*num2*num3)
}


######################################
######### Loading libraries ##########
######################################

library(fdasrvf)
library(geoR)
library(fields)
library(RandomFields)
library(plot3D)
library(fields)
library(rgl)
library(scatterplot3d)
library(scoringRules)

###############################################################
############ N is the number of points in each dimension ######
###############################################################

N=70
grdsize<-8
x<-seq(0,grdsize,length=N)
y<-seq(0,grdsize,length=N)


##################################
#### Creating Simulation grid ####
##################################

coord.x<-rep(x,times=N)
coord.y<-rep(y,each=N)

grid<-data.frame(x=coord.x,y=coord.y)


###########################################
#### Now assigning the parameter values####
###########################################

range.region1<-(0.1)^2
range.region2<-(0.45)^2
nu<-0.5
sigma<-1

param.nu<-rep(nu,times=(N^2))
param.sigma<-rep(sigma,times=(N^2))


#######################################################
####### setting range parameter for two regions #######
#######################################################
range<-numeric()

range[grid$x<(grdsize/2)]<-range.region1
range[grid$x>(grdsize/2)]<-range.region2


########################################################
########## Visualizing the parameters ##################
########################################################



par(mfrow=c(1,3))
quilt.plot(grid,param.nu,nx=N,ny=N,main="Smoothness")
quilt.plot(grid,param.sigma,nx=N,ny=N,main="Standard deviation")
quilt.plot(grid,range,nx=N,ny=N,main="Spatial Range")

################################################################################
############## Creating a dataframe for each parameter on each location ########
################################################################################

total.frame<-data.frame(x=coord.x,y=coord.y,param.range=range,param.nu=param.nu,param.sigma=param.sigma)




#########################################################################
######### Now we generate kernel matrices for each locations ############
#########################################################################

kernel<-array(dim=c(2,2,(N^2)))
phi<-0
P<-matrix(c(cos(phi),-sin(phi),sin(phi),cos(phi)),2,2,byrow = T)
P
for(i in 1:(N^2))
{
  kernel[,,i]<-P%*%matrix(c(range[i],0,0,range[i]),2,2,byrow=T)%*%t(P)
}

############################################################
############### Computing mahalanobis distances ############
############################################################

Q<-matrix(NA,(N^2),(N^2))
for(i in 1:(N^2))
{
  for(j in 1:(N^2))
  { si<-rbind(total.frame[i,1],total.frame[i,2])
  sj<-rbind(total.frame[j,1],total.frame[j,2])
  p1<-t(si-sj)
  p2<-(si-sj)
  mid<-(kernel[,,i]+kernel[,,j])/2
  A<-solve(mid)
  temp<-p1%*%A%*%p2
  Q[i,j]<-temp
  }
}
##########################################################################
############## Computing non-stationary matern covariance matrix #########
##########################################################################


Cns<-matrix(NA,(N^2),(N^2))
for(i in 1:(N^2))
{
  for(j in 1:(N^2))
  {
    part1<-param.sigma[i]*param.sigma[j]
    part2d<-gamma(nu)*(2^(nu-1))
    part2<-1/part2d
    part3a<-(det(kernel[,,i]))^(1/4)
    part3b<-(det(kernel[,,j]))^(1/4)
    part3<-part3a*part3b
    part4p<-(kernel[,,i]+kernel[,,j])/2
    part4d<-det(part4p)
    part4<-part4d^(-1/2)
    part5<-(2*sqrt(nu*Q[i,j]))^nu
    part6<-besselK(x=(2*sqrt(nu*Q[i,j])),nu=nu)
    Cns[i,j]<-part1*part2*part3*part4*part5*part6
  }
  
  part1<-param.sigma[i]*param.sigma[i]
  part3a<-(det(kernel[,,i]))^(1/4)
  part3b<-(det(kernel[,,i]))^(1/4)
  part3<-part3a*part3b
  part4p<-(kernel[,,i]+kernel[,,i])/2
  part4d<-det(part4p)
  part4<-part4d^(-1/2)
  Cns[i,i]<-part1*part3*part4
}



####### Computing 50 simulations from the specified covariance function ######
full_sim<-matrix(NA,nrow=N^2,ncol=52)
full_sim[,1]<-coord.x
full_sim[,2]<-coord.y
library(MASS)
set.seed(123)
sims60<-mvrnorm(n=60,mu=rep(0,times=N^2),Sigma=Cns)
for(i in 1:50)
{
  set.seed(i)
  full_sim[,i+2]<-sims60[i,]
}
par(mfrow=c(1,1))
par(cex.axis=2, cex.lab=3, cex.main=1, cex.sub=1,mar=c(4,4,1,4)+.1)
quilt.plot(total.frame[,1],total.frame[,2],full_sim[,52],ncol=N,nrow=N,pch=1,nx=N,ny=N,xlab="x",ylab="y")

abline(v=4)
text(x=2,y=4,"True Region 1",cex=3)
text(x=6,y=4,"True Region 2",cex=3)

###################################################################
####### Now we do multiple MLE for trueregion 1 in parallel #######
###################################################################
dist.mat.tr1<-rdist(treg1array[1,,c(1,2)])

init_tr1<-c(1,1/sqrt(range.region1))
tr1estimates<-foreach(i=1:50)%dopar%{
  z<-treg1array[i,,3]
  neg_loglikelihood_tr1<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.tr1,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_tr1,fn=neg_loglikelihood_tr1,control = list(maxit=3000,trace=6))
}


hs<-seq(0,sqrt(8),length=200)
varmattr1<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmattr1[,i]<-(my.matern(h=0,a=tr1estimates[[i]]$par[2],sigma=tr1estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=tr1estimates[[i]]$par[2],sigma=tr1estimates[[i]]$par[1],nu=0.5))/((tr1estimates[[i]]$par[1]^2))
  lines(hs,varmattr1[,i])
}




###############################################################
####### Now we do multiple MLE for region 2 in parallel #######
###############################################################
dist.mat.tr2<-rdist(treg2array[1,,c(1,2)])
init_tr2<-c(1,1/sqrt(range.region2))
tr2estimates<-foreach(i=1:50)%dopar%{
  z<-treg2array[i,,3]
  neg_loglikelihood_tr2<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.tr2,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_tr2,fn=neg_loglikelihood_tr2,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmattr2<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmattr2[,i]<-(my.matern(h=0,a=tr2estimates[[i]]$par[2],sigma=tr2estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=tr2estimates[[i]]$par[2],sigma=tr2estimates[[i]]$par[1],nu=0.5))/((tr2estimates[[i]]$par[1]^2))
  lines(hs,varmattr2[,i])
}

treg1_warp<-treg2_warp<-matrix(NA,nrow=200,ncol=50)

for(i in 1:50)
{
  time<-seq(0,sqrt(8),length=200)
  #par(mfrow=c(1,1))
  #plot(time,RFvariogram(fitr1,dist=time,dim=2)/fitr1@ml@globalvariance,type="l",col=1)
  #lines(time,RFvariogram(fitr2,dist=time,dim=2)/fitr2@ml@globalvariance,col=2)
  #lines(time,RFvariogram(fitr3,dist=time,dim=2)/fitr3@ml@globalvariance,col=3)
  #lines(time,RFvariogram(fitr4,dist=time,dim=2)/fitr4@ml@globalvariance,col=4)
  #legend("bottomright",lty=c(1,1,1,1),col=c(1,2,3,4),c("Region1","Region2","Region3","Region4"))
  
  my_data<-data.frame(f1=varmattr1[,i],f2=varmattr2[,i])
  my_out<-time_warping(f=as.matrix(my_data),time = time)
  #plot(my_out)
  
  inv1<-invertGamma(my_out$gam[,1])
  inv2<-invertGamma(my_out$gam[,2])
  
  
  scale_inv1<-max(time)*inv1
  scale_inv2<-max(time)*inv2
  
  treg1_warp[,i-2]<-scale_inv1
  treg2_warp[,i-2]<-scale_inv2
  
  
}

# Augmenting identity warping for the distances greater than sqrt(8) and less than equal to sqrt(128)
time2<-seq(sqrt(8)+0.05,sqrt(128),length = 600 )
par(mar=c(4,6,1,1))
par(mfrow=c(1,1))

par(cex.axis=2.5, cex.lab=2.5, cex.main=1, cex.sub=1,mar=c(4,6,0,0)+.1)
plot(NA,NA,xlim=c(0,(sqrt(128)/2)),ylim=c(0,(sqrt(128)/2)),xlab="||h||",ylab=expression(paste(phi['i'],"(||h||)")))
for(i in 1:50)
{
  lines(c(time,time2),c(treg1_warp[,i],time2),col=1)
  lines(c(time,time2),c(treg2_warp[,i],time2),col=2)
}
legend("bottomright",lty = c(1,1),col=c(1,2),c("True Region 1","True Region 2"),cex=1.3)









par(cex.axis=2, cex.lab=3, cex.main=1, cex.sub=1,mar=c(4,4,1,4)+.1)
quilt.plot(total.frame[,1],total.frame[,2],full_sim[,52],ncol=N,nrow=N,pch=1,nx=N,ny=N,xlab="x",ylab="y")

abline(v=grdsize/2)
abline(h=grdsize/2)
############# Drawing partitioning line ###########
#lines(x=c(1,1),y=c(0,grdsize/2),lwd=2)
#abline(h=1)
text(x=grdsize/4,y=grdsize*3/4,"G1",cex = 3)
text(x=grdsize*3/4,y=grdsize*3/4,"G2",cex = 3)
text(x=grdsize/4,y=grdsize/4,"G3",cex = 3)
text(x=grdsize*3/4,y=grdsize/4,"G4",cex = 3)


reg1_var<-reg2_var<-reg3_var<-reg4_var<-reg1_warp<-reg2_warp<-reg3_warp<-reg4_warp<-matrix(NA,nrow=200,ncol=50)
library(doParallel)
registerDoParallel(cores=30)
getDoParWorkers()
reg4array<-reg3array<-reg2array<-reg1array<-array(dim=c(50,((N^2)/4),3))
(N^2)/4
for(i in 1:50)
{
  reg1array[i,,]<-full_sim[coord.x<(grdsize/2)&coord.y>(grdsize/2),c(1,2,i+2)]
  reg2array[i,,]<-full_sim[coord.x>(grdsize/2)&coord.y>(grdsize/2),c(1,2,i+2)]
  reg3array[i,,]<-full_sim[coord.x<(grdsize/2)&coord.y<(grdsize/2),c(1,2,i+2)]
  reg4array[i,,]<-full_sim[coord.x>(grdsize/2)&coord.y<(grdsize/2),c(1,2,i+2)]
}
par(mfrow=c(2,2))#### Test plots
quilt.plot(reg1array[50,,],nx=35,ny=35,zlim=c(min(full_sim[,52]),max(full_sim[,52])))
quilt.plot(reg2array[50,,],nx=35,ny=35,zlim=c(min(full_sim[,52]),max(full_sim[,52])))
quilt.plot(reg3array[50,,],nx=35,ny=35,zlim=c(min(full_sim[,52]),max(full_sim[,52])))
quilt.plot(reg4array[50,,],nx=35,ny=35,zlim=c(min(full_sim[,52]),max(full_sim[,52])))




###############################################################
####### Now we do multiple MLE for region 1 in parallel #######
###############################################################
dist.mat.r1<-rdist(reg1array[1,,c(1,2)])

neg_loglikelihood_r1<-function(p)
{
  sigma1<-p[1]
  a1<-p[2]
  nu1<-0.5
  
  if(sigma1<=0||a1<=0||nu1<=0)
  {
    return(10000000)
  }
  else{
    
    C<-my.matern(h=dist.mat.r1,a=a1,sigma=sigma1,nu=nu1)
    eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
    eigen_values <- eigen_decomp$values # Eigenvalues
    
    eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
    C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
    nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
    if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
    return(nlog_likelihood)
  }
}
init_r1<-c(1,1/sqrt(range.region1))
r1estimates<-foreach(i=1:50)%dopar%{
  z<-reg1array[i,,3]
  neg_loglikelihood_r1<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r1,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r1,fn=neg_loglikelihood_r1,control = list(maxit=3000,trace=6))
}

r1estimates[[1]]$par[2]
hs<-seq(0,sqrt(8),length=200)
varmatr1<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr1[,i]<-(my.matern(h=0,a=r1estimates[[i]]$par[2],sigma=r1estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r1estimates[[i]]$par[2],sigma=r1estimates[[i]]$par[1],nu=0.5))/((r1estimates[[i]]$par[1]^2))
  lines(hs,varmatr1[,i])
}




###############################################################
####### Now we do multiple MLE for region 2 in parallel #######
###############################################################
dist.mat.r2<-rdist(reg2array[1,,c(1,2)])

neg_loglikelihood_r2<-function(p)
{
  sigma1<-p[1]
  a1<-p[2]
  nu1<-0.5
  
  if(sigma1<=0||a1<=0||nu1<=0)
  {
    return(10000000)
  }
  else{
    
    C<-my.matern(h=dist.mat.r2,a=a1,sigma=sigma1,nu=nu1)
    eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
    eigen_values <- eigen_decomp$values # Eigenvalues
    
    eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
    C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
    nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
    if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
    return(nlog_likelihood)
  }
}
init_r2<-c(1,1/sqrt(range.region2))
r2estimates<-foreach(i=1:50)%dopar%{
  z<-reg2array[i,,3]
  neg_loglikelihood_r2<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r2,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r2,fn=neg_loglikelihood_r2,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr2<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr2[,i]<-(my.matern(h=0,a=r2estimates[[i]]$par[2],sigma=r2estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r2estimates[[i]]$par[2],sigma=r2estimates[[i]]$par[1],nu=0.5))/((r2estimates[[i]]$par[1]^2))
  lines(hs,varmatr2[,i])
}

###############################################################
####### Now we do multiple MLE for region 3 in parallel #######
###############################################################
dist.mat.r3<-rdist(reg3array[1,,c(1,2)])

neg_loglikelihood_r3<-function(p)
{
  sigma1<-p[1]
  a1<-p[2]
  nu1<-0.5
  
  if(sigma1<=0||a1<=0||nu1<=0)
  {
    return(10000000)
  }
  else{
    
    C<-my.matern(h=dist.mat.r3,a=a1,sigma=sigma1,nu=nu1)
    eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
    eigen_values <- eigen_decomp$values # Eigenvalues
    
    eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
    C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
    nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
    if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
    return(nlog_likelihood)
  }
}
init_r3<-c(1,1/sqrt(range.region1))
r3estimates<-foreach(i=1:50)%dopar%{
  z<-reg3array[i,,3]
  neg_loglikelihood_r3<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r3,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  
  optim(par = init_r3,fn=neg_loglikelihood_r3,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr3<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr3[,i]<-(my.matern(h=0,a=r3estimates[[i]]$par[2],sigma=r3estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r3estimates[[i]]$par[2],sigma=r3estimates[[i]]$par[1],nu=0.5))/((r3estimates[[i]]$par[1]^2))
  lines(hs,varmatr3[,i])
}

###############################################################
####### Now we do multiple MLE for region 3 in parallel #######
###############################################################
dist.mat.r4<-rdist(reg4array[1,,c(1,2)])

neg_loglikelihood_r4<-function(p)
{
  sigma1<-p[1]
  a1<-p[2]
  nu1<-0.5
  
  if(sigma1<=0||a1<=0||nu1<=0)
  {
    return(10000000)
  }
  else{
    
    C<-my.matern(h=dist.mat.r4,a=a1,sigma=sigma1,nu=nu1)
    eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
    eigen_values <- eigen_decomp$values # Eigenvalues
    
    eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
    C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
    nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
    if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
    return(nlog_likelihood)
  }
}
init_r4<-c(1,1/sqrt(range.region2))
r34estimates<-foreach(i=1:50)%dopar%{
  z<-reg4array[i,,3]
  neg_loglikelihood_r4<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r4,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  
  optim(par = init_r4,fn=neg_loglikelihood_r4,control = list(maxit=3000,trace=6))
}
r4estimates<-r34estimates
hs<-seq(0,sqrt(8),length=200)
varmatr4<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr4[,i]<-(my.matern(h=0,a=r4estimates[[i]]$par[2],sigma=r4estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r4estimates[[i]]$par[2],sigma=r4estimates[[i]]$par[1],nu=0.5))/((r4estimates[[i]]$par[1]^2))
  lines(hs,varmatr4[,i])
}




for(i in 1:50)
{
  time<-seq(0,sqrt(8),length=200)
  #par(mfrow=c(1,1))
  #plot(time,RFvariogram(fitr1,dist=time,dim=2)/fitr1@ml@globalvariance,type="l",col=1)
  #lines(time,RFvariogram(fitr2,dist=time,dim=2)/fitr2@ml@globalvariance,col=2)
  #lines(time,RFvariogram(fitr3,dist=time,dim=2)/fitr3@ml@globalvariance,col=3)
  #lines(time,RFvariogram(fitr4,dist=time,dim=2)/fitr4@ml@globalvariance,col=4)
  #legend("bottomright",lty=c(1,1,1,1),col=c(1,2,3,4),c("Region1","Region2","Region3","Region4"))
  
  my_data<-data.frame(f1=varmatr1[,i],f2=varmatr2[,i],f3=varmatr3[,i],f4=varmatr4[,i])
  my_out<-time_warping(f=as.matrix(my_data),time = time)
  #plot(my_out)
  
  inv1<-invertGamma(my_out$gam[,1])
  inv2<-invertGamma(my_out$gam[,2])
  inv3<-invertGamma(my_out$gam[,3])
  inv4<-invertGamma(my_out$gam[,4])
  
  
  scale_inv1<-max(time)*inv1
  scale_inv2<-max(time)*inv2
  scale_inv3<-max(time)*inv3
  scale_inv4<-max(time)*inv4
  
  reg1_warp[,i-2]<-scale_inv1
  reg2_warp[,i-2]<-scale_inv2
  reg3_warp[,i-2]<-scale_inv3
  reg4_warp[,i-2]<-scale_inv4
  
}

par(mar=c(4,6,1,1))
par(mfrow=c(1,1))
par(cex.axis=2.5, cex.lab=2.5, cex.main=1, cex.sub=1,mar=c(4,6,0,0)+.1)

plot(NA,NA,xlim=c(0,(sqrt(128)/2)),ylim=c(0,(sqrt(128)/2)),xlab="||h||",ylab=expression(paste(phi['i'],"(||h||)")))
for(i in 1:50)
{
  lines(c(time,time2),c(reg1_warp[,i],time2),col=1)
  lines(c(time,time2),c(reg2_warp[,i],time2),col=2)
  lines(c(time,time2),c(reg3_warp[,i],time2),col=3)
  lines(c(time,time2),c(reg4_warp[,i],time2),col=4)
}
legend("bottomright",lty = c(1,1,1,1),col=c(1,2,3,4),c("G1","G2", "G3", "G4"),cex=2.0)




###### L2 distances between warping functions #####
d1<-sqrt((reg1_warp-reg2_warp)^2)
d2<-sqrt((reg1_warp-reg3_warp)^2)
d3<-sqrt((reg1_warp-reg4_warp)^2)
d4<-sqrt((reg2_warp-reg3_warp)^2)
d5<-sqrt((reg2_warp-reg4_warp)^2)
d6<-sqrt((reg3_warp-reg4_warp)^2)

m1<-rowMeans(d1)
m2<-rowMeans(d2)
m3<-rowMeans(d3)
m4<-rowMeans(d4)
m5<-rowMeans(d5)
m6<-rowMeans(d6)

plot(time,m1,type="l",ylim=c(0,2),xlab="||h||",ylab="L2 distance")
lines(time,m2,col=2)
lines(time,m3,col=3)
lines(time,m4,col=4)
lines(time,m5,col=5)
lines(time,m6,col=6)
legend("topright",lty = c(1,1,1,1,1,1),col=c(1,2,3,4,5,6),c("G1 and G2","G1 and G3","G1 and G4","G2 and G3","G2 and G4","G3 and G4"),cex=1.4)




#################################################################################
############# Boxplot for pairwise local distance warping functions  ############
#################################################################################
bp1<-colMeans(d1)
bp2<-colMeans(d2)
bp3<-colMeans(d3)
bp4<-colMeans(d4)
bp5<-colMeans(d5)
bp6<-colMeans(d6)

par(mfrow=c(2,3))
boxplot(bp1,main="Region1 and Region2")
boxplot(bp2,main="Region1 and Region3")
boxplot(bp3,main="Region1 and Region4")
boxplot(bp4,main="Region2 and Region3")
boxplot(bp5,main="Region2 and Region4")
boxplot(bp6,main="Region3 and Region4")


bpdata<-data.frame(RMSE=c(bp1,bp2,bp3,bp4,bp5,bp6),Pairs=c(rep("(G1,G2)",times=length(bp1)),rep("(G1,G3)",times=length(bp2)),rep("(G1,G4)",times=length(bp3)),rep("(G2,G3)",times=length(bp4)),rep("(G2,G4)",times=length(bp5)),rep("(G3,G4)",times=length(bp6))))
library(ggplot2)
ggplot(bpdata, aes(x=Pairs, y=RMSE)) +
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme_classic()
# Change box plot colors by groups
p<-ggplot(bpdata, aes(x=Pairs, y=RMSE, fill=Pairs)) +
  geom_boxplot()+ theme(legend.text=element_text(size=17),axis.text=element_text(size=16),
                        axis.title=element_text(size=16))
p




########## Now further we divide the subregion #########
par(mfrow=c(1,1))
par(cex.axis=2, cex.lab=3, cex.main=1, cex.sub=1,mar=c(4,4,1,4)+.1)
quilt.plot(total.frame[,1],total.frame[,2],full_sim[,52],ncol=N,nrow=N,pch=1,nx=N,ny=N,xlab="x",ylab="y")

abline(v=grdsize/2)
abline(h=grdsize/2)
abline(v=grdsize/4)
abline(h=grdsize/4)
abline(v=3*grdsize/4)
abline(h=3*grdsize/4)

############# Drawing partitioning line ###########
#lines(x=c(1,1),y=c(0,grdsize/2),lwd=2)
#abline(h=1)
text(x=1,y=7,"G11",cex = 3)
text(x=3,y=7,"G12",cex = 3)
text(x=1,y=5,"G13",cex = 3)
text(x=3,y=5,"G14",cex = 3)
text(x=5,y=7,"G21",cex = 3)
text(x=7,y=7,"G22",cex = 3)
text(x=5,y=5,"G23",cex = 3)
text(x=7,y=5,"G24",cex = 3)
text(x=1,y=3,"G31",cex = 3)
text(x=3,y=3,"G32",cex = 3)
text(x=1,y=1,"G33",cex = 3)
text(x=3,y=1,"G34",cex = 3)
text(x=5,y=3,"G41",cex = 3)
text(x=7,y=3,"G42",cex = 3)
text(x=5,y=1,"G43",cex = 3)
text(x=7,y=1,"G44",cex = 3)


full_sim[coord.x<(grdsize/4)&coord.y>(3*grdsize/4),c(1,2,i+2)]
reg11array<-reg22array<-reg33array<-reg44array<-array(dim=c(50,324,3))

reg14array<-reg41array<-reg32array<-reg23array<-array(dim=c(50,289,3))
reg12array<-reg13array<-reg21array<-reg24array<-reg31array<-reg34array<-reg42array<-reg43array<-array(dim=c(50,306,3))

for(i in 1:50)
{
  
  
  reg11array[i,,]<-reg1array[i,(reg1array[i,,1]<2&reg1array[i,,2]>6),]
  reg12array[i,,]<-reg1array[i,(reg1array[i,,1]>2&reg1array[i,,2]>6),]
  reg13array[i,,]<-reg1array[i,(reg1array[i,,1]<2&reg1array[i,,2]<6),]
  reg14array[i,,]<-reg1array[i,(reg1array[i,,1]>2&reg1array[i,,2]<6),]
  
  
  reg21array[i,,]<-reg2array[i,(reg2array[i,,1]<6&reg2array[i,,2]>6),]
  reg22array[i,,]<-reg2array[i,(reg2array[i,,1]>6&reg2array[i,,2]>6),]
  reg23array[i,,]<-reg2array[i,(reg2array[i,,1]<6&reg2array[i,,2]<6),]
  reg24array[i,,]<-reg2array[i,(reg2array[i,,1]>6&reg2array[i,,2]<6),]
  
  reg31array[i,,]<-reg3array[i,(reg3array[i,,1]<2&reg3array[i,,2]>2),]
  reg32array[i,,]<-reg3array[i,(reg3array[i,,1]>2&reg3array[i,,2]>2),]
  reg33array[i,,]<-reg3array[i,(reg3array[i,,1]<2&reg3array[i,,2]<2),]
  reg34array[i,,]<-reg3array[i,(reg3array[i,,1]>2&reg3array[i,,2]<2),]
  
  reg41array[i,,]<-reg4array[i,(reg4array[i,,1]<6&reg4array[i,,2]>2),]
  reg42array[i,,]<-reg4array[i,(reg4array[i,,1]>6&reg4array[i,,2]>2),]
  reg43array[i,,]<-reg4array[i,(reg4array[i,,1]<6&reg4array[i,,2]<2),]
  reg44array[i,,]<-reg4array[i,(reg4array[i,,1]>6&reg4array[i,,2]<2),]
  
}
par(mfrow=c(2,2))#### Test plots
quilt.plot(reg11array[50,,],nx=18,ny=18,zlim=c(min(full_sim[,52]),max(full_sim[,52])))
quilt.plot(reg12array[50,,],nx=18,ny=18,zlim=c(min(full_sim[,52]),max(full_sim[,52])))
quilt.plot(reg13array[50,,],nx=18,ny=18,zlim=c(min(full_sim[,52]),max(full_sim[,52])))
quilt.plot(reg14array[50,,],nx=18,ny=18,zlim=c(min(full_sim[,52]),max(full_sim[,52])))




############### Now we compute the likelihood routines ############
dist.mat.r11<-rdist(reg11array[1,,c(1,2)])
init_r11<-c(1,1/sqrt(range.region1))
r11estimates<-foreach(i=1:50)%dopar%{
  z<-reg11array[i,,3]
  neg_loglikelihood_r11<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r11,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r11,fn=neg_loglikelihood_r11,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr11<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr11[,i]<-(my.matern(h=0,a=r11estimates[[i]]$par[2],sigma=r11estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r11estimates[[i]]$par[2],sigma=r11estimates[[i]]$par[1],nu=0.5))/((r11estimates[[i]]$par[1]^2))
  lines(hs,varmatr11[,i])
}







dist.mat.r12<-rdist(reg12array[1,,c(1,2)])
init_r12<-c(1,1/sqrt(range.region1))
r12estimates<-foreach(i=1:50)%dopar%{
  z<-reg12array[i,,3]
  neg_loglikelihood_r12<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r12,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r12,fn=neg_loglikelihood_r12,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr12<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr12[,i]<-(my.matern(h=0,a=r12estimates[[i]]$par[2],sigma=r12estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r12estimates[[i]]$par[2],sigma=r12estimates[[i]]$par[1],nu=0.5))/((r12estimates[[i]]$par[1]^2))
  lines(hs,varmatr12[,i])
}





dist.mat.r13<-rdist(reg13array[1,,c(1,2)])
init_r13<-c(1,1/sqrt(range.region1))
r13estimates<-foreach(i=1:50)%dopar%{
  z<-reg13array[i,,3]
  neg_loglikelihood_r13<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r13,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r13,fn=neg_loglikelihood_r13,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr13<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr13[,i]<-(my.matern(h=0,a=r13estimates[[i]]$par[2],sigma=r13estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r13estimates[[i]]$par[2],sigma=r13estimates[[i]]$par[1],nu=0.5))/((r13estimates[[i]]$par[1]^2))
  lines(hs,varmatr13[,i])
}


dist.mat.r14<-rdist(reg14array[1,,c(1,2)])
init_r14<-c(1,1/sqrt(range.region1))
r14estimates<-foreach(i=1:50)%dopar%{
  z<-reg14array[i,,3]
  neg_loglikelihood_r14<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r14,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r14,fn=neg_loglikelihood_r14,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr14<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr14[,i]<-(my.matern(h=0,a=r14estimates[[i]]$par[2],sigma=r14estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r14estimates[[i]]$par[2],sigma=r14estimates[[i]]$par[1],nu=0.5))/((r14estimates[[i]]$par[1]^2))
  lines(hs,varmatr14[,i])
}








dist.mat.r21<-rdist(reg21array[1,,c(1,2)])
init_r21<-c(1,1/sqrt(range.region2))
r21estimates<-foreach(i=1:50)%dopar%{
  z<-reg21array[i,,3]
  neg_loglikelihood_r21<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r21,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r21,fn=neg_loglikelihood_r21,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr21<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr21[,i]<-(my.matern(h=0,a=r21estimates[[i]]$par[2],sigma=r21estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r21estimates[[i]]$par[2],sigma=r21estimates[[i]]$par[1],nu=0.5))/((r21estimates[[i]]$par[1]^2))
  lines(hs,varmatr21[,i])
}







dist.mat.r22<-rdist(reg22array[1,,c(1,2)])
init_r22<-c(1,1/sqrt(range.region2))
r22estimates<-foreach(i=1:50)%dopar%{
  z<-reg22array[i,,3]
  neg_loglikelihood_r22<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r22,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r22,fn=neg_loglikelihood_r22,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr22<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr22[,i]<-(my.matern(h=0,a=r22estimates[[i]]$par[2],sigma=r22estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r22estimates[[i]]$par[2],sigma=r22estimates[[i]]$par[1],nu=0.5))/((r22estimates[[i]]$par[1]^2))
  lines(hs,varmatr22[,i])
}





dist.mat.r23<-rdist(reg23array[1,,c(1,2)])
init_r23<-c(1,1/sqrt(range.region2))
r23estimates<-foreach(i=1:50)%dopar%{
  z<-reg23array[i,,3]
  neg_loglikelihood_r23<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r23,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r23,fn=neg_loglikelihood_r23,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr23<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr23[,i]<-(my.matern(h=0,a=r23estimates[[i]]$par[2],sigma=r23estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r23estimates[[i]]$par[2],sigma=r23estimates[[i]]$par[1],nu=0.5))/((r23estimates[[i]]$par[1]^2))
  lines(hs,varmatr23[,i])
}


dist.mat.r24<-rdist(reg24array[1,,c(1,2)])
init_r24<-c(1,1/sqrt(range.region2))
r24estimates<-foreach(i=1:50)%dopar%{
  z<-reg24array[i,,3]
  neg_loglikelihood_r24<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r24,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r24,fn=neg_loglikelihood_r24,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr24<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr24[,i]<-(my.matern(h=0,a=r24estimates[[i]]$par[2],sigma=r24estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r24estimates[[i]]$par[2],sigma=r24estimates[[i]]$par[1],nu=0.5))/((r24estimates[[i]]$par[1]^2))
  lines(hs,varmatr24[,i])
}








dist.mat.r31<-rdist(reg31array[1,,c(1,2)])
init_r31<-c(1,1/sqrt(range.region1))
r31estimates<-foreach(i=1:50)%dopar%{
  z<-reg31array[i,,3]
  neg_loglikelihood_r31<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r31,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r31,fn=neg_loglikelihood_r31,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr31<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr31[,i]<-(my.matern(h=0,a=r31estimates[[i]]$par[2],sigma=r31estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r31estimates[[i]]$par[2],sigma=r31estimates[[i]]$par[1],nu=0.5))/((r31estimates[[i]]$par[1]^2))
  lines(hs,varmatr31[,i])
}







dist.mat.r32<-rdist(reg32array[1,,c(1,2)])
init_r32<-c(1,1/sqrt(range.region1))
r32estimates<-foreach(i=1:50)%dopar%{
  z<-reg32array[i,,3]
  neg_loglikelihood_r32<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r32,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r32,fn=neg_loglikelihood_r32,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr32<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr32[,i]<-(my.matern(h=0,a=r32estimates[[i]]$par[2],sigma=r32estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r32estimates[[i]]$par[2],sigma=r32estimates[[i]]$par[1],nu=0.5))/((r32estimates[[i]]$par[1]^2))
  lines(hs,varmatr32[,i])
}





dist.mat.r33<-rdist(reg33array[1,,c(1,2)])
init_r33<-c(1,1/sqrt(range.region1))
r33estimates<-foreach(i=1:50)%dopar%{
  z<-reg33array[i,,3]
  neg_loglikelihood_r33<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r33,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r33,fn=neg_loglikelihood_r33,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr33<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr33[,i]<-(my.matern(h=0,a=r33estimates[[i]]$par[2],sigma=r33estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r33estimates[[i]]$par[2],sigma=r33estimates[[i]]$par[1],nu=0.5))/((r33estimates[[i]]$par[1]^2))
  lines(hs,varmatr33[,i])
}


dist.mat.r34<-rdist(reg34array[1,,c(1,2)])
init_r34<-c(1,1/sqrt(range.region1))
r34estimates<-foreach(i=1:50)%dopar%{
  z<-reg34array[i,,3]
  neg_loglikelihood_r34<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r34,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r34,fn=neg_loglikelihood_r34,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr34<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr34[,i]<-(my.matern(h=0,a=r34estimates[[i]]$par[2],sigma=r34estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r34estimates[[i]]$par[2],sigma=r34estimates[[i]]$par[1],nu=0.5))/((r34estimates[[i]]$par[1]^2))
  lines(hs,varmatr34[,i])
}










dist.mat.r41<-rdist(reg41array[1,,c(1,2)])
init_r41<-c(1,1/sqrt(range.region2))
r41estimates<-foreach(i=1:50)%dopar%{
  z<-reg41array[i,,3]
  neg_loglikelihood_r41<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r41,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r41,fn=neg_loglikelihood_r41,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr41<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr41[,i]<-(my.matern(h=0,a=r41estimates[[i]]$par[2],sigma=r41estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r41estimates[[i]]$par[2],sigma=r41estimates[[i]]$par[1],nu=0.5))/((r41estimates[[i]]$par[1]^2))
  lines(hs,varmatr41[,i])
}







dist.mat.r42<-rdist(reg42array[1,,c(1,2)])
init_r42<-c(1,1/sqrt(range.region2))
r42estimates<-foreach(i=1:50)%dopar%{
  z<-reg42array[i,,3]
  neg_loglikelihood_r42<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r42,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r42,fn=neg_loglikelihood_r42,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr42<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr42[,i]<-(my.matern(h=0,a=r42estimates[[i]]$par[2],sigma=r42estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r42estimates[[i]]$par[2],sigma=r42estimates[[i]]$par[1],nu=0.5))/((r42estimates[[i]]$par[1]^2))
  lines(hs,varmatr42[,i])
}





dist.mat.r43<-rdist(reg43array[1,,c(1,2)])
init_r43<-c(1,1/sqrt(range.region2))
r43estimates<-foreach(i=1:50)%dopar%{
  z<-reg43array[i,,3]
  neg_loglikelihood_r43<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r43,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r43,fn=neg_loglikelihood_r43,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr43<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr43[,i]<-(my.matern(h=0,a=r43estimates[[i]]$par[2],sigma=r43estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r43estimates[[i]]$par[2],sigma=r43estimates[[i]]$par[1],nu=0.5))/((r43estimates[[i]]$par[1]^2))
  lines(hs,varmatr43[,i])
}


dist.mat.r44<-rdist(reg44array[1,,c(1,2)])
init_r44<-c(1,1/sqrt(range.region2))
r44estimates<-foreach(i=1:50)%dopar%{
  z<-reg44array[i,,3]
  neg_loglikelihood_r44<-function(p)
  {
    sigma1<-p[1]
    a1<-p[2]
    nu1<-0.5
    
    if(sigma1<=0||a1<=0||nu1<=0)
    {
      return(10000000)
    }
    else{
      
      C<-my.matern(h=dist.mat.r44,a=a1,sigma=sigma1,nu=nu1)
      eigen_decomp <- eigen(C) # Eigen-decomposition of covariance matrix
      eigen_values <- eigen_decomp$values # Eigenvalues
      
      eigen_vector <- eigen_decomp$vectors # Matrix with eigenvectors
      C_inv <- eigen_vector %*% diag(1/eigen_values) %*% t(eigen_vector) #Inverse of covariance matrix
      nlog_likelihood <- 0.5*sum(log(eigen_values)) + 0.5*t(z) %*% C_inv %*% z #Negative log likelihood
      if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      return(nlog_likelihood)
    }
  }
  optim(par = init_r44,fn=neg_loglikelihood_r44,control = list(maxit=3000,trace=6))
}

hs<-seq(0,sqrt(8),length=200)
varmatr44<-matrix(NA,nrow=200,ncol=50)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,sqrt(8)),ylim=c(0,1.5))
for(i in 1:50)
{
  varmatr44[,i]<-(my.matern(h=0,a=r44estimates[[i]]$par[2],sigma=r44estimates[[i]]$par[1],nu=0.5)-my.matern(h=hs,a=r44estimates[[i]]$par[2],sigma=r44estimates[[i]]$par[1],nu=0.5))/((r44estimates[[i]]$par[1]^2))
  lines(hs,varmatr44[,i])
}


reg41_warp<-reg42_warp<-reg43_warp<-reg44_warp<-reg31_warp<-reg32_warp<-reg33_warp<-reg34_warp<-reg21_warp<-reg22_warp<-reg23_warp<-reg24_warp<-reg11_warp<-reg12_warp<-reg13_warp<-reg14_warp<-matrix(NA,nrow=200,ncol=50)


for(i in 1:50)
{
  time<-seq(0,sqrt(8),length=200)
  #par(mfrow=c(1,1))
  #plot(time,RFvariogram(fitr1,dist=time,dim=2)/fitr1@ml@globalvariance,type="l",col=1)
  #lines(time,RFvariogram(fitr2,dist=time,dim=2)/fitr2@ml@globalvariance,col=2)
  #lines(time,RFvariogram(fitr3,dist=time,dim=2)/fitr3@ml@globalvariance,col=3)
  #lines(time,RFvariogram(fitr4,dist=time,dim=2)/fitr4@ml@globalvariance,col=4)
  #legend("bottomright",lty=c(1,1,1,1),col=c(1,2,3,4),c("Region1","Region2","Region3","Region4"))
  
  my_data<-data.frame(f11=varmatr11[,i],f12=varmatr12[,i],f13=varmatr13[,i],f14=varmatr14[,i],
                      f21=varmatr21[,i],f12=varmatr22[,i],f23=varmatr23[,i],f14=varmatr24[,i],
                      f31=varmatr31[,i],f32=varmatr32[,i],f33=varmatr33[,i],f34=varmatr34[,i],
                      f41=varmatr41[,i],f42=varmatr42[,i],f43=varmatr43[,i],f44=varmatr44[,i])
  my_out<-time_warping(f=as.matrix(my_data),time = time)
  #plot(my_out)
  
  inv11<-invertGamma(my_out$gam[,1])
  inv12<-invertGamma(my_out$gam[,2])
  inv13<-invertGamma(my_out$gam[,3])
  inv14<-invertGamma(my_out$gam[,4])
  inv21<-invertGamma(my_out$gam[,5])
  inv22<-invertGamma(my_out$gam[,6])
  inv23<-invertGamma(my_out$gam[,7])
  inv24<-invertGamma(my_out$gam[,8])
  inv31<-invertGamma(my_out$gam[,9])
  inv32<-invertGamma(my_out$gam[,10])
  inv33<-invertGamma(my_out$gam[,11])
  inv34<-invertGamma(my_out$gam[,12])
  inv41<-invertGamma(my_out$gam[,13])
  inv42<-invertGamma(my_out$gam[,14])
  inv43<-invertGamma(my_out$gam[,15])
  inv44<-invertGamma(my_out$gam[,16])
  
  
  scale_inv11<-max(time)*inv11
  scale_inv12<-max(time)*inv12
  scale_inv13<-max(time)*inv13
  scale_inv14<-max(time)*inv14
  
  scale_inv21<-max(time)*inv21
  scale_inv22<-max(time)*inv22
  scale_inv23<-max(time)*inv23
  scale_inv24<-max(time)*inv24
  
  scale_inv31<-max(time)*inv31
  scale_inv32<-max(time)*inv32
  scale_inv33<-max(time)*inv33
  scale_inv34<-max(time)*inv34
  
  scale_inv41<-max(time)*inv41
  scale_inv42<-max(time)*inv42
  scale_inv43<-max(time)*inv43
  scale_inv44<-max(time)*inv44
  
  reg11_warp[,i-2]<-scale_inv11
  reg12_warp[,i-2]<-scale_inv12
  reg13_warp[,i-2]<-scale_inv13
  reg14_warp[,i-2]<-scale_inv14

  
  reg21_warp[,i-2]<-scale_inv21
  reg22_warp[,i-2]<-scale_inv22
  reg23_warp[,i-2]<-scale_inv23
  reg24_warp[,i-2]<-scale_inv24
  
  reg31_warp[,i-2]<-scale_inv31
  reg32_warp[,i-2]<-scale_inv32
  reg33_warp[,i-2]<-scale_inv33
  reg34_warp[,i-2]<-scale_inv34
  
  reg41_warp[,i-2]<-scale_inv41
  reg42_warp[,i-2]<-scale_inv42
  reg43_warp[,i-2]<-scale_inv43
  reg44_warp[,i-2]<-scale_inv44
  
  
    
}


par(mar=c(4,6,1,1))
plot(NA,NA,xlim=c(0,sqrt(128)/2),ylim=c(0,sqrt(128)/2),xlab="||h||",ylab=expression(paste(phi['i'],"(||h||)")))

for(i in 1:50)
{
  lines(c(time,time2),c(reg11_warp[,i],time2),col=1)
  lines(c(time,time2),c(reg12_warp[,i],time2),col=2)
  lines(c(time,time2),c(reg13_warp[,i],time2),col=3)
  lines(c(time,time2),c(reg14_warp[,i],time2),col=4)
  
  
  lines(c(time,time2),c(reg21_warp[,i],time2),col=5)
  lines(c(time,time2),c(reg22_warp[,i],time2),col=6)
  lines(c(time,time2),c(reg23_warp[,i],time2),col=7)
  lines(c(time,time2),c(reg24_warp[,i],time2),col=8)
  
  lines(c(time,time2),c(reg31_warp[,i],time2),col=9)
  lines(c(time,time2),c(reg32_warp[,i],time2),col=10)
  lines(c(time,time2),c(reg33_warp[,i],time2),col=11)
  lines(c(time,time2),c(reg34_warp[,i],time2),col=12)
  
  lines(c(time,time2),c(reg41_warp[,i],time2),col=13)
  lines(c(time,time2),c(reg42_warp[,i],time2),col=14)
  lines(c(time,time2),c(reg43_warp[,i],time2),col=15)
  lines(c(time,time2),c(reg44_warp[,i],time2),col=16)
  
  
}
legend("bottomright",lty = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),c("G11","G12", "G13", "G14","G21","G22","G23","G24","G31","G32","G33","G34","G41","G42","G43","G44"),cex=1.2)


