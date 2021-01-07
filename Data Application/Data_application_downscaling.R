##########################################################
###### Setting working directory and attaching data ######
##########################################################
setwd("/Users/qadirga/Documents/Project 1 (Deformation)/Technometrics submission/Revision/Revised codes/Data Application/Application Section (Section 4 of the main manuscript)/New Data Appliation code")
library( fields)
attach("RData.COmonthly.met")  ####### Attaching dataset ######
CO.id                          ####### Station ids
CO.loc<-CO.loc                 ####### Station locations
CO.ppt                         ####### precipitation data
library(fdasrvf)               ####### Loading libraries
library(geoR)
library(fields)
library(plot3D)
library(rgl)
library(scatterplot3d)
library(scoringRules)
library(mvtnorm)

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


#################################################################
#### Now we prepare a dataset for yearly total precipitation ####
#################################################################


yr<-rep(1895:1997,each=12)
ppt.data<-cbind(yr,CO.ppt)
annual.ppt<-matrix(NA,nrow=367,ncol = length(1895:1997))
for(j in 1895:1997)
{
  temp.data<-ppt.data[ppt.data[,1]==j,]
  temp.data<-temp.data[,-1]
  for(i in 1:367)
  {
    temp<-sum(!is.na(temp.data[,i]))
    if(temp==12)
    {
      annual.ppt[i,(j-1894)]<-sum(temp.data[,i])
    }
    else
      annual.ppt[i,(j-1894)]<-NA 
  }
}

colnames(annual.ppt)<-1895:1997
annual.ppt<-as.data.frame(annual.ppt)
annual.ppt$`1895`  #### Total annual precipitation for the year 1895 
annual.ppt$`1896`  #### Total annual precipitation for the year 1896

####################################################################
##### Lets check the number of observed stations in each year ######
####################################################################
yrs<-1895:1997
num.stations<-numeric(length = 103)
for(i in 1:103)
{
  num.stations[i]<-sum(!is.na(annual.ppt[,i]))
}

stn.dets<-cbind(yrs,num.stations)
###### Number of stations in a given year #############
stn.dets  

#################################################################################################
####### Now we do the analysis for year 1992 since it has the highest number of observation #####
#################################################################################################
my.year=which.max(stn.dets[,2]) #(year taken is 1992)#


##### creating a dataset in a simple format with log transformation of precipitation data ####
my.data<-data.frame(x=CO.loc[,1],y=CO.loc[,2],ppt=log(annual.ppt[,my.year]))
my.raw.data<-data.frame(x=CO.loc[,1],y=CO.loc[,2],ppt=log(annual.ppt[,my.year]))
quilt.plot(my.data) ###### Plotting data #####
#### We partitioning from the line x=-104.873 ###
mdline<--104.873
abline(v=mdline)


my.data<-na.omit(my.data) ###### Removing the NA values, or the unobserved station locations 
my.data$ppt<-(my.data$ppt-mean(my.data$ppt))/sd(my.data$ppt) ######### Standardizing the data
my.data<-my.data[order(my.data$y,my.data$x),]

quilt.plot(my.data,xlab="Longitude",ylab="Latitude")  ###### Plotting the standardized data
abline(v=mdline,lwd=1)

##### arranging data as reg1 data (x<partition line) and reg2 data (x>=partition line) ####
reg1<-my.data[my.data$x<mdline,]
reg2<-my.data[my.data$x>=mdline,]
full.data<-rbind(reg1,reg2)

reord.index<-order(full.data$y,full.data$x)


##### Checking normal assumptions ########
par(mfrow=c(1,2))
colnames(my.data)<-colnames(reg1)<-colnames(reg2)<-colnames(full.data)<-c("x","y","z")
hist(full.data$z)
qqnorm(full.data$z)
qqline(full.data$z)


#######################################################################
######### Creating a fine spatial grid for prediction #################
#######################################################################

xmin<-min(my.data$x)
xmax<-max(my.data$x)
ymin<-min(my.data$y)
ymax<-max(my.data$y)
grid.res<-50 ### grid resolution ###
x.grid<-seq(xmin,xmax,length.out = grid.res)
y.grid<-seq(ymin,ymax,length.out = grid.res)
f.grid<-expand.grid(x.grid,y.grid)

par(mfrow=c(1,1))
plot(f.grid,xlab="x",ylab="y")
quilt.plot(my.data,add = T)

##########################################################################
############### Now we estimate regional variograms ######################
##########################################################################

####################################################################
################ plotting regionwise empirical variograms ##########
####################################################################

vr1<-variog(coords = cbind(reg1$x,reg1$y),data=reg1$z,option = "bin",uvec=30,max.dist=(10)/2)
plot(vr1,main="Empirical Variogram for Region 1")

vr2<-variog(coords = cbind(reg2$x,reg2$y),data=reg2$z,uvec=30,max.dist = (10)/2,option = "bin")
plot(vr2,main="Empirical Variogram for Region 2")

##### First we fit the WLS variogram, and then use its estimates as initial values in MLE ####
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

temp.r1<-optim(my.var.loss,par = c(runif(n=1,min=0.1,max=30),0.5,sd(reg1$z)))
temp.r1

########## Now using MLE to estimate Matérn covariance parameters #######


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
  return(mle.comp.all(locs=cbind(reg1$x,reg1$y),z=reg1$z,p=par)$mlv
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
  return(mle.comp.all(locs=cbind(reg2$x,reg2$y),z=reg2$z,p=par)$mlv
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
temp.r2<-optim(my.var.loss,par = c(runif(n=1,min=0.1,max=30),0.5,sd(reg2$z)))
temp.r2

mle.r2<-optim_marg_comp.loglik.r2(par = temp.r2$par)
mle.r2


##################################################################################
######### Plotting regional variograms overlayed on empirical variograms #########
##################################################################################
lseq<-seq(0,max(vr1$u,vr2$u),length.out = 200)
plot(vr1$u,vr1$v,pch=19,col="black",ylim = c(0,max(vr1$v,vr2$v)),
   xlim = c(0,max(vr1$u,vr2$u)),
 main = "Empirical and Fitted variograms",xlab="||h||",ylab = bquote(gamma(h)))
points(vr2$u,vr2$v,pch=19,col="red")
lines(lseq,my.matern(h=lseq,a=mle.r1$par[1],nu=mle.r1$par[2],sigma = mle.r1$par[3],variog = T),
    col="black",lty=1)
lines(lseq,my.matern(h=lseq,a=mle.r2$par[1],nu=mle.r2$par[2],sigma = mle.r2$par[3],variog = T),
    col="red",lty=2)
legend("bottomright",lty=c(4,5),col=c("black","red"),c("Region1","Region2"))


#############################################
######### Aligning Variograms ###############
#############################################



#### variable time represents the distances at which the regional variograms are evaluated for registration ###

time<-seq(0,10,length=1000)

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

par(mfrow=c(1,3))
plot(vr1$u,vr1$v,pch=19,col="black",ylim = c(0,max(vr1$v,vr2$v)),
   xlim = c(0,max(vr1$u,vr2$u)),
   main = "Empirical and Fitted variograms",xlab="||h||",ylab = bquote(gamma(h)))
points(vr2$u,vr2$v,pch=19,col="red")
lines(time,my.matern(h=time,a=mle.r1$par[1],nu=mle.r1$par[2],sigma = mle.r1$par[3],variog = T),
    col="black",lty=1)
lines(time,my.matern(h=time,a=mle.r2$par[1],nu=mle.r2$par[2],sigma = mle.r2$par[3],variog = T),
    col="red",lty=2)
legend("bottomright",lty=c(4,5),col=c("black","red"),c("Region1","Region2"))
plot(time,my_out$fn[,1],type="l",ylim=c(0,1.4),main="Aligned Variograms",ylab=expression(paste(gamma,"(h)")),xlab="h")
lines(time,my_out$fn[,2],type="l",col="red",lty=2)
legend("bottomright",lty=c(1,2),col=c("black","red"),c("Region1","Region2"))

plot(time,scale_inv1,type="l",main="Local Distance Warping Functions",xlab="h",ylab=expression(paste(phi['i'],"(h)")))
lines(time,scale_inv2,lty=2,col="red")
legend("bottomright",lty=c(1,2),col=c("black","red"),c("Region1","Region2"))


#plot(vr1$u,vr1$v,pch=19,col="black",ylim = c(0,max(vr1$v,vr2$v)),
 #    xlim = c(0,max(vr1$u,vr2$u)),
  #   main = "Empirical and Fitted variograms",xlab="||h||",ylab = bquote(gamma(h)))
#points(vr2$u,vr2$v,pch=19,col="red")
#lines(time,my.matern(h=time,a=mle.r1$par[1],nu=mle.r1$par[2],sigma = mle.r1$par[3],variog = T)/mle.r1$par[3]^2,
#      col="black",lty=1)
#lines(time,my.matern(h=time,a=mle.r2$par[1],nu=mle.r2$par[2],sigma = mle.r2$par[3],variog = T)/mle.r2$par[3]^2,
 #     col="red",lty=2)
#legend("bottomright",lty=c(4,5),col=c("black","red"),c("Region1","Region2"))
#plot(time,my_out$fn[,1]/max(my_out$fn[,1]),type="l",ylim=c(0,1.4),main="Aligned Variograms",ylab=expression(paste(gamma,"(h)")),xlab="h")
#lines(time,my_out$fn[,2]/max(my_out$fn[,2]),type="l",col="red",lty=2)
#legend("bottomright",lty=c(1,2),col=c("black","red"),c("Region1","Region2"))

#plot(time,scale_inv1,type="l",main="Local Distance Warping Functions",xlab="h",ylab=expression(paste(phi['i'],"(h)")))
#lines(time,scale_inv2,lty=2,col="red")
#legend("bottomright",lty=c(1,2),col=c("black","red"),c("Region1","Region2"))

##############################################################################
########## b1 is the distance matrix for locations in subregion 1 ############
##############################################################################

### Dividing the f.grid into two regions #####
colnames(f.grid)<-c("x","y")
f.reg1<-f.grid[f.grid$x<mdline,]
f.reg2<-f.grid[f.grid$x>=mdline,]

w.reg1<-rbind(f.reg1,reg1[,-3])
w.reg2<-rbind(f.reg2,reg2[,-3])

b1<-rdist(cbind(w.reg1$x,w.reg1$y))
#max(b1)
#####################################################################################################
############ Creating local distance warping function for subregion 1 with kernel smoothing #########
############# Gaussian Kernel with bandwidth = 0.02 #################################################
#####################################################################################################



dwarp1<-function(h)
{
  
  h2<-matrix(NA,nrow=nrow(h),ncol=ncol(h))
  checkf<-data.frame(rows=c(row(h)), columns=c(col(h)),
                     values=c(h))
  checkf<-checkf[order(checkf$values),]
  h2[cbind(checkf$rows,checkf$columns)]<-ksmooth(x=time,y=scale_inv1,kernel = "normal",bandwidth = 0.02,x.points = checkf$values)$y
  return(h2)
  
}


####################################################################################################
############# R11 is the warped distance matrix for subregion 1 ####################################
####################################################################################################
#plot(time,scale_inv1)

R11<-dwarp1(h=b1)
diag(R11)<-0 #### forcing it to be zero, since kernel smoothing may have changed is neglibly different from zero


########## b2 is the distance matrix for locations in subregion 2 ############

b2<-rdist(cbind(w.reg2$x,w.reg2$y))

#####################################################################################################
############ Creating local distance warping function for subregion 2 with kernel smoothing #########
#####################################################################################################


dwarp2<-function(h)
{
  h2<-matrix(NA,nrow=nrow(h),ncol=ncol(h))
  checkf<-data.frame(rows=c(row(h)), columns=c(col(h)),
                     values=c(h))
  checkf<-checkf[order(checkf$values),]
  h2[cbind(checkf$rows,checkf$columns)]<-ksmooth(x=time,y=scale_inv2,kernel = "normal",bandwidth = 0.02,x.points = checkf$values)$y
  return(h2)
  
}

R22<-dwarp2(h=b2)
diag(R22)<-0 #### forcing it to be zero, since kernel smoothing may have changed is neglibly different from zero


#####################################################################################################################################
########## Now computing the matrix R12, which is the warped distances for the two points in region 1 and region 2 respectively #####
#####################################################################################################################################


midpoint.f<-function(i,j,xmid=mdline)
{
  x1<-w.reg1[i,]$x
  y1<-w.reg1[i,]$y
  x2<-w.reg2[j,]$x
  y2<-w.reg2[j,]$y
  ymid<-y1+((y2-y1)/(x2-x1))*(xmid-x1)
  return(ymid)
}

y.inter<-outer(1:length(w.reg1$x),1:length(w.reg2$x),midpoint.f)



dwarp12<-function(region1=w.reg1,region2=w.reg2,y.intersect=y.inter,xmid=mdline)
{
  lr1<-length(region1$x)
  lr2<-length(region2$x)
  h1<-t(apply(matrix(1:lr1,ncol=1),1,function(i) rdist(cbind(region1$x[i],region1$y[i]),cbind(xmid,y.intersect[i,]))))
  h2<-apply(matrix(1:lr2,ncol=1),1,function(i) rdist(cbind(xmid,y.intersect[,i]),cbind(region2$x[i],region2$y[i])))
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


psi<-which.max(-abs(nmse_wdist-1))

####################################################
########## Reordering deform coordinates ###########
####################################################
fr11.ln<-length(f.reg1$x)
fr12.ln<-length(reg1$x)
fr21.ln<-length(f.reg2$x)
fr22.ln<-length(reg2$x)

new.order<-c(1:fr11.ln,(fr11.ln+fr12.ln+1):(fr11.ln+fr12.ln+fr21.ln),(fr11.ln+1):(fr11.ln+fr12.ln),(fr11.ln+fr12.ln+fr21.ln+1):(fr11.ln+fr12.ln+fr21.ln+fr22.ln))

deform_ordered<-deform_coord_30d[new.order,]
deform_ordered.grid<-deform_ordered[1:length(f.grid$x),]

grid.reorder.index<-order(c(f.reg1$y,f.reg2$y),c(f.reg1$x,f.reg2$x))

deform_ordered.grid2<-deform_ordered.grid[grid.reorder.index,]
deform_ordered.grid2<-rbind(deform_ordered.grid2,deform_ordered[(length(f.grid$x)+1):(length(f.grid$x)+length(my.data$x)),])

#### Rotating for visualization purpose (around x axis) #####
pp=pi/1
deform_ordered3<-deform_ordered.grid2[,1:3]
deform_ordered3[,2]<-deform_ordered.grid2[,2]*cos(pp)-deform_ordered.grid2[,3]*sin(pp)
deform_ordered3[,3]<-deform_ordered.grid2[,2]*sin(pp)+deform_ordered.grid2[,3]*cos(pp)
#split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
s3<-scatterplot3d(x=deform_ordered3[,1],y=deform_ordered3[,2],z=deform_ordered3[,3],pch=19,cex.symbols = 0.8,color = "white",grid=T,box=FALSE,xlab="x'",ylab="y'",zlab="z'",angle = 10)
N=grid.res
for(i in 1:N)
{
 s3$points3d(x=deform_ordered3[((i-1)*N+1):(i*N),1],y=deform_ordered3[((i-1)*N+1):(i*N),2],z=deform_ordered3[((i-1)*N+1):(i*N),3],type = "l")

}

xvals<-numeric()
yvals<-numeric()
zvals<-numeric()
for(j in 1:N)
{
 for(i in 1:N)
{
 xvals[i]<-deform_ordered3[(i-1)*N+j,1]
yvals[i]<-deform_ordered3[(i-1)*N+j,2]
zvals[i]<-deform_ordered3[(i-1)*N+j,3]
}
s3$points3d(x=xvals,y=yvals,z=zvals,type = "l")

}
#(length(f.grid$x)+1):(length(f.grid$x)+length(my.data$x))
s3$points3d(x=deform_ordered3[(length(f.grid$x)+1):(length(f.grid$x)+length(my.data$x)),1],y=deform_ordered3[(length(f.grid$x)+1):(length(f.grid$x)+length(my.data$x)),2],z=deform_ordered3[(length(f.grid$x)+1):(length(f.grid$x)+length(my.data$x)),3],pch=19,cex = 0.5,col=color.scale(c(reg1$z,reg2$z),col=tim.colors(n=254),zlim=c(min(c(reg1$z,reg2$z))-0.001,max(c(reg1$z,reg2$z))+0.001)))

image.plot(legend.only = T,zlim=c(min(c(reg1$z,reg2$z))-0.001,max(c(reg1$z,reg2$z))+0.001),horizontal = T,legend.width = 0.5,legend.mar = c(1.0,5.1,5.1,2.1))


#######################################################################################################
########### Now fitting a stationary Matérn with nugget effect in the deformed space ##################
################## on entire data #####################################################################
#######################################################################################################


mle.comp_only_mlv.def<-function(par)
{
  return(mle.comp.all(locs=deform_ordered[(length(f.grid$x)+1):(length(f.grid$x)+length(my.data$x)),1:psi],z=c(reg1$z,reg2$z),p=par)$mlv
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
  return(mle.comp.all(locs=cbind(c(reg1$x,reg2$x),c(reg1$y,reg2$y)),z=c(reg1$z,reg2$z),p=par)$mlv
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


#### setting initial values as the average WLS estimates of the two regions #####
f.ini<-rbind(temp.r1$par,temp.r2$par)
f.ini<-colMeans(f.ini)

mle.def<-optim_marg_comp.loglik.def(par = f.ini)
mle.def

mle.geo<-optim_marg_comp.loglik.geo(par = f.ini)
mle.geo


####### Plotting nonstationary covariance function #######
h_geo<-rbind(f.reg1,f.reg2,reg1[,c(1,2)],reg2[,c(1,2)])



def_cov<-my.matern(h=rdist(deform_ordered[,1:psi]),a=mle.def$par[1],nu=mle.def$par[2],sigma = mle.def$par[3])
geo_cov<-my.matern(h=rdist(cbind(h_geo$x,h_geo$y)),a=mle.geo$par[1],nu=mle.geo$par[2],sigma = mle.geo$par[3])

par(mfrow=c(4,4))

rand.pl.index<-sample(1:(length(f.grid$x)+length(my.data$x)),16)

for(pl in 1:16){
quilt.plot(h_geo,def_cov[rand.pl.index[pl],]/def_cov[1,1],nx=N,ny=N)}


####### Now we do prediction in the geographic space ##########
rmiss<-1:length(f.grid$x)


geo.sig11<-geo_cov[rmiss,rmiss]
geo.sig22<-geo_cov[-rmiss,-rmiss]
geo.sig12<-geo_cov[rmiss,-rmiss]
rearr.data<-c(reg1$z,reg2$z)
geo_wts<-geo.sig12%*%solve(geo.sig22)
geo_pred<-geo_wts%*%rearr.data ### Predicted value
geo_pvar<-diag(geo.sig11-geo_wts%*%t(geo.sig12)) #### Prediction variance

##### Plotting predicted surface with prediction variances in Geographic space #####
par(mfrow=c(2,2))
quilt.plot(h_geo[rmiss,],geo_pred,main="Predicted (G)",nx=N,ny=N)

quilt.plot(h_geo[rmiss,],geo_pvar,main="Prediction variance (G)",nx=N,ny=N)


##############################################################################
########### Now we do prediction in the deformed space #######################
##############################################################################


def.sig11<-def_cov[rmiss,rmiss]
def.sig22<-def_cov[-rmiss,-rmiss]
def.sig12<-def_cov[rmiss,-rmiss]


def_wts<-def.sig12%*%solve(def.sig22)
def_pred<-def_wts%*%rearr.data ### Predicted values
def_pvar<-diag(def.sig11-def_wts%*%t(def.sig12)) ### Prediction variances

quilt.plot(h_geo[rmiss,],def_pred,main="Predicted (D)",nx=N,ny=N)
quilt.plot(h_geo[rmiss,],def_pvar,main="Prediction variance (D)",nx=N,ny=N)

###############################################################################
############### Now we compute plots for the manuscript #######################
###############################################################################

library(ggmap)
library(viridis)
###### Loading map for Colorado #######

##### Log transformed data plot #####

load("map.col.RData")
my.raw.data<-na.omit(my.raw.data)
par(mfrow=c(1,1))
ggmap(map.col)+
    scale_y_continuous(limits = c(ymin,ymax))+
    scale_x_continuous(limits = c(xmin,xmax))+
  geom_point(data = my.raw.data, aes(x, y, colour = ppt),pch=19, size=2,alpha=1) +
  scale_color_viridis("Log (precipitation (mm))", option = "D",guide = guide_colourbar(direction = "horizontal",barwidth = 20, barheight = 0.75,title.position="top")) +
  labs(x="Longitude (degrees)", y="Latitude (degrees)")+theme(legend.position = "bottom") 



############### Juxtaposition of Geographic space and deformed space ############

########## Geographic space ##########
par(mar=c(4,4,0.5,0.5))
quilt.plot(my.data,col = viridis(n=254),horizontal = T,xlim = c(xmin,xmax),ylim=c(ymin,ymax),ylab = "Latitude (degrees)",xlab = "Longitude (degrees)")

#plot(NA,NA,xlim = c(xmin,xmax),ylim=c(ymin,ymax),ylab = "Latitude (degrees)",xlab = "Longitude (degrees)")
abline(v=x.grid,col="grey")
abline(h=y.grid,col="grey")
abline(v=mdline)

#quilt.plot(my.data,col = viridis(n=254),add = T,horizontal = T)

######## Deformed space #####
par(mar=c(4,4,0.5,0.5))
s3<-scatterplot3d(x=deform_ordered3[,1],y=deform_ordered3[,2],z=deform_ordered3[,3],pch=19,cex.symbols = 0.8,color = "white",grid=T,box=FALSE,xlab="x'",ylab="y'",zlab="z'",angle = 10)
for(i in 1:N)
{
  s3$points3d(x=deform_ordered3[((i-1)*N+1):(i*N),1],y=deform_ordered3[((i-1)*N+1):(i*N),2],z=deform_ordered3[((i-1)*N+1):(i*N),3],type = "l",col="grey")
  
}

xvals<-numeric()
yvals<-numeric()
zvals<-numeric()
for(j in 1:N)
{
  for(i in 1:N)
  {
    xvals[i]<-deform_ordered3[(i-1)*N+j,1]
    yvals[i]<-deform_ordered3[(i-1)*N+j,2]
    zvals[i]<-deform_ordered3[(i-1)*N+j,3]
  }
  s3$points3d(x=xvals,y=yvals,z=zvals,type = "l",col="grey")
  
}
#(length(f.grid$x)+1):(length(f.grid$x)+length(my.data$x))
s3$points3d(x=deform_ordered3[(length(f.grid$x)+1):(length(f.grid$x)+length(my.data$x)),1],y=deform_ordered3[(length(f.grid$x)+1):(length(f.grid$x)+length(my.data$x)),2],z=deform_ordered3[(length(f.grid$x)+1):(length(f.grid$x)+length(my.data$x)),3],pch=15,cex = 0.5,col=color.scale(c(reg1$z,reg2$z),col=viridis(n=254),zlim=c(min(c(reg1$z,reg2$z))-0.001,max(c(reg1$z,reg2$z))+0.001)))

image.plot(legend.only = T,zlim=c(min(c(reg1$z,reg2$z))-0.001,max(c(reg1$z,reg2$z))+0.001),horizontal = T,legend.width = 0.5,col = viridis(n=254),legend.mar = 2 )



#################################################################################################
########################### Now plotting alignment results ######################################
#################################################################################################

par(mar=c(4,4,0.2,6))
plot(NA,NA,xlim = c(xmin,xmax),ylim=c(ymin,ymax),ylab = "Latitude (degrees)",xlab = "Longitude (degrees)")
abline(v=x.grid,col="grey")
abline(h=y.grid,col="grey")
quilt.plot(my.data,col = viridis(n=254),add = T,horizontal = F,legend.width = 0.5)
abline(v=mdline)

####### Standardised estimated variogram ########
std.est.var<-data.frame(h=c(time,time),
                        v=c(my.matern(h=time,a=mle.r1$par[1],nu=mle.r1$par[2],sigma = mle.r1$par[3],variog = T)/max(my.matern(h=time,a=mle.r1$par[1],nu=mle.r1$par[2],sigma = mle.r1$par[3],variog = T)),
                                           my.matern(h=time,a=mle.r2$par[1],nu=mle.r2$par[2],sigma = mle.r2$par[3],variog = T)/max(my.matern(h=time,a=mle.r2$par[1],nu=mle.r2$par[2],sigma = mle.r2$par[3],variog = T))),
                        Subregion=c(rep("Western subregion",times=length(time)),rep("Eastern subregion",times=length(time))))

ggplot(std.est.var, aes(x=h, y=v, fill=Subregion)) + 
  geom_line(aes(x=h,y=v,col=Subregion,lty=Subregion))+
  labs(x=expression(paste("||",bold(h),"||")),y=expression(paste(hat(gamma)[i],"(||",bold(h),"||)")))+
  theme(legend.position="bottom")

####### Standardised aligned variogram ########
align.var<-data.frame(h=c(time,time),
                        v=c(my_out$fn[,1]/max(my_out$fn[,1]),
                            my_out$fn[,2]/max(my_out$fn[,2])),
                        Subregion=c(rep("Western subregion",times=length(time)),rep("Eastern subregion",times=length(time))))

ggplot(align.var, aes(x=h, y=v, fill=Subregion)) + 
  geom_line(aes(x=h,y=v,col=Subregion,lty=Subregion))+
  labs(x=expression(paste("||",bold(h),"||")),y=expression(paste(hat(gamma)[i],"(",phi[i]^-1,"(||",bold(h),"||))")))+
  theme(legend.position="bottom")

####### Regional distance warping functions #######

regdist.war<-data.frame(h=c(time,time),
                      v=c(scale_inv1,
                          scale_inv2),
                      Subregion=c(rep("Western subregion",times=length(time)),rep("Eastern subregion",times=length(time))))

ggplot(regdist.war, aes(x=h, y=v, fill=Subregion)) + 
  geom_line(aes(x=h,y=v,col=Subregion,lty=Subregion))+
  labs(x=expression(paste("||",bold(h),"||")),y=expression(paste(phi['i'],"(||",bold(h),"||)")))+
  theme(legend.position="bottom")+geom_abline(intercept = 0, slope = 1, color="black",linetype="dotdash", size=0.4)



###################################################################################################
################# Plotting the estimated nonstationary correlation function #######################
###################################################################################################



pl_nons.corr<-def_cov[1:length(f.grid$x),1:length(f.grid$x)]/def_cov[1,1] ### considering only grid points for better visualization
pl_h_geo<-h_geo[1:length(f.grid$x),]

set.seed(10)
rand.pl.index<-sample(1:2500,10)
temp_cat<-numeric()
temp_locs<-NULL
for(i in 1:10)
{
  temp_cat<-c(temp_cat,pl_nons.corr[rand.pl.index[i],])
  temp_locs<-c(temp_locs,paste("Location ",i))
}


data_heat<-data.frame(x=rep(pl_h_geo$x,times=10),
                      y=rep(pl_h_geo$y,times=10),
                      z=temp_cat,
                      locs=rep(temp_locs,each=length(pl_h_geo$x)))
data_heat$locs<-factor(data_heat$locs,levels = c("Location  1","Location  2","Location  3","Location  4","Location  5","Location  6", 
                                                 "Location  7","Location  8","Location  9","Location  10" ))
ggplot(data_heat, aes(x=x, y=y))+geom_point(data = data_heat, aes(x, y, colour = z),pch=19, size=2) +
  facet_wrap(~locs,ncol=5)+scale_color_viridis("",option = "D")+labs(x="Longitude",y="Latitude")+theme(axis.text = element_text(size=9))



##############################################################################
############ Prediction surfaces and prediction standard deviations ##########
##############################################################################

krv.range<-range(c(geo_pred,def_pred))
krsd.range<-range(sqrt(c(geo_pvar,def_pvar)))
par(mar=c(4,4,0.5,0.5))
quilt.plot(h_geo[rmiss,],geo_pred,zlim=c(krv.range),nx=N,ny=N,col = viridis(n=3000),horizontal=T)

quilt.plot(h_geo[rmiss,],sqrt(geo_pvar),zlim=krsd.range,horizontal=T,nx=N,ny=N,col = viridis(n=3000))

quilt.plot(h_geo[rmiss,],def_pred,zlim=c(krv.range),horizontal=T,nx=N,ny=N,col = viridis(n=3000))
quilt.plot(h_geo[rmiss,],sqrt(def_pvar),zlim=krsd.range,horizontal=T,nx=N,ny=N,col = viridis(n=3000))



