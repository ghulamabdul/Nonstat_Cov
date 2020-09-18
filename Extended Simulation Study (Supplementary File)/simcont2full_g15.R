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


library(doParallel)
ncores<-detectCores()
##################################################
########## Setting working directory #############
##################################################
setwd("/ibex/scratch/qadirga/Project_1/Two_region_trial_Ibex")
registerDoParallel(cores = ncores-4)
sim_summary<-foreach(sim.index=1:100,.packages = c("fdasrvf","geoR","fields","plot3D","rgl","scatterplot3d","scoringRules","mvtnorm"))%dopar%{
  
##################################################
########## Setting working directory #############
##################################################

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
#plot(grid)

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
##### Seting nonstationary Matérn parameters #####
##################################################

nu<-0.8
sigma<-sqrt(5)

#### Anchor ranges ####
range.region1<-(0.1)^2
range.region2<-(0.9)^2

range.region3<-(0.1)^2
range.region4<-(0.9)^2  

#######################################################
####### Computing continuously varying ranges #########
#######################################################

weight_avg<-function(tloc=grid,anch.locs=cbind(c(0.5,1.5,0.5,1.5),c(1.5,1.5,0.5,0.5)),
                     anch.ranges=c(range.region1,range.region2,range.region3,range.region4),
                     bwidth=0.5)
{
  t.a.dist<-rdist(tloc,anch.locs)
  gauss.part<-exp(-(t.a.dist^2)/bwidth)
  denoms<-rowSums(gauss.part)
  
  num.part<-gauss.part*matrix(anch.ranges,nrow=nrow(gauss.part),ncol=ncol(gauss.part),byrow = T)
  num.part2<-rowSums(num.part)
  return(num.part2/denoms)
  
}

#### setting bandwidth lambda ####
cont.ranges<-weight_avg(bwidth = 0.15)

#quilt.plot(grid,cont.ranges,nx=N,ny=N)

cross.region<-outer(cont.ranges,cont.ranges,function(x,y) (x+y)/2)

norm.const1<-outer(cont.ranges,cont.ranges,function(x,y) sqrt(x*y))
norm.const2<-norm.const1/cross.region
fulldist<-rdist(grid)  ##### Computing distance matrix
nonst.dist<-2*sqrt(nu)*fulldist/sqrt(cross.region) ##### computing nonstationary matern covariance function
my_Cns<-norm.const2*my.matern(h=nonst.dist,a=1,sigma = sigma,nu=nu)
##################################################
#### Assigning subregions to the locations #######
##################################################

grid$subr<-ifelse(grid$x<1,"R1","R2")


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

#par(mfrow=c(1,1))
#quilt.plot(grid[-3],error,nx=N,ny=N)


############# Drawing partitioning line ###########
#abline(v=1)



##################################################################
######## Splitting into training and testing datasets ############
##################################################################

##### training size = 1200
##### testing size = 3700

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
#plot(vr1,main="Empirical Variogram for Region 1")

vr2<-variog(coords = cbind(reg2.tr$x,reg2.tr$y),data=reg2.tr$z,uvec=30,max.dist = (sqrt(8))/2,option = "bin")
#plot(vr2,main="Empirical Variogram for Region 2")

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

#par(mfrow=c(3,2))
#plot(my_out)


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
midpoint.f<-function(i,j,xmid=1)
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
  h1<-t(apply(matrix(1:lr1,ncol=1),1,function(i) rdist(cbind(region1$x[i],region1$y[i]),cbind(1,y.intersect[i,]))))
  h2<-apply(matrix(1:lr2,ncol=1),1,function(i) rdist(cbind(1,y.intersect[,i]),cbind(region2$x[i],region2$y[i])))
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
#plot(1:30,nmse_wdist)


###### choosing the psi for which the NMSE is closest to 1 ######
psi<-which.max(-abs(nmse_wdist-1))
#nmse_wdist[psi]




##### first re_ordering deformed coordinates ######
deform_ordered<-deform_coord_30d[ordering.index,]



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

#par(mfrow=c(2,2))
#quilt.plot(grid[,-3],def_cov[2890,]/def_cov[1,1],nx=N,ny=N)
#quilt.plot(grid[,-3],my_Cns[2890,],nx=N,ny=N)
#quilt.plot(grid[,-3],def_cov[3128,]/def_cov[1,1],nx=N,ny=N)
#quilt.plot(grid[,-3],my_Cns[3128,],nx=N,ny=N)


####### Now we do prediction in the geographic space ##########

geo.sig11<-geo_cov[rmiss,rmiss]
geo.sig22<-geo_cov[-rmiss,-rmiss]
geo.sig12<-geo_cov[rmiss,-rmiss]

geo_wts<-geo.sig12%*%solve(geo.sig22)
geo_pred<-geo_wts%*%error[-rmiss] #### Prediction values
geo_pvar<-diag(geo.sig11-geo_wts%*%t(geo.sig12)) #### Prediction variances


##### scores_function compute all the prediction scores ######

scores_func<-function(true=error[rmiss],pred,predvar)
{
  rmse<-sqrt(mean((true-pred)^2))
  mae<-mean(abs(true-pred))
  nmse<-1-(sum((pred-true)^2))/(sum((true-mean(true))^2))
  mycrps<-mean(crps(true, "norm", mean = c(pred), sd = c(sqrt(predvar))))
  mylogs<-mean(logs(true, "norm", mean = c(pred), sd = c(sqrt(predvar))))
  return(list(RMSE=rmse,NMSE=nmse,MAE=mae,mCRPS=mycrps,mLogS=mylogs))
}

geo_scores<-scores_func(pred = geo_pred,predvar = geo_pvar) #### scores for the stationary model


##############################################################################
########### Now we do prediction in the deformed space #######################
##############################################################################


def.sig11<-def_cov[rmiss,rmiss]
def.sig22<-def_cov[-rmiss,-rmiss]
def.sig12<-def_cov[rmiss,-rmiss]


def_wts<-def.sig12%*%solve(def.sig22)
def_pred<-def_wts%*%error[-rmiss] ### prediction values
def_pvar<-diag(def.sig11-def_wts%*%t(def.sig12)) ### Prediction variances
def_scores<-scores_func(pred = def_pred,predvar = def_pvar) #prediction scores for the deformed space


#####################################################################
######## Exploring the accuracy of prediction uncertainty ###########
#####################################################################

######## Computing G-statistics and Wp, K ########
prob_seq<-seq(0.01,0.99,length.out = 4*99)
uncert_p<-function(pr,truev=error[rmiss],predval,predvar)
{
  li<-qnorm(p=(1-pr)/2,mean =predval,sd=sqrt(predvar) )
  ui<-qnorm(p=(1+pr)/2,mean =predval,sd=sqrt(predvar) )
  kj_p<-as.numeric(truev<ui & truev>li)
  k_p<-mean(kj_p)
  temp<-kj_p*(ui-li)
  w_p<-mean(temp)/k_p
  return(list(kjp=kj_p,kp=k_p,wp=w_p))
}

#### Goodness statistics ##
G_stat<-function(p,kp)
{
  ap<-as.numeric(kp>p)
  s1<-3*ap-2
  s2<-kp-p
  s<-s1*s2
  dp<-p[2]-p[1]
  return(1-dp*sum(s))
}


geoun<-defun<-list()
geowp<-geokp<-defwp<-defkp<-numeric(length = length(prob_seq))

for(i in 1:length(prob_seq))
{
  geoun[[i]]<-uncert_p(pr=prob_seq[i],predval =geo_pred,predvar = geo_pvar )
  geowp[i]<-geoun[[i]]$wp
  geokp[i]<-geoun[[i]]$kp
  
  defun[[i]]<-uncert_p(pr=prob_seq[i],predval =def_pred,predvar = def_pvar )
  defwp[i]<-defun[[i]]$wp
  defkp[i]<-defun[[i]]$kp
  
}
geo_G<-G_stat(p=prob_seq,kp=geokp)
def_G<-G_stat(p=prob_seq,kp=defkp)

list(geo_scores=geo_scores,def_scores=def_scores,geoG=geo_G,defG=def_G,
     geokp=geokp,geowp=geowp,defkp=defkp,defwp=defwp,
     def.coord=deform_ordered,mle.geo=mle.geo,mle.def=mle.def,
     data=error,
     grid=grid,
     psi=psi,range.region1=range.region1,range.region2=range.region2,range.region3=range.region3,
     range.region4=range.region4,
     geo_pred=geo_pred,
     geo_pvar=geo_pvar,
     def_pred=def_pred,
     def_pvar=def_pvar,
     scale_inv1=scale_inv1,
     scale_inv2=scale_inv2,
     valid.set=error[rmiss],
     r1.v=my.matern(h=time,a=mle.r1$par[1],nu=mle.r1$par[2],sigma = mle.r1$par[3],variog = T),
     r2.v=my.matern(h=time,a=mle.r2$par[1],nu=mle.r2$par[2],sigma = mle.r2$par[3],variog = T),
     time=time,
     my_out=my_out
)

}


save.image("conti2rsim1_010_090g_lamda_015.RData")

time<-seq(0,sqrt(8),length=1000)
prob_seq<-seq(0.01,0.99,length.out = 4*99)
geo_G<-def_G<-geo_nmse<-def_nmse<-geo_rmse<-def_rmse<-geo_mae<-def_mae<-geo_crps<-def_crps<-geo_logs<-def_logs<-numeric(length = 100)
geo_kp<-def_kp<-geo_wp<-def_wp<-matrix(NA,nrow=100,ncol=length(prob_seq))
geo_pred<-geo_pvar<-def_pred<-def_pvar<-valid.set<-list()
my_out<-r1.v<-r2.v<-list()

for(i in 1:100)
{
  geo_G[i]<-sim_summary[[i]]$geoG
  def_G[i]<-sim_summary[[i]]$defG
  geo_nmse[i]<-sim_summary[[i]]$geo_scores$NMSE
  def_nmse[i]<-sim_summary[[i]]$def_scores$NMSE
  geo_rmse[i]<-sim_summary[[i]]$geo_scores$RMSE
  def_rmse[i]<-sim_summary[[i]]$def_scores$RMSE
  geo_mae[i]<-sim_summary[[i]]$geo_scores$MAE
  def_mae[i]<-sim_summary[[i]]$def_scores$MAE
  geo_crps[i]<-sim_summary[[i]]$geo_scores$mCRPS
  def_crps[i]<-sim_summary[[i]]$def_scores$mCRPS
  geo_logs[i]<-sim_summary[[i]]$geo_scores$mLogS
  def_logs[i]<-sim_summary[[i]]$def_scores$mLogS
  geo_kp[i,]<-sim_summary[[i]]$geokp
  geo_wp[i,]<-sim_summary[[i]]$geowp
  def_kp[i,]<-sim_summary[[i]]$defkp
  def_wp[i,]<-sim_summary[[i]]$defwp
  geo_pred[[i]]<-sim_summary[[i]]$geo_pred
  def_pred[[i]]<-sim_summary[[i]]$def_pred
  geo_pvar[[i]]<-sim_summary[[i]]$geo_pvar
  def_pvar[[i]]<-sim_summary[[i]]$def_pvar
  valid.set[[i]]<-sim_summary[[i]]$valid.set
  r1.v[[i]]<-sim_summary[[i]]$r1.v
  r2.v[[i]]<-sim_summary[[i]]$r2.v
  my_out[[i]]<-sim_summary[[i]]$my_out
}


###### Plots for prediction scores ###
par(mfrow=c(2,3))


boxplot(geo_rmse,def_rmse,names = c("Stat","Nstat"),col = "gray",ylim=c(range(geo_rmse,def_rmse)+c(0,+0.01)),main="RMSE")
points(c(1,2),c(mean(geo_rmse),mean(def_rmse)),pch=19,col="red")
text(x=c(1,2),y=c(max(geo_rmse)+0.005,max(def_rmse)+0.005),round(c(mean(geo_rmse),mean(def_rmse)),2),col="blue")

boxplot(geo_nmse,def_nmse,names = c("Stat","Nstat"),col = "gray",ylim=c(range(geo_nmse,def_nmse)+c(0,+0.1)),main="NMSE")
points(c(1,2),c(mean(geo_nmse),mean(def_nmse)),pch=19,col="red")
text(x=c(1,2),y=c(max(geo_nmse)+0.05,max(def_nmse)+0.05),round(c(mean(geo_nmse),mean(def_nmse)),2),col="blue")

boxplot(geo_mae,def_mae,names = c("Stat","Nstat"),col = "gray",ylim=c(range(geo_mae,def_mae)+c(0,+0.01)),main="MAE")
points(c(1,2),c(mean(geo_mae),mean(def_mae)),pch=19,col="red")
text(x=c(1,2),y=c(max(geo_mae)+0.005,max(def_mae)+0.005),round(c(mean(geo_mae),mean(def_mae)),2),col="blue")


boxplot(geo_crps,def_crps,names = c("Stat","Nstat"),col = "gray",ylim=c(range(geo_crps,def_crps)+c(0,+0.01)),main="mCRPS")
points(c(1,2),c(mean(geo_crps),mean(def_crps)),pch=19,col="red")
text(x=c(1,2),y=c(max(geo_crps)+0.005,max(def_crps)+0.005),round(c(mean(geo_crps),mean(def_crps)),2),col="blue")

boxplot(geo_logs,def_logs,names = c("Stat","Nstat"),col = "gray",ylim=c(range(geo_logs,def_logs)+c(0,+0.05)),main="mLogS")
points(c(1,2),c(mean(geo_logs),mean(def_logs)),pch=19,col="red")
text(x=c(1,2),y=c(max(geo_logs)+0.017,max(def_logs)+0.017),round(c(mean(geo_logs),mean(def_logs)),2),col="blue")

boxplot(geo_G,def_G,names = c("Stat","Nstat"),col = "gray",ylim=c(range(geo_G,def_G)+c(0,+0.01)),main="G-stat")
points(c(1,2),c(mean(geo_G),mean(def_G)),pch=19,col="red")
text(x=c(1,2),y=c(max(geo_G)+0.005,max(def_G)+0.005),round(c(mean(geo_G),mean(def_G)),2),col="blue")

par(mfrow=c(1,2))
plot(prob_seq,colMeans(geo_kp),col="red",type="l",lty=3,xlab = "p-Probability interval",ylab="Proportion within probability interval",main="Accuracy plot")
lines(prob_seq,colMeans(def_kp),col="blue",lty=2)
lines(prob_seq,prob_seq,lty=1,col="black")
legend("bottomright",lty=c(3,2),col=c("red","blue"),c("Stationary","Nonstationary"))

plot(prob_seq,colMeans(geo_wp),col="red",type="l",lty=3,xlab = "p-Probability interval",ylab="Probability interval width",main="Average width plot")
lines(prob_seq,colMeans(def_wp),col="blue",lty=2)
legend("bottomright",lty=c(3,2),col=c("red","blue"),c("Stationary","Nonstationary"))


######################################################################
####### Computing accuracy and width plot over combined runs #########
######################################################################


test.set<-def.krig<-def.krig.var<-geo.krig<-geo.krig.var<-NULL

for(i in 1:100)
{
  test.set<-c(c(test.set),valid.set[[i]])
  def.krig<-c(c(def.krig),def_pred[[i]])
  def.krig.var<-c(c(def.krig.var),def_pvar[[i]])
  geo.krig<-c(c(geo.krig),geo_pred[[i]])
  geo.krig.var<-c(c(geo.krig.var),geo_pvar[[i]])
}
#prob_seq<-seq(0.01,0.99,length.out = 4*99)

uncert_p<-function(pr,truev=test.set,predval,predvar)
{
  li<-qnorm(p=(1-pr)/2,mean =predval,sd=sqrt(predvar) )
  ui<-qnorm(p=(1+pr)/2,mean =predval,sd=sqrt(predvar) )
  kj_p<-as.numeric(truev<ui & truev>li)
  k_p<-mean(kj_p)
  temp<-kj_p*(ui-li)
  w_p<-mean(temp)/k_p
  return(list(kjp=kj_p,kp=k_p,wp=w_p))
}

#### Goodness statistics ##
G_stat<-function(p,kp)
{
  ap<-as.numeric(kp>p)
  s1<-3*ap-2
  s2<-kp-p
  s<-s1*s2
  dp<-p[2]-p[1]
  return(1-dp*sum(s))
}


geoun<-defun<-list()
geowp<-geokp<-defwp<-defkp<-numeric(length = length(prob_seq))

for(i in 1:length(prob_seq))
{
  geoun[[i]]<-uncert_p(pr=prob_seq[i],predval =geo.krig,predvar = geo.krig.var )
  geowp[i]<-geoun[[i]]$wp
  geokp[i]<-geoun[[i]]$kp
  
  defun[[i]]<-uncert_p(pr=prob_seq[i],predval =def.krig,predvar = def.krig.var )
  defwp[i]<-defun[[i]]$wp
  defkp[i]<-defun[[i]]$kp
  
}


par(mfrow=c(1,2))
plot(prob_seq,geokp,col="red",type="l",lty=3,xlab = "p-Probability interval",ylab="Proportion within probability interval",main="Accuracy plot")
lines(prob_seq,defkp,col="blue",lty=2)
lines(prob_seq,prob_seq,lty=1,col="black")
legend("bottomright",lty=c(3,2),col=c("red","blue"),c("Stationary","Nonstationary"))

plot(prob_seq,geowp,col="red",type="l",lty=3,xlab = "p-Probability interval",ylab="Probability interval width",main="Average width plot")
lines(prob_seq,defwp,col="blue",lty=2)
legend("bottomright",lty=c(3,2),col=c("red","blue"),c("Stationary","Nonstationary"))



###################################################################
##################### Plots for the manuscript ####################
###################################################################
library(RColorBrewer)
library(ggplot2)
rounding=3
rmse.data<-data.frame(RMSE=c(geo_rmse,def_rmse),Model=c(rep("Stationary",times=100),rep("Nonstationary",times=100)))
rmse.mean.geo<-round(mean(geo_rmse),rounding)
rmse.sd.geo<-round(sd(geo_rmse),rounding)
rmse.mean.def<-round(mean(def_rmse),rounding)
rmse.sd.def<-round(sd(def_rmse),rounding)

ggplot(rmse.data, aes(x=Model, y=RMSE, fill=Model)) + 
  geom_boxplot(alpha=1) +
  scale_y_continuous(limits = c(min(rmse.data$RMSE),max(rmse.data$RMSE)+0.1))+
  labs(title = "Root mean squared error")+
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")+
  geom_label(aes(y = max(RMSE)+0.1, x = c(1), label = paste("mean=",rmse.mean.def)))+
  geom_label(aes(y = max(RMSE)+0.1, x = c(2), label = paste("mean=",rmse.mean.geo)))+
  geom_label(aes(y = max(RMSE)+0.05, x = c(1), label = paste("sd=",rmse.sd.def)))+
  geom_label(aes(y = max(RMSE)+0.05, x = c(2), label = paste("sd=",rmse.sd.geo)))

nmse.data<-data.frame(NMSE=c(geo_nmse,def_nmse),Model=c(rep("Stationary",times=100),rep("Nonstationary",times=100)))
nmse.mean.geo<-round(mean(geo_nmse),rounding)
nmse.sd.geo<-round(sd(geo_nmse),rounding)
nmse.mean.def<-round(mean(def_nmse),rounding)
nmse.sd.def<-round(sd(def_nmse),rounding)

ggplot(nmse.data, aes(x=Model, y=NMSE, fill=Model)) + 
  geom_boxplot(alpha=1) +
  scale_y_continuous(limits = c(min(nmse.data$NMSE),max(nmse.data$NMSE)+0.1))+
  labs(title = "Normalized mean squared error")+
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")+
  geom_label(aes(y = max(NMSE)+0.1, x = c(1), label = paste("mean=",nmse.mean.def)))+
  geom_label(aes(y = max(NMSE)+0.1, x = c(2), label = paste("mean=",nmse.mean.geo)))+
  geom_label(aes(y = max(NMSE)+0.040, x = c(1), label = paste("sd=",nmse.sd.def)))+
  geom_label(aes(y = max(NMSE)+0.040, x = c(2), label = paste("sd=",nmse.sd.geo)))

mae.data<-data.frame(MAE=c(geo_mae,def_mae),Model=c(rep("Stationary",times=100),rep("Nonstationary",times=100)))
mae.mean.geo<-round(mean(geo_mae),rounding)
mae.sd.geo<-round(sd(geo_mae),rounding)
mae.mean.def<-round(mean(def_mae),rounding)
mae.sd.def<-round(sd(def_mae),rounding)

ggplot(mae.data, aes(x=Model, y=MAE, fill=Model)) + 
  geom_boxplot(alpha=1) +
  scale_y_continuous(limits = c(min(mae.data$MAE),max(mae.data$MAE)+0.1))+
  labs(title = "Mean absolute error")+
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")+
  geom_label(aes(y = max(MAE)+0.1, x = c(1), label = paste("mean=",mae.mean.def)))+
  geom_label(aes(y = max(MAE)+0.1, x = c(2), label = paste("mean=",mae.mean.geo)))+
  geom_label(aes(y = max(MAE)+0.06, x = c(1), label = paste("sd=",mae.sd.def)))+
  geom_label(aes(y = max(MAE)+0.06, x = c(2), label = paste("sd=",mae.sd.geo)))

mcrps.data<-data.frame(mCRPS=c(geo_crps,def_crps),Model=c(rep("Stationary",times=100),rep("Nonstationary",times=100)))
mcrps.mean.geo<-round(mean(geo_crps),rounding)
mcrps.sd.geo<-round(sd(geo_crps),rounding)
mcrps.mean.def<-round(mean(def_crps),rounding)
mcrps.sd.def<-round(sd(def_crps),rounding)

ggplot(mcrps.data, aes(x=Model, y=mCRPS, fill=Model)) + 
  geom_boxplot(alpha=1) +
  scale_y_continuous(limits = c(min(mcrps.data$mCRPS),max(mcrps.data$mCRPS)+0.1))+
  labs(title = "Mean continuous ranked probability score")+
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")+
  geom_label(aes(y = max(mCRPS)+0.1, x = c(1), label = paste("mean=",mcrps.mean.def)))+
  geom_label(aes(y = max(mCRPS)+0.1, x = c(2), label = paste("mean=",mcrps.mean.geo)))+
  geom_label(aes(y = max(mCRPS)+0.06, x = c(1), label = paste("sd=",mcrps.sd.def)))+
  geom_label(aes(y = max(mCRPS)+0.06, x = c(2), label = paste("sd=",mcrps.sd.geo)))

mlogs.data<-data.frame(mLogS=c(geo_logs,def_logs),Model=c(rep("Stationary",times=100),rep("Nonstationary",times=100)))
mlogs.mean.geo<-round(mean(geo_logs),rounding)
mlogs.sd.geo<-round(sd(geo_logs),rounding)
mlogs.mean.def<-round(mean(def_logs),rounding)
mlogs.sd.def<-round(sd(def_logs),rounding)

ggplot(mlogs.data, aes(x=Model, y=mLogS, fill=Model)) + 
  geom_boxplot(alpha=1) +
  scale_y_continuous(limits = c(min(mlogs.data$mLogS),max(mlogs.data$mLogS)+0.1))+
  labs(title = "Mean logarithmic score")+
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")+
  geom_label(aes(y = max(mLogS)+0.1, x = c(1), label = paste("mean=",mlogs.mean.def)))+
  geom_label(aes(y = max(mLogS)+0.1, x = c(2), label = paste("mean=",mlogs.mean.geo)))+
  geom_label(aes(y = max(mLogS)+0.05, x = c(1), label = paste("sd=",mlogs.sd.def)))+
  geom_label(aes(y = max(mLogS)+0.05, x = c(2), label = paste("sd=",mlogs.sd.geo)))


G.data<-data.frame(Gstat=c(geo_G,def_G),Model=c(rep("Stationary",times=100),rep("Nonstationary",times=100)))
G.mean.geo<-round(mean(geo_G),rounding)
G.sd.geo<-round(sd(geo_G),rounding)
G.mean.def<-round(mean(def_G),rounding)
G.sd.def<-round(sd(def_G),rounding)

ggplot(G.data, aes(x=Model, y=Gstat, fill=Model)) + 
  geom_boxplot(alpha=1) +
  scale_y_continuous(limits = c(min(G.data$Gstat),max(G.data$Gstat)+0.1))+
  labs(title = "Goodness statistic")+
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")+
  geom_label(aes(y = max(Gstat)+0.1, x = c(1), label = paste("mean=",G.mean.def)))+
  geom_label(aes(y = max(Gstat)+0.1, x = c(2), label = paste("mean=",G.mean.geo)))+
  geom_label(aes(y = max(Gstat)+0.05, x = c(1), label = paste("sd=",G.sd.def)))+
  geom_label(aes(y = max(Gstat)+0.05, x = c(2), label = paste("sd=",G.sd.geo)))

####################################################################################
############### Now we do the accuracy plots and the average width plots ###########
####################################################################################

#kp.data<-data.frame(kp=c(colMeans(geo_kp),colMeans(def_kp)),p=c(prob_seq,prob_seq),Model=c(rep("Stationary",times=length(prob_seq)),rep("Nonstationary",times=length(prob_seq))))
#ggplot(kp.data, aes(x=p, y=kp, fill=Model)) + 
#  geom_line(aes(x=p,y=kp,col=Model,lty=Model))+geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.4)+
#  labs(title = "Accuracy plot",x="p-Probability interval",y="Proportion within probability interval")+
#  theme(legend.position="bottom")


#wp.data<-data.frame(wp=c(colMeans(geo_wp),colMeans(def_wp)),p=c(prob_seq,prob_seq),Model=c(rep("Stationary",times=length(prob_seq)),rep("Nonstationary",times=length(prob_seq))))
#ggplot(wp.data, aes(x=p, y=wp, fill=Model)) + 
#  geom_line(aes(x=p,y=wp,col=Model,lty=Model))+
#  labs(title = "Average width plot",x="p-Probability interval",y="Probability interval width")+
#  theme(legend.position="bottom") 

#####################################################################################################
################# Now we plot accracy and average width plot for combined runs ######################
#####################################################################################################

kp.data2<-data.frame(kp=c(geokp,defkp),p=c(prob_seq,prob_seq),Model=c(rep("Stationary",times=length(prob_seq)),rep("Nonstationary",times=length(prob_seq))))
ggplot(kp.data2, aes(x=p, y=kp, fill=Model)) + 
  geom_line(aes(x=p,y=kp,col=Model,lty=Model))+geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.4)+
  labs(title = "Accuracy plot",x="p-Probability interval",y="Proportion within probability interval")+
  theme(legend.position="bottom")


wp.data2<-data.frame(wp=c(geowp,defwp),p=c(prob_seq,prob_seq),Model=c(rep("Stationary",times=length(prob_seq)),rep("Nonstationary",times=length(prob_seq))))
ggplot(wp.data2, aes(x=p, y=wp, fill=Model)) + 
  geom_line(aes(x=p,y=wp,col=Model,lty=Model))+
  labs(title = "Average width plot",x="p-Probability interval",y="Probability interval width")+
  theme(legend.position="bottom") 


############################################################################
######### defining stationary Matern covariance/variogram function #########
############################################################################
N=70

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

################## recomputing my_Cns for comparison purpose ############
x<-seq(0,2,length=N)
y<-seq(0,2,length=N)

grid<-expand.grid(x,y)
colnames(grid)=c("x","y")


##################################################
##### Seting nonstationary Matérn parameters #####
##################################################

nu<-0.8
sigma<-sqrt(5)

#### Anchor ranges ####
range.region1<-(0.1)^2
range.region2<-(0.9)^2

range.region3<-(0.1)^2
range.region4<-(0.9)^2  

#######################################################
####### Computing continuously varying ranges #########
#######################################################

weight_avg<-function(tloc=grid,anch.locs=cbind(c(0.5,1.5,0.5,1.5),c(1.5,1.5,0.5,0.5)),
                     anch.ranges=c(range.region1,range.region2,range.region3,range.region4),
                     bwidth=0.5)
{
  t.a.dist<-rdist(tloc,anch.locs)
  gauss.part<-exp(-(t.a.dist^2)/bwidth)
  denoms<-rowSums(gauss.part)
  
  num.part<-gauss.part*matrix(anch.ranges,nrow=nrow(gauss.part),ncol=ncol(gauss.part),byrow = T)
  num.part2<-rowSums(num.part)
  return(num.part2/denoms)
  
}

#### Setting bandwidth lambda ####
mylambda=0.15
cont.ranges<-weight_avg(bwidth = mylambda)
par(mfrow=c(1,1))
quilt.plot(grid,sqrt(cont.ranges),nx=N,ny=N,col = viridis(n=5000))
par(mar=c(4,4,0.5,0.5))

quilt.plot(grid,sqrt(cont.ranges),nx=N,ny=N,col = viridis(n=5000),horizontal=T) #### 453 by 487 saving resolution
quilt.plot(grid,sqrt(cont.ranges),nx=N,ny=N,col = viridis(n=5000),horizontal=T,zlim=range(sqrt(cont.ranges)),xlab="x",ylab="y",
           axis.args = list(at = seq(from = round(min(sqrt(cont.ranges)),2),to=round(max(sqrt(cont.ranges)),2),length.out = 5), labels=seq(from = round(min(sqrt(cont.ranges)),2),to=round(max(sqrt(cont.ranges)),2),length.out = 5))) #### 453 by 487 saving resolution




