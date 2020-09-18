##########################################################
###### Setting working directory and attaching data ######
##########################################################

setwd("/ibex/scratch/qadirga/Project_1/Data_Application_revised")
library( fields)
attach("RData.COmonthly.met")  ####### Attaching dataset ######
CO.id                          ####### Station ids
CO.loc<-CO.loc                 ####### Station locations
CO.ppt                         ####### precipitation data
library(fdasrvf)               ####### Loading libraries
library(doParallel)            ####### Loading libraries
ncores<-detectCores()
registerDoParallel(cores = ncores-4) ####### Setting parallel environment
getDoParWorkers()
library(rgl)

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

quilt.plot(my.data) ###### Plotting data #####
max(my.data$y)


#### We partitioning from the line x=-104.873 ###
mdline<--104.873
abline(v=mdline)


my.data<-na.omit(my.data) ###### Removing the NA values, or the unobserved station locations 
my.data$ppt<-(my.data$ppt-mean(my.data$ppt))/sd(my.data$ppt) ######### Standardizing the data
my.data<-my.data[order(my.data$y,my.data$x),]

quilt.plot(my.data,xlab="Longitude",ylab="Latitude")  ###### Plotting the standardized data
abline(v=mdline,lwd=3)

##### arranging data as reg1 data (x<partition line) and reg2 data (x>=partition line) ####
reg1<-my.data[my.data$x<mdline,]
reg2<-my.data[my.data$x>=mdline,]
full.data<-rbind(reg1,reg2)

reord.index<-order(full.data$y,full.data$x)
#full.data[reord.index,]==my.data

##### Checking normal assumptions ########
par(mfrow=c(1,2))
colnames(my.data)<-colnames(full.data)<-c("x","y","z")
hist(full.data$z)
qqnorm(full.data$z)
qqline(full.data$z)

#########################################################
######## Performing random fold cross-validation ########
#########################################################
library(fdasrvf)
library(geoR)
library(fields)
library(plot3D)
library(rgl)
library(scatterplot3d)
library(scoringRules)
library(mvtnorm)

cross_validation_summary<-foreach(myseed=1:100,.packages = c("fdasrvf","geoR","fields","plot3D","rgl","scatterplot3d","scoringRules","mvtnorm"))%dopar%{
  set.seed(myseed) ### Setting seed
  rmiss<-sample(1:length(my.data$x),size=30) #### Generating test location indices
  
  test.data<-my.data[rmiss,]
  train.data<-my.data[-rmiss,]
  
  ###########################################################################
  ############ Dividing dataset regionwise from the training data ###########
  ###########################################################################
  
  reg1.tr<-train.data[train.data$x<mdline,]
  reg2.tr<-train.data[train.data$x>=mdline,]
  
  
  ###########################################################################
  ########## Fitting WLS Matérn variogram for the MLE initial values ########
  ###########################################################################
  
  #####################################################################
  ################ Computing regionwise empirical variograms ##########
  #####################################################################
  
  vr1<-variog(coords = cbind(reg1.tr$x,reg1.tr$y),data=reg1.tr$z,option = "bin",uvec=30,max.dist=(10)/2)
  vr2<-variog(coords = cbind(reg2.tr$x,reg2.tr$y),data=reg2.tr$z,uvec=30,max.dist = (10)/2,option = "bin")
  
  
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
  
  #########################################################################
  ########## Now using MLE to estimate Matérn covariance parameters #######
  #########################################################################
  
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
  
  myvario<-vr2
  temp.r2<-optim(my.var.loss,par = c(runif(n=1,min=0.1,max=30),0.5,sd(reg2.tr$z)))
  
  mle.r2<-optim_marg_comp.loglik.r2(par = temp.r2$par)
  
  
  
  
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
  
  #########################################################
  ######### Now we compute the warped distances ###########
  #########################################################
  
  
  ##############################################################################
  ########## b1 is the distance matrix for locations in subregion 1 ############
  ##############################################################################
  
  b1<-rdist(cbind(reg1$x,reg1$y))
  
  #####################################################################################################
  ############ Creating local distance warping function for subregion 1 with kernel smoothing #########
  ############# Gaussian kernel with bandwidth 0.02 ###################################################
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
  
  R11<-dwarp1(h=b1)
  diag(R11)<-0 ### forcing diagonals to be zero, as the smoothing might have made it neglibly different from zero 
  
  
  
  ########## b2 is the distance matrix for locations in subregion 2 ############
  
  b2<-rdist(cbind(reg2$x,reg2$y))
  
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
  diag(R22)<-0 ### forcing diagonals to be zero, as the smoothing might have made it neglibly different from zero 
  
  
  ####################################################################################################
  ############# R22 is the warped distance matrix for subregion 2 ####################################
  ####################################################################################################
  
  #####################################################################################################################################
  ########## Now computing the matrix R12, which is the warped distances for the two points in region 1 and region 2 respectively #####
  #####################################################################################################################################
  
  midpoint.f<-function(i,j,xmid=mdline)
  {
    x1<-reg1[i,]$x
    y1<-reg1[i,]$y
    x2<-reg2[j,]$x
    y2<-reg2[j,]$y
    ymid<-y1+((y2-y1)/(x2-x1))*(xmid-x1)
    return(ymid)
  }
  
  y.inter<-outer(1:length(reg1$x),1:length(reg2$x),midpoint.f)
  
  
  
  
  dwarp12<-function(region1=reg1,region2=reg2,y.intersect=y.inter,xmid=mdline)
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
  
  ################################################################################################################################
  ########### setting the diagonals of the distance matrix to be zero, as due to smoothing they might not exactly be zero ########
  ################################################################################################################################
  
  diag(wdist)<-0 ### forcing diagonals to be zero, as the smoothing might have made it neglibly different from zero 
  
  
  ##################################################################
  ####Constructing Deformed space coordinates in 30d ###############
  ##################################################################
  
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
  #plot(1:30,nmse_wdist)
  
  
  
  psi<-which.max(-abs(nmse_wdist-1))
  #nmse_wdist[psi]
  
  
  
  ###################################################
  ##### first re_ordering deformed coordinates ######
  ###################################################
  deform_ordered<-deform_coord_30d[reord.index,]
  
  
  ####################################################################################
  ########### Now fitting a stationary Matérn in the deformed space ##################
  ########### on training data #######################################################
  ####################################################################################
  mle.comp_only_mlv.def<-function(par)
  {
    return(mle.comp.all(locs=deform_ordered[-rmiss,1:psi],z=my.data$z[-rmiss],p=par)$mlv
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
    return(mle.comp.all(locs=my.data[-rmiss,-3],z=my.data$z[-rmiss],p=par)$mlv
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
  
  
  ##### setting the initial values by averaging the region wise estimates ######
  f.ini<-rbind(temp.r1$par,temp.r2$par)
  f.ini<-colMeans(f.ini)
  mle.def<-optim_marg_comp.loglik.def(par = f.ini)
  mle.geo<-optim_marg_comp.loglik.geo(par = f.ini)
  
  
  ####### Plotting nonstationary covariance function #######
  
  def_cov<-my.matern(h=rdist(deform_ordered[,1:psi]),a=mle.def$par[1],nu=mle.def$par[2],sigma = mle.def$par[3])
  geo_cov<-my.matern(h=rdist(cbind(my.data$x,my.data$y)),a=mle.geo$par[1],nu=mle.geo$par[2],sigma = mle.geo$par[3])
  
  
  ####### Now we do prediction in the geographic space ##########
  
  geo.sig11<-geo_cov[rmiss,rmiss]
  geo.sig22<-geo_cov[-rmiss,-rmiss]
  geo.sig12<-geo_cov[rmiss,-rmiss]
  
  geo_wts<-geo.sig12%*%solve(geo.sig22)
  geo_pred<-geo_wts%*%my.data$z[-rmiss] ###### predicted value
  
  geo_pvar<-diag(geo.sig11-geo_wts%*%t(geo.sig12)) ##### prediction variance
  
  ####### scores_func computes all the scores for the given predictions #####
  scores_func<-function(true=my.data$z[rmiss],pred,predvar)
  {
    rmse<-sqrt(mean((true-pred)^2))
    mae<-mean(abs(true-pred))
    nmse<-1-(sum((pred-true)^2))/(sum((true-mean(true))^2))
    mycrps<-mean(crps(true, "norm", mean = c(pred), sd = c(sqrt(predvar))))
    mylogs<-mean(logs(true, "norm", mean = c(pred), sd = c(sqrt(predvar))))
    return(list(RMSE=rmse,NMSE=nmse,MAE=mae,mCRPS=mycrps,mLogS=mylogs))
  }
  geo_scores<-scores_func(pred = geo_pred,predvar = geo_pvar)
  
  
  ##############################################################################
  ########### Now we do prediction in the deformed space #######################
  ##############################################################################
  
  
  def.sig11<-def_cov[rmiss,rmiss]
  def.sig22<-def_cov[-rmiss,-rmiss]
  def.sig12<-def_cov[rmiss,-rmiss]
  
  
  def_wts<-def.sig12%*%solve(def.sig22)
  def_pred<-def_wts%*%my.data$z[-rmiss] ### predicted value
  
  def_pvar<-diag(def.sig11-def_wts%*%t(def.sig12)) #### prediction variance
  
  def_scores<-scores_func(pred = def_pred,predvar = def_pvar)
  
  
  #####################################################################
  ######## Exploring the accuracy of prediction uncertainty ###########
  #####################################################################
  
  ######## Computing G-statistics and Wp, K ########
  prob_seq<-seq(0.01,0.99,length.out = 4*99)
  uncert_p<-function(pr,truev=my.data$z[rmiss],predval,predvar)
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
       data=my.data,
       psi=psi,nmse_psi=nmse_wdist[psi],
       def.predval=def_pred,
       def.predvariance=def_pvar,
       geo_predval=geo_pred,
       geo_predvariance=geo_pvar,
       valid_set=my.data[rmiss,],
       rmiss=rmiss
  )
  
}

save.image("Data_crossvals9.RData")

sim_summary<-cross_validation_summary
prob_seq<-seq(0.01,0.99,length.out = 4*99)
geo_G<-def_G<-geo_nmse<-def_nmse<-geo_rmse<-def_rmse<-geo_mae<-def_mae<-geo_crps<-def_crps<-geo_logs<-def_logs<-numeric(length = 100)
geo_kp<-def_kp<-geo_wp<-def_wp<-matrix(NA,nrow=100,ncol=length(prob_seq))



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
  
  
}


###### Plots for prediction scores ###
par(mfrow=c(2,3))
rounding=4

boxplot(geo_rmse,def_rmse,names = c("Stat","Nstat"),col = "gray",ylim=c(range(geo_rmse,def_rmse)+c(0,+0.01)),main="RMSE")
points(c(1,2),c(mean(geo_rmse),mean(def_rmse)),pch=19,col="red")
text(x=c(1,2),y=c(max(geo_rmse)+0.005,max(def_rmse)+0.005),round(c(mean(geo_rmse),mean(def_rmse)),rounding),col="blue")

boxplot(geo_nmse,def_nmse,names = c("Stat","Nstat"),col = "gray",ylim=c(range(geo_nmse,def_nmse)+c(0,+0.1)),main="NMSE")
points(c(1,2),c(mean(geo_nmse),mean(def_nmse)),pch=19,col="red")
text(x=c(1,2),y=c(max(geo_nmse)+0.05,max(def_nmse)+0.05),round(c(mean(geo_nmse),mean(def_nmse)),rounding),col="blue")

boxplot(geo_mae,def_mae,names = c("Stat","Nstat"),col = "gray",ylim=c(range(geo_mae,def_mae)+c(0,+0.01)),main="MAE")
points(c(1,2),c(mean(geo_mae),mean(def_mae)),pch=19,col="red")
text(x=c(1,2),y=c(max(geo_mae)+0.005,max(def_mae)+0.005),round(c(mean(geo_mae),mean(def_mae)),rounding),col="blue")


boxplot(geo_crps,def_crps,names = c("Stat","Nstat"),col = "gray",ylim=c(range(geo_crps,def_crps)+c(0,+0.01)),main="mCRPS")
points(c(1,2),c(mean(geo_crps),mean(def_crps)),pch=19,col="red")
text(x=c(1,2),y=c(max(geo_crps)+0.005,max(def_crps)+0.005),round(c(mean(geo_crps),mean(def_crps)),rounding),col="blue")

boxplot(geo_logs,def_logs,names = c("Stat","Nstat"),col = "gray",ylim=c(range(geo_logs,def_logs)+c(0,+0.05)),main="mLogS")
points(c(1,2),c(mean(geo_logs),mean(def_logs)),pch=19,col="red")
text(x=c(1,2),y=c(max(geo_logs)+0.017,max(def_logs)+0.017),round(c(mean(geo_logs),mean(def_logs)),rounding),col="blue")

boxplot(geo_G,def_G,names = c("Stat","Nstat"),col = "gray",ylim=c(range(geo_G,def_G)+c(0,+0.01)),main="G-stat")
points(c(1,2),c(mean(geo_G),mean(def_G)),pch=19,col="red")
text(x=c(1,2),y=c(max(geo_G)+0.005,max(def_G)+0.005),round(c(mean(geo_G),mean(def_G)),rounding),col="blue")


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
  test.set<-c(c(test.set),sim_summary[[i]]$valid_set$z)
  def.krig<-c(c(def.krig),sim_summary[[i]]$def.predval)
  def.krig.var<-c(c(def.krig.var),sim_summary[[i]]$def.predvariance)
  geo.krig<-c(c(geo.krig),sim_summary[[i]]$geo_predval)
  geo.krig.var<-c(c(geo.krig.var),sim_summary[[i]]$geo_predvariance)
}
prob_seq<-seq(0.01,0.99,length.out = 4*99)

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


#geo_unc<-uncert_p(pr=prob_seq[1],predval =geo_pred,predvar = geo_pvar )
###
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
  geom_label(aes(y = max(NMSE)+0.034, x = c(1), label = paste("sd=",nmse.sd.def)))+
  geom_label(aes(y = max(NMSE)+0.034, x = c(2), label = paste("sd=",nmse.sd.geo)))

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
  geom_label(aes(y = max(mLogS)+0.03, x = c(1), label = paste("sd=",mlogs.sd.def)))+
  geom_label(aes(y = max(mLogS)+0.03, x = c(2), label = paste("sd=",mlogs.sd.geo)))


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
  geom_label(aes(y = max(Gstat)+0.07, x = c(1), label = paste("sd=",G.sd.def)))+
  geom_label(aes(y = max(Gstat)+0.07, x = c(2), label = paste("sd=",G.sd.geo)))

####################################################################################
############### Now we do the accuracy plots and the average width plots ###########
####################################################################################

kp.data<-data.frame(kp=c(geokp,defkp),p=c(prob_seq,prob_seq),Model=c(rep("Stationary",times=length(prob_seq)),rep("Nonstationary",times=length(prob_seq))))
ggplot(kp.data, aes(x=p, y=kp, fill=Model)) + 
  geom_line(aes(x=p,y=kp,col=Model,lty=Model))+geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.4)+
  labs(title = "Accuracy plot",x="p-Probability interval",y="Proportion within probability interval")+
  theme(legend.position="bottom")


wp.data<-data.frame(wp=c(geowp,defwp),p=c(prob_seq,prob_seq),Model=c(rep("Stationary",times=length(prob_seq)),rep("Nonstationary",times=length(prob_seq))))
ggplot(wp.data, aes(x=p, y=wp, fill=Model)) + 
  geom_line(aes(x=p,y=wp,col=Model,lty=Model))+
  labs(title = "Average width plot",x="p-Probability interval",y="Probability interval width")+
  theme(legend.position="bottom") 






