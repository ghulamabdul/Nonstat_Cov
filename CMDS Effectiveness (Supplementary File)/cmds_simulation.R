#############################################
#### simulations on cmds effectiveness ######
#############################################
library(fields)
#### warping function ####
warp<-function(t,a)
{ if(a!=0){
  num1<-(exp(a*t/sqrt(8))-1)
  num2<-(exp(a)-1)
  return(sqrt(8)*(num1/num2))
}
  else
  {
    return(t)
  }
}

h<-seq(0,sqrt(8),length=100)
plot(h,warp(t=h,a=1))
points(h,warp(t=h,a=-1))
warp_param=1  #### Warping parameter (-ve leads to stretching, postive leads to squeezing)
  

  

########### Now creating a set of locations in [0,2]^[0,2] ##########
N<-30
x<-y<-seq(0,2,length=N)
  
##################################
#### Creating Simulation grid ####
##################################

coord.x<-rep(x,times=N)
coord.y<-rep(y,each=N)

data<-data.frame(x=coord.x,y=coord.y)
plot(data,pch=19,cex=0.5)



################################################################
############ Extracting points for subregion 1 #################
################################################################


reg1<-data[data$x<1,]
reg2<-data[data$x>1,]
par(mfrow=c(1,2))
plot(reg1,pch=19,cex=0.5)
plot(reg2,pch=19,cex=0.5)

#########################################################
##### Now we first compute the marginal distances #######
#########################################################
ordered_data<-data.frame(x=c(reg1$x,reg2$x),y=c(reg1$y,reg2$y))


############################################
###### Matrix decomposition ################
############################################

HR1<-HR2<-W1<-W2<-matrix(NA,nrow=length(reg1[,1]),ncol=length(reg2[,1]))
for(i in 1:nrow(HR1))
{
  for(j in 1:ncol(HR1))
  {
    x1<-reg1[i,]$x
    y1<-reg1[i,]$y
    x2<-reg2[j,]$x
    y2<-reg2[j,]$y
    ymid<-y1+((y2-y1)/(x2-x1))*(1-x1)
    xmid=1
    h1<-sqrt( ((x1-xmid)^2) +  ((y1-ymid)^2) )
    h2<-sqrt( ((x2-xmid)^2) +  ((y2-ymid)^2) )
    h<-h1+h2
    W1[i,j]<-h1/h
    W2[i,j]<-h2/h
    HR1[i,j]<-h
    HR2[i,j]<-h
    
  }
}



####### set warp_parameter #######
warp_param<-2.5


########## b1 is the distance matrix for locations in subregion 1 ############
b1<-rdist(cbind(reg1$x,reg1$y))
#####################################################################################################
############ Creating local distance warping function for subregion 1 with kernel smoothing #########
#####################################################################################################

dwarp1<-function(h,warp_param)
{
  return(warp(t=h,a=warp_param))
}

####################################################################################################
############# R11 is the warped distance matrix for subregion 1 ####################################
####################################################################################################

R11<-dwarp1(h=b1,warp_param = warp_param)
rm(b1) ####### Removing b1 matrix #######


########## b2 is the distance matrix for locations in subregion 2 ############

b2<-rdist(cbind(reg2$x,reg2$y))

#####################################################################################################
############ Creating local distance warping function for subregion 2 with kernel smoothing #########
#####################################################################################################


dwarp2<-function(h,warp_param)
{
  return(warp(t=h,a=-warp_param))
}

####################################################################################################
############# R22 is the warped distance matrix for subregion 2 ####################################
####################################################################################################

R22<-matrix(NA,nrow=length(reg2[,1]),ncol=length(reg2[,1]))
R22<-dwarp2(h=b2,warp_param = warp_param)


##########################################################
########### Now we compute all the cross distances #######
##########################################################

############################################################################################################
############# R12 is the warped cross distance matrix for subregion 1,2 ####################################
############################################################################################################


####### Matrix decomposition ######
R12<-W1*dwarp1(h=HR1,warp_param = warp_param)+W2*dwarp2(h=HR2,warp_param = warp_param)







################################################################################################
############# Constructing the global distance matrix wdist (Delta in the manuscript) ##########
################################################################################################

row1<-cbind(R11,R12)
row2<-cbind(t(R12),R22)
wdist<-rbind(row1,row2)
rm(R11,R22,R12,row1,row2)

################################################################################################################################
########### setting the diagonals of the distance matrix to be zero, as due to smoothing they might not exactly be zero ########
################################################################################################################################

diag(wdist)<-0

###################################################################################
####Constructing Deformed space coordinates in 2d and 3d using cmds ###############
###################################################################################

deform_coord<-cmdscale(wdist,k=2)
deform_coord_3d<-cmdscale(wdist,k=3)
plot(deform_coord)
deform_coord_30d<-cmdscale(wdist,k=30)

#Visualization part

###############################################
###### Plotting Geographic Space Gridded ###### 
###############################################
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,2),ylim=c(0,2),xlab="x",ylab="y",main="Original space")
for(i in 1:N)
{
  lines(x=c(x[i],x[i]),y=c(0,2))
  lines(x=c(0,2),y=c(y[i],y[i]))
}
points(coord.x,coord.y,pch=15,cex=0.8)
abline(v=1,col="red",lwd=2)
text(0.5,1,"Subregion1",col="Red",cex=3)
text(1.5,1,"Subregion2",col="Red",cex=3)

##################################################################
################ Setting data in the original ordering ###########
##################################################################
ord_geo<-data.frame(x=coord.x, y=coord.y)
unord_def<-data.frame(x=deform_coord_3d[,1],y=deform_coord_3d[,2],z=deform_coord_3d[,3])
unord_geo<-rbind(reg1,reg2)
co.order<-with(unord_geo,order(y,x))
ord_def<-unord_def[co.order,]


##### Gridded Fields ######

par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=1,mar=c(5,6,4,1)+.1)

#######################################################################################################
###### Now we need to sort the x and y axis of deformed space according to the original ordering ######
#######################################################################################################

dummy1<-data.frame(x=c(reg1$x,reg2$x),y=c(reg1$y,reg2$y))
coord_ordering<-with(dummy1,order(y,x))
deform_ordered<-deform_coord_3d[coord_ordering,]
#quilt.plot(deform_ordered[,c(1,2)],error,nx=70,ny=70,cex=5)
#plot(NA,NA,xlim=c(min(deform_ordered[,1]),max(deform_ordered[,1])),ylim=c(min(deform_ordered[,2]),max(deform_ordered[,2])))
#for(i in 1:N)
#{
#  lines(x=deform_ordered[((i-1)*N+1):(i*N),1],y=deform_ordered[((i-1)*N+1):(i*N),2])
  
#}
#xvals<-numeric()
#yvals<-numeric()
#for(j in 1:N)
#{
#  for(i in 1:N)
 # {
 #   xvals[i]<-deform_ordered[(i-1)*N+j,1]
 #   yvals[i]<-deform_ordered[(i-1)*N+j,2]
 #   
 # }
 # lines(x=xvals,y=yvals)
#}
#points(deform_ordered[,1],deform_ordered[,2],pch=19,cex=0.5,col=color.scale(error,col=tim.colors()))



#####################################################
####### Plotting Deformed gridded space in 2D #######
#####################################################

###### rotating 90 degrees ######

deform_ordered1<-deform_ordered
deform_ordered1[,1]<--deform_ordered[,2]
deform_ordered1[,2]<-deform_ordered[,1]


plot(NA,NA,xlim=c(min(deform_ordered1[,1]),max(deform_ordered1[,1])),ylim=c(min(deform_ordered1[,2]),max(deform_ordered1[,2])),xlab="x'",ylab="y'",main="Deformed space in 2d")
for(i in 1:N)
{
  lines(x=deform_ordered1[((i-1)*N+1):(i*N),1],y=deform_ordered1[((i-1)*N+1):(i*N),2])
  
}
xvals<-numeric()
yvals<-numeric()
for(j in 1:N)
{
  for(i in 1:N)
  {
    xvals[i]<-deform_ordered1[(i-1)*N+j,1]
    yvals[i]<-deform_ordered1[(i-1)*N+j,2]
    
  }
  lines(x=xvals,y=yvals)
}
points(deform_ordered1[,1],deform_ordered1[,2],pch=15,cex=0.8)

#####################################################
####### Plotting Deformed gridded space in 3D #######
#####################################################
library(scatterplot3d)
#### Rotating for visualization purpose (around x axis at 45 degree) #####
par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=1,mar=c(8,6,4,1)+.1)
par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=1,mar=c(8,5,1,1)+.1)
pp=pi/4
deform_ordered3<-deform_ordered
deform_ordered3[,2]<-deform_ordered[,2]*cos(pp)-deform_ordered[,3]*sin(pp)
deform_ordered3[,3]<-deform_ordered[,2]*sin(pp)+deform_ordered[,3]*cos(pp)
#split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
s3<-scatterplot3d(x=deform_ordered3[,1],y=deform_ordered3[,2],z=deform_ordered3[,3],pch=19,cex.symbols = 0.8,color = "white",grid=T,box=FALSE,xlab="x'",ylab="y'",zlab="z'",angle = 10,main="Deformed space in 3d")
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
s3$points3d(x=deform_ordered3[,1],y=deform_ordered3[,2],z=deform_ordered3[,3],pch=15,cex = 0.8)




############################
#### computing NMSE ########
############################
dimensions_k<-2:30
nmse_k<-numeric(length=29)
for(i in 2:30)
{
  est_dist<-rdist(deform_coord_30d[,1:i])
  true_dist<-c(wdist)
  cmds_dist<-c(est_dist)
  nmse_k[i-1]<-1-(sum((cmds_dist-true_dist)^2))/(sum((true_dist-mean(true_dist))^2))

}
par(mfrow=c(1,1))
par(mar=c(4,4,4,4))
h<-seq(0,sqrt(8),length=100)
plot(h,warp(t=h,a=warp_param),type = "l",col=1,main = "Warping functions",ylab = " ",xlab=" ")
lines(h,warp(t=h,a=-warp_param),col=2,lty=2)
legend("bottomright",lty=c(1,2),col=c(1,2),c("Subregion 1 ", "Subregion 2"),cex=1.3)
par(mar=c(4,6,4,4))
plot(dimensions_k,nmse_k,xlim = c(2,30),xlab="dimensions (starting from 2 to going till 30)",ylab = "Normalised mean squared error",cex=1.5,pch=19,ylim=c(0,1))
mxn<-max(nmse_k)
dim_max<-which.max(nmse_k)+1
legend("bottomright",c(paste("max(NMSE)=",round(mxn,4)),paste("at dimension=",dim_max)))



################################################################################################################
############# From here we will do everything again but for the warping functions created by beta cdfs #########
################################################################################################################




betadwarp1<-function(h,a1,b1)
{
  return(sqrt(8)*pbeta(h/sqrt(8),a1,b1))
}

betadwarp2<-function(h,a2,b2)
{
  return(sqrt(8)*pbeta(h/sqrt(8),a2,1/b2))
}

####### set beta cdfs parameters a1,a2,b1 ####
t.h<-seq(0,sqrt(8),length=200)
#a1=1/4  #### case 1: a1=0.7,1,1/4,1/4
#a2=2 #### case 1:a2 =1,1,8,2
#b1=1.4  ##### case 1: b1=1.5,1.4,1.4,1.4
#b2=1    ##### case 1:b2=2.5,2.5,0.7/1.4,1
a1=0.25
b1=1.4
a2=2
b2=1

par(mfrow=c(1,1))
plot(t.h,betadwarp1(h=t.h,a1=a1,b1=b1),type="l",col=1)
lines(t.h,betadwarp2(h=t.h,a2=a2,b2=b2),col=2)
lines(t.h,t.h,col=3)
########## b1 is the distance matrix for locations in subregion 1 ############
b1dist<-rdist(cbind(reg1$x,reg1$y))

####################################################################################################
############# R11 is the warped distance matrix for subregion 1 ####################################
####################################################################################################

R11<-betadwarp1(h=b1dist,a1=a1,b1=b1)
####### Removing b1 matrix #######


########## b2 is the distance matrix for locations in subregion 2 ############

b2dist<-rdist(cbind(reg2$x,reg2$y))


####################################################################################################
############# R22 is the warped distance matrix for subregion 2 ####################################
####################################################################################################

R22<-matrix(NA,nrow=length(reg2[,1]),ncol=length(reg2[,1]))
R22<-betadwarp2(h=b2dist,a2=a2,b2=b2)


##########################################################
########### Now we compute all the cross distances #######
##########################################################

############################################################################################################
############# R12 is the warped cross distance matrix for subregion 1,2 ####################################
############################################################################################################


####### Matrix decomposition ######
R12<-W1*betadwarp1(h=HR1,a1=a1,b1=b1)+W2*betadwarp2(h=HR2,a2=a2,b2=b2)







################################################################################################
############# Constructing the global distance matrix wdist (Delta in the manuscript) ##########
################################################################################################

row1<-cbind(R11,R12)
row2<-cbind(t(R12),R22)
wdist<-rbind(row1,row2)
rm(R11,R22,R12,row1,row2)

################################################################################################################################
########### setting the diagonals of the distance matrix to be zero, as due to smoothing they might not exactly be zero ########
################################################################################################################################

diag(wdist)<-0

###################################################################################
####Constructing Deformed space coordinates in 2d and 3d using cmds ###############
###################################################################################

deform_coord<-cmdscale(wdist,k=2)
deform_coord_3d<-cmdscale(wdist,k=3)
plot(deform_coord)
deform_coord_30d<-cmdscale(wdist,k=30)

#Visualization part

###############################################
###### Plotting Geographic Space Gridded ###### 
###############################################
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,2),ylim=c(0,2),xlab="x",ylab="y",main="Original space")
for(i in 1:N)
{
  lines(x=c(x[i],x[i]),y=c(0,2))
  lines(x=c(0,2),y=c(y[i],y[i]))
}
points(coord.x,coord.y,pch=15,cex=0.8)
abline(v=1,col="red",lwd=2)
text(0.5,1,"Subregion1",col="Red",cex=3)
text(1.5,1,"Subregion2",col="Red",cex=3)

##################################################################
################ Setting data in the original ordering ###########
##################################################################
ord_geo<-data.frame(x=coord.x, y=coord.y)
unord_def<-data.frame(x=deform_coord_3d[,1],y=deform_coord_3d[,2],z=deform_coord_3d[,3])
unord_geo<-rbind(reg1,reg2)
co.order<-with(unord_geo,order(y,x))
ord_def<-unord_def[co.order,]


##### Gridded Fields ######

par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=1,mar=c(5,6,4,1)+.1)

#######################################################################################################
###### Now we need to sort the x and y axis of deformed space according to the original ordering ######
#######################################################################################################

dummy1<-data.frame(x=c(reg1$x,reg2$x),y=c(reg1$y,reg2$y))
coord_ordering<-with(dummy1,order(y,x))
deform_ordered<-deform_coord_3d[coord_ordering,]
#quilt.plot(deform_ordered[,c(1,2)],error,nx=70,ny=70,cex=5)
#plot(NA,NA,xlim=c(min(deform_ordered[,1]),max(deform_ordered[,1])),ylim=c(min(deform_ordered[,2]),max(deform_ordered[,2])))
#for(i in 1:N)
#{
#  lines(x=deform_ordered[((i-1)*N+1):(i*N),1],y=deform_ordered[((i-1)*N+1):(i*N),2])

#}
#xvals<-numeric()
#yvals<-numeric()
#for(j in 1:N)
#{
#  for(i in 1:N)
# {
#   xvals[i]<-deform_ordered[(i-1)*N+j,1]
#   yvals[i]<-deform_ordered[(i-1)*N+j,2]
#   
# }
# lines(x=xvals,y=yvals)
#}
#points(deform_ordered[,1],deform_ordered[,2],pch=19,cex=0.5,col=color.scale(error,col=tim.colors()))



#####################################################
####### Plotting Deformed gridded space in 2D #######
#####################################################

###### rotating 90 degrees ######

deform_ordered1<-deform_ordered
deform_ordered1[,1]<--deform_ordered[,2]
deform_ordered1[,2]<-deform_ordered[,1]


plot(NA,NA,xlim=c(min(deform_ordered1[,1]),max(deform_ordered1[,1])),ylim=c(min(deform_ordered1[,2]),max(deform_ordered1[,2])),xlab="x'",ylab="y'",main="Deformed space in 2d")
for(i in 1:N)
{
  lines(x=deform_ordered1[((i-1)*N+1):(i*N),1],y=deform_ordered1[((i-1)*N+1):(i*N),2])
  
}
xvals<-numeric()
yvals<-numeric()
for(j in 1:N)
{
  for(i in 1:N)
  {
    xvals[i]<-deform_ordered1[(i-1)*N+j,1]
    yvals[i]<-deform_ordered1[(i-1)*N+j,2]
    
  }
  lines(x=xvals,y=yvals)
}
points(deform_ordered1[,1],deform_ordered1[,2],pch=15,cex=0.8)

#####################################################
####### Plotting Deformed gridded space in 3D #######
#####################################################
library(scatterplot3d)
#### Rotating for visualization purpose (around x axis at 45 degree) #####
par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=1,mar=c(8,6,4,1)+.1)
par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=1,mar=c(8,5,1,1)+.1)
pp=pi/4
deform_ordered3<-deform_ordered
deform_ordered3[,2]<-deform_ordered[,2]*cos(pp)-deform_ordered[,3]*sin(pp)
deform_ordered3[,3]<-deform_ordered[,2]*sin(pp)+deform_ordered[,3]*cos(pp)
#split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
s3<-scatterplot3d(x=deform_ordered3[,1],y=deform_ordered3[,2],z=deform_ordered3[,3],pch=19,cex.symbols = 0.8,color = "white",grid=T,box=FALSE,xlab="x'",ylab="y'",zlab="z'",angle = 10,main="Deformed space in 3d")
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
s3$points3d(x=deform_ordered3[,1],y=deform_ordered3[,2],z=deform_ordered3[,3],pch=15,cex = 0.8)




############################
#### computing NMSE ########
############################
dimensions_k<-2:30
nmse_k<-numeric(length=29)
for(i in 2:30)
{
  est_dist<-rdist(deform_coord_30d[,1:i])
  true_dist<-c(wdist)
  cmds_dist<-c(est_dist)
  nmse_k[i-1]<-1-(sum((cmds_dist-true_dist)^2))/(sum((true_dist-mean(true_dist))^2))
  
}
par(mfrow=c(1,1))
par(mar=c(4,4,4,4))
h<-seq(0,sqrt(8),length=100)
plot(h,betadwarp1(h=h,a1=a1,b1=b1),type = "l",col=1,main = "Warping functions",ylab = " ",xlab=" ")
#lines(h,h,col=3)
lines(h,betadwarp2(h=h,a2=a2,b2=b2),col=2,lty=2)
legend("bottomright",lty=c(1,2),col=c(1,2),c("Subregion 1 ", "Subregion 2"),cex=1.3)
par(mar=c(4,6,4,4))
plot(dimensions_k,nmse_k,xlim = c(2,30),xlab="dimensions (starting from 2 to going till 30)",ylab = "Normalised mean squared error",cex=1.5,pch=19,ylim=c(0,1))
mxn<-max(nmse_k)
dim_max<-which.max(nmse_k)+1
legend("bottomright",c(paste("max(NMSE)=",round(mxn,4)),paste("at dimension=",dim_max)))



