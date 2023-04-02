#################################
library(BayesLogit)#  for Polya-Gamma
library(MCMCpack)  #  for riwish
library(msm)       #  for rtnorm
library(mvnfast)   #  for rmvn
library(Matrix)    #  For matrix manipulation
library(kableExtra)# For viewing dataframes
######################
####################### 
#Data simulation begins
######################

######################
######## Spatial #####
######################
n<-500
##############################
#### Sample spatial effect 
#############################
# Download adjacency matrix 
A <- read.csv("https://raw.githubusercontent.com/eosafu/Model/master/Nig_adj_Matrix.txt",sep=",",header = F) # Nigeria adjacency matrix
colnames(A) <- NULL
A <- as(as.matrix(A[,-1]),"sparseMatrix")
image(A)
# Compute spatial effect precision matrix
P <- as(diag(rowSums(A)),"sparseMatrix")
q <- nrow(P)
rho <- 0.99
D <- P-rho*A
Q <- D
sig1 <- 2.0
sig2 <- 1.0
rhocorr <- 0.8
sig21<-sig12<-sqrt(sig1*sig2)*rhocorr

# Sample spatial effect
phi1 <- rnorm(q)
phi2 <- rnorm(q,0,1)
mu1_2   <-(sig12/sig2)*phi2
Sigma1_2<-((sig1*sig2-sig12*sig21)/sig2)*solve(Q)
mu2_1   <-(sig21/sig1)*phi1
Sigma2_1<-((sig1*sig2-sig12*sig21)/sig1)*solve(Q)
phi1 <- t(rmvn (1,mu1_2,Sigma1_2))
phi2 <- t(rmvn (1,mu2_1,Sigma2_1))
phi1 <- phi1 - mean(phi1)   # Constrain on theta1
phi2 <- phi2 - mean(phi2)   # Constrain on theta2

# Sample spatial id
V1id <- sample(1:q,n,replace = TRUE)
V1 <- matrix(0,nrow = n,ncol = q)
for(i in 1:n) V1[i,V1id[i]]<-1
V2 <-V1
####################
# Fixed effect
####################
beta1 <- matrix(c(2.30,-0.5,3.0))
beta2 <- matrix(c(1.10,-1.5,0.01))

# Sample design matrix
p <- nrow(beta1)
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
X1 <- X2 <- cbind(x1,x2,x3)


######################
#### Random walk model 
######################
# Precision matrix
r<-10
sigmatheta1 <- 1
sigmatheta2 <- 2
BB <- matrix(0,nrow=r,ncol=r)
for(i in 1:r){
  BB[i,]<-c(rep(1,i),rep(0,r-i))
}
Qtheta2 <-  BB %*% t(BB)
Qtheta2 <-Qtheta1 <-   chol2inv(chol(Qtheta2))

# Desing matrix
Z1id <- sample(1:r,n,replace = TRUE)
Z1 <- matrix(0,nrow = n,ncol = r)
for(i in 1:n) Z1[i,Z1id[i]]<-1
Z2 <-Z1

theta1 <- t(rmvn (1,matrix(0,r),sigmatheta1*solve(Qtheta1)))
theta2 <- t(rmvn (1,matrix(0,r),sigmatheta2*solve(Qtheta2)))
###########################
# Sample response
###########################
# linear predictor
eta1 <- X1%*%beta1+Z1%*%theta1+V1%*%phi1
eta2 <- X2%*%beta2+Z2%*%theta2 +V2%*%phi2

# probability of unstructured zero
pi <- exp(eta1)/(1+exp(eta1))

# Negative binomial probability
pneg <- exp(eta2)/(1+exp(eta2))

# draw response
xi<-1
auxid<-y<-NULL
for (i in 1:n) {
  auxid[i] <-  sample(c(0,1),1,prob =c(1-pi[i],pi[i]))
  if(auxid[i]==0){
    y[i]<-0
  }else{
    y[i]<-rnbinom(1,size=xi,prob = (1-pneg[i]))
  }
}
# Check proportion of zero
mean(y==0)
####################### 
#Data simulation Ends
######################

#####################################################
# Estimation using mcmc with polya-Gamma augmentation
#####################################################

# starting values
sigmatheta1 <- 1 # fixed at true value
sigmatheta2 <- 2 # fixed at true value
lambda<-6
ext_xi    <- 0.5
s<-0.005
H <- I<- diag(c(1,1))
Q_beta1<-Q_beta2<-0.01*Diagonal(p)
ext_phi1   <- matrix(rnorm(q))
ext_phi2   <- matrix(rnorm(q))
ext_theta1 <- matrix(rnorm(r))
ext_theta2 <- matrix(rnorm(r))
ext_beta1  <- matrix(rnorm(p))
ext_beta2  <- matrix(rnorm(p))

Burn<-2000
thinning<-2
Num<-Burn+5000
Accept <- 0

# Storage
BETA1 <-matrix(0,p,(Num-Burn)/thinning)
BETA2 <-matrix(0,p,(Num-Burn)/thinning)
THETA1<-matrix(0,r,(Num-Burn)/thinning)
THETA2<-matrix(0,r,(Num-Burn)/thinning)
PHI1  <-matrix(0,q,(Num-Burn)/thinning)
PHI2  <-matrix(0,q,(Num-Burn)/thinning)
XI    <-vector("numeric",(Num-Burn)/thinning)
HH    <-array(0,dim=c(2,2,(Num-Burn)/thinning))

## Estimation Begins

for (iter in 1:Num){
  # Using M-H
  # Update all the sigma using inv-Wishart prior and proposal
  auxphi<-  rbind(ext_phi1,ext_phi2)
  Hnew <-  riwish(lambda,H)
  ratio <- dmvn(t(auxphi),matrix(0,2*q),kronecker(Hnew,solve(Q)),log=T)+log(diwish(Hnew,lambda,I))+log(diwish(H,lambda,Hnew))-
  dmvn(t(auxphi),matrix(0,2*q),kronecker(H,solve(Q)),log=T)-log(diwish(H,lambda,I))-log(diwish(Hnew,lambda,H))
  
  U <- runif(1)
  if(U<exp(ratio)){
    H <-  Hnew
    Accept <- Accept+1
  }else{
    H <- H
  }
  
  sig1<-H[1,1]
  sig2<-H[2,2]
  sig12<-sig21<-H[1,2]
  
  # Evaluate pi and pneg
  eta2 <- X2%*%ext_beta2+Z2%*%ext_theta2+V2%*%ext_phi2
  eta1 <- X1%*%ext_beta1+Z1%*%ext_theta1+V1%*%ext_phi1
  pneg <- pmax(0.01,pmin(0.99, exp(eta2)/(1+exp(eta2))))
  pi <-   pmax(0.01,pmin(0.99,exp(eta1)/(1+exp(eta1))))
  prob <- ((1-pneg)^ext_xi*pi)/(1+((1-pneg)^ext_xi-1)*pi)
  
  # Compute latent variable
  w <-NULL
  w <- rbinom(n,1,prob = prob)
  w[y>0] <- 1
  
  nu1 <- rpg(n,1,eta1)
  z1 <- (w-0.5)/nu1
  Omega1 <- diag(nu1)
  
  # Update fixed, nonlinear and spatial effects for binary component
  SigmaBeta1 <- solve(Q_beta1+t(X1)%*%Omega1%*%X1)
  Mubeta1   <- SigmaBeta1%*%(t(X1)%*%Omega1%*%(z1-V1%*%ext_phi1-Z1%*%ext_theta1))
  ext_beta1 <- t(rmvn(1,Mubeta1,SigmaBeta1))
  #
  Q_theta1 <- (1/sigmatheta1)*Qtheta1  # assign prior to sigmatheta1
  Sigmatheta1 <- solve(Q_theta1+t(Z1)%*%Omega1%*%Z1)
  Mutheta1   <- Sigmatheta1%*%(t(Z1)%*%Omega1%*%(z1-V1%*%ext_phi1-X1%*%ext_beta1))
  ext_theta1<-t(rmvn(1,Mutheta1,Sigmatheta1))
  ext_theta1 <- ext_theta1-mean(ext_theta1)
  #
  mu1_2   <-(sig12/sig2)*ext_phi2 # mu1|2
  prec1_2<-(sig2/(sig1*sig2-sig12*sig21))*Q
  Sigmaphi1 <- solve(prec1_2+t(V1)%*%Omega1%*%V1)
  muphi12    <- Sigmaphi1%*%(prec1_2%*%mu1_2+t(V1)%*%Omega1%*%(z1-X1%*%ext_beta1-Z1%*%ext_theta1))
  ext_phi1 <- t(rmvn(1,muphi12,Sigmaphi1))
  ext_phi1 <- ext_phi1-mean(ext_phi1)
  
  ### Update xi (over dispersion parameter) using M-H
  nstar <- sum(w)
  
  # Update eta2
  
  eta2 <- X2[w==1,]%*%ext_beta2+Z2[w==1,]%*%ext_theta2+V2[w==1,]%*%ext_phi2
  
  nu2 <- rpg(nstar,ext_xi+y[w==1],eta2)
  Omega2 <- diag(nu2)
  z2 <- (y[w==1]-ext_xi)/(2*nu2)
  pneg <- pmax(0.01,pmin(0.99, exp(eta2)/(1+exp(eta2))))
  
  xinew<-rtnorm(1,ext_xi,sqrt(s),lower=0)
  ratio<-sum(dnbinom(y[w==1],xinew,(1-pneg),log=T))+dtnorm(ext_xi,xinew,sqrt(s),0,log=T)-
    sum(dnbinom(y[w==1],ext_xi,(1-pneg),log=T))-dtnorm(xinew,ext_xi,sqrt(s),0,log=T) # Uniform prior
  U <- runif(1)
  
  if(U<exp(ratio)){
    ext_xi <-  xinew
    Accept <- Accept+1
  }else{
    ext_xi <-  ext_xi
  }
  
  # Update fixed, nonlinear, and spatial effects of the count component
  
  SigmaBeta2 <- solve(Q_beta2+t(X2[w==1,])%*%Omega2%*%X2[w==1,])
  Mubeta2 <- SigmaBeta2 %*% (t(X2[w==1,])%*%Omega2%*%(z2-V2[w==1,]%*%ext_phi2-Z2[w==1,]%*%ext_theta2))
  ext_beta2 <- t(rmvn(1,Mubeta2,SigmaBeta2))
  #### 
  Q_theta2 <- (1/sigmatheta2)*Qtheta2  # assign prior to sigmatheta1
  Sigmatheta2 <- solve(Q_theta2+t(Z2[w==1,])%*%Omega2%*%Z2[w==1,])
  Mutheta2   <- Sigmatheta2%*%(t(Z2[w==1,])%*%Omega2%*%(z2-V2[w==1,]%*%ext_phi2-X2[w==1,]%*%ext_beta2))
  ext_theta2<-t(rmvn(1,Mutheta2,Sigmatheta2))
  ext_theta2 <- ext_theta2-mean(ext_theta2)
  ###
  ###
  mu2_1   <-(sig12/sig1)*ext_phi1 # mu1|2
  prec2_1<-(sig1/(sig1*sig2-sig12*sig21))*Q
  Sigmaphi2 <- solve(prec2_1+t(V2[w==1,])%*%Omega2%*%V2[w==1,])
  muphi21    <- Sigmaphi2%*%(prec2_1%*%mu2_1+t(V2[w==1,])%*%Omega2%*%(z2-X2[w==1,]%*%ext_beta2-Z2[w==1,]%*%ext_theta2))
  ext_phi2 <- t(rmvn(1,muphi21,Sigmaphi2))
  ext_phi2 <- ext_phi2-mean(ext_phi2)
  
  if(iter%%500==0) cat(iter,"\n")
  
  if(iter>Burn & iter%% thinning==0){
    iter2<-(iter-Burn)/thinning
    BETA1[, iter2] <-ext_beta1
    BETA2[,iter2] <-ext_beta2
    THETA1[,iter2]<-ext_theta1
    THETA2[,iter2]<-ext_theta2
    PHI1[,iter2]  <-ext_phi1
    PHI2[,iter2]  <-ext_phi2
    XI[iter2]     <-ext_xi
    HH[,,iter2]   <-H
  }
  
}

###########################

Est <- rbind(matrix(rowMeans(BETA1)),
             matrix(rowMeans(BETA2)),
             matrix(rowMeans(THETA1)),
             matrix(rowMeans(THETA2)),
             matrix(rowMeans(PHI1)),
             matrix(rowMeans(PHI2)),
             matrix(mean(XI)),
             matrix(mean(HH[1,1,])),
             matrix(mean(HH[2,2,])),
             matrix(mean(HH[1,2,]))
)
Interval_Func <- function(y,c=0.025)apply(y, 1, function(x)quantile(x,probs=c(c)))

Estlower <- rbind(matrix(Interval_Func(BETA1)),
             matrix(Interval_Func(BETA2)),
             matrix(Interval_Func(THETA1)),
             matrix(Interval_Func(THETA2)),
             matrix(Interval_Func(PHI1)),
             matrix(Interval_Func(PHI2)),
             matrix(Interval_Func(t(as.matrix(XI)))),
             matrix(Interval_Func(t(as.matrix(HH[1,1,])))),
             matrix(Interval_Func(t(as.matrix(HH[2,2,])))),
             matrix(Interval_Func(t(as.matrix(HH[1,2,]))))
)
Estupper <- rbind(matrix(Interval_Func(BETA1,c=0.975)),
                  matrix(Interval_Func(BETA2,c=0.975)),
                  matrix(Interval_Func(THETA1,c=0.975)),
                  matrix(Interval_Func(THETA2,c=0.975)),
                  matrix(Interval_Func(PHI1,c=0.975)),
                  matrix(Interval_Func(PHI2,c=0.975)),
                  matrix(Interval_Func(t(as.matrix(XI)),c=0.975)),
                  matrix(Interval_Func(t(as.matrix(HH[1,1,])),c=0.975)),
                  matrix(Interval_Func(t(as.matrix(HH[2,2,])),c=0.975)),
                  matrix(Interval_Func(t(as.matrix(HH[1,2,])),c=0.975))
)

TRU <- rbind(beta1,beta2,theta1,theta2,phi1,phi2,xi,sig1,sig2,sig12)

Result <- data.frame(True = TRU, Post.Mean=Est, CI2.5 =Estlower,CI9.75 = Estupper)
#Visualize estimates
Result %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = F, 
                font_size = 20)
#rerun medium and low with these iterations
############################
### residual Diagnostic ####
############################

## Dropped

################
# Plot estimate
################ For plots
require(rgdal)
require(ggplot2)
library(ggpubr)
library(tidyverse)
###############
# download shp file
path ='https://africaopendata.org/dataset/83582021-d8f7-4bfd-928f-e9c6c8cb1247/resource/372a616a-66cc-41f7-ac91-d8af8f23bc2b/download/nigeria-lgas.zip'
shp <- readOGR(dsn ="C:\\Users\\OsafuEgbon\\OneDrive\\Documents\\R\\Rwork\\ZNIB-20220626T053915Z-001\\nigeria-lgas\\new_lga_nigeria_2003.shp", stringsAsFactors = F)

shp.map <- fortify(shp, region = "STATE")

shp.map$count <- 0
s <- function(region,est){
  shp.map$count[which(shp.map$group==as.character(region))] <- est
  return(shp.map)
}
############ Particular for the Nigerian graph
ss  <- read.csv("https://raw.githubusercontent.com/eosafu/Poly--GammaInference/5a230ab3e24457bfa18ff682d80d669401995f50/NIGBND.csv")
state2 <-  ss$state2
cloc2 <-   ss$cloc2
############
pest <- list(phi1,phi2,   # True
             rowMeans(PHI1),rowMeans(PHI2)) # Estimated
titl <- c("True-binary comp.","True-count comp.","Est-binary comp.","Est-count comp.")
plt<-list()

for(k in 1: 4){
# True
est_data <- data.frame(region=state2,
                       est= pest[[k]][cloc2])

for (i in 1:37) {
  shp.map <-s(est_data[i,1],est_data[i,2]) 
}
shp.map$count[which(shp.map$group=="Lagos.2")] <- est_data[25,2]
shp.map$count[which(shp.map$group=="Lake.1")] <- est_data[25,2]

plt[[k]] <- ggplot(shp.map, aes(x = long, y = lat, group = group, fill = count)) +
        geom_polygon(colour = "black", size = 0.5, aes(group = group)) +
        theme_bw()+ scale_fill_gradient(name="Post. mean",high='white', low='orange')+
  labs(x= "longitude",y="latitude",title = titl[k])
}

ggarrange(plt[[1]],plt[[2]],
          plt[[3]],plt[[4]],nrow = 2,ncol = 2)


#### Plot non-linear Effect
par(mar = c(1, 2, 0.3, 0.3),mfrow=c(1,2))
plot(scale(rowMeans(THETA1)),type="l",lty=2,ylim=c(-1.5,3),xlab="id",ylab="Effect",xaxt="n")
points(scale(rowMeans(THETA1)))
lines(scale(theta1))
points(scale(theta1),pch=19)

plot(scale(rowMeans(THETA2)),type="l",lty=2,ylim=c(-2,3),xlab="id",ylab="",xaxt="n")
points(scale((rowMeans(THETA2))))
lines(scale(theta2))
points(scale((theta2)),pch=19)

### Trace plot
f <- function(j){
  
  eta2 <- X2%*%BETA2[,j]+Z2%*%THETA2[,j]+V2%*%PHI2[,j]
  eta1 <- X1%*%BETA1[,j]+Z1%*%THETA1[,j]+V1%*%PHI1[,j]
  
  ext_pi <- exp(eta1)/(1+exp(eta1))
  ext_pneg <- exp(eta2)/(1+exp(eta2))
  ext_xi=(XI[j])
  
  pm <- pi*dnbinom(yobs,size=xi,prob = (1-pneg))
  auxpm <- rep(0,length(pm))
  auxpm[yobs==0] <- (1-pi[yobs==0])
  pm <- log(auxpm+pm+min((auxpm+pm)[auxpm+pm>0]))
  logL <- sum(pm)
  
  return(logL)
}
par(mar = c(2, 2, 0.3, 0.3),mfrow=c(2,3))
plot(BETA1[2,],type="l")
plot(THETA1[2,],type="l")
plot(PHI1[2,],type="l")
plot(BETA2[2,],type="l")
plot(THETA2[2,],type="l")
plot(PHI2[2,],type="l")

