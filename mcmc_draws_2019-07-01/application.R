# Application of model to Nurture data
# Carter Allen
# April 05, 2019

# Packages
library(magrittr)
library(sn) # for rsn
library(truncnorm) # for rtruncnorm
library(mvtnorm) # for rmvtnorm
library(ggplot2)
library(patchwork)
library(mvnmle)
library(Hmisc) # For rMultinom function (easy sampling of class indicators C)
library(bayesplot)
library(coda)
library(MCMCpack)
library(matrixsampling)
library(QRM)        # For rGumbel function
library(nnet)       # To fit multinom function
library(BayesLogit) # For rpg function
library(tidyverse)

setwd("~/Documents/School/Summer_2019/Research/MVSN-FMM/mcmc_draws_2019-07-01")
X <- load(X)
Y <- load(Y)

k <- ncol(Y)
p <- ncol(X)
h <- 3
n <- nrow(Y)

w3 <- nurt_quart$sex
w4 <- nurt_quart$baby_race
W <- cbind(1,w3,w4) # design matrix
v <- ncol(W) # number of multinomial predictors

c <- rep_len(seq(1:h),length.out = n)

# Generate class-specific parameters and data 
# Could also define priors here if you wanted priors to be class specific
beta.mle.list <- list(0)
beta.init.list <- list(0)
sig2.mle.list <- list(0)
sig2.init.list <- list(0)
psi.init.list <- list(0)
V0.list <- list(0)
L0.list <- list(0)
B0.list <- list(0)

# Generate class-specific parameters and data 
# Could also define priors here if you wanted priors to be class specific
beta.true.list <- list(0)
beta.mle.list <- list(0)
beta.init.list <- list(0)
betastar.mle.list <- list(0)
sig2.true.list <- list(0)
sig2.mle.list <- list(0)
sig2.init.list <- list(0)
psi.true.list <- list(0)
psi.init.list <- list(0)
V0.list <- list(0)
L0.list <- list(0)
B0.list <- list(0)

for(l in 1:h)
{
    X.l <- as.matrix(X[c == l,])
    Y.l <- as.matrix(Y[c == l,])
    n.l <- sum(c == l) # number of observations truly assigned to cluster l

    # Cluster specific MLEs
    sn_fit <- sn::msn.mle(X.l,Y.l)
    alpha.l <- sn_fit$dp$alpha
    Omega.l <- sn_fit$dp$Omega
    psi.mle.l <- (1/sqrt(1+t(alpha.l) %*% Omega.l %*% alpha.l)) %*% t(Omega.l %*% alpha.l)
    beta.mle.l <- sn_fit$dp$beta
    beta.mle.list[[l]] <- beta.mle.l
    betastar.mle.l <- rbind(beta.mle.l,psi.mle.l) # MLE of Beta_star matrix
    betastar.mle.list[[l]] <- betastar.mle.l
    sig2.mle.l <- Omega.l# MLE of v-cov matrix
    sig2.mle.list[[l]] <- sig2.mle.l
    sig2.init.list[[l]] <- sig2.mle.l # convinient to just set the inits to sig2.mle.l here
    
    # Cluster specific inits
    psi.init.l <- psi.mle.l
    psi.init.list[[l]] <- as.vector(psi.init.l)
    beta.init.l <- beta.mle.l
    beta.init.list[[l]] <- beta.init.l
    B0.list[[l]] <- matrix(rep(0,(p+1)*k),nrow = p+1,ncol = k) # prior mean of beta
    V0.list[[l]] <- sig2.mle.l # kxk hyperparam for Sigma
    L0.list[[l]] <- diag(rep(1,p+1))  # appears in Kronecker of beta0.vec prior
}

# Define prior structure
print("Defining prior structure")
a <- rep(2,h) # prior parameter vector for pi1,...,pih
nu0 <- rep(2,h)  # scalar df hyperparam for Sigma (same across classes)
V0 <- diag(rep(1,k)) # kxk hyperparam for Sigma (same across classes)
B0 <- rbind(beta.init.list[[1]],rep(0,k)) # zero matrix
L0 <- diag(rep(1,p+1)) # appears in Kronecker of beta0.vec prior
delta0 <- rep(0,v) # prior mean for delta coefficients (multinomial regression)
S0 <- diag(1,v) # prior covariance for delta coefficients (multinomial regression)

# Sample storage
nsim <- 1000 # number of iterations
burn <- 300 # number of iterations to save
n.iter <- nsim - burn # number of saved iterations
Z <- matrix(0,nrow = n.iter,ncol = n) # large matrix where each row is the value of z at a specific iteration
PI <- matrix(0,nrow = n.iter,ncol = h) # matrix w/ each row as pi vector at each iteration
BETA.list <- vector("list",h) # storage for Beta (vectorized by row)
SIGMA.list <- vector("list",h) # storage for S matrix (vectorized by row)
PSI.list <- vector("list",h) # storage for psi vector
DELTA <- matrix(0,nrow = n.iter,ncol = v*(h-1))

# Initialized values
# pi <- rep(1/h,h) # mixing probabilities
pi <- table(c)/n
print(pi)
z <- c
Bn.list <- B0.list
Ln.list <- L0.list
nun <- nu0
Vn.list <- V0.list
beta.list <- beta.init.list
sig2.list <- sig2.init.list
psi.list <- psi.init.list

delta <- matrix(0,nrow = v, ncol = h-1)
t.init <- rtruncnorm(nrow(X),0,Inf)
Xstar <- cbind(X,t.init)
Y <- as.matrix(Y)

# MCMC
start.time<-proc.time()
print(paste("Started MCMC of",nsim,"iterations with a burn in of",burn))
pb <- txtProgressBar(min = 0, max = nsim, style = 3)
for(i in 1:nsim)
{
    # Step 0:
    # perform conditional mv-normal imputation
    for(j in 1:n)
    {
        # perform conditional normal imputation
        na.ind <- na.mat[j,]
        y.j <- Y[j,] # outcome for observation j
        if(any(na.ind) & !all(na.ind))
        {
            a.j <- y.j[!na.ind] # observed components of y.j
            xstar.j <- Xstar[j,] # covariates for observation j
            z.j <- z[j] # current cluster id of y.j
            A.j <- sig2.list[[z.j]] # covmatrix of y.j
            betastar.j <- rbind(beta.list[[z.j]],psi.list[[z.j]]) # betstar for y.j
            mu.j <- xstar.j %*% betastar.j # mean of observation j if it belonged to each class
            
            mu.1 <- mu.j[na.ind]
            mu.2 <- mu.j[!na.ind]
            
            sig.11 <- A.j[na.ind,na.ind]
            sig.12 <- A.j[na.ind,!na.ind]
            sig.21 <- A.j[!na.ind,na.ind]
            sig.22 <- A.j[!na.ind,!na.ind]
            
            mu.cond <- mu.1 + sig.12 %*% solve(sig.22) %*% (a.j - mu.2)
            sig.cond <- sig.11 - sig.12 %*% solve(sig.22) %*% sig.21
            
            y.imp <- rmvnorm(1,mu.cond,sig.cond)
            y.j[na.ind] <- y.imp
            Y[j,] <- y.j
        }
    }
    Y <- as.matrix(Y)
    # End step 0
    
    # Step 3:
    # Update pi1,...,piK 
    # pi <- rdirichlet(1,a + n.z)
    eta<-cbind(rep(0,n),W%*%delta)
    pi<-exp(eta)/(1+apply(as.matrix(exp(eta[,-1])),1,sum))
    # End Step 3
    
    # Step 1A:
    # Update z_j from multinomial conditional posterior
    # for each observed y, compute the probability of belonging in each class
    # this is done iteratively in a loop but there is probably a better way
    for(j in 1:n)
    {
        y.j <- Y[j,] # outcome for observation j
        xstar.j <- Xstar[j,] # covariates for observation j
        # mu.j <- matrix(rep(0,k*h),nrow = h,ncol = k) # empty storage for class specific X %*% betas
        p.j <- rep(0,h)
        for(l in 1:h)
        {
            betastar.l <- rbind(beta.list[[l]],psi.list[[l]])
            mu.j <- xstar.j %*% betastar.l # mean of observation j if it belonged to each class
            sig2.l <- sig2.list[[l]]
            p.j[l] <- dmvnorm(y.j,mu.j,sig2.l)
        }
        pi.j <- pi[j,]
        p.z <- pi.j*p.j/sum(pi.j*p.j)
        # print(p.z)
        z[j] <- which(rmultinom(1,1,p.z) == 1) # this is very slow 
        # z <- c # for testing purposes hold z to true value
    }
    # convert z to factor so that empty classes are included as 0 counts below and not dropped
    z <- factor(z,levels = 1:h)
    # # compute the sample size within each class
    # # n should be a vector of length h
    # # if there is an empty class then the program should stop
    n.z <- as.vector(unname(table(z))) # gives the number of members currently in each class
    # print(n.z/n)
    if(any(n.z == 0)) # check for empty clusters, if so stop
    {
        print("MCMC terminated due to an empty class.")
        break
    }
    if(any(n.z == 1)) # check for singleton clusters, if so stop
    {
        print("MCMC terminated due to a singleton class.")
        break
    }
    # End Step 1A
    
    # Step 1B:
    # Update multinomial regression parameters
    for(l in 1:(h-1))
    {
        delta.l <- delta[,l]
        delta.notl <- delta[,-l]
        u.l <- 1*(c == l+1)
        c.l <- log(1 + exp(W %*% delta.notl))
        eta <- W %*% delta.l - c.l
        w <- rpg(n, 1, eta)
        z.l <- (u.l - 1/2)/w + c.l 
        
        V <- solve(S0 + crossprod(W*sqrt(w)))  
        M <- V %*% (S0 %*% delta0 + t(w*W) %*% z.l)
        delta.l <- c(rmvnorm(1,M,V))
        delta[,l] <- delta.l
    }
    # End Step 1B
    
    # Step 2: 
    # Update class-specific regression parameters
    for(l in 1:h) # loop through each cluster
    {
        X.l <- X[z == l,] # all covariates for those in class l
        Y.l <- Y[z == l,] # all outcomes for those in class l
        n.l <- sum(z == l) # current class specific sample size
        psi.l <- psi.list[[l]] # current class specific psis 
        B0.l <- B0.list[[l]]
        V0.l <- V0.list[[l]]
        L0.l <- L0.list[[l]]
        nu0.l <- nu0[l]
        nun.l <- nun[l]
        Bn.l <- Bn.list[[l]]
        Vn.l <- Vn.list[[l]]
        Ln.l <- Ln.list[[l]]
        sig2.l <- sig2.list[[l]]
        beta.l <- beta.list[[l]]
        betastar.l <- rbind(beta.l,psi.l)
        
        # Update t
        A <- solve(1 + t(psi.l) %*% solve(sig2.l) %*% psi.l)
        t.l <- rep(0,n.l)
        for(j in 1:n.l)
        {
            ai <- A %*% t(psi.l) %*% solve(sig2.l) %*% t((t(Y.l[j,]) - t(X.l[j,]) %*% beta.l))
            t.l[j] <- rtruncnorm(n = 1,a = 0, b = Inf,mean = ai, sd = sqrt(A))
        }
        Xstar.l <- cbind(X.l,t.l) # add updated t's back to Xstar matrix
        
        # Update sigma
        # Same as matrix normal regression update
        Vn.l <- V0.l + t(Y.l - Xstar.l %*% Bn.l) %*% (Y.l - Xstar.l %*% Bn.l) + t(Bn.l - B0.l) %*% L0.l %*% (Bn.l - B0.l)
        nun.l <- nu0.l + n.l
        sig2.l <- riwish(nun.l,Vn.l)
        sig2.list[[l]] <- sig2.l
        # sig2.list <- sig2.true.list # used for testing purposes
        
        # Update beta
        # Same as matrix normal regression update with added psi
        Bn.l <- solve(t(Xstar.l) %*% Xstar.l + L0.l) %*% (t(Xstar.l) %*% Y.l + L0.l %*% B0.l)
        Bn.list[[l]] <- Bn.l
        Ln.l <- t(Xstar.l) %*% Xstar.l + L0.l
        Ln.list[[l]] <- Ln.l
        betastar.l <- rmatrixnormal(1,Bn.l,solve(Ln.l),sig2.l, checkSymmetry = FALSE)
        betastar.l <- matrix(betastar.l,nrow = p+1, ncol = k)
        beta.l <- betastar.l[1:p,]
        psi.l <- betastar.l[p+1,]
        beta.list[[l]] <- beta.l
        psi.list[[l]] <- psi.l
        # beta.list <- beta.true.list # used for testing purposes
        # psi.list <- psi.true.list # used for testing purposes
    }
    # End step 2
    
    
    # Store results
    if (i > burn)
    {
        j <- i-burn
        Z[j,] <- z
        PI[j,] <- table(z)/n
        DELTA[j,] <- c(delta)
        for(l in 1:h)
        {
            BETA.list[[l]] <- rbind(BETA.list[[l]],c(beta.list[[l]]))
            SIGMA.list[[l]] <- rbind(SIGMA.list[[l]],c(sig2.list[[l]]))
            PSI.list[[l]] <- rbind(PSI.list[[l]],psi.list[[l]])
        }
    }
    
    setTxtProgressBar(pb, i)
}
close(pb)
run.time<-proc.time()-start.time
print(paste("Finished MCMC after",run.time[1],"seconds"))

run_date <- Sys.Date()
store <- paste("APPLICATION_",run_date,sep = "")
if(!dir.exists(store))
{
    dir.create(store)
}
meta_file <- paste(store,"/META.txt",sep = "")
file.create(meta_file)
write(paste("MCMC Run Date:",Sys.time()),file = meta_file, append = TRUE)
#write(paste("Seed:",seed),file = meta_file, append = TRUE)
write(paste("Number of subjects (n):",n),file = meta_file, append = TRUE)
write(paste("Number of outcomes (k):",k),file = meta_file, append = TRUE)
write(paste("Number of predictors (p):",p),file = meta_file, append = TRUE)
write(paste("Number of multinomial predictors (v):",v),file = meta_file, append = TRUE)
write(paste("Number of clusters (h):",h),file = meta_file, append = TRUE)
# write(paste("Number of observations (N):",),file = meta_file, append = TRUE)
# write(paste("Missing Mechanism:",mech.miss),file = meta_file, append = TRUE)
# write(paste("Missing Proportion:",prop.miss),file = meta_file, append = TRUE)
write(paste("Number of simulations (nsim):",nsim),file = meta_file, append = TRUE)
write(paste("Number of burn-in simulations (burn):",burn),file = meta_file, append = TRUE)

#save(Y,file = paste(store,"/Y",sep = ""))
#save(X,file = paste(store,"/X",sep = ""))
save(Z,file = paste(store,"/Z_noskew",sep = ""))
save(PI,file = paste(store,"/PI_noskew",sep = ""))
save(DELTA,file = paste(store,"/DELTA_noskew",sep = ""))
save(BETA.list,file = paste(store,"/BETA_noskew",sep = ""))
save(SIGMA.list,file = paste(store,"/SIGMA_noskew",sep = ""))
save(PSI.list,file = paste(store,"/PSI_noskew",sep = ""))

