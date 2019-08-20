# Multivariate skew normal regression with clustering
# Includes built-in mvn conditional imputation
# Includes regression on multinomial regression coefficients
# Carter Allen
# March 19, 2019

# Packages
library(sn) # for rsn
library(truncnorm) # for rtruncnorm
library(mvtnorm) # for rmvtnorm
library(mvnmle)
library(Hmisc) # For rMultinom function (easy sampling of class indicators C)
library(coda)
library(MCMCpack)
library(matrixsampling)
library(QRM)        # For rGumbel function
library(nnet)       # To fit multinom function
library(BayesLogit) # For rpg function

#seed <- rnorm(1)
#set.seed(seed)
# No scientific notation
options(scipen=999)
options(warn=2)
# Generate data
k <- 4 # number of outcomes
p <- 2 # number of predictors 
h <- 2 # number of clusters
n <- 1000 # number of observations (takes alot of data to estimate this many parameters)
print(paste("Generating data with",n,"observations of dimension",k,
            "with",p,"covariates and",h,"clusters"))

N <- n*k # number of measurments 
P <- matrix(rep(rep(1/h,h),each=n),n,h) # probability matrix governing cluster membership

# Generate covariates
X <- matrix(rnorm(n*p),nrow = n,ncol = p) # random predictors
X[,1] <- 1
t <- rtruncnorm(n,0,Inf,0,1) # trunc norm random effects
Xstar <- cbind(X,t) # Combine X and t

w1 <- rbinom(n,1,0.5) # binary predictor
W <- cbind(1,w1) # design matrix
v <- ncol(W) # number of multinomial predictors
delta.true <- matrix(rtruncnorm(n = v*(h-1),
                                a = -1, 
                                b = 1,
                                mean = 0),
                     nrow = v,
                     ncol = h-1) # true multinomial regression coefficients
eta.true <- W %*% delta.true # true eta term for each observation
U <- matrix(0,nrow = n,ncol = h) # make empty U matrix
U[,1] <- rGumbel(n, mu = 0, sigma = 1)
for(i in 2:h) # populate U matrix with gumbel values
{
    U[,i] <- rGumbel(n, mu = eta.true[,i-1], sigma = 1)
}
c <- c(apply(U,1,which.max))
W <- W[order(c),]
eta.true <- eta.true[order(c),]
U <- U[order(c),]
c <- sort(c)
pi.true <- table(c)/n

# Generate class-specific parameters and data 
# Could also define priors here if you wanted priors to be class specific
beta.true.list <- list(0)
beta.mle.list <- list(0)
beta.init.list <- list(0)
betastar.mle.list <- list(0)
betastar.true.list <- list(0)
sig2.true.list <- list(0)
sig2.mle.list <- list(0)
sig2.init.list <- list(0)
psi.true.list <- list(0)
psi.init.list <- list(0)
V0.list <- list(0)
L0.list <- list(0)
B0.list <- list(0)
Y <- NULL

sig2.true.ar1 <- 0.5^(abs(outer(1:k,1:k,"-")))
beta.true.list <- list(
    matrix(c(110,115,120,125,
             1,1.5,2,2.5),
           nrow = p,
           ncol = k,
           byrow = TRUE),
    matrix(c(90,85,80,75,
             -1,-1.5,-2,-2.5),
           nrow = p,
           ncol = k,
           byrow = TRUE),
    matrix(c(100,100,100,100,
             -1,1,-1,1),
           nrow = p,
           ncol = k,
           byrow = TRUE),
    matrix(c(100,100,100,100,
             -1,1,-1,1),
           nrow = p,
           ncol = k,
           byrow = TRUE)
)
psi.true.list <- list(
    c(2,3,4,5),
    c(-2,-3,-4,-5),
    c(0,1,0,1),
    c(0,1,0,1)
)
for(l in 1:h)
{
    X.l <- as.matrix(X[c == l,])
    Xstar.l <- Xstar[c == l,] # subset of covariates belonging to cluster l
    n.l <- sum(c == l) # number of observations truly assigned to cluster l
    
    # Cluster specific betas and psis
    beta.true.l <- beta.true.list[[l]]
    psi.true.l <- psi.true.list[[l]]
    betastar.true.l <- rbind(beta.true.l,psi.true.l)
    betastar.true.list[[l]] <- betastar.true.l
    
    sig2.true.list[[l]] <- sig2.true.l <- sig2.true.ar1 # Add to list of true class-specific sig2 matrices
    
    # Cluster specific Ys
    E.l <- rmvnorm(n.l,rep(0,k),sig2.true.l) # zero centered MVN error matrix
    Y.l <- Xstar.l %*% betastar.true.l + E.l # generated data
    Y <- rbind(Y,Y.l) # add class specific observations to total outcome matrix
    
    # Cluster specific MLEs
    beta.mle.l <- solve(t(X.l) %*% X.l) %*% t(X.l) %*% Y.l # MLE of Beta matrix
    beta.mle.list[[l]] <- beta.mle.l
    betastar.mle.l <- solve(t(Xstar.l) %*% Xstar.l) %*% t(Xstar.l) %*% Y.l # MLE of Beta_star matrix
    betastar.mle.list[[l]] <- betastar.mle.l
    sig2.mle.l <- mlest(Y.l - Xstar.l %*% betastar.true.l)$sigmahat # MLE of v-cov matrix
    sig2.mle.list[[l]] <- sig2.mle.l
    sig2.init.list[[l]] <- sig2.true.ar1 # convinient to set the inits to sig2.mle.l here
    
    # Cluster specific inits/priors
    # psi.init.l <- betastar.mle.l[p+1,]
    psi.init.l <- psi.true.l
    psi.init.list[[l]] <- psi.init.l
    beta.init.l <- beta.true.l
    beta.init.list[[l]] <- beta.init.l
    # B0.list[[l]] <- matrix(rep(0,(p+1)*k),nrow = p+1,ncol = k) # prior mean of beta
    B0.list[[l]] <- betastar.true.l
    # V0.list[[l]] <- sig2.mle.l # kxk hyperparam for Sigma
    V0.list[[l]] <- sig2.true.l
    L0.list[[l]] <- diag(rep(1,p+1))  # appears in Kronecker of beta0.vec prior
}

# Define prior structure
print("Defining prior structure")
nu0 <- rep(6,h)  # scalar df hyperparam for Sigma (same across classes)
delta0 <- rep(0,v) # prior mean for delta coefficients (multinomial regression)
S0 <- diag(1,v) # prior covariance for delta coefficients (multinomial regression)

# Sample storage
nsim <- 1000 # number of iterations
burn <- 0 # number of iterations to save
n.iter <- nsim - burn # number of saved iterations
Z <- matrix(0,nrow = n.iter,ncol = n) # large matrix where each row is the value of z at a specific iteration
PI <- matrix(0,nrow = n.iter,ncol = h) # matrix w/ each row as pi vector at each iteration
BETA.list <- vector("list",h) # storage for Beta (vectorized by row)
SIGMA.list <- vector("list",h) # storage for S matrix (vectorized by row)
PSI.list <- vector("list",h) # storage for psi vector
DELTA <- matrix(0,nrow = n.iter,ncol = v*(h-1))
OMEGA.list <- vector("list",h) # storage for omega matrix (vectorized by row)
ALPHA.list <- vector("list",h)

# Initialized values
pi <- table(c)/n
print(pi)
# z <- sample(1:h,n,replace = TRUE,prob = pi) # just generate some random class labels to start with
z <- c
Bn.list <- B0.list
Ln.list <- L0.list
nun <- nu0
Vn.list <- V0.list
beta.list <- beta.init.list
sig2.list <- sig2.init.list
omega.list <- sig2.init.list
psi.list <- psi.init.list
alpha.list <- psi.init.list
delta <- matrix(0,nrow = v, ncol = h-1)

load("~/Documents/School/Summer_2019/Research/MVSN-FMM/simulation_3/MSN_K4_DATA/Y")

# MCMC
start.time<-proc.time()
print(paste("Started MCMC of",nsim,"iterations with a burn in of",burn))
pb <- txtProgressBar(min = 0, max = nsim, style = 3)
for(i in 1:nsim)
{
    Y <- as.matrix(Y)
    # End step 0
    
    # Step 3:
    # Update pi1,...,piK 
    eta <- cbind(rep(0,n),W%*%delta)
    pi <- exp(eta)/(1+apply(as.matrix(exp(eta[,-1])),1,sum))
    # End Step 3
    
    # Step 1A:
    # Update z_j from multinomial conditional posterior
    # for each observed y, compute the probability of belonging in each class
    # this is done iteratively in a loop but there is probably a better way
    for(j in 1:n)
    {
        y.j <- Y[j,] # outcome for observation j
        xstar.j <- Xstar[j,] # covariates for observation j
        p.j <- rep(0,h)
        for(l in 1:h)
        {
            betastar.l <- rbind(beta.list[[l]],psi.list[[l]])
            mu.j <- xstar.j %*% betastar.l # mean of observation j if it belonged to each class
            sig2.l <- sig2.list[[l]]
            p.j[l] <- dmvnorm(y.j,mu.j,sig2.l) 
            if(any(p.j == 0))
            {
                p.j[l] <- dmvnorm(y.j,mu.j,sig2.l)
            }
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
    W <- as.matrix(W)
    for(l in 1:(h-1))
    {
        delta.l <- delta[,l]
        delta.notl <- delta[,-l]
        u.l <- 1*(c == l+1)
        c.l <- log(1 + exp(rowSums(W %*% delta.notl)))
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
        X.l <- as.matrix(X[z == l,]) # all covariates for those in class l
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
        omega.l <- sig2.l + outer(psi.l,psi.l)
        omega.list[[l]] <- omega.l
        
        # Update beta
        # Same as matrix normal regression update with added psi
        Bn.l <- solve(t(Xstar.l) %*% Xstar.l + L0.l) %*% (t(Xstar.l) %*% Y.l + L0.l %*% B0.l)
        Bn.list[[l]] <- Bn.l
        Ln.l <- t(Xstar.l) %*% Xstar.l + L0.l
        Ln.list[[l]] <- Ln.l
        betastar.l <- rmatrixnormal(1,Bn.l,solve(Ln.l), sig2.l, checkSymmetry = FALSE)
        betastar.l <- matrix(betastar.l,nrow = p+1, ncol = k)
        beta.l <- betastar.l[1:p,]
        psi.l <- betastar.l[p+1,]
        beta.list[[l]] <- beta.l
        psi.list[[l]] <- psi.l
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
            # ALPHA.list[[l]] <- rbind(ALPHA.list[[l]],alpha.list[[l]])
            OMEGA.list[[l]] <- rbind(OMEGA.list[[l]],c(omega.list[[l]]))
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
store <- paste("mcmc_draws_",run_date,sep = "")
if(!dir.exists(store))
{
    dir.create(store)
} 
meta_file <- paste(store,"/META.txt",sep = "")
file.create(meta_file)
write(paste("MCMC Run Date:",Sys.time()),file = meta_file, append = TRUE)
write(paste("Number of subjects (n):",n),file = meta_file, append = TRUE)
write(paste("Number of outcomes (k):",k),file = meta_file, append = TRUE)
write(paste("Number of predictors (p):",p),file = meta_file, append = TRUE)
write(paste("Number of multinomial predictors (v):",v),file = meta_file, append = TRUE)
write(paste("Number of clusters (h):",h),file = meta_file, append = TRUE)
write(paste("Number of observations (N):",N),file = meta_file, append = TRUE)
write(paste("Number of simulations (nsim):",nsim),file = meta_file, append = TRUE)
write(paste("Number of burn-in simulations (burn):",burn),file = meta_file, append = TRUE)

save(Y,file = paste(store,"/Y",sep = ""))
save(X,file = paste(store,"/X",sep = ""))
save(Xstar,file = paste(store,"/Xstar",sep = ""))
save(Z,file = paste(store,"/Z",sep = ""))
save(PI,file = paste(store,"/PI",sep = ""))
save(DELTA,file = paste(store,"/DELTA",sep = ""))
save(BETA.list,file = paste(store,"/BETA",sep = ""))
save(SIGMA.list,file = paste(store,"/SIGMA",sep = ""))
save(PSI.list,file = paste(store,"/PSI",sep = ""))
save(beta.true.list,file = paste(store,"/beta_true",sep = ""))
save(delta.true,file = paste(store,"/delta_true",sep = ""))
save(sig2.true.list,file = paste(store,"/sig2_true",sep = ""))
save(psi.true.list,file = paste(store,"/psi_true",sep = ""))
save(pi.true,file = paste(store,"/pi_true",sep = ""))
save(c,file = paste(store,"/c_true",sep = ""))
save(W,file = paste(store,"/W",sep = ""))

