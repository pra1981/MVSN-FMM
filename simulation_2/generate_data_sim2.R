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
library(mice)

#seed <- rnorm(1)
#set.seed(seed)
# No scientific notation
options(scipen=999)
options(warn=2)
# Generate data
k <- 4 # number of outcomes
p <- 2 # number of predictors 
h <- 3 # number of clusters
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
           byrow = TRUE)
)
psi.true.list <- list(
    c(2,3,4,5),
    c(-2,-3,-4,-5),
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

store = "simulation_2/Data"
save(Y,file = paste(store,"/Y",sep = ""))
save(X,file = paste(store,"/X",sep = ""))
save(Xstar,file = paste(store,"/Xstar",sep = ""))
save(beta.true.list,file = paste(store,"/beta_true",sep = ""))
save(delta.true,file = paste(store,"/delta_true",sep = ""))
save(sig2.true.list,file = paste(store,"/sig2_true",sep = ""))
save(psi.true.list,file = paste(store,"/psi_true",sep = ""))
save(pi.true,file = paste(store,"/pi_true",sep = ""))
save(c,file = paste(store,"/c_true",sep = ""))
save(W,file = paste(store,"/W",sep = ""))

# ampute and store data
Y_amp <- ampute(Y,
                mech = "MAR",
                prop = 0.3)$amp
save(Y_amp,file = paste(store,"/Y_amp",sep = ""))

# Generate M imputations a-la Z&R 
M <- 1000
Y_multi <- list(0)
mu_est <- colMeans(Y_amp,na.rm = TRUE)
sig2_est <- var(Y_amp,na.rm = TRUE)
pb <- txtProgressBar(min = 0, max = M, style = 3)
for(m in 1:M)
{
    Y_multi[[m]] <- t(apply(Y_amp, 
                          1, 
                          impute_conditional_mvn,
                          mu_est,
                          sig2_est))
    setTxtProgressBar(pb, m)
}
close(pb)

save(Y_multi,file = paste(store,"/Y_multi_imputed",sep = ""))
