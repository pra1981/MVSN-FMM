library(sn)

load("Y")
load("X")
load("Z")
load("SIGMA")
load("BETA")

setwd("~/Documents/School/Summer_2019/Research/MVSN-FMM/mcmc_draws_2019-07-09/mcmc_draws_NOSKEW_2019-07-09/")

n <- nrow(Y)
K <- length(BETA.list)
J <- ncol(Y)
p <- ncol(X)
S <- nrow(Z)

lppd <- 0
pwaic <- 0

pb <- txtProgressBar(min = 0, max = n, style = 3)
for(i in 1:n)
{
    y_i <- Y[i,]
    x_i <- X[i,]
    z_i <- Z[,i]
    p_yi <- rep(0,S)
    for(s in 1:S)
    {
        k <- z_i[s]
        
        beta_k <- BETA.list[[k]]
        beta_k_s <- beta_k[s,]
        beta_ks <- matrix(beta_k_s,
                          nrow = p,
                          ncol = J,
                          byrow = TRUE)
        
        sigma_k <- SIGMA.list[[k]]
        sigma_k_s <- sigma_k[s,]
        sigma_ks <- matrix(sigma_k_s,
                           nrow = J,
                           ncol = J,
                           byrow = TRUE)
        
        zeta_iks <- x_i %*% beta_ks
        p_yi[s] <- dmvnorm(y_i,zeta_iks,sigma_ks)
    }
    setTxtProgressBar(pb, i)
    
    print(lppd <- lppd + log(mean(p_yi, na.rm = TRUE)))
    print(pwaic <- pwaic + var(log(p_yi)[is.finite(log(p_yi))], na.rm = TRUE))
}
WAIC <- -2*(lppd - pwaic)

close(pb)