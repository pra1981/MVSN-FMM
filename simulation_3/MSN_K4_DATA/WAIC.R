library(sn)
setwd("~/Documents/School/Summer_2019/Research/MVSN-FMM/simulation_3/MSN_K4_DATA/mcmc_draws_2019-08-05/")

load("Y")
load("X")
load("Xstar")
load("W")
load("Z")
load("SIGMA")
load("PSI")
load("BETA")
load("DELTA")

n <- nrow(Y)
K <- length(BETA.list)
J <- ncol(Y)
p <- ncol(X)
r <- ncol(W)
S <- nrow(Z)

lppd <- 0
pwaic <- 0

pb <- txtProgressBar(min = 0, max = n, style = 3)
for(i in 1:n)
{
    y_i <- Y[i,]
    x_i <- X[i,]
    z_i <- Z[,i]
    t_i <- Xstar[i,p+1]
    p_yi <- rep(0,S)
    for(s in 1:S)
    {
        k <- z_i[s]
        
        beta_k <- BETA.list[[k]]
        beta_k_s <- beta_k[s,]
        beta_ks <- matrix(beta_k_s,
                          nrow = p,
                          ncol = J,
                          byrow = FALSE)
        
        sigma_k <- SIGMA.list[[k]]
        sigma_k_s <- sigma_k[s,]
        sigma_ks <- matrix(sigma_k_s,
                           nrow = J,
                           ncol = J,
                           byrow = TRUE)
        
        psi_k <- PSI.list[[k]]
        psi_ks <- psi_k[s,]
        
        zeta_iks <- x_i %*% beta_ks
        
        mu_iks <- zeta_iks + t_i %*% psi_ks
        
        p_yi[s] <- dmvnorm(y_i,mu_iks,sigma_ks)
    }
    setTxtProgressBar(pb, i)
    
    lppd <- lppd + log(mean(p_yi, na.rm = TRUE))
    pwaic <- pwaic + var(log(p_yi), na.rm = TRUE)
}
WAIC <- -2*(lppd - pwaic)
save(WAIC,file = "WAIC")

close(pb)