impute_conditional_mvn <- function(y,mu,sig2)
{
    na.ind = is.na(y)
    if(any(na.ind))
    {
        a = y[!na.ind]
        A = sig2
        mu.1 = mu[na.ind]
        mu.2 = mu[!na.ind]
        
        sig.11 = A[na.ind,na.ind]
        sig.12 = A[na.ind,!na.ind]
        sig.21 = A[!na.ind,na.ind]
        sig.22 = A[!na.ind,!na.ind]
        
        mu.cond = mu.1 + sig.12 %*% solve(sig.22) %*% (a - mu.2)
        sig.cond = sig.11 - sig.12 %*% solve(sig.22) %*% sig.21
        
        y.imp = rmvnorm(1,mu.cond,sig.cond)
        y[na.ind] = y.imp
    }
    return(y)
}