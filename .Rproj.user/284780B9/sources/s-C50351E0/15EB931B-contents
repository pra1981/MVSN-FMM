# Plot MCMC draws

library(bayesplot)
library(tidyverse)
library(vapoRwave)

vaporwave <- vapoRwave:::jazzCup_palette


load_dir <- "mcmc_draws_2019-03-28"
load(paste(load_dir,"/BETA",sep = ""))
load(paste(load_dir,"/beta_true",sep = ""))
load(paste(load_dir,"/DELTA",sep = ""))
load(paste(load_dir,"/delta_true",sep = ""))
load(paste(load_dir,"/PI",sep = ""))
load(paste(load_dir,"/PSI",sep = ""))
load(paste(load_dir,"/psi_true",sep = ""))
load(paste(load_dir,"/SIGMA",sep = ""))
load(paste(load_dir,"/sig2_true",sep = ""))
load(paste(load_dir,"/Z",sep = ""))
load(paste(load_dir,"/c_true",sep = ""))

h <- length(BETA.list)
k <- ncol(beta.true.list[[1]])
p <- nrow(beta.true.list[[1]])
v <- nrow(delta.true)

burn <- 100

colnames(Z) <- paste("Z",seq(1,ncol(Z)),sep = "")
Z <- Z[-(1:burn),]
colnames(DELTA) <- paste("D",seq(1,v*(h-1)),sep = "")
DELTA <- DELTA[-(1:burn),]
colnames(PI) <- paste("PI",seq(1,ncol(PI)),sep = "")
PI <- PI[-(1:burn),]

for(l in 1:h)
{
    colnames(BETA.list[[l]]) <- paste("B",seq(1,p*k),sep = "")
    BETA.list[[l]] <- BETA.list[[l]][-(1:burn),]
    colnames(PSI.list[[l]]) <- paste("PSI",seq(1,length(psi.true.list[[1]])),sep = "")
    PSI.list[[l]] <- PSI.list[[l]][-(1:burn),]
    colnames(SIGMA.list[[l]]) <- paste("SIG",seq(1,ncol(SIGMA.list[[l]])),sep="")
    SIGMA.list[[l]] <- SIGMA.list[[l]][-(1:burn),]
}

## BETA plots

BETA_trace <- mcmc_trace(BETA.list) + 
    scale_color_manual(name = "Cluster",
                         values = vaporwave) + 
    ggtitle("Trace plots of Beta Coefficients")
ggsave(BETA_trace,filename = paste(load_dir,"beta_trace.pdf",sep = "/"))

## SIGMA plots

SIGMA_trace <- mcmc_trace(SIGMA.list) + 
    scale_color_manual(name = "Cluster",
                       values = vaporwave) + 
    ggtitle("Trace plots of Sigma Coefficients")
ggsave(SIGMA_trace,filename = paste(load_dir,"sigma_trace.pdf",sep = "/"))

## PSI plots

PSI_trace <- mcmc_trace(PSI.list) + 
    scale_color_manual(name = "Cluster",
                       values = vaporwave) + 
    ggtitle("Trace plots of Psi Coefficients")
ggsave(PSI_trace,filename = paste(load_dir,"psi_trace.pdf",sep = "/"))

## DELTA plots

DELTA_trace <- mcmc_trace(DELTA) + 
    scale_color_manual(name = "Cluster",
                       values = vaporwave) + 
    ggtitle("Trace plots of Psi Coefficients")
ggsave(DELTA_trace,filename = paste(load_dir,"delta_trace.pdf",sep = "/"))

## PI plots

PI_trace <- mcmc_trace(PI) + 
    scale_color_manual(name = "Cluster",
                       values = vaporwave)+ 
    ggtitle("Trace plots of Pi Coefficients") + 
    scale_y_continuous(limits = c(0,1))
ggsave(PI_trace,filename = paste(load_dir,"pi_trace.pdf",sep = "/"))

## Z plots 

z_means <- colMeans(Z[-(1:burn),])
ggplot() + 
    geom_line(aes(x = seq_along(z_means),y = c),color = vaporwave[1]) + 
    geom_line(aes(x = seq_along(z_means),y = z_means),color = vaporwave[2])
