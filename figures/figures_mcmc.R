# Plot MCMC draws

library(bayesplot)
library(tidyverse)
library(vapoRwave)

vaporwave <- vapoRwave:::jazzCup_palette

setwd("~/Documents/School/Summer_2019/Research/MVSN-FMM")
load_dir <- "mcmc_draws_2019-07-12"
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
load(paste(load_dir,"/Y",sep = ""))

h <- length(BETA.list)
k <- ncol(beta.true.list[[1]])
p <- nrow(beta.true.list[[1]])
v <- nrow(delta.true)

burn <- 0

colnames(Y) <- paste("Y",seq(1,ncol(Y)),sep = "")
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
    ggtitle("Trace plots of Delta Coefficients")
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

## Y densities
y_density_plot <- data.frame(Y,c) %>%
    rename(clust = c) %>%
    gather(variable,value,-clust) %>%
    ggplot(.) +
    stat_density(aes(x = value),geom = "line") +
    facet_grid(cols = vars(variable), rows = vars(clust)) +
    theme_minimal() + 
    xlab("Y") + 
    ylab("Density") 
ggsave(y_density_plot,
       filename = paste(load_dir,"y_densities.pdf",
                        sep = "/"))

## Y means
y_mean_plot <- data.frame(Y,c) %>%
    rename(clust = c) %>%
    gather(variable,value,-clust) %>%
    group_by(variable,clust) %>%
    dplyr::summarize(mean_y = mean(value)) %>%
    ungroup() %>% 
    mutate(clust = as.factor(clust)) %>%
    ggplot(.,aes(x = variable, 
                 y = mean_y, 
                 color = clust,
                 group = clust)) + 
    geom_line() +
    theme_minimal()
ggsave(y_mean_plot,
       filename = paste(load_dir,"y_means.pdf",
                        sep = "/"))

