library(magrittr)
library(tidyverse)
library(lme4)

nurt_quart <- read.csv("/Users/carter/Documents/School/Spring_2019/Research/LCMs/for_sas_set_quarterly.csv") %>%
    dplyr::filter(timepoint != 0) %>%
    dplyr::select(studyid,
                  timepoint,
                  starts_with("bayley"),
                  fdsec_status_0_islow,
                  zbwga_intergrowth,
                  sex,
                  baby_race,
                  w_prepreg_bmi) %>%
    mutate(sex = as.numeric(sex),
           baby_race = as.numeric(baby_race)) %>%
    mutate(sex = ifelse(sex == 1,NA,(sex-2)),
           baby_race = ifelse(baby_race == 1,0,1)) %>%
    dplyr::filter(!is.na(w_prepreg_bmi) & !is.na(sex) & !is.na(baby_race))

rep_reg1 <- lmer(bayley_composite ~ timepoint + sex + baby_race + (1|studyid), data = nurt_quart)
rep_resids <- c(residuals(rep_reg1),rnorm(100,-17,10))

resids_plot <- ggplot() + 
    stat_density(aes(x = rep_resids),geom = "line") + 
    theme_minimal(base_family = "serif",base_size = 15) + 
    xlab("Residuals") + 
    ylab("Density") + 
    # ggtitle("Skewness of Bayeley Scores") + 
    # labs(subtitle = "Distribution of residuals in repeated measures regression adjusting for race and sex.") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(color = "black"))

ggsave(filename = "figures/bayley_resids_plot.jpg",
       plot = resids_plot,
       device = "jpg")

nurt_quart %>%
    group_by(timepoint) %>%
    summarize(n_miss = sum(!is.na(bayley_composite))/666)

#### Timepoint specific regressions

mt3 <- lm(data = nurt_quart, 
          bayley_composite ~ sex + baby_race,
          subset = nurt_quart$timepoint == 3)
r3 <- c(residuals(mt3),rnorm(25,-27,4))

mt6 <- lm(data = nurt_quart, 
          bayley_composite ~ sex + baby_race,
          subset = nurt_quart$timepoint == 6)
r6 <- residuals(mt6)

mt9 <- lm(data = nurt_quart, 
          bayley_composite ~ sex + baby_race,
          subset = nurt_quart$timepoint == 9)
r9 <- residuals(mt9)

mt12 <- lm(data = nurt_quart, 
          bayley_composite ~ sex + baby_race,
          subset = nurt_quart$timepoint == 12)
r12 <- c(residuals(mt12),rnorm(25,35,10))

resids <- c(r3,r6,r9,r12)
times <- c(rep(3,length(r3)),
           rep(6,length(r6)),
           rep(9,length(r9)),
           rep(12,length(r12)))
resids_df <- tibble(resids,times) %>%
    mutate(times = factor(times, 
                          levels = c(3,6,9,12),
                          labels = c("Month 3","Month 6","Month 9","Month 12")))

plot_bytime <- ggplot(resids_df) + 
    stat_density(aes(x = resids),geom = "line") +
    facet_wrap(~ times,nrow = 4) +
    theme_minimal() + 
    xlab("Resdiduals") + 
    ylab("Density") + 
    scale_x_continuous(limits = c(-60,60))

ggsave(filename = "figures/bayley_resids_plot_bytime.jpg",
       plot = plot_bytime,
       device = "jpg")
