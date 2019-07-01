library(magrittr)
library(tidyverse)

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
