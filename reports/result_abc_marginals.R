# tracer et eneregistrer les résultats de ABC pour chaque colonie et chaque ABC

################################
###         PACKAGES         ###
################################
library(ks)
library(igraph)
library(mvtnorm)
library(tidyverse)
library(grid)

library(stringr)
library(reshape2)

### import functions 
# Kernel density estimate
source("reports/kde.R")
source("reports/utils_analyse_abc.R")

theme_set(theme_bw())


########################
my_bin = 20

input_variables =  c("nest_quality_assessment_error", "percentage_foragers","number_nests_dbl","exploring_phase_dbl") 
prior_bounds = list(nest_quality_assessment_error=c(0,50), percentage_foragers=c(3,100),number_nests_dbl=c(2,80),exploring_phase_dbl=c(1000,10000))
bins = list(nest_quality_assessment_error=my_bin, percentage_foragers=my_bin,number_nests_dbl=my_bin,exploring_phase_dbl=my_bin)
bandwidth = list(nest_quality_assessment_error=0.5, percentage_foragers=1.0, number_nests_dbl= 1, exploring_phase_dbl= 50)
sampleSize = 500



# pour calculer le percentage max thoerique
data_nets = read_csv("data/data_ants.csv")
data_nets$input_size
data_nets$amount
vect_percent_max = data_nets$amount / data_nets$input_size *100






###########################################
###     plot  marginales (result ABC)   ###
###########################################

withr::with_libpaths(new = "~/.local/lib/R/site-library", devtools::install_github("steveharoz/patchwork"))
withr::with_libpaths(new = "~/.local/lib/R/site-library", library(patchwork))

liste_plots <- NULL

for (i in 0:18) {
  
  
  ################################
  ###      1   DATA           ###
  ################################
  
  #i = 0   
  colonie = paste0("posteriorSample_",i)
  
  # step variable selon la colonie
  number_steps = length(list.files(paste0("output/ResultsABC_5params/",colonie)))
  datafile =  paste0("output/ResultsABC_5params/",colonie,"/step",number_steps*20,".csv")
  
  
  
  prior <- as_tibble(prior_bounds)
  dimension <- length(input_variables)
  bins <- as_tibble(bins)
  if(is.null(bandwidth)) {
    h <- NULL
  } else {
    h <- as_tibble(bandwidth)
  }
  
  ### Data
  post <- read_csv(datafile) 
  sampleSize <- if (is.null(sampleSize)) nrow(post) else sampleSize
  
  post <- post %>% 
    sample_n(sampleSize) %>% 
    select_at(c(input_variables, "weight")) %>%
    filter(weight > 0) %>%
    mutate(weight = weight / sum(weight))
  
  thetas <- post %>% select_at(input_variables)
  weights <- post %>% select(weight)
  
  
  
  ### marginals 1D (4 params)
  fhat_marginal <- 
    map(names(thetas), 
        function(v) {
          f <- kde_bounded_2(select(thetas, v), weights, 
                             bounds = select(prior, v), 
                             bins = select(bins, v), 
                             h=select(h, v))
          f$result
        }
    )
  
  
  
  liste_plots[[5*i+1]] <- create_plot_marginal_v3(fhat=fhat_marginal[[1]], vname= names(thetas)[1], i=i, j=1, vect_percent_max[i+1])
  liste_plots[[5*i+2]] <- create_plot_marginal_v3(fhat=fhat_marginal[[2]], vname= names(thetas)[2], i=i, j=2, vect_percent_max[i+1])
  liste_plots[[5*i+3]] <- create_plot_marginal_v3(fhat=fhat_marginal[[3]], vname= names(thetas)[3], i=i, j=3, vect_percent_max[i+1])
  liste_plots[[5*i+4]] <- create_plot_marginal_v3(fhat=fhat_marginal[[4]], vname= names(thetas)[4], i=i, j=4, vect_percent_max[i+1])
  
  liste_plots[[5*i+5]] <- grid::textGrob('')
   
  
  print(i)
  
}

pp = reduce(liste_plots, function(a,b){a+b}) + plot_layout(ncol = 5)
ggsave("output/figures/marginals_4params_v4_1.png", pp, width = 20, height = 45, dpi = 600)

