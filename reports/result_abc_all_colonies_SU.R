# tracer et eneregistrer les résultats de ABC pour chaque colonie et chaque ABC

################################
###         PACKAGES         ###
################################
library(ks)
library(igraph)
library(mvtnorm)
library(tidyverse)
library(patchwork)
library(grid)

library(stringr)
library(reshape2)

### import functions 
# Kernel density estimate
source("kde.R")
source("utils_analyse_abc.R")

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

devtools::install_github("steveharoz/patchwork")
library(patchwork)

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
















################################
###      2    RESAMPLER      ###
################################

df = matrix(, ncol = 11, nrow = sampleSize*19)
dim(df)
REsampleSize = 2000
resample = matrix(, ncol = 6, nrow = REsampleSize*19)



for (i in 0:18) {
  
  #i = 0   
  colonie = paste0("posteriorSample_",i)
  
  number_steps = length(list.files(paste0("output/ResultsABC_5params/",colonie)))
  datafile =  paste0("output/ResultsABC_5params/",colonie,"/step",number_steps*20,".csv")
  
  post <- read_csv(datafile)
  
  post2 <- post %>% 
    sample_n(sampleSize) %>% 
    select_at(c(input_variables, "weight")) %>%
    filter(weight > 0) %>%
    mutate(weight = weight / sum(weight))
  
  
  temp_resample <- sample_n(post2, REsampleSize, replace=TRUE, weight=weight) %>% select(nest_quality_assessment_error, percentage_foragers, number_nests_dbl, exploring_phase_dbl, weight)
  
  post = post %>% filter(weight > 0) %>% mutate(weight = weight / sum(weight))
  post = as.matrix(post)
  df[ ((i*sampleSize)+1) : ((i+1)*sampleSize),1:10] <- post
  temp = rep(i,sampleSize)
  
  df[ ((i*sampleSize)+1) : ((i+1)*sampleSize), 11] <- temp
  
  resample[ ((i*REsampleSize)+1) : ((i+1)*REsampleSize), 1:5] <- as.matrix(temp_resample)
  
  resample[ ((i*REsampleSize)+1) : ((i+1)*REsampleSize), 6] <- rep(i,REsampleSize)
  
}




post <- read_csv(datafile)
colnames(df) <- c( colnames(post), "colony")
df= as.data.frame(df)
df = as_tibble(df)
colnames(resample) <- c( "nest_quality_assessment_error", "percentage_foragers", "number_nests_dbl", "exploring_phase_dbl", "weight", "colony")
resample= as.data.frame(resample)
resample = as_tibble(resample)




data_nets = read_csv("data/data_ants.csv")
data_nets$colony <- 0:18
resample2 = inner_join(resample, data_nets %>% select(colony,input_size, amount)) 
write_csv(resample2, "output/resampleABCs_4_params_v3.csv")


















###################################
###     4   calcul des erreurs   ###
###################################


simu_ReSample = read_csv("output/ResultsReSample_4_params_v3.csv")
colnames(simu_ReSample)

colnames(simu_ReSample)[8:17] <- paste0("nest_",0:9)


simu_ReSample.m <- melt(simu_ReSample,
                        measure.vars = 8:17,
                        variable.name = "nest_number")



data_nets = read_csv("data/data_ants_NA.csv")
data_nets$colony <- 0:18
colnames(data_nets)[5:14] <- paste0("nest_",0:9)
data_nets.m <- melt(data_nets,
                    measure.vars = 5:14,
                    variable.name = "nest_number")
colnames(data_nets.m)[6] <- "value_data"


# join data and simu
df_data_simu = inner_join(simu_ReSample.m,data_nets.m)
df_data_simu = as_tibble(df_data_simu)
colnames(df_data_simu)





# tracer barplot
temp_df = df_data_simu %>%  filter( colony == 1)

pp <- ggplot(temp_df, aes(x = nest_number, y= value) ) +
  #geom_point( size = 0.5) +
  geom_boxplot() +
  geom_point( aes(x = nest_number, y=value_data), colour = "red")
pp




# avec facet 
pp <- ggplot(df_data_simu, aes(x = nest_number, y= value) ) +
  geom_boxplot() +
  geom_point( aes(x = nest_number, y=value_data), colour = "red") +
  facet_wrap( ~colony, labeller = label_both) 
pp

 ggsave("erreur_all_colonies_4_params_v3_1.png" ,pp, width = 15, height = 15, dpi = 400)








# #########################################
# #########################################
# ###   juste les erreurs
# 
# library(rlist)
# library(patchwork)
# liste_plots <- NULL
# 
# 
# # for (i in 0:18) {
# for (i in 0:4) {
#   for (j in 0:3){
#   #i=0
#   
#   temp_df = df_data_simu %>%  filter( colony == i+j )
#   temp_df = temp_df %>% mutate(nest_number =  str_replace(nest_number, "nest_", "") ) 
#   size_x_title = if(i %%5==0  ){18}else{0}
#   size_y_title = if( i %in% c(18,20,21,22) ){20}else{0}
#   
#   
#   theme_set(theme_bw())  ## added TM
#   
#   pp <- ggplot(temp_df, aes(x = nest_number, y= value) ) +
#     geom_boxplot() +
#     geom_point( aes(x = nest_number, y=value_data), colour = "red")  +
#     scale_x_discrete(breaks= paste0(0:9),
#                      labels= c("1","","","4","","","7","","","10"))+
#     theme( axis.text.y = element_text(size=18), 
#            axis.text.x = element_text(size= 18),
#            axis.title.y = element_text(size=size_x_title),
#            axis.title.x = element_text(size=size_y_title) )
#   list.append(liste_plots,pp + ggtitle(''))
#   
#   if (i %in% c(4,9,14,19)){  list.append(liste_plots, grid::textGrob(' ')) }else{}
#   
#   print(i)
#   }
# }
# 
# 
# pp = reduce(liste_plots, function(a,b){a+b}) + plot_layout(ncol = 4)
# #ggsave("erreur_all_colonies_4_params_v3_3.png" ,pp, width = 15, height = 15, dpi = 400)










#########################################
#########################################
###   juste les erreurs    2
 
library(patchwork)
liste_plots <- NULL

for (i in 0:18) {
  #i=0
  ### erreur 

  grob = grobTree(textGrob(i + 1, x = 0.8,  y = 0.8, hjust=0, gp=gpar(fontsize=20, fontface="bold")))
  
  temp_df = df_data_simu %>%  filter( colony == i )
  temp_df = temp_df %>% mutate(nest_number =  str_replace(nest_number, "nest_", "") ) 
  #unique(temp_df$nest_number2)
  
  size_x_title = if(i %%4==0  ){18}else{0}
  size_y_title = if( i %in% c(15,16,17,18) ){20}else{0}
  
  pp <- ggplot(temp_df, aes(x = nest_number, y = value) ) +
    #geom_point( size = 0.5) +
    geom_boxplot() +
    geom_point( aes(x = nest_number, y = value_data), colour = "red", size = 3)  +
    #theme( axis.text.x=element_blank() )
    scale_x_discrete(breaks= paste0(0:9),
                     labels= c("1","","","4","","","7","","","10")) +
    xlab("New colonies") +
    ylab("Resources") +
    theme( axis.text.y = element_text(size=18), 
           axis.text.x = element_text(size= 18),
           axis.title.y = element_text(size=size_x_title),
           axis.title.x = element_text(size=size_y_title) ) +
    annotation_custom(grob)
    # theme(axis.title.y = element_text(), axis.title.y.left = element_text(angle = 90)) +
    # labs(title=i+1) #, face = bold, size = 30, hjust = "right", vjust = "top")
    # annotate(geom="text", x = 10, y = 250, label=i + 1, face = "bold")
    
  
  # axis.title.y = element_text(size=20), 
  # axis.title.x = element_text(size=18) )
  #pp
  #  liste_plots[[5*i+5]] <- grid::textGrob(' ')
  liste_plots[[i+1]] = pp 
  
  print(i)
  
}


pp = reduce(liste_plots, function(a,b){a+b}) + plot_layout(ncol = 4)
ggsave("erreur_all_colonies_4_params_v3_2.png" ,pp, width = 15, height = 15, dpi = 400)












