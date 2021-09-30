library(tidyverse)
library(reshape2)

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

ggsave("output/figures/erreur_all_colonies_4_params_v3_1.png" ,pp, width = 15, height = 15, dpi = 400)





# #########################################
# #########################################
# ###   juste les erreurs
# 
# library(rlist)
# withr::with_libpaths(new = "~/.local/lib/R/site-library", devtools::install_github("steveharoz/patchwork"))
# withr::with_libpaths(new = "~/.local/lib/R/site-library", library(patchwork))
# 
# plots <- ggplot()
# 
# # for (i in 0:18) {
# for (i in 0:4) {
#   for (j in 0:3){
#       #i=0
# 
#       temp_df = df_data_simu %>%  filter( colony == i+j )
#       temp_df = temp_df %>% mutate(nest_number =  str_replace(nest_number, "nest_", "") ) 
#       size_x_title = if(i %%5==0  ){18}else{0}
#       size_y_title = if( i %in% c(18,20,21,22) ){20}else{0}
# 
#       theme_set(theme_bw())  ## added TM
# 
#       pp <- ggplot(temp_df, aes(x = nest_number, y= value) ) +
#         geom_boxplot() +
#         geom_point( aes(x = nest_number, y=value_data), colour = "red")  +
#         scale_x_discrete(breaks= paste0(0:9),
#                          labels= c("1","","","4","","","7","","","10"))+
#         theme( axis.text.y = element_text(size=18), 
#                axis.text.x = element_text(size= 18),
#                axis.title.y = element_text(size=size_x_title),
#                axis.title.x = element_text(size=size_y_title) ) +
#         ggtitle("")
# 
# 
#       plots <- plots + pp
# 
#       if (i %in% c(4,9,14,19)){  
#           plots <- list.append(liste_plots, grid::textGrob(' ')) 
#       } else {}
#   }
# }
# 
# 
# pp <- plots + plot_layout(ncol = 4) 
# ggsave("output/figures/erreur_all_colonies_4_params_v3_3.png" ,pp, width = 15, height = 15, dpi = 400)



