# tracer et eneregistrer les résultats de ABC pour chaque colonie et chaque ABC

# COMPARER L'EFFICACITÉ ENTRE LES DIFFÉRENTS ABCs (appelé expe):

################################
###         PACKAGES         ###
################################
library(ks)
library(igraph)
library(mvtnorm)
library(tidyverse)
library(patchwork)

library(stringr)
library(reshape2)

### import functions 
# Kernel density estimate
source("kde.R")
source("utils_analyse_abc.R")



################################
###         1                ###
################################
# pour créer un dossier par expe, puis un sous dossier par colonie, et enregistrer
# la densité 2D, les marginales 1D, et des indicateurs sur les peaks
# + sanity_check (vecteur) pour vérifier que le volume sous la densité estimée à partir des points en 2D est proche de 1

liste_expes = list.files("data")
#sanity_check = matrix(, ncol = 3, nrow = length(liste_expes)*19)
sanity_check = c()

### boucle sur les fichiers d'une même expé
i = 0
j = 1
#j = 1   #j : expe
for (j in 1: length(liste_expes)) {
  expe = liste_expes[j]
  
  for (i in 0:18) {
  
  ################################
  ###           DATA           ###
  ################################
  # expérience:
  #expe = "ResultsABC_N=40_explo_phase=2500"
  #expe = "ResultsABC_N=20_explo_phase=2500"
  #expe = "ResultsABC_N=60_explo_phase=2500"
  #expe = "ResultsABC_N=40_explo_phase=1000"
  #expe = "ResultsABC_N=20_explo_phase=1000"
  #expe = "ResultsABC_N=60_explo_phase=1000"
  #expe = "ResultsABC_N=40_explo_phase=4000"
  #expe = "ResultsABC_N=20_explo_phase=4000"
  #expe = "ResultsABC_N=60_explo_phase=4000"
  #expe = "ResultsABC_N=40_explo_phase=6000"
  #expe = "ResultsABC_N=40_explo_phase=8000"
    
  colonie = paste0("posteriorSample_",i)
  
  # step variable selon la colonie
  number_steps = length(list.files(paste0("data/",expe,"/",colonie)))
  datafile =  paste0("data/",expe,"/",colonie,"/step",number_steps*20,".csv")
  
  
  
  input_variables =  c("nest_quality_assessment_error", "percentage_foragers") 
  prior_bounds = list(nest_quality_assessment_error=c(0,50), percentage_foragers=c(3,100))
  bins = list(nest_quality_assessment_error=20, percentage_foragers=20)
  bandwidth = list(nest_quality_assessment_error=0.5, percentage_foragers=0.5)
  sampleSize = 500
  
  
  
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
  
  
  
  ### KDE
  fhat <- kde_bounded_2(thetas, weights, bounds=prior, bins = bins,h = h)
  
  
  
  dim(fhat[[1]])
  
  
  # Sanity check of the density estimates:
  (fhat$result %>% summarize(sum(density*vcell)) )[[1]]
  
  sanity_check = rbind(sanity_check, c(expe,i,(fhat$result %>% summarize(sum(density*vcell)) )[[1]]) )
  
  
  ### results
  # dir to save results:
  dir.create("results")
  dir_results_parent = paste0("results/plots_",expe)
  dir.create(dir_results_parent)
  dir_results = paste0(dir_results_parent,"/results_colonie_",i)
  dir.create(dir_results)
  
  ### Marginales
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
  
  
  # tracer
  # pmap(list(fhat=fhat_marginal, vname=names(thetas)),
  #      plot_marginal)
  
  # que enregistrer
  pmap(list(fhat=fhat_marginal, vname=names(thetas), dir = dir_results),
       save_marginal)
  
  
  
  
  ### Marginales 2D
  variable_pairs <- t(combn(names(thetas), 2))
  marginals <- tibble(v1=variable_pairs[,1], v2=variable_pairs[,2])
  
  fhat_marginal_2 <- pmap(marginals, function(v1, v2) {
    f <- kde_bounded_2(thetas[c(v1, v2)], weights, bins=select(bins, v1, v2), 
                       bounds=select(prior, v1, v2), h = select(h, v1, v2))
    f$result
  })
  
  
  #dim(fhat_marginal_2[[1]])
  
  # tracer
  # pmap(list(fhat=fhat_marginal_2, v1=marginals$v1, v2=marginals$v2),
  #      plot_marginal_2)
  
  # save
  pmap(list(fhat=fhat_marginal_2, v1=marginals$v1, v2=marginals$v2, dir = dir_results),
       save_marginal_2)
  
  
  
  
  ### Is the full joint posterior distribution flat or peaked?
  # Volume of the set GTau where the density is > tau, as a function of tau.
  taus <- seq(0, max(fhat$result$density), length.out=50)
  pp <- ggplot(data.frame(tau=taus, 
                    volume=map_dbl(taus, function(x) vgtau(fhat, x)),
                    cdf=map_dbl(taus, function(x) pgtau(fhat, x)))) +
    geom_line(aes(x=tau, y=volume)) +
    scale_x_reverse() +
    ylab("Volume")
  
  
  # save
  ggsave( paste0(dir_results,"/volume_density_over_tau",".png") ,pp, width = 15, height = 15, dpi = 400)
  
  
  
  
  ### Where are the peaks?
  # Disjoint components of the subset where the density is > tau, as a function of tau.
  taus <- seq(0, max(fhat$result$density), length.out=100)
  cgtaus <- bind_rows(lapply(taus, function(x) cgtau(fhat, x)))
  
  # Highest position of each peak.
  cgtaus %>% group_by_at(input_variables) %>% summarize(density=max(density)) %>% ungroup() %>% arrange(desc(density)) 
  
  ggplot(cgtaus %>% group_by(tau) %>% summarize(clusters = n())) +
    geom_line(aes(x=tau, y=clusters)) +
    scale_x_reverse()
  
  
  
  
  ### How close to the highest peak is the distribution?
  # Posterior probability of the parameters being within a sphere centered on the highest density point, as a function of its radius. 
  
  # most likely
  ml <- fhat$result[which.max(fhat$result$density),]
  
  selected_points <- fhat$result %>% filter(density > 0)
  
  # Distance of each point from the highest density
  sp <- selected_points %>% select(input_variables)
  repml <- matrix(rep(as.numeric(ml[,input_variables]), each=nrow(sp)), ncol=dimension)
  dhd <- sqrt(rowSums((sp - repml) ** 2))
  
  # Cumulative probability as a function of the radius
  ppoint <- selected_points$density * selected_points$vcell
  sortind <- order(dhd)
  prad<- cumsum(ppoint[sortind])
  
  
  pp <- ggplot(data.frame(radius = dhd[sortind], posterior = prad)) +
    geom_line(aes(x=radius, y=posterior))
  
  # save
  ggsave( paste0(dir_results,"/volume_distance_hughest_peak",".png") ,pp, width = 15, height = 15, dpi = 400)
  
  
  print(c(j,i))
  
  }
}  


# write sanity check
colnames(sanity_check) <- c("params","colony","sanity_check")
sanity_check = as.data.frame(sanity_check)
# write_csv(sanity_check, "sanity_check_ABCs.csv")
sanity_check$sanity_check = as.numeric(as.character(sanity_check$sanity_check))
min(sanity_check$sanity_check) ; max(sanity_check$sanity_check)












################################
###         2                ###
################################
# créer dataframe (df) avec tous les résultats des abcs (expe et colonies)
# créer aussi un csv pour resampler (toutes les colonies et les expes)
# rq: on renoarmalise les poids (somme à 1)


input_variables =  c("nest_quality_assessment_error", "percentage_foragers") 
prior_bounds = list(nest_quality_assessment_error=c(0,50), percentage_foragers=c(3,100))
bins = list(nest_quality_assessment_error=20, percentage_foragers=20)
bandwidth = list(nest_quality_assessment_error=0.5, percentage_foragers=0.5)

sampleSize = 500
liste_expes = list.files("data")

#df = c()
df = matrix(, ncol = 12, nrow = sampleSize*length(liste_expes)*19)
dim(df)
# les colonnes du df (8)
# epsilon, pAcc, t, ts,	rhos,	weight,	nest_quality_assessment_error, percentage_foragers
# ajouter number_nest, exploring_phase, le numéro de la colonie un id pour l'expérience (équivalent à nest_quality_assessment_error et percentage_foragers)
# pour resampler selon les lois a posteriori
REsampleSize = 800
resample = matrix(, ncol = 6, nrow = REsampleSize*length(liste_expes)*19)
# les colonnes du df (6)
# nest_quality_assessment_error, percentage_foragers, weight
# ajouter number_nest, exploring_phase, le numéro de la colonie


#j = 1   #j : expe
for (j in 1: length(liste_expes)) {
  expe = liste_expes[j]
  extract_numbers = expe %>% str_match_all("[0-9]+") %>% unlist %>% as.numeric
  number_nest = extract_numbers[1]
  exploring_phase = extract_numbers[2]
  
  #i = 0 #i: la colonie (dans l'expe j)
  for (i in 0:18) {
    colonie = paste0("posteriorSample_",i)
    # step variable selon la colonie
    number_steps = length(list.files(paste0("data/",expe,"/",colonie)))
    datafile =  paste0("data/",expe,"/",colonie,"/step",number_steps*20,".csv")
    
    post <- read_csv(datafile)
    
    post2 <- post %>% 
      sample_n(sampleSize) %>% 
      select_at(c(input_variables, "weight")) %>%
      filter(weight > 0) %>%
      mutate(weight = weight / sum(weight))
    
    
    temp_resample <- sample_n(post2, REsampleSize, replace=TRUE, weight=weight) %>% select(nest_quality_assessment_error, percentage_foragers, weight)

    post = post %>% filter(weight > 0) %>% mutate(weight = weight / sum(weight))
    post = as.matrix(post)
    df[((19*(j-1)+i)*sampleSize+1) : ((19*(j-1)+(i+1))*sampleSize),1:8] <- post
    #((19*(j-1)+i)*sampleSize+1) : ((19*(j-1)+(i+1))*sampleSize) 
    temp = cbind( rep(number_nest,sampleSize), rep(exploring_phase,sampleSize), rep(i,sampleSize), rep(j,sampleSize))
    # dans l'odre: number_nest, exploring_phase, le numéro de la colonie un id pour l'expérience (équivalent à nest_quality_assessment_error et percentage_foragers)
    df[((19*(j-1)+i)*sampleSize+1) : ((19*(j-1)+(i+1))*sampleSize),9:12] <- temp
    
    resample[((19*(j-1)+i)*REsampleSize+1) : ((19*(j-1)+(i+1))*REsampleSize),1:3] <- as.matrix(temp_resample)
    #resample[((19*(j-1)+i)*REsampleSize+1) : ((19*(j-1)+(i+1))*REsampleSize),4:6] <- temp[1:REsampleSize,1:3]
      
    resample[((19*(j-1)+i)*REsampleSize+1) : ((19*(j-1)+(i+1))*REsampleSize),4:6] <- cbind(rep(temp[1,1],REsampleSize), rep(temp[1,2],REsampleSize), rep(temp[1,3],REsampleSize) ) 
    
   }
}





df[dim(df)[1],]

post <- read_csv(datafile)
colnames(df) <- c( colnames(post), c("number_nest", "exploring_phase", "colony", "id") )
df= as.data.frame(df)
df = as_tibble(df)
colnames(resample) <- c( c("nest_quality_assessment_error", "percentage_foragers", "weight"), c("number_nests", "exploring_phase", "colony") )
resample= as.data.frame(resample)
resample = as_tibble(resample)



dim(df %>% group_by(id, colony))
# epsilon en fonction de l'experience (number nest et exploring phase) pour toutes les colonies
pp <- ggplot(df %>% group_by(id, colony), aes(x = id, y = epsilon), group = colony ) +
  geom_line( aes(colour = factor(colony)) ) #+
  #scale_color_identity()

pp


# faire le plot dans la densité 2D pour les différentes expe d'une même colony
# i.e facet_wrap( ~ nest_quality_assessment_error, labeller = label_both, ncol = 5)  sur number nest et exploring phase


# ajouter la taille de pop au df (prise dans le df de données), avec join
data_nets = read_csv("data2/data_ants.csv")
data_nets$colony <- 0:18
resample2 = inner_join(resample, data_nets %>% select(colony,input_size)) 
#write_csv(resample2, "resampleABCs_3.csv")



dim(resample2)
dim(resample2 %>% filter(percentage_foragers >= 70) )









################################
###         3                ###
################################
# créer un dataframe  avec tous les densités calculés à partir des résultats de ABC

# possibilité de charger le df (si on l'a enregistré)
# charger le df
#df_densities = read.csv("df_densities_all_colonies_all_exp.csv")
df_densities = read.csv("df_densities_all_colonies_all_exp_2.csv")
df_densities = df_densities[,-1]

# le même pour tous
input_variables =  c("nest_quality_assessment_error", "percentage_foragers") 
prior_bounds = list(nest_quality_assessment_error=c(0,50), percentage_foragers=c(3,100))
bins = list(nest_quality_assessment_error=20, percentage_foragers=20)
bandwidth = list(nest_quality_assessment_error=0.5, percentage_foragers=0.5)
sampleSize = 500


prior <- as_tibble(prior_bounds)
dimension <- length(input_variables)
bins <- as_tibble(bins)
if(is.null(bandwidth)) {
  h <- NULL
} else {
  h <- as_tibble(bandwidth)
}



liste_expes = list.files("data")

# 441 vient de 21*21 (bins)
taille_grille = 441
df_densities = matrix(, ncol = 8, nrow = taille_grille *length(liste_expes)*19)
dim(df_densities)
# les colonnes du df: celles de fhat$result (4)
# nest_quality_assessment_error percentage_foragers  density vcell
# ajouter (4) number_nest, exploring_phase, le numéro de la colonie un id pour l'expérience (équivalent à number_nest, exploring_phase)


# fin même pour tous


#j = 1   #j : expe
for (j in 1: length(liste_expes)) {
  expe = liste_expes[j]
  extract_numbers = expe %>% str_match_all("[0-9]+") %>% unlist %>% as.numeric
  number_nest = extract_numbers[1]
  exploring_phase = extract_numbers[2]
  
  #i = 0 #i: la colonie (dans l'expe j)
  for (i in 0:18) {
    colonie = paste0("posteriorSample_",i)
    # step variable selon la colonie
    number_steps = length(list.files(paste0("data/",expe,"/",colonie)))
    datafile =  paste0("data/",expe,"/",colonie,"/step",number_steps*20,".csv")
    
    post <- read_csv(datafile)
    
    post <- post %>% 
      sample_n(sampleSize) %>% 
      select_at(c(input_variables, "weight")) %>%
      filter(weight > 0) %>%
      mutate(weight = weight / sum(weight))
    
    thetas <- post %>% select_at(input_variables)
    weights <- post %>% select(weight)
    
    
    ### KDE
    fhat <- kde_bounded_2(thetas, weights, bounds=prior, bins = bins,h = h)
    
    fhat_res = as.matrix(fhat$result)
    df_densities[((19*(j-1)+i)*taille_grille+1) : ((19*(j-1)+(i+1))*taille_grille),1:4] <- fhat_res
    temp = cbind( rep(number_nest,taille_grille), rep(exploring_phase,taille_grille), rep(i,taille_grille), rep(j,taille_grille))
    # dans l'odre: number_nest, exploring_phase, le numéro de la colonie un id pour l'expérience (number_nest, exploring_phase)
    df_densities[((19*(j-1)+i)*taille_grille+1) : ((19*(j-1)+(i+1))*taille_grille),5:8] <- temp
    print(c(j,i))
  }
}



df_densities

# ajouter nom aux colonnes
colnames(df_densities) <- c( c("nest_quality_assessment_error", "percentage_foragers", "density", "vcell"), c("number_nest", "exploring_phase", "colony", "id") )
#c("nest_quality_assessment_error", "percentage_foragers", "density", "vcell")
# colnames(fhat[[1]])
df_densities = as.data.frame(df_densities)
df_densities = as_tibble(df_densities)


# enregistrer le df pour le charger ensuite
#write.csv(df_densities, "df_densities_all_colonies_all_exp_2.csv", row.names = T)


# filtré sur colonie
ggplot(df_densities %>% filter(colony == 0)) + geom_raster(aes_string(x="nest_quality_assessment_error" , y="percentage_foragers", fill="density")) +
  scale_fill_viridis_c(limits = c(0,NA)) + 
  facet_wrap( number_nest ~ exploring_phase, labeller = label_both)




(df_densities %>% filter(colony == 0, number_nest == 20, exploring_phase == 1000))$nest_quality_assessment_error
(df_densities %>% filter(colony == 0, number_nest == 40, exploring_phase == 1000))$nest_quality_assessment_error










################################
###         4                ###
################################

# calcul des erreur : on utilise le resample des paramètres selon la loi a posteriori (fait au point 2) )
# et déjà resimulé avec OpenMOLE
# on a simulé le nombre de fourmi dans chaque nid après fission, et on compare avec le terrain

#simu_ReSample = read_csv("data2/ResultsReSample_2.csv")
simu_ReSample = read_csv("data2/ResultsReSample_3.csv")
colnames(simu_ReSample)

# ajout d'un indice sur les param uniques:
simu_ReSample_Params = simu_ReSample %>% dplyr:: select(number_nests,exploring_phase,colony)
uniqueParam <- unique(simu_ReSample_Params)
dim(uniqueParam)
uniqueParam$paramID <-  1:nrow(uniqueParam)  # ajoute une colonne qui numerote les exps
uniqueParam <-  as.data.frame(uniqueParam)
#JOIN uniqueParam with df then split into human history records and zombies history records
simu_ReSample <- inner_join(simu_ReSample, uniqueParam)
colnames(simu_ReSample)

# changer nom colonne
colnames(simu_ReSample)[8:17] <- paste0("nest_",0:9)


simu_ReSample.m <- melt(simu_ReSample,
                           measure.vars = 8:17,
                           variable.name = "nest_number")



data_nets = read_csv("data2/data_ants.csv")
data_nets$colony <- 0:18
colnames(data_nets)[4:13] <- paste0("nest_",0:9)
data_nets.m <- melt(data_nets,
                        measure.vars = 4:13,
                        variable.name = "nest_number")
colnames(data_nets.m)[5] <- "value_data"

###### ATTENTION: number_nest: dans les données, le nombre de nif effectif après fission, 
### alors que number_nests: nombre de nids potentiels dans simu

# join data and simu
df_data_simu = inner_join(simu_ReSample.m,data_nets.m)
df_data_simu = as_tibble(df_data_simu)
colnames(df_data_simu)

# tracer barplot
temp_df = df_data_simu %>%  filter( paramID == 1)

pp <- ggplot(temp_df, aes(x = nest_number, y= value) ) +
  #geom_point( size = 0.5) +
  geom_boxplot() +
  geom_point( aes(x = nest_number, y=value_data), colour = "red")
pp


temp_df %>%  group_by(nest_number) %>% summarise( mean = mean(value))


# avec facet ?
pp <- ggplot(df_data_simu %>%  filter( number_nests == 40, exploring_phase == 2500), aes(x = nest_number, y= value) ) +
  geom_boxplot() +
  geom_point( aes(x = nest_number, y=value_data), colour = "red") +
  facet_wrap( ~colony, labeller = label_both) 
pp

# ggsave("erreur_all_colonies_N=40_explo=2500.png" ,pp, width = 15, height = 15, dpi = 400)











################################
###         5                ###
################################

##########################################
##########################################
### à paramètres fixés, comparer les colonies
### enregistredes graphe différents: marginales 2D, 1D et erreur sur une même ligne pour une colonie
## autant de lignes de que de colonies (pour une expe: valeur de nest_numbers et explo phase fixées)

#install.packages("devtools")
#devtools::install_github("thomasp85/patchwork")
devtools::install_github("steveharoz/patchwork")
library(patchwork)

liste_expes = list.files("data")
### boucle sur les fichiers d'une même expé
i = 0
j = 1
#j = 1   #j : expe

for (j in 1:length(liste_expes)){
  liste_plots = NULL
  expe = liste_expes[j]

  ### pour erreur 
  #expe = liste_expes[j]
  extract_numbers = expe %>% str_match_all("[0-9]+") %>% unlist %>% as.numeric
  number_nests_1 = extract_numbers[1]
  exploring_phase_1 = extract_numbers[2]
  
  
  for (i in 0:18) {
    
    ################################
    ###           DATA           ###
    ################################
    # expérience:
    colonie = paste0("posteriorSample_",i)
    
    # step variable selon la colonie
    number_steps = length(list.files(paste0("data/",expe,"/",colonie)))
    datafile =  paste0("data/",expe,"/",colonie,"/step",number_steps*20,".csv")
    
    input_variables =  c("nest_quality_assessment_error", "percentage_foragers") 
    prior_bounds = list(nest_quality_assessment_error=c(0,50), percentage_foragers=c(3,100))
    bins = list(nest_quality_assessment_error=20, percentage_foragers=20)
    bandwidth = list(nest_quality_assessment_error=0.5, percentage_foragers=0.5)
    sampleSize = 500
    
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
    
    
    ### KDE
    #fhat <- kde_bounded_2(thetas, weights, bounds=prior, bins = bins,h = h)
    
    
    ### Marginales 2D
    variable_pairs <- t(combn(names(thetas), 2))
    marginals <- tibble(v1=variable_pairs[,1], v2=variable_pairs[,2])
    
    fhat_marginal_2 <- pmap(marginals, function(v1, v2) {
      f <- kde_bounded_2(thetas[c(v1, v2)], weights, bins=select(bins, v1, v2), 
                         bounds=select(prior, v1, v2), h = select(h, v1, v2))
      f$result
    })
    
    
    liste_plots[[4*i+1]] = create_plot_marginal_2(fhat=fhat_marginal_2[[1]], v1=marginals$v1, v2=marginals$v2)
    

    ### Marginales
    ### Marginale 1
    kk=1
    v= names(thetas)[kk]
    fhat_marginal <- kde_bounded_2(select(thetas, v), weights, 
                               bounds = select(prior, v), 
                               bins = select(bins, v), 
                               h=select(h, v))$result

    
    liste_plots[[4*i+2]] = create_plot_marginal(fhat=fhat_marginal, vname=v)

    ### Marginale 2
    kk=2
    v= names(thetas)[kk]
    
    fhat_marginal <- kde_bounded_2(select(thetas, v), weights, 
                       bounds = select(prior, v), 
                       bins = select(bins, v), 
                       h=select(h, v))$result
    
    
    liste_plots[[4*i+3]] = create_plot_marginal(fhat=fhat_marginal, vname=v)
    
    ### erreur 
  
    temp_df = df_data_simu %>%  filter( exploring_phase == exploring_phase_1, number_nests == number_nests_1, colony == i )
    
    pp <- ggplot(temp_df, aes(x = nest_number, y= value) ) +
      #geom_point( size = 0.5) +
      geom_boxplot() +
      geom_point( aes(x = nest_number, y=value_data), colour = "red")
    #pp
    
    liste_plots[[4*i+4]] = pp 
      
    print(i)
    
  }  
  pp = reduce(liste_plots, function(a,b){a+b}) + plot_layout(ncol = 4)
  
  #ggsave("marginals_all_colonies_2.png" ,pp, width = 15, height = 45, dpi = 400)
  ggsave(paste0("marginals_all_colonies_3_N=",number_nests_1,"_explo_phase",exploring_phase_1,".png") ,pp, width = 15, height = 45, dpi = 400)
  
  
}







################################
###         6                ###
################################

##########################################
# comparer les paramètres (expe) pour une même colonie
# enregistrer les graphes de l'erreur et les densité 2D


dir_res = "result_compar_param"
dir.create(dir_res)

for (i in 0:18){  # i : colony
  
  sous_dir_res = paste0("colonie_",i)
  sous_dir_res_full = paste0(dir_res,"/",sous_dir_res)
  dir.create(sous_dir_res_full)
  
  # densité (2D)
  pp <- ggplot(df_densities %>% filter(colony == i)) + geom_raster(aes_string(x="nest_quality_assessment_error" , y="percentage_foragers", fill="density")) +
    scale_fill_viridis_c(limits = c(0,NA)) + 
    xlim(0,50) + ylim(0,50) + 
    facet_wrap( number_nest ~ exploring_phase, labeller = label_both)
  
  # save
  ggsave( paste0(sous_dir_res_full,"/compar_densities_colonie_",i,".png") ,pp, width = 15, height = 15, dpi = 400)
  
  
  ### erreur 
  pp <- ggplot(df_data_simu %>% filter(colony == i), aes(x = nest_number, y= value) ) +
    #geom_point( size = 0.5) +
    geom_boxplot() +
    geom_point( aes(x = nest_number, y=value_data), colour = "red") +
    facet_wrap( number_nests ~ exploring_phase, labeller = label_both)

  # save
  ggsave( paste0(sous_dir_res_full,"/compar_erreur_colonie_",i,".png") ,pp, width = 15, height = 15, dpi = 400)
  print(i)  
}










################################
###         7                ###
################################

##########################################
# comparer les paramètres pour une même colonie
# construire une matrice avec un rho "moyen" (tient compte des poids)
# avec en ligne les colonies (19) et en colonnes les expériences (11)
# utilise df (créé au point 2: concaténation des résultat brut des ABCs) (rq poids déjà renormalisés à somme =1)
# la colonne id de df renvoit au choix de l'expe (hors colonies)

colnames(df)
temp_df = df %>% filter(colony == 2, id == 3)
max(temp_df$rhos)
temp_df$epsilon[1]
temp_df$weight
# epsilon c'est le max des rho



mat_compare_rho = matrix(, nrow = 19, ncol = 11)
for (i in 1:19){     # i: colonie
  for (j in 1:11){    # j: expe
    temp_df = df %>% filter(colony == (i-1), id == j)
    mean(temp_df$rhos)
    rho_agg = sum(temp_df$rhos * temp_df$weight)
    mat_compare_rho[i,j] = rho_agg
  }
}


mat_compare_rho

# par ligne, renvoyer le numéro de la colonne qui fait le rho minimal 
ind_min_rho = apply(mat_compare_rho, 1, function(x){which.min(x)}) 

table(ind_min_rho)
liste_expes

# pour résumer et n'avoir qu'un indicateur (comparer les colonies), on divise par la taille de la colonie ? 
data_nets = read_csv("data2/data_ants.csv")
data_nets$input_size

res_rho_div = apply(mat_compare_rho, 2, function(x){x /data_nets$input_size})


mat_heatmap = melt(mat_compare_rho)
ggplot(data = mat_heatmap, aes(y=Var1, x=Var2, fill=value)) + 
  geom_tile()



# somme des erreur par expe (après avoir renormalisé par taille depart)
colSums(res_rho_div)







#####################

colnames(df_data_simu)
unique(df_data_simu$nest_number)


df_data_simu
df_data_simu$colony
# renormalisation par la taille
df_data_simu_err = df_data_simu %>% group_by(number_nests,exploring_phase) %>% mutate(erreur_abs = abs(value_data-  value) )  #/ input_size
# sans renormalisation par la taille
# df_data_simu_err = df_data_simu %>% group_by(number_nests,exploring_phase, nest_number) %>% mutate(erreur_abs = abs(value_data-  value) )
# colnames(df_data_simu_err)

df_data_simu_err$erreur_abs

df_data_simu_err %>% summarise( erreur_moyenne = sum(erreur_abs)/(200*19*10) )

view(df_data_simu[51,])









