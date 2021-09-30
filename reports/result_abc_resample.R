library(tidyverse)

sampleSize = 500
input_variables =  c("nest_quality_assessment_error", "percentage_foragers","number_nests_dbl","exploring_phase_dbl") 




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
















