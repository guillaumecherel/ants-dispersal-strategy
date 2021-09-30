# forumis:fonctions utiles pour tracer courbes de résultats abc

theme_set(theme_bw())


# Marginales
plot_marginal <- function(fhat, vname) {
  ggplot(fhat) + geom_line(aes_string(x=vname, y="density")) +
    ylim(0,NA)
}


save_marginal <- function(fhat, vname, dir) {
  pp <- ggplot(fhat) + geom_line(aes_string(x=vname, y="density")) +
    ylim(0,NA)
  ggsave( paste0(dir,"/marginale_",vname,".png") ,pp, width = 15, height = 15, dpi = 400)
}



create_plot_marginal <- function(fhat, vname) {
  #input_variables =  c("nest_quality_assessment_error", "percentage_foragers","number_nests_dbl","exploring_phase_dbl") 
  
  size_x_text = 22
  size_x_title = 22
  
  if(vname == "nest_quality_assessment_error") {
    borne_inf_x = 0
    borne_sup_x = 50
    borne_sup_y = 0.08 #0.14 
    size_x_title = 18
    } else if (vname == "percentage_foragers") {
      borne_inf_x = 3
      borne_sup_x = 100
      borne_sup_y = 0.05 #0.07 
      } else if (vname == "number_nests_dbl") {
        borne_inf_x = 15
        borne_sup_x = 80
        borne_sup_y = 0.05 #0.07 
        } else if (vname == "exploring_phase_dbl"){
          borne_inf_x = 1000
          borne_sup_x = 10000
          borne_sup_y = 0.0006 #0.0004 
          size_x_text = 15
          }
    
    
    pp <- ggplot(fhat) + geom_line(aes_string(x=vname, y="density")) +
      ylim(0,borne_sup_y) +
      xlim(borne_inf_x,borne_sup_x)+
      theme( axis.text.y = element_text(size=20), 
             axis.text.x = element_text(size= size_x_text),
             axis.title.y = element_text(size=0), 
             axis.title.x = element_text(size=0) )
  pp
}





create_plot_marginal_v3 <- function(fhat,vname,i,j, percent_max) {
  #input_variables =  c("nest_quality_assessment_error", "percentage_foragers","number_nests_dbl","exploring_phase_dbl") 
  
  size_x_text = 22
  size_x_title = 22
  
  if(vname == "nest_quality_assessment_error") {
    borne_inf_x = 0
    borne_sup_x = 50
    borne_sup_y = 0.06 #0.14 
    size_x_title = 18
    vname2 = "nest_quality_assessment_error"
  } else if (vname == "percentage_foragers") {
    borne_inf_x = 3
    borne_sup_x = percent_max
    borne_sup_y = 0.1 #0.035 #0.07 
    vname2 = "percentage_foragers"
  } else if (vname == "number_nests_dbl") {
    borne_inf_x = 2  ###  avant 15
    borne_sup_x = 80
    borne_sup_y = 0.05 #0.07 
    vname2 = "number_nests"
  } else if (vname == "exploring_phase_dbl"){
    borne_inf_x = 1000
    borne_sup_x = 10000
    borne_sup_y = 0.002 #0.0003   à modifier pour ne pas tronquer en hauteur
    #size_x_text = 15
    vname2 = "exploring_phase"
  }

  size_x = if(i==18){size_x_title}else{0}
  size_y = if(j==1){20}else{0}
  
  pp <- ggplot(fhat) + geom_line(aes_string(x=vname, y="density")) +
    ylim(0,borne_sup_y) +
    xlim(borne_inf_x,borne_sup_x)+
    labs(x = vname2) +
    theme( axis.text.y = element_text(size=20), 
           axis.text.x = element_text(size= size_x_text),
           axis.title.y = element_text(size = size_y), 
           axis.title.x = element_text(size = size_x) )
  pp
}






# Marginales 2D
plot_marginal_2 <- function(fhat, v1, v2) {
  ggplot(fhat) + geom_raster(aes_string(x=v1, y=v2, fill="density")) +
    scale_fill_viridis_c(limits = c(0,NA))
}


save_marginal_2 <- function(fhat, v1, v2, dir) {
  pp <- ggplot(fhat) + geom_raster(aes_string(x=v1, y=v2, fill="density")) +
    scale_fill_viridis_c(limits = c(0,NA))
  ggsave( paste0(dir,"/marginale2_",v1,"_",v2,".png") ,pp, width = 15, height = 15, dpi = 400)
}


create_plot_marginal_2 <- function(fhat, v1, v2) {
  pp <- ggplot(fhat) + geom_raster(aes_string(x=v1, y=v2, fill="density")) +
    xlim(0,50) +
    ylim(0,50) +
    scale_fill_viridis_c(limits = c(0,NA))
  pp
}


# Is the full joint posterior distribution flat or peaked?
vgtau <- function(fhat, tau) {
  fhat$result %>% 
    filter(density >= tau) %>% 
    summarize(sum(vcell)) %>% 
    deframe()
}

pgtau <- function(fhat, tau) {
  fhat$result %>% 
    filter(density >= tau) %>% 
    summarize(sum(vcell * density)) %>% 
    deframe()
}


# Where are the peaks?

cgtau <- function(fhat, tau) {
  if(tau <= 0) {
    c <- rep(1, nrow(fhat$result))
  } else {
    nodes <- (1:nrow(fhat$result))[fhat$result$density >= tau]
    g <- induced_subgraph(make_lattice(sapply(fhat$grid_borders, length)), nodes)
    c <- components(g)$membership
  }
  fhat$result %>% 
    filter(density >=tau) %>% 
    group_by(cluster = c) %>%
    filter(density == max(density)) %>%
    ungroup() %>%
    select_at(c(input_variables, "density")) %>%
    mutate(tau=tau)
}









