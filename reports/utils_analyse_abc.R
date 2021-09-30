# forumis:fonctions utiles pour tracer courbes de résultats abc


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
  pp <- ggplot(fhat) + geom_line(aes_string(x=vname, y="density")) +
    ylim(0,NA) +
    xlim(0,50)
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









