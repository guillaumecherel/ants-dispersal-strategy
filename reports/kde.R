library(tidyverse)

kde <- function(thetas, weights, bins) {
  if(ncol(thetas) == 1) {
    H <- hpi(thetas[[1]], deriv.order = 4)
    Hmax <- H
  } else {
    H <- Hpi(x=thetas, deriv.order = 3)
    Hmax <- max(diag(H))
  }
  
  supp <- 3.7
  grid_min <- gather(summarize_all(thetas, min))$value - Hmax*supp
  grid_max <- gather(summarize_all(thetas, max))$value + Hmax*supp
  grid_def <- tibble(from=grid_min,
                     to=grid_max,
                     by=(to - from) / bins)
  
  grid_borders <- pmap(grid_def, seq)

  points <- as_tibble(expand.grid(grid_borders))
  names(points) <- names(thetas)
  
  # La longeur de la cellule j pour la dimension d est donnée par 
  # min(grid_def$by[d] / 2, points[j,d] - grid_def$from[d]) + 
  #   min(grid_def$by[d] / 2, grid_def$to[d] - points[j,d]) 
  half_grid_def <- grid_def$by / 2
  limit_below <- points - grid_def$from
  limit_above <- grid_def$to - points
  cond_below <- (half_grid_def < limit_below)
  cond_above <- (half_grid_def < limit_above)
  length_below <- (half_grid_def) * cond_below + limit_below * (!cond_below)
  length_above <- (half_grid_def) * cond_above + limit_above * (!cond_above)
  dimcell <- length_below + length_above
  vcell <- reduce(dimcell, `*`)
  
  if(ncol(thetas) == 1) {
    fhat <- ks::kde(x=deframe(thetas), h=H, supp=supp, w=deframe(weights) * nrow(thetas), eval.points=deframe(points))
  } else {
    fhat <- ks::kde(x=thetas, H=H, supp=supp, w=deframe(weights) * nrow(thetas), eval.points=points)
  }

  list(result=points %>% mutate(density=fhat$estimate, vcell=vcell), grid_def=grid_def, grid_borders=grid_borders, dim=length(grid_def))
}

# kde_bounded <- function(thetas, weights, bounds, bins) {
#   # For bounded support, transform the space through the logit function, see
#   # http://thirdorderscientist.org/homoclinic-orbit/2013/10/24/kernel-density-estimation-for-random-variables-with-bounded-support-mdash-the-transformation-trick
#   C <- rep(1, ncol(thetas))
#   scale_unit <- function(x) {
#     (x - bounds$min) / (bounds$max - bounds$min)
#   }
#   
#   transform <- function(x) {
#     u <- scale_unit(x)
#     #log(u / (1 - u))
#     C * log(u / (1 - u))
#   }
# 
#   inverse <- function(x) {
#     (bounds$max * exp(x/C) + bounds$min) / (1 + exp(x/C))
#   }
# 
#   if(ncol(thetas) == 1) {
#     H <- hpi(transform(thetas[[1]]), deriv.order = 4)
#     Hmax <- H
#   } else {
#     H <- Hpi(transform(thetas), deriv.order = 3)
#     Hmax <- max(diag(H))
#   }
#   
#   supp <- 3.7
#   grid_min_t <- transform(gather(summarize_all(thetas, min))$value) - Hmax*supp
#   grid_max_t <- transform(gather(summarize_all(thetas, max))$value) + Hmax*supp
#   
#   grid_min <- inverse(grid_min_t)
#   grid_max <- inverse(grid_max_t)
#   
#   #grid_def <- tibble(from=bounds$min + margin,
#   #                   to=bounds$max - margin,
#   #                   by=(to - from) / bins)
#   grid_def <- tibble(from=grid_min,
#                      to=grid_max,
#                      by=(to - from) / bins)
#   
#   grid_borders <- pmap(grid_def, seq)
# 
#   points <- as_tibble(expand.grid(grid_borders))
#   names(points) <- names(thetas)
#   
#   # La longeur de la cellule j pour la dimension d est donnée par 
#   # min(grid_def$by[d] / 2, points[j,d] - grid_def$from[d]) + 
#   #   min(grid_def$by[d] / 2, grid_def$to[d] - points[j,d]) 
#   half_grid_def <- grid_def$by / 2
#   limit_below <- points - bounds$min
#   limit_above <- bounds$max - points
#   cond_below <- (half_grid_def < limit_below)
#   cond_above <- (half_grid_def < limit_above)
#   length_below <- (half_grid_def) * cond_below + limit_below * (!cond_below)
#   length_above <- (half_grid_def) * cond_above + limit_above * (!cond_above)
#   dimcell <- length_below + length_above
#   vcell <- reduce(dimcell, `*`)
#   
#   if(ncol(thetas) == 1) {
#     fhat <- ks::kde(x=deframe(transform(thetas)), h=H, supp=supp, w=deframe(weights) * nrow(thetas), eval.points=deframe(transform(points)))
#   } else {
#     fhat <- ks::kde(x=transform(thetas), H=H, w=deframe(weights) * nrow(thetas), eval.points=transform(points))
#   }
#   #density <- fhat$estimate / (apply(points * (1 - points), 1, prod) * prod(bounds$max - bounds$min))
#   scaled_points <- scale_unit(points)
#   #browser()
#   density <- fhat$estimate / (apply((bounds$max - bounds$min) * (1 / C) * scaled_points / ((1 + scaled_points) ** 2), 1, prod))
#   list(result=tibble(points, density=density, vcell=vcell), grid_def=grid_def, grid_borders=grid_borders, dim=ncol(thetas))
# }

kde_bounded_2 <- function(thetas, weights, bounds, bins=NULL, h=NULL, evalp = NULL) {
  # KDE for compact support as described in 
  # Bouezmarni and Rombouts 2010 NONPARAMETRIC DENSITY ESTIMATION FOR MULTIVARIATE BOUNDED DATA

  if(!is.null(evalp) && is.null(h))
    stop("kde_bounded_2: argument h cannot be null when evalp is provided.")
  if(!is.null(evalp) && !is.null(bins))
    stop("kde_bounded_2: argument bins is unused when evalp is provided.")


  bounds_min <- slice(bounds, 1)
  bounds_max <- slice(bounds, 2)

  scale_unit <- function(x) {
    as_tibble(pmap(list(xx = x, bmi = bounds_min, bma = bounds_max), 
      function(xx, bmi, bma) (xx - bmi) / (bma - bmi)))
  }

  unscale <- function(x) {
    as_tibble(pmap(list(xx = x, bmi = bounds_min, bma = bounds_max),
      function(xx, bmi, bma) xx * (bma - bmi) + bmi))
  }
  
  thetas_u <- scale_unit(thetas)
   
  if (is.vector(weights) && length(weights) == 1) 
    weights <- tibble(weight=1) %>% slice(rep(1, nrow(thetas)))

  weights <- deframe(weights / sum(weights))
 
  thetas_min_u <- as_vector(thetas_u %>% summarize_all(min))
  thetas_max_u <- as_vector(thetas_u %>% summarize_all(max))
  thetas_range_u <- thetas_max_u - thetas_min_u

  if(is.null(h)){
    if(ncol(thetas) == 1){
      h_u <- hpi(as_vector(thetas_u))
    } else {
      h_u <- diag(Hpi.diag(thetas_u))
    }
  } else {
    h_u <- h / (bounds_max - bounds_min)
  }

  if (is.null(evalp)) {
    grid_min_u <- as_tibble(pmap(list(t = thetas_min_u, h = h_u), function(t, h)
      qbeta(0.01, t / h + 1, (1 - t) / h + 1)))
    grid_max_u <- as_tibble(pmap(list(t = thetas_max_u, h = h_u), function(t, h)
      qbeta(0.99, t / h + 1, (1 - t) / h + 1)))
    grid_step_u <- (grid_max_u - grid_min_u) / bins

    grid_borders_u <- pmap(list(from=grid_min_u,
                                to=grid_max_u,
                                by=grid_step_u), 
                           seq)
    
    evalp_u <- as_tibble(expand.grid(grid_borders_u))
  } else {
    evalp_u <- scale_unit(evalp)
  }
  names(evalp_u) <- names(thetas)
  
  long_thetas_u <- thetas_u %>% mutate(sample_index=row_number()) %>% gather("variable", "sample_value", -sample_index)
  long_evalp_u <- evalp_u %>% mutate(evalp_index=row_number()) %>% gather("variable", "evalp_value", -evalp_index)
  long_h_u <- gather(h_u, "variable", "h_value")
  fhat_u <- long_thetas_u %>% 
    left_join(long_evalp_u, by="variable") %>% 
    left_join(long_h_u, by="variable") %>%
    mutate(param_alpha = evalp_value / h_value + 1,
           param_beta = (1 - evalp_value) / h_value + 1,
           density = dbeta(sample_value, param_alpha, param_beta)) %>%
    group_by(evalp_index, sample_index) %>%
    summarize(density=prod(density)) %>%
    summarize(density=sum(density * weights))

  # Volume autour de chaque point (dans l'espace d'origine). La longeur de la cellule j pour la dimension d est donnée par 
  # min(grid_def$by[d] / 2, evalp[j,d] - grid_def$from[d]) + 
  #   min(grid_def$by[d] / 2, grid_def$to[d] - evalp[j,d])
  if(is.null(evalp)) {
    grid_min <- unscale(grid_min_u)
    grid_max <- unscale(grid_max_u)
    grid_step <- (grid_max - grid_min) / bins
    evalp <- unscale(evalp_u)
    half_grid_step <- slice(grid_step / 2, rep(1,nrow(evalp)))
    limit_below <- as_tibble(map2(evalp, grid_min, `-`))
    limit_above <- as_tibble(map2(grid_max, evalp, `-`)) 
    cond_below <- half_grid_step < limit_below
    cond_above <- half_grid_step < limit_above
    length_below <- (half_grid_step) * cond_below + limit_below * (!cond_below)
    length_above <- (half_grid_step) * cond_above + limit_above * (!cond_above)
    dimcell <- length_below + length_above
    vcell <- reduce(dimcell, `*`)
   } else {
    vcell <- NA
    grid_borders_u = NULL
  }

  # Scale the density to adapt to the original space, see "transformation of random variables".
  fhat <- evalp %>% mutate(density = fhat_u$density / prod(bounds_max - bounds_min), vcell = vcell)

  list(result=fhat, 
       grid_borders=unscale(grid_borders_u),
       dim=ncol(thetas),
       h=h_u * (bounds_max - bounds_min))
}
