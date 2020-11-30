library(spatstat)
library(pbmcapply)
library(GET)
library(ggplot2)
library(gridExtra)

# Functions ---------------------------------------------------------------

Kcyl <- function(X, base_radius, heights, dir) {
# This function was kindly shared with us by Andreas Dyreborg Christoffersen.
  ## This function computes the empirical cylindrical K-function
  ## X -- point pattern of class pp3
  ## base_radius -- numeric vector constituting the base_radii of the cylinder
  ## heights -- numeric vector constituting the heights of the cylinder
  ## dir -- the direction of the cylinder. One of "x", "y" or "z"
  
  stopifnot(verifyclass(X, "pp3"))
  dom <- X$domain
  switch(dir,
         "x" = {
           X_mat <- as.matrix(X)[, c(2, 3, 1)]
           dom <- box3(xrange = dom$yrange, yrange = dom$zrange, zrange = dom$xrange)
         },
         "y" = {
           X_mat <- as.matrix(X)[, c(1, 3, 2)]
           dom <- box3(xrange = dom$xrange, yrange = dom$zrange, zrange = dom$yrange)
         },
         "z" = {
           X_mat <- as.matrix(X)
         })
  
  ## Find pairs of points within distance max(r) in the (x, y)-plane
  cl <- closepairs.ppp(as.ppp(X_mat[,-3], W = bounding.box.xy(X_mat[,-3])), rmax = max(base_radius))
  cl <- as.data.frame(cl)
  
  ## Find the pairs of cl satisfying that the z-coordinate is within distance max(h)
  z <- X_mat[,3]
  cl$dz <- z[cl$i]-z[cl$j]
  cl <- cl[cl$dz<=max(heights) & cl$dz >= 0,]
  
  ## Compute the weights for the Osher-Stoyan taranslation edge correction
  cl$w <- 1 / ((diff(dom$xrange) - abs(cl$dx))
               * (diff(dom$yrange) - abs(cl$dy))
               * (diff(dom$zrange) - abs(cl$dz)))
  
  ## Determine between which r and h values the pairs are located
  cell_r <- findInterval(cl$d, base_radius)
  cell_h <- findInterval(cl$dz, heights)
  
  ## Create grid with info on relevant points (points satisfying that they are cylinder close) (index and weights)
  nr <- length(base_radius)
  nh <- length(heights)
  grid <- expand.grid(r = 0:nr, h = 0:nh)
  grid$both <- paste(grid$r, grid$h, sep = "-")
  cell <- factor(paste(cell_r, cell_h, sep = "-"), levels = grid$both)
  grid$w <- unlist(lapply(split(cl$w, cell), sum, na.rm = TRUE))
  
  ## Create matrix where each entry is the sum of weights of the points satisfying they are cylinder-close
  mat <- matrix(grid$w, nrow = nr+1, ncol = nh+1)
  mat <- t(apply(mat, 1, cumsum))
  mat <- apply(mat, 2, cumsum)
  mat <- mat[-nrow(mat),-ncol(mat), drop = FALSE]
  
  colnames(mat) <- paste0("height = ", heights)
  rownames(mat) <- paste0("base radius = ", base_radius)
  mat <- 2 * mat
  
  ## Compute the estimated K-function
  np <- npoints(X)
  K2 <- volume.box3(dom)^2 / (np * (np - 1)) * mat
  
  ## Return result
  list(Kest = K2,
       base_radii= base_radius,
       heights = heights,
       direction = dir,
       edge_correction = "Osher-Stoyan")
}

ge_Kcyl <- function(pp3, br, t, directions = c("x", "y", "z"), nsim = 4000, ncores = 1){
  ## This function computes p-values and global envelopes for an extreme rank length global envelope test (GET) 
  ## based on the cylindrical K-function where the null hypothesis is that the observed point pattern is a 
  ## realisation of complete spatial randomness (CSR)
  ## pp3 -- the observed 3D point pattern
  ## br -- numeric vector constituting the base radii of the cylinders of the cylindrical K-function
  ## t -- numeric vector constituting the heights of the cylinders of the cylindrical K-function
  ## directions -- the directions of the cylinders of the cylindrical K-function. A GET will be calculated
  ## for each direction.
  ## nsim -- number of simulations under CSR used for the GET
  ## ncores -- if greater than 1 some computatins will be run in parallel on ncores cores
  
  sims <- rpoispp3(intensity(pp3), pp3$domain, nsim = nsim)
  
  ge <- as.list(rep(NA, length(directions)))
  names(ge) <- directions
  
  for (d in directions) {
    print(paste("direction", d, ":"))
    
    K <- Kcyl(pp3, br, t, dir = d)
    K_sim <- pbmclapply(sims, Kcyl, base_radius = br, heights = t, dir = d, mc.cores = ncores)
    K_sim_M <- lapply(K_sim, function(x){x$Kest})
    K_sim_A <- array(dim = c(length(br), length(t), nsim))
    for (i in 1:nsim) {
      K_sim_A[,,i] <- K_sim_M[[i]]
    }
    im_set <- create_image_set(list(r = list(br, t), obs = K$Kest, sim_m = K_sim_A))
    ge[[d]] <- global_envelope_test(im_set, type = "erl")
  }
  
  return(ge)
}

ge_Kcyl_fixed_t <- function(pp3, br, t, directions = c("x", "y", "z"), nsim = 2000, ncores = 1){
  ## This function computes p-values and global envelopes for an extreme rank length global envelope test (GET) 
  ## based on the cylindrical K-function where the null hypothesis is that the observed point pattern is a 
  ## realisation of complete spatial randomness (CSR)
  ## pp3 -- the observed 3D point pattern.
  ## br -- numeric vector constituting the base radii of the cylinders of the cylindrical K-function.
  ## t -- numeric value giving the fixed height of the cylinders of the cylindrical K-function.
  ## directions -- the directions of the cylinders of the cylindrical K-function. A GET will be calculated
  ## for each direction.
  ## nsim -- number of simulations under CSR used for the GET.
  ## ncores -- if greater than 1 some computatins will be run in parallel on this number of cores.
  sims <- rpoispp3(intensity(pp3), pp3$domain, nsim = nsim)
  
  sum_fun <- function(pp3){
    c(Kcyl(pp3, br, t, dir = "x")$Kest) - 2 * pi * br ^ 2 * t
  }
  
  K_sims <- pbmclapply(sims, sum_fun, mc.cores = ncores)
  K_sims <- simplify2array(K_sims)
  
  K_obs <- ge_list <- as.list(rep(NA, length(directions)))
  names(K_obs) <- names(ge_list) <- directions
  for (d in directions) {
    K_obs[[d]] <- c(Kcyl(pp3, br, t, dir = d)$Kest) - 2 * pi * br ^ 2 * t
    curves <- create_curve_set(list(r = br, obs = K_obs[[d]], sim_m = K_sims))
    ge_list[[d]] <- global_envelope_test(curves)
  }
  
  return(ge_list)
}

plot_ge_Kcyl <- function(genv_Kcyl){
  ## This function plots the returned objekt from the function ge_Kcyl. It returns a list
  ## with a plot for each direction.
  ## genv_Kcyl -- an object returned by ge_Kcyl.
  p <- as.list(rep(NA, length(genv_Kcyl)))
  names(p) <- names(genv_Kcyl)
  
  for (d in names(genv_Kcyl)){
    ge <- genv_Kcyl[[d]]
    data <- data.frame(Radius = ge$x, Height = ge$y, Z = "within")
    data$Z[ge$obs < ge$lo] <- "below"
    data$Z[ge$obs > ge$hi] <- "above"
    data$Z <- factor(data$Z, levels = c("above", "below", "within"))
    
    p[[d]] <- ggplot(data = data, mapping = aes(x = Radius, y = Height, fill = Z)) + 
      geom_tile() + 
      scale_fill_manual(values = c("above" = "black", "below" = "gray95", "within" = "grey")) +
      theme_bw() +
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = 22),
            axis.title = element_text(size = 20),
            axis.text = element_text(size = 20)) +
      labs(fill = "") +
      ggtitle(paste("Direction", d))
  }
  
  if (length(p) == 1){
    print(p[[1]])
  } 
  
  if (length(p) == 2){
    print(grid.arrange(p[[1]], p[[2]], nrow = 1))
  }
  
  if (length(p) == 3){
    print(grid.arrange(p[[1]], p[[2]], p[[3]], nrow = 1))
  }
  return(p)
}

plot_ge_Kcyl_fixed_t <- function(genv_Kcyl_fixed_t){
  ## This function creates a plot of the returned objekt from the function ge_Kcyl_fixed_t. 
  ## genv_Kcyl_fixed_t -- an object returned by ge_Kcyl_fixed_t.
  
  obs <- lapply(genv_Kcyl_fixed_t, function(l){l$obs})
  br <- genv_Kcyl_fixed_t[[1]]$r
  dat_obs <- data.frame(r = rep(br, length(obs)),
                        obs = unlist(obs),
                        dir = rep(names(obs), each = length(br)),
                        lo = rep(genv_Kcyl_fixed_t[[1]]$lo, length(obs)),
                        hi = rep(genv_Kcyl_fixed_t[[1]]$hi, length(obs))) 
  p <- ggplot(dat_obs) + 
    geom_ribbon(aes(x = r, ymin = lo, ymax = hi), fill = "grey") +
    geom_line(aes(x = r, y = obs, linetype = dir)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(angle = 90, size = 20, hjust = 0.5),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          plot.title = element_text(size = 22)) +
    xlab("Radius")
  
  return(p)
}


# Example of use ----------------------------------------------------------
# Load point pattern data. (Such a point pattern can be created with the function pp3 
# from the package spatstat. see the help page for this function for help on how).
# Note that the file path in the below load command will depend on where the pp_data file is saved.
load(file = "pp_data")

#The following uses the data from subject 2 (pp_2) as an example.
# Make global envelopes.
ge <- ge_Kcyl(pp_2, br = seq(0, 25, length.out = 64), t = seq(0, 80, length.out = 64))
ge_fixed_t <- ge_Kcyl_fixed_t(pp_2, br = seq(0, 25, length.out = 64), t = 80)

# Plot global envelopes like in Figure 6.
p <- plot_ge_Kcyl(ge)
plot_ge_Kcyl_fixed_t(ge_fixed_t)
