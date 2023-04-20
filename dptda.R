
library(TDA)
library(dplyr)
library(purrr)

generate_simul_data <- function(n,seed){
  set.seed(seed)
  X = circleUnif(n,0.7)
  par(mfrow=c(1,2))
  plot(X)
  list(X=X)
}

generate_simul_data2 <- function(n,seed){
  set.seed(seed)
  n1 <- n2 <- n %/% 2
  
  X1=circleUnif(n1,r=1)-1.5
  X2=circleUnif(n2,r=1.5)+1.5
  X = rbind(X1, X2)
  plot(X)
  list(X=X)
}

get_birth_death_pair <- function(min, max){
  b <- 0
  d <- 0
  while(d-b <= 0){
    b <- runif(1, min=min, max=max)
    d <- runif(1, min=min, max=max)  
  }
  return(c(b,d))
}

evaluate_utility <- function(true_diag, current_diag){
  dists <- map_dbl(0:max_dim, function(dim) {
    bottleneck(true_diag[["diagram"]], current_diag[["diagram"]], dimension = dim) + bottleneck(true_diag[["diagram"]], current_diag[["diagram"]], dimension = dim)
  })
  return(sum(dists))
}

get_private_diag <- function(x, noise_sd){
  
  pdi = x$diagram
  ix <- sample(1:nrow(pdi), size = 1)
  dim <- pdi[ix,][1]
  type = dim
  
  bd <- pdi[ix,c(2,3)]
  
  gen <- mvtnorm::rmvnorm(1, mean = bd, sigma = diag(2) * noise_sd)
  
  if((gen[1,1] >0) & (gen[1,2] >0) & (gen[1,2] > gen[1,1])){
    pdi[ix,c(2,3)] <- gen[1,c(1,2)]
  }
  
  return(list(diagram=pdi, max_dim = type))
}


run_example <- function(n, noise_sd, r, max_dim, dp_epsilon, true_m0, m0, by, seed, type, remove_ldt=FALSE){
  
  Xlim = c(-2.5, 2.5)
  lim = cbind(Xlim, Xlim)
  r = 1

  if(type == 1){
    C = 2*2 / m0 / (n^(1/r))
    data <- generate_simul_data(n,seed)  
    X <- data$X

  }else if(type == 2){
    C = 7.78 /m0 / (n^(1/r))
    data <- generate_simul_data2(n,seed)  
    X <- data$X
  }
  
  n <- nrow(X)
  
  A <- X %>% as.matrix
  true_Diag <-  gridDiag(X = A, FUN = dtm, r = r, m0 = true_m0, lim = lim, by = by, maxdimension = max_dim, sublevel = T)
  
  # ============ Remove largest Death time ======== # 
  if(remove_ldt == TRUE){
    temp_true_diag <- true_Diag$diagram
    class(temp_true_diag) <- NULL
    temp_true_diag  <- data.frame(temp_true_diag)
    
    max_death <- data.frame(temp_true_diag) %>% pull(Death) %>% max
    temp_true_diag_gr1 <- temp_true_diag %>% filter(dimension != 0)
    temp_true_diag_dim0 <- temp_true_diag %>% filter(dimension == 0, Death < max_death)
    
    temp_true_diag <- rbind(temp_true_diag_dim0, temp_true_diag_gr1) %>% as.matrix()
    class(temp_true_diag) <- 'diagram'
    true_Diag$diagram <- temp_true_diag
  }

  # ============ Generate Initial Diagram ============ # 
  features_true <- true_Diag$diagram[,1]
  true_pd <- true_Diag$diagram[,c(2,3)]
  x_min <- min(true_pd)
  x_max <- max(true_pd)
  
  #  ======== Type 3 ========
  init_diag <- map(0:max_dim, function(x){
    bd_matrix <- map(1:5, ~ get_birth_death_pair(x_min, x_max)) %>% 
      flatten_dbl() %>% 
      matrix(ncol = 2, byrow=T)
    cbind(rep(x,5),bd_matrix)
  })
  init_diag <- do.call(rbind, init_diag)  
  colnames(init_diag) <- colnames(true_Diag$diagram)
  class(init_diag) <- 'diagram'
  
  current_Diag <- list()
  current_Diag$diagram <- init_diag
  attributes(current_Diag$diagram)['maxdimension'] = max_dim
  attr(current_Diag$diagram, 'scale') <- c(min(init_diag), max(init_diag))
  
  p_distance2 <- evaluate_utility(true_Diag, current_Diag)
  
  p1 <- c()
  Diag_list <- list()
  accept_trace <- c()
  d0 <- c()
  d1 <- c()
  
  for(t in 1:10000){
    new_Diag <- get_private_diag(current_Diag, noise_sd)
    dim <- new_Diag[["max_dim"]]
    
    # True vs New 
    p_distance1 <- evaluate_utility(true_Diag, new_Diag)
    
    p1[t] <- p_distance1
    accept_prob = min(0, dp_epsilon/(2*C) * (-p_distance1) - dp_epsilon/(2*C) * (-p_distance2)) # We may need other prior
    u = log(runif(1))
    
    if(u<=accept_prob){
      current_Diag <- new_Diag
      p_distance2 <- evaluate_utility(true_Diag, current_Diag)
      accept_trace[t] <- 1
    }else{
      accept_trace[t] <- 0
    }
    
    Diag_list[[t]] <- current_Diag
    d0[t] <- bottleneck(true_Diag[["diagram"]], current_Diag[["diagram"]], dimension = 0)
    d1[t] <- bottleneck(true_Diag[["diagram"]], current_Diag[["diagram"]], dimension = 1)
    
  }
  return(list(PerD = Diag_list[-9800:-1], d0 = d0, d1 = d1, TrueD = true_Diag))
}
