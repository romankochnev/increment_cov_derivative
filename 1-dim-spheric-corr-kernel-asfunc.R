


#####  Define dimensionality k,l ##### 
# For generating purposes only
# Those are not used in functions: functions define k/l
# based on the length of xi_1 and xi_2, that were
# but into functions
k = 2
l = 9

#####  1. Define xi_1 (k x 1) and xi_2 (l x 1) ##### 
xi_1 = matrix(c(2.5, 7.3), nrow = k, ncol = 1)
xi_2 = matrix(seq(0.1, 14, length.out = l), nrow = l, ncol = 1)

# You can choose any xi_1 and xi_2. Select/remove any of them.
xi_2 = c(0.1, 0.2, 4.4, 4.6, 8.4, 8.2, 15,16,17)
xi_2 = matrix(seq(0.1, 20, length.out = l), nrow = l, ncol = 1)

xi_2 = matrix(c(1, 5,13, 18), nrow = 4, ncol = 1)
xi_2 = matrix(seq(1, 14, length.out = 6), nrow = 6, ncol = 1)

#####  2. Define sigma2 and phi ##### 
sigma2 = 1.5
phi = 2

compute_1dim_Sigma_2 = function(xi_1, xi_2, sigma2, phi, nugget=1, r_small = 1e-23){
  k = length(xi_1)
  l = length(xi_2)
  
  partial_sill = sigma2-nugget
  
  #####  3. Define A (k x k), B (1 x k), C_xi_1 (k x k, symmetric) ##### 
  D1 = as.matrix(dist(xi_1))
  D2 = as.matrix(dist(xi_2))
  
  # r_small - rounding parameter. For  small values like 0.00000000000000001 or like that
  C_xi_1 = sigma2 * (D1==0 | D1<=r_small) + (partial_sill) * exp(-D1 / phi) * (D1>r_small)
  C2 = sigma2 * (D2==0 | D2<=1e-23) + (partial_sill) * exp(-D2 / phi) * (D2>r_small)
  
  F1 = xi_1
  F2 = xi_2
  
  B = solve(t(F1) %*% solve(C_xi_1) %*% F1) %*% t(F1) %*% solve(C_xi_1)
  A = solve(C_xi_1) %*% (diag(1,k) - F1 %*% B)
  K = C_xi_1
  
  #####  4. Compute D(xi1, xi2) (k x l matrix of absolute differences) ##### 
  D_xi1_xi2 = matrix(0, nrow = k, ncol = l)
  for (i in 1:k) {
    for (j in 1:l) {
      D_xi1_xi2[i, j] = abs(xi_1[i] - xi_2[j])
    }
  }
  
  #####  5. Compute C(xi1, xi2) = sigma2 * exp(-D(xi1, xi2)/phi) ##### 
  C_xi1_xi2 = (
    sigma2 * (D_xi1_xi2==0 | D_xi1_xi2<=r_small)
    + sigma2 * exp(-D_xi1_xi2 / phi) * (D_xi1_xi2>r_small)
               )
  
  BKBT = B %*% K %*% t(B)  # 1 x 1
  AKBT = A %*% K %*% t(B)  # k x 1
  AKAT = A %*% K %*% t(A)  # k x k
  
  ##### Calculate Σ_2 #####
  
  Sigma_2 = (xi_2 %*% BKBT %*% t(xi_2) 
             + xi_2 %*% B%*%K%*%t(A) %*% C_xi1_xi2 
             + t(C_xi1_xi2) %*% AKBT %*% t(xi_2) 
             + t(C_xi1_xi2) %*% AKAT %*% C_xi1_xi2
             - t(C_xi1_xi2) %*% t(A) %*% C_xi1_xi2
             - t(C_xi1_xi2) %*% t(B) %*% t(xi_2)
             - xi_2 %*% B %*% C_xi1_xi2
             - t(C_xi1_xi2) %*% A %*% C_xi1_xi2
             + C2)
  return(Sigma_2)
}

compute_1dim_der = function(xi_1, xi_2, sigma2, phi, nugget=1, r_small = 1e-23, time_control=F){
  k = length(xi_1)
  l = length(xi_2)
  
  partial_sill = sigma2-nugget
  
  #####  3. Define A (k x k), B (1 x k), C_xi_1 (k x k, symmetric) ##### 
  D1 = as.matrix(dist(xi_1))
  D2 = as.matrix(dist(xi_2))
  
  C_xi_1 = sigma2 * (D1==0 | D1<=r_small) + (partial_sill) * exp(-D1 / phi) * (D1>r_small)
  C2 = sigma2 * (D2==0 | D2<=1e-23) + (partial_sill) * exp(-D2 / phi) * (D2>r_small)
  
  F1 = xi_1
  F2 = xi_2
  
  B = solve(t(F1) %*% solve(C_xi_1) %*% F1) %*% t(F1) %*% solve(C_xi_1)
  A = solve(C_xi_1) %*% (diag(1,k) - F1 %*% B)
  K = C_xi_1
  
  #####  4. Compute D(xi1, xi2) (k x l matrix of absolute differences) ##### 
  D_xi1_xi2 = matrix(0, nrow = k, ncol = l)
  for (i in 1:k) {
    for (j in 1:l) {
      D_xi1_xi2[i, j] = abs(xi_1[i] - xi_2[j])
    }
  }
  
  #####  5. Compute C(xi1, xi2) = sigma2 * exp(-D(xi1, xi2)/phi) ##### 
  C_xi1_xi2 = (
    sigma2 * (D_xi1_xi2==0 | D_xi1_xi2<=r_small)
    + sigma2 * exp(-D_xi1_xi2 / phi) * (D_xi1_xi2>r_small)
  )
  
  #####  6. Compute C(xi2) (l x l self-similarity matrix) ##### 
  D_xi2 = as.matrix(dist(xi_2))  # Distance matrix (l x l)
  C_xi2 = sigma2 * exp(-D_xi2 / phi)
  
  
  ##### --- Core Functions --- ##### 
  compute_dC_dxi2 = function(C_xi1_xi2, xi1, xi2, phi) {
    k = nrow(C_xi1_xi2)
    l = ncol(C_xi1_xi2)
    dC_dxi2 = array(0, dim = c(k, l, l))
    
    for (i in 1:k) {
      for (j in 1:l) {
        for (c in 1:l) {
          if (j == c) {
            dC_dxi2[i, j, c] = -(C_xi1_xi2[i, j] / phi) * sign(xi1[i] - xi2[j])
          }
        }
      }
    }
    return(dC_dxi2)
  }
  
  compute_dCxi2_dxi2 = function(C_xi2, xi2, phi) {
    l = nrow(C_xi2)
    dCxi2_dxi2 = array(0, dim = c(l, l, l))
    
    for (i in 1:l) {
      for (j in 1:l) {
        for (c in 1:l) {
          if (i == c) {
            dCxi2_dxi2[i, j, c] = -(C_xi2[i, j] / phi) * sign(xi2[i] - xi2[j])
          } else if (j == c) {
            dCxi2_dxi2[i, j, c] = (C_xi2[i, j] / phi) * sign(xi2[i] - xi2[j])
          }
        }
      }
    }
    return(dCxi2_dxi2)
  }
  
  compute_dSigma_dxi2 = function(xi2, A, B, K, C_xi1_xi2, C_xi2, dC_dxi2) {
    l = length(xi2)
    dSigma_dxi2 = array(0, dim = c(l, l, l))
    
    # Precompute repeated terms
    BKBT = B %*% K %*% t(B)  # 1 x 1
    AKBT = A %*% K %*% t(B)  # k x 1
    AKAT = A %*% K %*% t(A)  # k x k
    
    BKAT = B %*% K %*% t(A)  # 1 x k
    
    for (i in 1:l) {
      for (j in 1:l) {
        for (c in 1:l) {
          # Term 1: ξ2 BKB' ξ2'
          term1 = BKBT * xi2[j] * (i == c) + BKBT * xi2[i] * (j == c)
          
          # Term 2: ξ2 BKA' C(ξ1, ξ2)'
          term2 = (BKAT %*% C_xi1_xi2[, j]) * (i == c) + 
            sum(xi2 %*% (BKAT) * dC_dxi2[, j, c])
          
          # Term 3: C(ξ1, ξ2) AKB' ξ2'
          term3 = sum(dC_dxi2[, i, c] * (AKBT %*% t(xi2))) + 
            (C_xi1_xi2[, i] %*% AKBT) * (j == c)
          
          # Term 4: C(ξ1, ξ2)' AKA' C(ξ1, ξ2)
          term4 = 2 * sum(dC_dxi2[, i, c] %*% (AKAT %*% C_xi1_xi2)[, j])
          
          # Term 5: -C(ξ1, ξ2)' A' C(ξ1, ξ2)
          term5 = -sum(dC_dxi2[, i, c] %*% (t(A) %*% C_xi1_xi2)[, j]) - 
            sum((t(C_xi1_xi2) %*% t(A))[i, ] %*% dC_dxi2[, j, c])
          
          # Term 6: -C(ξ1, ξ2)' B' ξ2'
          term6 = -sum(dC_dxi2[, i, c] %*% (t(B) %*% t(xi2))) - 
            (t(C_xi1_xi2) %*% t(B))[i] * (j == c)
          
          # Term 7: -ξ2 B C(ξ1, ξ2)
          term7 = -(B %*% C_xi1_xi2)[j] * (i == c) - 
            sum(xi2 %*% B * dC_dxi2[, j, c])
          
          # Term 8: -C(ξ1, ξ2)' A C(ξ1, ξ2)
          term8 = -sum(dC_dxi2[, i, c] * (A %*% C_xi1_xi2)[, j]) - 
            sum((t(C_xi1_xi2) %*% A)[i, ] * dC_dxi2[, j, c])
          
          # Term 9: C(ξ2) (derivative computed separately)
          dCxi2_dxi2 = compute_dCxi2_dxi2(C_xi2, xi2, phi)
          term9 = dCxi2_dxi2[i, j, c]
          
          dSigma_dxi2[i, j, c] = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9
        }
      }
    }
    return(dSigma_dxi2)
  }
  
  ##### Compute the derivative #####
  start_time = Sys.time()
  
  dC_dxi2 = compute_dC_dxi2(C_xi1_xi2, xi_1, xi_2, phi)
  dSigma_dxi2 = compute_dSigma_dxi2(xi_2, A, B, K, C_xi1_xi2, C_xi2, dC_dxi2)
  
  end_time = Sys.time()
  time_taken = end_time - start_time
  
  if (time_control){print(time_taken)}
  
  
  ##### --- ∂ det(Σ)/∂ξ) --- #####
  
  BKBT = B %*% K %*% t(B)  # 1 x 1
  AKBT = A %*% K %*% t(B)  # k x 1
  AKAT = A %*% K %*% t(A)  # k x k
  
  ##### Calculate Σ_2 #####
  
  Sigma_2 = (xi_2 %*% BKBT %*% t(xi_2) 
             + xi_2 %*% B%*%K%*%t(A) %*% C_xi1_xi2 
             + t(C_xi1_xi2) %*% AKBT %*% t(xi_2) 
             + t(C_xi1_xi2) %*% AKAT %*% C_xi1_xi2
             - t(C_xi1_xi2) %*% t(A) %*% C_xi1_xi2
             - t(C_xi1_xi2) %*% t(B) %*% t(xi_2)
             - xi_2 %*% B %*% C_xi1_xi2
             - t(C_xi1_xi2) %*% A %*% C_xi1_xi2
             + C2)
  
  result = array(0, dim = l)
  
  start_time_2 = Sys.time()
  
  for (c in 1:l) {
    result[c] = det(Sigma_2) * sum(diag(solve(Sigma_2) %*% dSigma_dxi2[, , c]))
  }
  
  end_time_2 = Sys.time()
  time_taken_2 = end_time_2 - start_time_2
  if (time_control){print(time_taken_2)}
  
  result_list = list(result, dSigma_dxi2, Sigma_2)
  
  return(result_list)
}

##### Final result #####

result_list = compute_1dim_der(xi_1, xi_2, sigma2, phi)




##### Gradient optimization function #####

do_optimisation = function(xi_1, xi_2, sigma2, phi, 
                           l_rate = 0.05,
                           bounds = c(0,20),
                           iteration_stop = ((bounds[2]-bounds[1])/l_rate)*1.5, 
                           avoid_double=l_rate, repetitions=2, remove_last = 0,
                           draw_plots = T, 
                           debug_mode_global = c(F,F,F),
                           show_fails_global = c(F,F,F)){
  l = length(xi_2)
  start_time = Sys.time()
  
  xi_2_updated = xi_2
  iters = rep(0,l)
  codes = rep(0,l)
  
  determinants = det(compute_1dim_der(xi_1, xi_2, sigma2, phi)[[3]])
  
  # Get distances from xi_2 locations to borders
  # Needed for optimisation order
  D_bounds = matrix(0, nrow = 2, ncol = l)
  for (b in 1:2) {
    for (j in 1:l) {
      D_bounds[b, j] = abs(bounds[b] - xi_2[j])
    }
  }
  optimisation_order = rep(0,2)
  
  # Decide which loc to start first with.
  # We take a loc that is the closest to the border side.
  # Thus we get 2 points, because we have 2 sides.
  # Then we look, which of those 2 points is further from the border.
  # The one further from the border is preferable.
  if (min(D_bounds[2,]) > min(D_bounds[1,])){
    optimisation_order[2] = match(min(D_bounds[1,]),D_bounds)
    optimisation_order[1] = match(min(D_bounds[2,]),D_bounds[2,])
  } else {
    optimisation_order[1] = match(min(D_bounds[1,]),D_bounds)
    optimisation_order[2] = match(min(D_bounds[2,]),D_bounds[2,])
  }
  
  # Optimisation start (locs closest to the bounds)
  
  do_optim_loc1_1dim = function(xi_1, xi_2, sigma2, phi, i, 
                                l_rate = 0.05,
                                bounds = c(0,20),
                                iteration_stop = ((bounds[2]-bounds[1])/l_rate)*1.5,
                                debug_mode=F, 
                                avoid_double=l_rate, 
                                show_fails=debug_mode){
    xi_2_updated = xi_2
    xi_2_prev = xi_2
    xi_2_curr = xi_2
    
    det_start = det(compute_1dim_der(xi_1, xi_2_prev, sigma2, phi)[[3]])
    
    dochange = T
    fail_counter = 0
    
    direction = 1
    iteration = 1
    last_records = rep(0,8)
    code = -1
    double_alert = 0
    
    while (dochange == T) {
      
      # Indicator to track last 8 changes
      iteration_o_prev = iteration_o
      iteration_o = iteration %% 8
      if(iteration_o==0){iteration_o = 8}
      
      if (iteration == 1){
        xi_2_curr[i] = xi_2_curr[i] + direction * l_rate
      }
      
      if (avoid_double!=0){
        while (round(xi_2_curr[i],10) 
               %in%
               round(xi_2_curr[-i],10)){
          xi_2_curr[i] = xi_2_curr[i] + direction * avoid_double
          double_alert = 300
        }
      }
      
      result_prev = compute_1dim_der(xi_1, xi_2_prev, sigma2, phi)
      result_curr = compute_1dim_der(xi_1, xi_2_curr, sigma2, phi)
      
      # Debug messages
      if (debug_mode){
        print(paste("Debug ON. dochange:", dochange, 
                    "|iteration_o", iteration_o,
                    "|iteration", iteration
        )
        )
        print(paste("i:",i, "|Det prev:", round(det(result_prev[[3]]),4),
                    "|Det curr:", round(det(result_curr[[3]]),4)
        )
        )
        print("Points (prev) are:")
        print(t(xi_2_prev))
        print("Points (curr) are:")
        print(t(xi_2_curr))
        print("Gradient (prev) is:")
        print(t(result_prev[[1]]))
        print("Gradient (curr) is:")
        print(t(result_curr[[1]]))
      }
      
      # Change values
      if(abs(result_curr[[1]][i]) < abs(result_prev[[1]][i])){
        direction = direction*1
        fail_counter = 0
      } else {
        fail_counter = fail_counter + 1
        
        if (show_fails){
          print(paste("fail_counter:",fail_counter))
        }
        
        if (fail_counter>=2){
          dochange = F
          code = 5
        }
        
        if (dochange == T){
          direction = direction*(-1)
          xi_2_curr[i] = xi_2_curr[i] + direction * l_rate # Abort change. Back to starting point
        }
      }
      
      if (dochange == T){
        # Locations update
        xi_2_prev[i] = xi_2_curr[i]
        xi_2_curr[i] = xi_2_curr[i] + direction * l_rate
        
        # Control block (stopping)
        last_records[iteration_o] = direction
        
        # Adjust locations if they go out of bounds
        if(xi_2_curr[i] <= bounds[1] | xi_2_curr[i] >= bounds[2] ){
          if(xi_2_curr[i] <= bounds[1]){
            xi_2_curr[i] = bounds[1]
            code = 10
          } else {
            xi_2_curr[i] = bounds[2]
            code = 11
          }
          dochange = F
          result_curr = compute_1dim_der(xi_1, xi_2_curr, sigma2, phi)
        }
        
        iteration = iteration+1
        
      }
      
      # Stop because of long travel
      if (iteration >= iteration_stop){
        dochange = F
        code = 1
      }
      
      # Stop because of direction change (suitable for 1-dim case)
      if (code!=4 & all(tail(last_records,2)!=0) & sum(tail(last_records,2)) == 0){
        dochange = F
        code = 2
      }
      
      ## Stop because of marking time ←→←→←→← / →←→←→←→←
      #if (iteration >= iteration_stop | all(last_records == rep(c(1,-1),4)) | all(last_records == rep(c(-1,1),4)))
      #  {
      #dochange = F
      #code = 3
      #}
      
      if (dochange == F){
        if(abs(result_curr[[1]][i]) < abs(result_prev[[1]][i])){
          updated_loc_i = xi_2_curr[i]
        } else {
          updated_loc_i = xi_2_prev[i]
        }
        
        # Do we need it to avoid values like 1.387779e-17 instead of 0?
        if (code==10){
          updated_loc_i = xi_2_curr[i]
        }
        
        xi_2_updated[i] = updated_loc_i
      }
    }
    
    # Debug messages
    if (debug_mode & dochange == F){
      print("=== Success! The final result is below ===")
      print(paste("Debug ON. dochange:", dochange, 
                  "|iteration_o", iteration_o,
                  "|iteration", iteration
      )
      )
      print(paste("i:",i, 
                  "|Det prev:", round(det(result_prev[[3]]),4),
                  "|Det curr:", round(det(result_curr[[3]]),4),
                  "|Det at start:", round(det_start,4)
      )
      )
      print("Points (prev) are:")
      print(t(xi_2_prev))
      print("Points (curr) are:")
      print(t(xi_2_curr))
      print("Gradient (prev) is:")
      print(t(result_prev[[1]]))
      print("Gradient (curr) is:")
      print(t(result_curr[[1]]))
    }
    
    if (double_alert != 0){
      code = code + double_alert
    }
    
    
    result = list(locs= xi_2_updated, 
                  iterations = iteration, 
                  code = code, 
                  last_records = last_records)
    return(result)
  }
  
  for (i in optimisation_order){
    
    optim_single = do_optim_loc1_1dim(xi_1, xi_2_updated, sigma2, phi, i, 
                                      l_rate = l_rate,
                                      bounds = bounds,
                                      iteration_stop = iteration_stop,
                                      debug_mode=debug_mode_global[1],
                                      show_fails = show_fails_global[1])
    xi_2_updated = optim_single$locs
    iters[i] = optim_single$iterations
    codes[i] = optim_single$code
  }
  
  stamp_sides = Sys.time()
  
  xi_2_history = cbind(xi_2,xi_2_updated)
  
  determinants = c(determinants,
                   det(compute_1dim_der(xi_1, xi_2_updated, sigma2, phi)[[3]]))
  
  if(draw_plots){
    
    plot(x = 1,
         type = "n",
         xlim = c(-1, 21), 
         ylim = c(-0, 2),
         pch = 20,
         xlab = "x", 
         ylab = "y",
         main = paste0("Optimisation step: sides\nDet old: ",
                       determinants[1],
                       " | Det new: ", determinants[2])
    )
    
    points(x = xi_2,
           y = rep(1, length(xi_2)),
           pch = 20,
           col = "green4")
    
    points(x = xi_2_updated,
           y = rep(1, length(xi_2)),
           pch = 20,
           col = "red")
    
    points(x = xi_1,
           y = rep(1, length(xi_1)),
           pch = 20,
           col = "blue")
    abline(v=bounds[1], col="gold3")
    abline(v=bounds[2], col="gold3")
    legend(0, 0.6, 
           legend=c("xi_2 before", 
                    "xi_2 after",
                    "xi_1"),
           col=c('green4',"red", "blue"), 
           pch=16, cex=1.5)
  }
  
  
  do_optim_loc1_1dim_rest = function(xi_1, xi_2_updated, sigma2, phi, exclude,
                                     l_rate = 0.05,
                                     bounds = c(0,20),
                                     iteration_stop = ((bounds[2]-bounds[1])/l_rate)*1.5,
                                     debug_mode=F, remove_last = 0, 
                                     show_fails=debug_mode){
    
    l = length(xi_2_updated)
    
    iters = rep(0,l)
    codes = rep(0,l)
    
    xi_2_rest_points = 1:l
    xi_2_rest_points = xi_2_rest_points[!xi_2_rest_points %in% exclude]
    
    rest_points_amount = length(xi_2_rest_points)
    finals = 1
    
    debug_locs_order = c()
    
    while (finals <= rest_points_amount) {
      loc_to_update = 0
      benchmark_distance = 0
      
      D2_updated = as.matrix(dist(xi_2_updated))
      D_bounds_updated = matrix(0, nrow = 2, ncol = l)
      
      for (b in 1:2) {
        for (j in 1:l) {
          D_bounds_updated[b, j] = abs(bounds[b] - xi_2_updated[j])
        }
      }
      
      ### Choose the location to optimize next
      
      if (length(xi_2_rest_points)>1){
        for(i in xi_2_rest_points){
          for (j in xi_2_rest_points[xi_2_rest_points!=i]){
            ## We are looking for locations with the biggest distance to the nearest neighbors
            #nearest_by_i = sum(sort(D2_updated[i,][D2_updated[i,]!=0])[1:2])
            #nearest_by_j = sum(sort(D2_updated[j,][D2_updated[i,]!=0])[1:2])
            #
            #if (nearest_by_i > benchmark_distance | nearest_by_i > benchmark_distance){
            #  if (nearest_by_i>nearest_by_j){
            #    loc_to_update = i
            #    benchmark_distance = nearest_by_i
            #  } else {
            #    loc_to_update = j
            #    benchmark_distance = nearest_by_j
            #  }
            ## Removed, beacause it does not take left/right into account
            #nearest_by_i = sum(sort(D2_updated[i,][D2_updated[i,]!=0])[1:2])
            #nearest_by_j = sum(sort(D2_updated[j,][D2_updated[i,]!=0])[1:2])
            
            # Acceptable for 1-dim ONLY
            nearest_by_i = sum(D2_updated[i, (i-1):(i+1)])
            nearest_by_j = sum(D2_updated[j, (j-1):(j+1)])
            if (nearest_by_i > benchmark_distance | nearest_by_i > benchmark_distance){
              if (nearest_by_i>nearest_by_j){
                loc_to_update = i
                benchmark_distance = nearest_by_i
              } else if (nearest_by_i<nearest_by_j) {
                loc_to_update = j
                benchmark_distance = nearest_by_j
              } else {
                temp = apply(D_bounds_updated,2,min)
                if (temp[i] < temp[j]){
                  loc_to_update = i
                } else if (temp[i] > temp[j]) {
                  loc_to_update = j
                } else {
                  if (i>j){loc_to_update = i} else {loc_to_update = j}
                }
              }
            }
            
          }
        }
      } else {
        loc_to_update = xi_2_rest_points
      }
      
      debug_locs_order = c(debug_locs_order, loc_to_update)
      
      ### Do optimization
      optim_single = do_optim_loc1_1dim(xi_1, xi_2_updated, sigma2, phi, 
                                        loc_to_update, 
                                        l_rate = l_rate,
                                        bounds = bounds,
                                        iteration_stop = iteration_stop,
                                        debug_mode=debug_mode)
      xi_2_updated = optim_single$locs
      iters[loc_to_update] = optim_single$iterations
      codes[loc_to_update] = optim_single$code
      
      ### Remove the location we have already optimized
      xi_2_rest_points = xi_2_rest_points[xi_2_rest_points!=loc_to_update]
      
      # New cycle begins
      finals = finals + 1
      if (remove_last > 0 & length(xi_2_rest_points)<=remove_last){
        finals = finals + remove_last
      }
    }
    
    result = list(locs = xi_2_updated, 
                  iterations = iters,
                  code = codes,
                  debug_locs_order)
    return(result)
    
  }
  
  optim_rest = do_optim_loc1_1dim_rest(xi_1, xi_2_updated, sigma2, phi, 
                                       optimisation_order, 
                                       l_rate = l_rate,
                                       bounds = bounds,
                                       iteration_stop = iteration_stop,
                                       debug_mode = debug_mode_global[2], 
                                       remove_last = remove_last,
                                       show_fails = show_fails_global[2])
  
  stamp_rest = Sys.time()
  
  xi_2_history = cbind(xi_2_history, optim_rest$locs)
  
  determinants = c(determinants,
                   det(compute_1dim_der(xi_1, optim_rest$locs, sigma2, phi)[[3]])
                   )
  
  iters_history = iters + optim_rest$iterations
  codes_history = codes + optim_rest$code
  
  if (draw_plots){
    plot(x = 1,
         type = "n",
         xlim = c(-1, 21), 
         ylim = c(-0, 2),
         pch = 20,
         xlab = "x", 
         ylab = "y",
         main = paste0("Optimisation step: rest\nDet old: ",
                       determinants[2],
                       " | Det new: ", determinants[3])
    )
    
    points(x = xi_2_history[,2],
           y = rep(1, length(xi_2)),
           pch = 20,
           col = "green4")
    
    points(x = optim_rest$locs,
           y = rep(1, length(xi_2)),
           pch = 20,
           col = "red")
    
    points(x = xi_1,
           y = rep(1, length(xi_1)),
           pch = 20,
           col = "blue")
    abline(v=bounds[1], col="gold3")
    abline(v=bounds[2], col="gold3")
    legend(0, 0.6, 
           legend=c("xi_2 before", 
                    "xi_2 after",
                    "xi_1"),
           col=c('green4',"red", "blue"), 
           pch=16, cex=1.5)
  }
  
  was_repeated = 0
  while (repetitions>0) {
    repetitions = repetitions-1
    was_repeated = was_repeated + 1
    
    optim_rest = do_optim_loc1_1dim_rest(xi_1, optim_rest$locs, sigma2, phi, 
                                         optimisation_order, 
                                         l_rate = l_rate,
                                         bounds = bounds,
                                         iteration_stop = iteration_stop,
                                         debug_mode = debug_mode_global[3], 
                                         remove_last = remove_last,
                                         show_fails = show_fails_global[3])
  }
  
  stamp_repetitions = Sys.time()
  
  if (was_repeated>0){
    
    xi_2_history = cbind(xi_2_history, optim_rest$locs)
    
    iters_history = rbind(iters_history, iters + optim_rest$iterations)
    codes_history = rbind(codes_history, codes + optim_rest$code)
    
    determinants = c(determinants,
                     det(compute_1dim_der(xi_1, optim_rest$locs, sigma2, phi)[[3]])
    )
    
    if (draw_plots){
      plot(x = 1,
           type = "n",
           xlim = c(-1, 21), 
           ylim = c(-0, 2),
           pch = 20,
           xlab = "x", 
           ylab = "y",
           main = paste0("Optimisation step: repeats\nDet old: ",
                         determinants[3],
                         " | Det new: ", determinants[4],
                         " | Rep: +", was_repeated)
      )
      
      points(x = xi_2_history[,3],
             y = rep(1, length(xi_2)),
             pch = 20,
             col = "green4")
      
      points(x = optim_rest$locs,
             y = rep(1, length(xi_2)),
             pch = 20,
             col = "red")
      
      points(x = xi_1,
             y = rep(1, length(xi_1)),
             pch = 20,
             col = "blue")
      abline(v=bounds[1], col="gold3")
      abline(v=bounds[2], col="gold3")
      legend(0, 0.6, 
             legend=c("xi_2 before", 
                      "xi_2 after",
                      "xi_1"),
             col=c('green4',"red", "blue"), 
             pch=16, cex=1.5)
    }
  }
  
  end_time = Sys.time()
  
  time_taken = list(total = end_time - start_time,
                    sides = stamp_sides - start_time,
                    rest = stamp_rest - stamp_sides,
                    rep = stamp_repetitions - stamp_rest)
  
  result = list(locs = optim_rest$locs,
                history_locs = xi_2_history,
                history_det = determinants,
                iters = iters_history,
                codes = codes_history,
                time=time_taken)
  return(result)
}

### The function returns:
# $locs - locations *in the end* of algorithm
# $history_locs - locations at every step: 
#                                          side optimisation, 
#                                          rest locs optimisation, 
#                                          repeated rest locs optimisation
# $history_det - Sigma_2 determinants after every step
# $iters - amount of iterations after first try (1. line) and after repetitions (2. line)
# $codes - codes indicating why algorithm stopped
#                10 - reached left bound.
#                11 - reached right bound.
#                 5 - marking time (←→ / →←)
#                 2 - same as 5, but different tracking method
#               3.. - location has overlaped another point 
#                     (not working, because R cheats: it gives you 2 points
#                      like |loc_1 - loc_2| = 1.386457e-15). Still try to fix it.


##### Optimization try #####

## Obligatory to define.
# xi_1, xi_2, sigma2, phi and the following:
l_rate = 0.05
bounds = c(0,20)

# Iterations per every location.
iteration_stop = ((bounds[2]-bounds[1])/l_rate)*1.5


## Not obligatory to define.
# Amount of algorithm repeats.
repetitions=2

# Debug-modes for every step: side locations, rest locations, repeated rest locations.
# Obligated to have 3 elements.
debug_mode_global=c(F,F,F)

# Extra debug mode (specific message) for the same steps.
# Obligated to have 3 elements.
show_fails = c(F,F,F)

testing = do_optimisation(xi_1, xi_2, sigma2, phi, 
                          l_rate = l_rate,
                          bounds = c(0,20),
                          iteration_stop = iteration_stop, 
                          avoid_double=l_rate, 
                          repetitions=repetitions,
                          draw_plots = T, 
                          debug_mode_global=debug_mode_global)