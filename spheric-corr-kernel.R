

#####  Define dimensionality k,l ##### 
k = 2
l = 9

#####  1. Define xi_1 (k x 1) and xi_2 (l x 1) ##### 
xi_1 = matrix(c(2.5, 7.3), nrow = k, ncol = 1)
xi_2 = matrix(seq(0.1, 14, length.out = l), nrow = l, ncol = 1)

#####  2. Define sigma2 and phi ##### 
sigma2 = 4
phi = 3

#####  3. Define A (k x k), B (1 x k), C_xi_1 (k x k, symmetric) ##### 
D1 = as.matrix(dist(xi_1))+1e-23  # Avoid division by zero
D2 = as.matrix(dist(xi_2))+1e-23  # Avoid division by zero

C_xi_1 = sigma2 * exp(-D1 / phi)
C2 = sigma2 * exp(-D2 / phi)

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
C_xi1_xi2 = sigma2 * exp(-D_xi1_xi2 / phi)

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
  
  for (i in 1:l) {
    for (j in 1:l) {
      for (c in 1:l) {
        # Term 1: ξ2 BKB' ξ2'
        term1 = BKBT * xi2[j] * (i == c) + BKBT * xi2[i] * (j == c)
        
        # Term 2: ξ2 BKA' C(ξ1, ξ2)'
        term2 = (B %*% K %*% t(A) %*% C_xi1_xi2[, j]) * (i == c) + 
          sum(xi2 %*% (B %*% K %*% t(A)) * dC_dxi2[, j, c])
        
        # Term 3: C(ξ1, ξ2) AKB' ξ2'
        term3 = sum(dC_dxi2[, i, c] * (A %*% K %*% t(B) %*% t(xi2))) + 
          (C_xi1_xi2[, i] %*% A %*% K %*% t(B)) * (j == c)
        
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
print(time_taken)

##### A slice (e.g., ∂Σ/∂ξ_{2,1}) #####
dSigma_dxi2[, , 1]

##### --- ∂ det(Σ)/∂ξ) --- #####

BKBT = B %*% K %*% t(B)  # 1 x 1
AKBT = A %*% K %*% t(B)  # k x 1
AKAT = A %*% K %*% t(A)  # k x k

##### Calculate Σ_2 #####

Sigma_2 = (xi_2 %*% BKBT %*% t(xi_2) 
+ xi_2 %*% B%*%K%*%t(A) %*% C_xi1_xi2 
+ t(C_xi1_xi2) %*% AKBT %*% t(xi_2) 
+ t(C_xi1_xi2) %*% AKAT %*% C_xi1_xi2
- xi_2 %*% B %*% C_xi1_xi2
- t(C_xi1_xi2) %*% A %*% C_xi1_xi2
- t(C_xi1_xi2) %*% t(B) %*% t(xi_2)
- t(C_xi1_xi2) %*% t(A) %*% C_xi1_xi2
+ C2)

##### Final result #####

result = array(0, dim = l)

start_time_2 = Sys.time()

for (c in 1:l) {
      result[c] = det(Sigma_2) * sum(diag(solve(Sigma_2) %*% dSigma_dxi2[, , c]))
}

end_time_2 = Sys.time()
time_taken_2 = end_time_2 - start_time_2
print(time_taken_2)