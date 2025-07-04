library(matrixcalc) # For matrix operations and kronecker products
set.seed(1337)

##### Data input #####
# ARBITRARY POINTS START ----
# comment this section in, if you do not want to use those points
k = 3  # rows in xi_1
l = 5  # rows in xi_2
p = 1  # columns in xi_1 and xi_2
xi_1 = matrix(c(20,22,23,
                21,21,25), nrow = k, ncol = p)


xi_2_col1 = seq(-8, 30, by=7)
#xi_2_col1 = c(1,6,11,16,21,26,31,36) # step is = 5
xi_2_col2 = xi_2_col1
#xi_2_col2 = rep(1,length(xi_2_col1))

xi_2 = cbind(xi_2_col1,xi_2_col2)

k = nrow(xi_1)  # rows in xi_1
l = nrow(xi_2)  # rows in xi_2
p = ncol(xi_1)  # columns in xi_1 and xi_2

# ARBITRARY POINTS END ----

# Synthetic data START ----
# Set dimensions ----
#k = 3  # rows in xi_1
#l = 5  # rows in xi_2
#p = 2  # columns in xi_1 and xi_2

#xi_1 = matrix(rnorm(k * p), nrow = k, ncol = p)  # xi_1 ∈ R^{k×p}
#xi_2 = matrix(rnorm(l * p), nrow = l, ncol = p)  # xi_2 ∈ R^{l×p}
# Synthetic data END ----

## Define COV-parameters: sigma2, phi, kappa and total set of points ----
xi_1_2 = rbind(xi_1,xi_2)

sigma2 = 1.5
phi = 2
kappa = 1.5

# Time control
# time vars are: time_o (overall), time_1 (term 1), time_2 (rest terms)
time_o = Sys.time()

## Define distances ----
D1 = as.matrix(dist(xi_1))
D2 = as.matrix(dist(xi_2))
D12 = as.matrix(dist(xi_1_2))
D12 = D12[1:k,-(1:k)]

F1 = xi_1
F2 = xi_2

# Helper functions ----
# Rho function: ρ(D, ϕ, κ)
rho = function(D, phi, kappa) {
  (1 / (2^(kappa - 1)) * gamma(kappa)) * (D / phi)^kappa * bessel_k(D / phi, kappa)
}


##### Doing C1, C2, C12 ----

C1 = sigma2 * rho(D1 / phi, phi, kappa)
C2 = sigma2 * rho(D2 / phi, phi, kappa)
C12 = sigma2 * rho(D12 / phi, phi, kappa)

diag(C1) = sigma2 * rho(1e-24, phi, kappa)  # Avoid division by zero
diag(C2) = sigma2 * rho(1e-24, phi, kappa)  # Avoid division by zero

# Define symmetric matrices A, B, K ----

B = solve(t(F1) %*% solve(C1) %*% F1) %*% t(F1) %*% solve(C1)
A = solve(C1) %*% (diag(1,k) - F1 %*% B)
K = C1

##### Derivatives: rho, C12, C2 ====
##### Derivative of ρ w.r.t. D
d_rho_dD = function(D, phi, kappa) {
  term1 = (kappa / phi) * (D / phi)^(kappa - 1) * bessel_k(D / phi, kappa)
  term2 = (D / phi)^kappa * (-bessel_k(D / phi, kappa - 1) - (kappa / (D / phi)) * bessel_k(D / phi, kappa)) / phi
  (1 / (2^(kappa - 1) * gamma(kappa))) * (term1 + term2)
}

# Compute ∂C_{1,2}/∂ξ₂ (tensor of shape (k×l) × (l×p))
d_C12_dxi2 = array(0, dim = c(k, l, l, p))
for (m in 1:k) {
  for (i in 1:l) {
    D_mi = D12[m, i]
    s_mi = (sigma2 / (2^(kappa - 1) * gamma(kappa))) * (D_mi^(kappa - 1) / (phi^(kappa + 1))) * bessel_k(D_mi / phi, kappa - 1)
    d_C12_dxi2[m, i, i, ] = -s_mi * (xi_1[m, ] - xi_2[i, ])
  }
}

# Compute ∂C_2/∂ξ₂ (tensor of shape (l×l) × (l×p))
d_C2_dxi2 = array(0, dim = c(l, l, l, p))
for (i in 1:l) {
  for (j in 1:l) {
    D_ij = D2[i, j]
    if (i != j) {
      s_ij = (sigma2 / (2^(kappa - 1) * gamma(kappa))) * (D_ij^(kappa - 1) / (phi^(kappa + 1))) * bessel_k(D_ij / phi, kappa - 1)
      d_C2_dxi2[i, j, i, ] = s_ij * (xi_2[i, ] - xi_2[j, ])
      d_C2_dxi2[i, j, j, ] = s_ij * (xi_2[j, ] - xi_2[i, ])
    }
  }
}

##### Terms precomputation: Sigma_2, T1-T9 ====
# Compute Σ(ξ₂)
M = A %*% K %*% t(A)
G = B %*% K %*% t(A)
H = A %*% K %*% t(B)

T1 = xi_2 %*% B %*% K %*% t(B) %*% t(xi_2) ###  ξ₂BKB'ξ₂'
T2 = xi_2 %*% G %*% C12                    ###  ξ₂BKA'C₁₂
T3 = t(C12) %*% H %*% t(xi_2)              ###  C'₁₂ AKB'C₁₂
T4 = t(C12) %*% M %*% C12                  ###  C'₁₂ AKA'C₁₂
T5 = -t(C12) %*% t(A) %*% C12              ### -C'₁₂ A'C₁₂
T6 = -t(C12) %*% t(B) %*% t(xi_2)          ### -C'₁₂ B'ξ'₂
T7 = -xi_2 %*% B %*% C12                   ### -ξ₂ B C₁₂
T8 = -t(C12) %*% A %*% C12                 ### -C'₁₂ A C₁₂
T9 = C2                                    ###  C2

Sigma = T1 + T2 + T3 + T4 + T5 + T6 + T7 + T8 + T9

dT1_xi_2_prepack = 2 * xi_2 %*% B%*%K%*%t(B)
dT2_xi_2_prepack_1 = G %*% C12          #see below
dT2_xi_2_prepack_2 = xi_2 %*% G         #dC12 is omitted
dT3_xi_2_prepack_1 = t(H) %*% C12       #see below
dT3_xi_2_prepack_2 = xi_2 %*% t(H)      #dC12 is omitted
dT4_xi_2_prepack = 2*M %*% C12          #dC12 is omitted
dT5_xi_2_prepack = -2 * t(A) %*% C12    #dC12 is omitted
dT6_xi_2_prepack_1 = -t(B) %*% t(xi_2)  #dC12 is omitted
dT6_xi_2_prepack_2 = -t(C12) %*% t(B)   #see above
dT7_xi_2_prepack_1 = -B%*%C12           #see below
dT7_xi_2_prepack_2 = -xi_2%*%B          #dC12 is omitted
dT8_xi_2_prepack = -2*A %*% C12         #dC12 is omitted
dT9_xi_2_prepack = sigma2               #dC2 is omitted

# Dimensionality control (just in case)
rbind(dim(dT1_xi_2_prepack),
     dim(dT2_xi_2_prepack_1),
     dim(dT2_xi_2_prepack_2),
     dim(dT3_xi_2_prepack_1),
     dim(dT3_xi_2_prepack_2),
     dim(dT4_xi_2_prepack),
     dim(dT5_xi_2_prepack),
     dim(dT6_xi_2_prepack_1),
     dim(dT6_xi_2_prepack_2),
     dim(dT7_xi_2_prepack_1),
     dim(dT7_xi_2_prepack_2),
     dim(dT8_xi_2_prepack),
     dim(dT9_xi_2_prepack)
)
matrix(c("k", "l", "p", 
         k,l,p),
       2,3, byrow=T)
# Debug-panel (just in case)
print_result = c(0, # for d term 1
                 0, # for d term 2
                 0, # for d term 3
                 0, # for d term 4
                 0, # for d term 5
                 0, # for d term 6
                 0, # for d term 7
                 0, # for d term 8
                 0, # for d term 9
                 0) # for indexes i,j,m,n

##### Compute ∂Σ/∂ξ₂ term-by-term #####
d_Sigma_dxi2 = array(0, dim = c(l, l, l, p))

time_1 = Sys.time()
for (j in 1:l){
    for (n in 1:p){
      
      # Term 1: ∂(ξ₂ B K B' ξ₂')/∂ξ₂
      d_Sigma_dxi2[,j,,n] =  d_Sigma_dxi2[,j,,n] + dT1_xi_2_prepack[j,p]
      dT1 = dT1_xi_2_prepack[j,p]
      
      
    }
}
time_1 = Sys.time() - time_1

time_2 = Sys.time()
for (i in 1:l){
  for (j in 1:l){
    for (m in 1:l){
      for (n in 1:p){
        # Term 2: ∂(ξ₂ G C_{1,2}')/∂ξ₂
        dT2 = (
          (i==m) * dT2_xi_2_prepack_1[n,j] + 
            dT2_xi_2_prepack_2[i,] %*% d_C12_dxi2[,j,m,n]
          )
        
        d_Sigma_dxi2[i,j,m,n] = d_Sigma_dxi2[i,j,m,n] + 2*dT2 # do term 2 and term 3
        
        # Term 3: ∂(C_{1,2} H ξ₂')/∂ξ₂ (symmetric to Term 2)
        
        # Term 4: ∂(C_{1,2}' M C_{1,2})/∂ξ₂
        dT4 = dT4_xi_2_prepack[,i] %*% d_C12_dxi2[,j,m,n]
        
        d_Sigma_dxi2[i,j,m,n] = d_Sigma_dxi2[i,j,m,n] + dT4
        
        # Term 5: ∂(-C_{1,2}' A' C_{1,2})/∂ξ₂
        dT5 = dT5_xi_2_prepack[,i] %*% d_C12_dxi2[,j,m,n]
        
        d_Sigma_dxi2[i,j,m,n] = d_Sigma_dxi2[i,j,m,n] + dT5
        
        # Term 6: ∂(-C_{1,2}' B' ξ₂')/∂ξ₂
        dT6 = (
          (i==m)*dT6_xi_2_prepack_2[j,n] + 
            dT6_xi_2_prepack_1[,i] %*% d_C12_dxi2[,j,m,n]
        )
        
        d_Sigma_dxi2[i,j,m,n] = d_Sigma_dxi2[i,j,m,n] + dT6
        
        # Term 7: ∂(-ξ₂ B C_{1,2})/∂ξ₂
        dT7 = (
          (i==m) * dT7_xi_2_prepack_1[n,j] + 
            dT7_xi_2_prepack_2[i,] %*% d_C12_dxi2[,j,m,n]
        )
        
        d_Sigma_dxi2[i,j,m,n] = d_Sigma_dxi2[i,j,m,n] + 2*dT7
        
        # Term 8: ∂(-C_{1,2}' A C_{1,2})/∂ξ₂
        dT8 = dT8_xi_2_prepack[,i] %*% d_C12_dxi2[,j,m,n]
        
        d_Sigma_dxi2[i,j,m,n] = d_Sigma_dxi2[i,j,m,n] + dT8
        
        # Term 9: ∂(C_2)/∂ξ₂
        dT9 = d_C2_dxi2[i,j,m,n]
        
        d_Sigma_dxi2[i,j,m,n] = d_Sigma_dxi2[i,j,m,n] + dT9
        
        # Debug printing: controlled via "print_result"
        if (print_result[1]==1){cat("dT1 is", dT1, "\n")}
        if (print_result[2]==1){cat("dT2 is", dT2, "\n")}
        if (print_result[3]==1){cat("dT3 is", dT3, "\n")}
        if (print_result[4]==1){cat("dT4 is", dT4, "\n")}
        if (print_result[5]==1){cat("dT5 is", dT5, "\n")}
        if (print_result[6]==1){cat("dT6 is", dT6, "\n")}
        if (print_result[7]==1){cat("dT7 is", dT7, "\n")}
        if (print_result[8]==1){cat("dT8 is", dT8, "\n")}
        if (print_result[9]==1){cat("dT9 is", dT9, "\n")}
        if (print_result[10]==1){cat("i=",i, "j=",j, "m=",m, "n=",n, "\n\n")}
      }
    }
  }
}

time_2 = Sys.time() - time_2

# Compute ∂det(Σ)/∂ξ₂ using Jacobi's formula
det_Sigma = det(Sigma)
Sigma_inv = solve(Sigma)
grad_det = matrix(0, nrow = l, ncol = p)

for (i in 1:l) {
  for (j in 1:p) {
    grad_det[i, j] = det_Sigma * sum(Sigma_inv * d_Sigma_dxi2[, , i, j])
  }
}

time_o = Sys.time() - time_o

#### Output ####
print("Gradient of det(Σ) w.r.t. ξ₂:")
print(grad_det)
print(det_Sigma)
cat(" Time taken overall:","\t", time_o, "\n",
    "Time taken at T1:","\t", time_1, "\n",
    "Time taken at T2-T9:","\t", time_2, "\n")
