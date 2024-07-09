
KO_Test <- function(Y, X){
  
  # This function implements the KO test of Hansen and Kazemi (2024). 
  
  # Inputs:
  # Y: A TxN matrix of asset excess returns
  # X: A TxK matrix of factor realizations
  
  # NOTE: No NAs are allowed. Panels must be balanced.
  
  # Output:
  # p: Estimated rank of the beta matrix
  # KO: A Kx1 vector whose elements are either ID or NoID, indicating whether a
  #     given factor risk-premium is identified according to the KO test
  # pvals: A Kx1 vector of p-values for the KO test
  
  # Needed packages: dplyr, matrixcalc, Rsolnp, MASS
  
  
  # Number of periods  
  Time <- dim(Y)[1]  
  
  # Number of factors 
  K <- dim(X)[2]
  
  # Number of assets
  N <- dim(Y)[2]
  
  # De-mean Y and X vars
  Y <- Y - t(replicate(Time,colMeans(Y)))
  X <- X - t(replicate(Time,colMeans(X)))
  
  
  # Betas
  Betas <- t(solve(t(X)%*%X) %*% t(X) %*% Y)
  
  ################################### Cragg and Donald #############
  # Vectorized beta
  B <- vec(Betas)
  
  # Cov mat of factors ((1/T)*X'X)
  CovF <- cov(X)
  
  # Regression residuals
  Eps <- (Y - t(Betas %*% t(X)))
  
  # Resid. cov mat
  CovMatE <- (Time - K)^(-1) * t(Eps) %*% Eps
  
  # Cov mat vec(Beta) (under iid assumption)
  Omega_R <- CovMatE %x% solve(CovF)
  
  # Matrix from which we extract eigenvalues
  PsiMat <- CovF %*% t(Betas) %*% solve(CovMatE) %*% Betas
  
  # Sorted eigenvalues
  Ei <- sort(eigen(PsiMat)$values)
  
  # Loop over eigenvalues
  for (p in 0:K){
    
    # Objective function (Cragg and donald (1997) / Ahn, Horenstein, and Wang (2018))
    
    if (p < K){
      Obj <- Time * sum(Ei[1:(K-p)])
    } else {
      Obj <- 0
      
    }
    
    # Set p to first instance where chi-square is not rejected
    if (pchisq(Obj,(N-p)*(K-p), lower.tail = F) > 0.05){
      break
    }
    
  }
  
  
  if (p<K){
  ###################### Main Test #################
  
  
  # Estimated rank
  r <- p
  # Deficiency
  d <- K - r
  
  # Initialize
  KO <- matrix(NA, nrow = K)
  pval <- matrix(NA, nrow = K)
  
  # Loop over each factor
  for (pos in 1:K){
    
    # Equality constraints
    eqcon_fun <- function(u){
      
      # Make vectorized u into matrix
      U <- matrix(u, nrow = K-1, ncol = d)
      
      # If testing the first position
      if (pos == 1){
        # Place row of zeros on top (where the t(U) * e =0)
        U <- rbind(rep(0,d),U)
        # If testing last position
      } else if (pos == K){
        # Put row of 0s at the end
        U <- rbind(U, rep(0, d))
        # If not first or last position in test
      } else {
        # Place row of zeros where position would be
        U <- rbind(U[1:(pos-1),],
                   rep(0, d),
                   U[(pos):(K-1),])
      }
      
      
      
      # U'U 
      UpU <- t(U) %*% U
      # UpU is symmetric, so top right is redundant. Keep only lower triangular part
      LowerTriIndex <- lower.tri(UpU, TRUE)
      # Select lower triangular part and vectorize
      Constraint <-UpU[LowerTriIndex]
      
      return(c(Constraint))
    }
    
    
    # Optim function
    optfun <- function(u){
      # Make vectorized u in matrix
      U <- matrix(u, nrow = K-1, ncol = d)
      # If testing the first position
      if (pos == 1){
        # Place row of zeros on top (where the t(U) * e =0)
        U <- rbind(rep(0,d),U)
        # If testing last position
      } else if (pos == K){
        # Put row of 0s at the end
        U <- rbind(U, rep(0, d))
        # If not first or last position in test
      } else {
        # Place row of zeros where position would be
        U <- rbind(U[1:(pos-1),],
                   rep(0, d),
                   U[(pos):(K-1),])
      }
      
      
      # Diag matrix NxN
      IN <- diag(N)
      # Impled CovMat of U
      OmegaU<- (IN %x% t(U)) %*% Omega_R %*% (IN %x% U)
      # Test statistic
      Out <-Time* t(vec(t(Betas %*% U))) %*% ginv(OmegaU) %*% vec(t(Betas %*% U))
      return(Out)
    }
    #
    
    # Create identity matrix with 0 row. Lower triangular (including diag)
    # should be the set of constraints
    Z <- matrix(0, K-1, d)
    diag(Z) <- 1
    # This is first equality constraint requirement
    # If testing the first position
    if (pos == 1){
      # Place row of zeros on top (where the t(U) * e =0)
      Z <- rbind(rep(0,d),Z)
      # If testing last position
    } else if (pos == K){
      # Put row of 0s at the end
      Z <- rbind(Z, rep(0, d))
      # If not first or last position in test
    } else {
      # Place row of zeros where position would be
      Z <- rbind(Z[1:(pos-1),],
                 rep(0, d),
                 Z[(pos):(K-1),])
    }
    
    
    # Z'Z
    ZZ <- t(Z) %*% Z
    
    # Lower triangular of Z'Z is what lower triangular of U'U should equal 
    Z1Sol <- ZZ[lower.tri(ZZ, TRUE)]
    
    # Pick random starting point
    Start <- rnorm((K-1)*d)
    
    # Solve
    Sol <- solnp((Start), fun = optfun, eqfun = eqcon_fun, eqB = c(Z1Sol))
    
    # Null: Is ID
    if (pchisq(Sol$values[length(Sol$values)],N*d, lower.tail = F) > 0.05){
      KO[pos] <- "ID"
    } else {
      KO[pos] <- "NoID"
    }
    pval[pos] <- pchisq(Sol$values[length(Sol$values)],N*d, lower.tail = F)
    
  }
  } else {
   KO <- "AllID"
   pval <- NA
  }
  return(list(p,KO,pval))
}