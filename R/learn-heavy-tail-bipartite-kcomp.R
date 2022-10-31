library(spectralGraphTopology)
library(CVXR)
library(quadprog)

#' @export
learn_heavy_tail_kcomp_bipartite_graph <- function(X,
                                                   r,
                                                   q,
                                                   k,
                                                   nu = 2.001,
                                                   rho = 1,
                                                   learning_rate = 1e-4,
                                                   maxiter = 1000,
                                                   reltol = 1e-5,
                                                   init = "default",
                                                   verbose = TRUE,
                                                   record_objective = FALSE) {
  X <- as.matrix(X)
  # number of samples
  n <- nrow(X)
  # number of nodes
  p <- r + q
  # auxiliary constants
  ones_r <- rep(1, r)
  G <- vector(mode = "list", length = n)
  h <- vector(mode = "list", length = n)
  for (i in 1:n) {
    XX_t <- X[i, ] %*% t(X[i, ])
    XX_t12 <- XX_t[1:r, (r+1):p]
    G[[i]] <- diag(XX_t)[(r+1):p] %*% t(ones_r) - 2 * t(XX_t12)
    h[[i]] <- sum(diag(XX_t)[1:r])
  }
  Ones <- matrix(1, r, r)
  Iq <- diag(q)
  # Laplacian initialization
  Sinv <- MASS::ginv(cor(X))
  L_ <- L(spectralGraphTopology:::w_init("naive", Sinv))
  if (init == "default") {
    # B initialization
    B <- as.matrix(-L_[1:r, (r+1):p] + 1e-5)
    B <- B / rowSums(B)
  } else {
    B <- init / rowSums(init)
  }
  # Y initilization
  Y <- matrix(0, p, p)
  # ADMM constants
  mu <- 2
  tau <- 2
  # diagnose
  dual_residual_seq <- c()
  primal_residual_seq <- c()
  aug_lag_seq <- c()
  rho_seq <- c()
  aug_lag_seq <- c(compute_lagrangian_kcomp(B, L_, Y, G, h, p, n, nu, rho))
  elapsed_time <- c()
  start_time <- proc.time()[3]
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta",
                                     total = maxiter, clear = FALSE, width = 80)
  for (i in 1:maxiter) {
    # update L_
    eig <- eigen(rho * from_B_to_laplacian(B) - Y, symmetric = TRUE)
    R <- eig$vectors[, 1:(p-k)]
    gamma <- eig$values[1:(p-k)]
    L_ <- R %*% diag((gamma + sqrt(gamma^2 + 4 * rho)) / (2 * rho)) %*% t(R)
    # update B
    M <- compute_M(G, B, h, nu, r, q, n)
    H <- compute_H(Y, L_, rho, ones_r, r, p)
    grad_B <- t(H) + M + rho * Ones %*% B + 2 * rho * B
    B_update <- project_onto_simplex(B - learning_rate * grad_B)
    # update Y
    primal_residual <- L_ - from_B_to_laplacian(B_update)
    Y <- Y + rho * primal_residual
    # update rho
    dual_residual <- rho * (B_update - B)
    norm_dual_residual <- rho * norm(from_B_to_laplacian(dual_residual), "F")
    norm_primal_residual <- norm(primal_residual, "F")
    dual_residual_seq <- c(dual_residual_seq, norm_dual_residual)
    primal_residual_seq <- c(primal_residual_seq, norm_primal_residual)
    rho_seq <- c(rho_seq, rho)
    aug_lag_seq <- c(aug_lag_seq, compute_lagrangian_kcomp(B_update, L_, Y, G, h, p, n, nu, rho))
    if (norm_primal_residual > mu * norm_dual_residual)
      rho <- rho * tau
    else if (norm_dual_residual > mu * norm_primal_residual)
      rho <- rho / tau
    if (verbose)
      pb$tick()
    elapsed_time <- c(elapsed_time, proc.time()[3] - start_time)
    has_converged = norm(B - B_update, 'F')/norm(B, 'F') < reltol
    B <- B_update
    if (has_converged)
      break
  }
  B[B < 1e-5] <- 0
  results <- list(laplacian = from_B_to_laplacian(B),
                  adjacency = from_B_to_adjacency(B),
                  B = B,
                  maxiter = i,
                  convergence = has_converged,
                  dual_residual = dual_residual_seq,
                  primal_residual = primal_residual_seq,
                  aug_lag = aug_lag_seq,
                  rho = rho_seq,
                  elapsed_time = elapsed_time)
  return(results)
}


compute_H <- function(Y, L_, rho, ones_r, r, p) {
  Y_rq <- Y[1:r, (r+1):p]
  L_rq <- L_[1:r, (r+1):p]
  y_q <- diag(Y)[(r+1):p]
  l_q <- diag(L_)[(r+1):p]
  return(2*t(Y_rq) - y_q %*% t(ones_r) - rho * (l_q %*% t(ones_r) - 2 * t(L_rq)))
}

compute_lagrangian_kcomp <- function(B, L_, Y, G, h, p, n, nu, rho) {
  log_term <- 0
  for (i in 1:n) {
    log_term <- log_term + log(1 + (h[[i]] + sum(B * t(G[[i]]))))
  }
  eigvals <- eigen(L_, symmetric = TRUE)$values
  eigvals <- eigvals[eigvals > 1e-8]
  L_B <- from_B_to_laplacian(B)
  L_diff <- L_ - L_B
  return(((p+nu)/n) * log_term - log(sum(eigvals)) + sum(L_diff * Y) + (rho * 0.5) * norm(L_diff, "F")^2)
}
