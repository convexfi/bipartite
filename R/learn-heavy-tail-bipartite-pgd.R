library(spectralGraphTopology)
library(CVXR)
library(quadprog)

#' @export
learn_heavy_tail_bipartite_graph_pgd <- function(X,
                                                 r,
                                                 q,
                                                 nu = 2.001,
                                                 learning_rate = 1e-4,
                                                 maxiter = 1000,
                                                 reltol = 1e-5,
                                                 init = "default",
                                                 verbose = TRUE,
                                                 record_objective = FALSE,
                                                 backtrack = TRUE) {
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
  J_rr <- matrix(1, r, r) / p
  J_rq <- matrix(1, r, q) / p
  J_qq <- matrix(1, q, q) / p
  if (init == "default") {
    # Laplacian initialization
    Sinv <- MASS::ginv(cor(X))
    L_ <- L(spectralGraphTopology:::w_init("naive", Sinv))
    # B initialization
    B <- as.matrix(-L_[1:r, (r+1):p] + 1e-5)
    B <- B / rowSums(B)
  } else {
    B <- init / rowSums(init)
  }
  invI_J <- solve(diag(rep(1, r)) + J_rr)
  M <- compute_M(G, B, h, nu, r, q, n)
  # projected gradient descent
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta",
                                     total = maxiter, clear = FALSE, width = 80)
  lr_seq <- c(learning_rate)
  obj_seq <- c(compute_obj_function(B, compute_g_B(B, J_qq, invI_J, J_rq), M))
  elapsed_time <- c()
  start_time <- proc.time()[3]
  for (i in 1:maxiter) {
    B_shift <- B - J_rq
    g_B <- compute_g_B(B, J_qq, invI_J, J_rq)
    g_inv <- solve(g_B)
    grad_B <- ones_r %*% t(diag(-g_inv)) - 2 * invI_J %*% B_shift %*% t(-g_inv) + M
    if (backtrack) {
      while (TRUE) {
        B_update <- project_onto_simplex(B - learning_rate * grad_B)
        has_converged = norm(B - B_update, 'F')/norm(B, 'F') < reltol
        if (has_converged)
          break
        g_B_update <- compute_g_B(B_update, J_qq, invI_J, J_rq)
        success <- assert_backtracking(B_update, B, g_B_update, g_B, grad_B, M, learning_rate)
        if (success[1]) {
          obj_seq <- c(obj_seq, success[2])
          learning_rate <- 2 * learning_rate
          B <- B_update
          break
        } else {
          learning_rate <- 0.5 * learning_rate
        }
      }
      lr_seq <- c(lr_seq, learning_rate)
    } else {
      B_update <- project_onto_simplex(B - learning_rate * grad_B)
      has_converged = norm(B - B_update, 'F')/norm(B, 'F') < reltol
      if (has_converged)
        break
      else
        B <- B_update
    }
    M <- compute_M(G, B_update, h, nu, r, q, n)
    if (verbose)
      pb$tick()
    elapsed_time <- c(elapsed_time, proc.time()[3] - start_time)
    if (has_converged)
      break
  }
  results <- list(laplacian = from_B_to_laplacian(B),
                  adjacency = from_B_to_adjacency(B),
                  B = B,
                  maxiter = i,
                  convergence = has_converged,
                  lr_seq = lr_seq,
                  obj_seq = obj_seq,
                  elapsed_time = elapsed_time)
  return(results)
}


compute_M <- function(G, B, h, nu, r, q, n) {
  M <- matrix(0, q, r)
  for (i in 1:n) {
    M <- M + G[[i]] / (nu + h[[i]] + sum(B * t(G[[i]])))
  }
  return(t(M) * (r + q + nu) / n)
}
