library(spectralGraphTopology)
library(CVXR)
library(quadprog)

#' @export
learn_connected_bipartite_graph_pgd <- function(S,
                                                r,
                                                q,
                                                init = "naive",
                                                learning_rate = 1e-4,
                                                maxiter = 1000,
                                                reltol = 1e-5,
                                                verbose = TRUE,
                                                record_objective = FALSE,
                                                backtrack = TRUE) {
  # number of nodes
  p <- r + q
  ones_r <- rep(1, r)
  J_rr <- matrix(1, r, r) / p
  J_rq <- matrix(1, r, q) / p
  J_qq <- matrix(1, q, q) / p
  # Laplacian initialization
  L_ <- L(spectralGraphTopology:::w_init(init, MASS::ginv(S)))
  # B initialization
  B <- as.matrix(-L_[1:r, (r+1):p] + 1e-5)
  B <- B / rowSums(B)
  # data cropping
  S_rq <- S[1:r, (r+1):p]
  diag_S <- diag(S)[(r+1):p]
  invI_J <- solve(diag(rep(1, r)) + J_rr)
  lin_term <- ones_r %*% t(diag_S) - 2 * S_rq
  # projected gradient descent
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta",
                                     total = maxiter, clear = FALSE, width = 80)
  lr_seq <- c(learning_rate)
  obj_seq <- c(compute_obj_function(B, compute_g_B(B, J_qq, invI_J, J_rq), lin_term))
  elapsed_time <- c()
  start_time <- proc.time()[3]
  for (i in 1:maxiter) {
    B_shift <- B - J_rq
    g_B <- compute_g_B(B, J_qq, invI_J, J_rq)
    g_inv <- solve(g_B)
    grad_B <- ones_r %*% t(diag(-g_inv)) - 2 * invI_J %*% B_shift %*% t(-g_inv) + lin_term
    if (backtrack) {
      while (TRUE) {
        B_update <- project_onto_simplex(B - learning_rate * grad_B)
        has_converged = norm(B - B_update, 'F')/norm(B, 'F') < reltol
        if (has_converged)
          break
        g_B_update <- compute_g_B(B_update, J_qq, invI_J, J_rq)
        success <- assert_backtracking(B_update, B, g_B_update, g_B, grad_B, lin_term, learning_rate)
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

compute_g_B <- function(B, J_qq, invI_J, J_rq) {
  B_shift <- B - J_rq
  return(diag(colSums(B)) + J_qq - t(B_shift) %*% invI_J %*% B_shift)
}

compute_obj_function <- function(B, g_B, lin_term) {
  eigvals <- eigen(g_B)$values
  return(-sum(log(eigvals)) + sum(B * lin_term))
}

assert_backtracking <- function(B_update, B, g_B_update, g_B, grad_B, lin_term, learning_rate) {
  if (sum(eigen(g_B_update)$values < 1e-7) >= 1)
    return(FALSE)
  obj <- compute_obj_function(B, g_B, lin_term)
  obj_update <- compute_obj_function(B_update, g_B_update, lin_term)
  if (obj_update < (obj + sum(grad_B * (B_update - B)) + 0.5 * (1/learning_rate) * norm(B_update - B, 'F')^2)) {
    return(c(TRUE, obj_update))
  } else {
    return(c(FALSE, NaN))
  }
}

project_onto_simplex <- function(B) {
  return(solve_subproblem_B_quadprog(t(B)))
}
