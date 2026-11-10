optim_procedure <- function(W,
                            svd_obj,
                            x,
                            tol,
                            max_iter,
                            print_every,
                            D_list,
                            alpha = 0.01,
                            beta1 = 0.9,
                            beta2 = 0.999,
                            epsilon = 1e-8) {
  n <- length(svd_obj$d)
  p <- length(x)/n
  K <- length(D_list)
  triangle_num <- p*(p-1)/2
  large_triangle <- (p+1)*p/2

  # initialise parameters
  estim_vals <- get_estim_inits(p = p, K = K)
  likelihood_vals <- c(likelihood_from_estim(estim = estim_vals, z = x, D_list = D_list, svd_obj = svd_obj, new = TRUE))
  m_vec <- rep(0, p*K*2+K*triangle_num)
  v_vec <- rep(0, p*K*2+K*triangle_num)

  for (current_iter in 1:max_iter) {
    if (current_iter %% print_every == 0) {print(current_iter)}
    current_vals_long <- estim_vals[current_iter,]
    # extract P_k matrices
    P_list <- list()
    for (current_k in 1:K) {
      vals_to_take <- c(((current_k-1)*p + 1):(current_k*p))
      P_list <- P_list %>% append(list(diag(current_vals_long[vals_to_take])))
    }
    # Extract Sigma_k matrices
    Sigma_list <- list()
    for (current_k in 1:K) {
      vals_to_take <- c(((current_k-1)*large_triangle + 1):(current_k*large_triangle)) + p*K
      current_Sigmak_entries <- current_vals_long[vals_to_take]
      current_Sigmak <- matrix(0, nrow = p, ncol = p)
      current_Sigmak[lower.tri(current_Sigmak, diag = TRUE)] <- current_Sigmak_entries # need to set upper triangle as well
      diag_entries <- diag(current_Sigmak)
      current_Sigmak <- current_Sigmak + t(current_Sigmak) - diag(diag_entries)
      Sigma_list <- Sigma_list %>% append(list(current_Sigmak))
    }

    Ainv_list <- get_ainv(P_list = P_list,
                          d = svd_obj$d,
                          U_inv = svd_obj$u_inv,
                          V_inv = svd_obj$v_inv,
                          tUxV = svd_obj$tUxV)
    D_Ainv_list <- lapply(c(1:K), function(k){
      return(Ainv_list[[k]] * diag(D_list[[k]]))
    })
    Ainv_D_list <- lapply(c(1:K), function(k){
      return(t(t(Ainv_list[[k]]) * diag(D_list[[k]])))
    })
    current_V <- calcul_V_P(current_Ainv_list = Ainv_list,
                            current_Sigma_list = Sigma_list,
                            D_list = D_list)

    grad_out <- calcul_grad_P(current_P_list = P_list,
                              current_Sigma_list = Sigma_list,
                              D_list = D_list,
                              W = W,
                              svd_obj = svd_obj,
                              x = x,
                              Ainv_list = Ainv_list,
                              D_Ainv_list = D_Ainv_list,
                              Ainv_D_list = Ainv_D_list,
                              current_V = current_V)

    P_vec <- lapply(grad_out[c(1:K)], function(grad_mat) {
      return(diag(grad_mat))
    }) %>% unlist()
    Sigma_vec <- lapply(grad_out[(K+1):(2*K)], function(grad_mat) {
      return(grad_mat[lower.tri(grad_mat, diag = TRUE)])
    }) %>% unlist()
    grads_long <- c(P_vec, Sigma_vec) * (-1) # grads wrt neg log-likelihood

    m_vec <- m_vec * beta1 + (1 - beta1) * grads_long
    v_vec <- v_vec * beta2 + (1 - beta2) * grads_long^2

    m_hat <- m_vec/(1 - beta1^current_iter)
    v_hat <- v_vec/(1 - beta2^current_iter)

    new_vals_long <- estim_vals[current_iter,] - alpha * m_hat/(sqrt(v_hat) + epsilon)

    estim_vals <- estim_vals %>% rbind(new_vals_long)

    likelihood_vals[current_iter + 1] <- likelihood_from_estim(estim = estim_vals, z = x, D_list = D_list, svd_obj = svd_obj, new = TRUE)

    # check if all converged
    likelihood_converged <- (likelihood_vals[current_iter] - likelihood_vals[current_iter + 1])/likelihood_vals[current_iter] < tol
    all_conv <- all(abs(c(new_vals_long - current_vals_long) / current_vals_long) < tol) & likelihood_converged
    if (all_conv) {break}
  }
  return(estim_vals)
}




