
WprevSeSp <- function(x, n, w, cn, mn, cp, mp,
                      method = c("binomial", "poisson"),
                      conf.level = 0.95,
                      nmc = 1e5, seed = 49201) {
  method <- match.arg(method)
  
  if (length(unique(sapply(list(x, n, w), length))) > 1) {
    stop("x, n, and w must be of equal length")
  }
  
  if (!isTRUE(all.equal(sum(w), 1))) {
    stop("Weights must sum to 1.")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  g <- function(theta_1, phi_n, phi_p) {
    result <- numeric(length(theta_1))
    
    ind_for_computed <- which((phi_p >= theta_1) & (theta_1 >= phi_n))
    computed <- (theta_1 - phi_n) / (phi_p - phi_n)
    
    result[(phi_n < phi_p) & (phi_p < theta_1)] <- 1
    result[ind_for_computed] <- computed[ind_for_computed]
    result
  }
  
  alpha <- 1 - conf.level
  theta_hat <- x / n
  s_w_theta_hat <- sum(w * theta_hat)
  
  samples_B_phi_n_L <- rbeta(n = nmc, shape1 = cn, shape2 = mn - cn + 1)
  samples_B_phi_n_U <- rbeta(n = nmc, shape1 = cn + 1, shape2 = mn - cn)
  
  samples_B_phi_p_L <- rbeta(n = nmc, shape1 = cp, shape2 = mp - cp + 1)
  samples_B_phi_p_U <- rbeta(n = nmc, shape1 = cp + 1, shape2 = mp - cp)
  
  if (method == "binomial") {
    n_eff <-
      if (s_w_theta_hat == 0) {
        sum(n)
      } else {
        (s_w_theta_hat * (1 - s_w_theta_hat)) / (sum(w^2 / n * theta_hat))
      }
    x_eff <- n_eff * s_w_theta_hat
    
    samples_B_KG_L <- rbeta(n = nmc, shape1 = x_eff, shape2 = n_eff - x_eff + 1)
    samples_B_KG_U <- rbeta(n = nmc, shape1 = x_eff + 1, shape2 = n_eff - x_eff)
    
    samples_g_L <- g(samples_B_KG_L, samples_B_phi_n_U, samples_B_phi_p_U)
    samples_g_U <- g(samples_B_KG_U, samples_B_phi_n_L, samples_B_phi_p_L)
  } else if (method == "poisson") {
    y <- sum(w / n * x)
    v <- sum((w / n)^2 * x)
    wm <- max(w / n)
    y_star <- y + wm
    v_star <- v + wm^2
    
    samples_G_beta_star_L <-
      if (y == 0) {
        numeric(nmc)
      } else {
        rgamma(n = nmc, y^2 / v, scale = v / y)
      }
    
    samples_G_beta_star_U <- rgamma(n = nmc, y_star^2 / v_star, scale = v_star / y_star)
    
    samples_g_L <- g(samples_G_beta_star_L, samples_B_phi_n_U, samples_B_phi_p_U)
    samples_g_U <- g(samples_G_beta_star_U, samples_B_phi_n_L, samples_B_phi_p_L)
  }
  
  confidence_limit_L <- quantile(samples_g_L, prob = alpha / 2)
  confidence_limit_U <- quantile(samples_g_U, prob = 1 - alpha / 2)
  
  ci <- c(confidence_limit_L, confidence_limit_U)
  AP <- g(s_w_theta_hat, cn / mn, cp / mp)
  Sp <- 1 - cn / mn
  Se <- cp / mp
  
  estimate <- c("adjusted prevalence" = g(s_w_theta_hat, cn / mn, cp / mp))
  data <- paste0("Unadjusted prevalence=", s_w_theta_hat)
  statistic <- Se
  names(statistic) <- paste0("Sensitivity (using nSe=", mp, ")")
  parameter <- Sp
  names(parameter) <- paste0("Specificity (using nSp=", mn, ")")
  method_text <- paste0(
    "Prevalence Adjusted for Sensitivity and Specificity (CI by ",
    ifelse(method == "binomial", "Korn-Graubard", "Fay and Feuer"),
    " with melding)"
  )
  output <- list(
    estimate = estimate,
    statistic = statistic,
    parameter = parameter,
    conf.int = ci,
    data.name = data,
    method = method_text
  )
  class(output) <- "htest"
  output
}


# Example -----------------------------------------------------------------
