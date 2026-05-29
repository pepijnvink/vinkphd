#' Compute AIC for a `mHMM_list` object
#'
#' @param models mHMM_list opject containing multiple chains of a multilevel hidden markov model.
#'
#' @returns Named vector containing loglikelihood, AIC, and AICc
#'
#' @export
mHMM_AIC <- function(models) {
  if (!inherits(models, 'mHMM_list')) {
    cli::cli_abort(c(
      "{.var models} must be an object of class {.cls mHMM_list}",
      "x" = "You provided an object of class {.cls {class(models)}}"
    ))
  }
  input <- models[[1]]$input
  n_subj <- input$n_subj
  burn_in <- input$burn_in
  J <- input$J
  m <- input$m
  n_vary <- input$n_vary
  data_distr <- input$data_distr
  n_dep <- input$n_dep

  LL <- numeric(n_subj)
  for (s in 1:n_subj) {
    LL[s] <- lapply(models, \(x) {
      x$PD_subj[[s]]$log_likl[((burn_in + 1):J), 1]
    }) |>
      unlist() |>
      stats::median()
  }
  if (data_distr == 'categorical') {
    q_emiss <- input$q_emiss
    n_par <- sum((q_emiss - 1) * m) + (m - 1) * m
    AIC <- 2 * n_par - (2 * LL)
    AICc <- ((2 * n_vary * n_par) / (n_vary - n_par - 1)) - (2 * LL)
  } else if (data_distr == 'continuous') {
    n_par <- m * n_dep * 2 + (m - 1) * m
    AIC <- 2 * n_par - (2 * LL)
    AICc <- ((2 * n_vary * n_par) / (n_vary - n_par - 1)) - (2 * LL)
  } else if (data_distr == 'count') {
    n_par <- m * n_dep + (m - 1) * m
    AIC <- 2 * n_par - (2 * LL)
    AICc <- ((2 * n_vary * n_par) / (n_vary - n_par - 1)) - (2 * LL)
  }
  output <- c(LL = mean(LL), AIC = mean(AIC), AICc = mean(AICc))
  return(output)
}
