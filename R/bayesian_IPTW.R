#' Functions for Constructing Inverse Probability of Treatment Weights for Bayesian Panel Models
#'
#' Functions for estimating the posterior distribution of stabilized inverse probability of
#' treatment weights for marginal structural models applied to time series cross-sectional
#' data.
#'
#' @aliases IPTW_Weights
#' @importFrom Rdpack reprompt
#'
#' @param num_model A required argument supplying either a model of class \code{brmsfit}
#' or an \eqn{m \times n} matrix of posterior expectations where \eqn{m} denotes draws
#' from the posterior predictive distribution and \eqn{n} is the number of observations
#' in the data. This argument corresponds to the numerator of the stabilized weights.
#' @param denom_model A required argument supplying either a model of class \code{brmsfit}
#' or an \eqn{m \times n} matrix of posterior expectations where \eqn{m} denotes draws
#' from the posterior predictive distribution and \eqn{n} is the number of observations
#' in the data. This argument corresponds to the denominator of the stabilized weights
#' and must be of same class as \code{num_model}.
#' @param method An argument of type character indicating whether \code{num_model} and
#' \code{denom_model} are \code{matrix} or \code{brmsfit} objects. Defaults to \code{"epred"}
#' which corresponds to \code{brmsfit} objects but can be set to \code{"average"} when
#' passing matrices of posterior predictions which allows for the use of model averaged
#' or stacked posterior predictive distributions.
#' @param response A required argument specifying a numeric vector of length \eqn{n}
#' containing the data for the response used in the treatment propensity model.
#' @param id A required argument of length \eqn{n} containing the grouping structure
#' of the cross-sectional dimension of the data (i.e., states, countries, etc.)
#' @param ts A required argument of length \eqn{n} containing the temporal structure
#' of the time series dimension of the data (i.e., years, months, etc.)
#' @param trim A numeric value by which to adjust the predictions before
#' calculating the weights. Defaults to \code{FALSE} which does not make any
#' adjustments to the posterior predictions. This can be useful in the case of
#' predicted treatment probabilities near 0 or 1.
#' @param log_weights Whether to return the summarized weights in raw or log transformed
#' form. This can be useful when the stabilized weights exhibit a high degree of variance.
#' Defaults to \code{FALSE} which just returns the weights in their raw form. If
#' \code{TRUE}, summary estimates are transformed by calling \code{log1p}.
#' @param ... Additional arguments to be passed to \code{brms::posterior_epred}.
#' Only used if method = "epred" and ignored otherwise.
#'
#' @return The function returns a list containing a matrix of the full posterior
#' distribution of the inverse probability of treatment weights and by-observation
#' summary statistics.
#'
#' @details Following \insertCite{Blackwell2018}{PoliSciTools}, the stabilized inverse probability weights
#' for the average treatment effect (ATE) of a binary treatment can be expressed as
#' \deqn{\hat{SW}_{i[t]} = \prod^{t}_{t = 1} \frac{\Pr[X_{i[t]} | \bar{X}_{i[t-1]}, V_i]}{\Pr[X_{i[t]} | \bar{X}_{i[t-1]}, Y_{i[t-1]}, C_{i[t]}, V_{[i]}]}}
#'
#' @references
#'     \insertAllCited{}
#'
#' @importFrom brms is.brmsfit posterior_epred
#' @importFrom dplyr mutate group_by across starts_with ungroup
#' @importFrom tidyr replace_na
#' @importFrom tibble as_tibble
#' @importFrom stats median nobs sd
#'
#'
#' @export IPTW_Weights
#' @export
IPTW_Weights <- function(num_model,
                         denom_model,
                         method = c("epred", "average"),
                         response = NULL,
                         id,
                         ts,
                         trim = FALSE,
                         log_weights = FALSE,
                         ...) {

  ## Get the method for generating predictions
  .method <- method[1] ### Defaults to epred

  ## Check that method is correctly specified
  stopifnot(isTRUE(.method == "epred" || .method == "average"))

  ## Initiate a list object to store the results in
  weights <- list()

  ## Check if both arguments are matrices if method is "average"
  if (isTRUE(.method == "average")) {

    # Stop if method is "average" and both models are not matricies
    stopifnot(exprs = {
      is.matrix(num_model) && is.matrix(denom_model)
      is.vector(response) && is.vector(id) && is.vector(ts)
      all.equal(
        length(response),
        ncol(num_model),
        ncol(denom_model),
        length(id),
        length(ts)
      )
    })

    # Transpose the matrix of posterior expectations for the denominator model
    preds_denom <- denom_model

    # Transpose the matrix of posterior expectations for the numerator model
    preds_num <- num_model

  }

  ## Check if both arguments are brmsfit if method is "epred"
  else if (isTRUE(.method == "epred")) {

    # Stop if method is "epred" and both models are not of class brmsfit
    stopifnot(exprs = {
      brms::is.brmsfit(num_model) && brms::is.brmsfit(denom_model)
      is.vector(response) && is.vector(id) && is.vector(ts)
      all.equal(
        length(response),
        nobs(num_model),
        nobs(denom_model),
        length(id),
        length(ts)
      )
    })

    # Generate posterior expectations denominator model
    preds_denom <- brms::posterior_epred(denom_model, ...)

    # Generate posterior expectations for the numerator model
    preds_num <- brms::posterior_epred(num_model, ...)

  }

  ## Check that the lengths are the same
  stopifnot(all.equal(dim(preds_denom), dim(preds_num)))

  ## Extract the number of observations
  n_draws <- 1:nrow(preds_num)

  # Fill in the matrix with the adjusted predictions
  if(isTRUE(is.numeric(trim))) {
    for (i in seq_along(n_draws)) {
      preds_num[i, ] <- pmax(pmin(preds_num[i, ], 1 - trim), trim)
      preds_denom[i, ] <- pmax(pmin(preds_denom[i, ], 1 - trim), trim)
    }
  }

  # Calculate a vector of weights for each posterior draw from the models
  for (i in seq_along(n_draws)) {
    preds_num[i, ] <- ifelse(response[i] == 1, preds_num[i, ], 1 - preds_num[i, ])
    preds_denom[i, ] <- ifelse(response[i] == 1, preds_denom[i, ], 1 - preds_denom[i, ])
  }

  # Generate a matrix to store the IPT weights in
  weights[["Weights_Matrix"]] <- matrix(data = NA_real_, ncol = ncol(preds_denom), nrow =  nrow(preds_denom))

  # Calculate a vector of weights for each posterior draw
  for (i in seq_along(n_draws)) {
    weights[["Weights_Matrix"]][i, ] <- preds_num[i, ]/preds_denom[i, ]
  }

  ## Coerce the output to a tibble
  weights[["Weights_Matrix"]] <- tibble::as_tibble(t(weights[["Weights_Matrix"]]), .name_repair = "minimal")

  ## Assign column names to each draw
  colnames(weights[["Weights_Matrix"]]) <- paste("draw_", 1:ncol(weights[["Weights_Matrix"]]), sep = "")

  weights[["Weights_Matrix"]] <- weights[["Weights_Matrix"]] |>
    ## Add group and time columns
    dplyr::mutate(ts = ts, id = id, .before = 1) |>
    ## Group the data by identifier
    dplyr::group_by(id) |>
    ## Calculate the cumulative product of the weights by country
    dplyr::mutate(dplyr::across(
      dplyr::starts_with("draw"),
      ~ cumprod(tidyr::replace_na(.x, 1)))
    ) |>
    ## Ungroup the Data
    ungroup()

  ## If log weights is true, take the logarithm of the stabilized weights
  if(isTRUE(log_weights)) {
    ## Store the location  and scale of the weights matrix
    weights[["Summary"]] <- list(
      weights_log_mean = log1p(rowMeans(weights[["Weights_Matrix"]][, -c(1:2)])),
      weights_log_median = log1p(apply(weights[["Weights_Matrix"]][, -c(1:2)], 1, median)),
      weights_log_sd = log1p(apply(weights[["Weights_Matrix"]][, -c(1:2)], 1, sd))
    )
  }

  ## Otherwise, return the raw summary statistics
  else {
    ## Store the location  and scale of the weights matrix
    weights[["Summary"]] <- list(
      weights_mean = rowMeans(weights[["Weights_Matrix"]][, -c(1:2)]),
      weights_median = apply(weights[["Weights_Matrix"]][, -c(1:2)], 1, median),
      weights_sd = apply(weights[["Weights_Matrix"]][, -c(1:2)], 1, sd)
    )
  }

  # Return the weights list object
  return(weights)
}
