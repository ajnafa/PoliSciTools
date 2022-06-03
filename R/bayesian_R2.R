#' Marginal and Conditional Bayes R-Square for Bayesian Models
#'
#' This function returns estimates for both the marginal and
#' conditional Bayes \eqn{R^2} for multilevel models fit with
#' either \pkg{brms} or \pkg{rstanarm} and stores them
#' in the original model object.
#'
#' @aliases Bayesian_R2
#' @importFrom Rdpack reprompt
#'
#' @param model A required argument supplying the \code{brmsfit}
#' or \code{stanreg} model object for which to calculate
#' Bayes \eqn{R^2}
#' @param marginal Logical argument indicating whether to calculate
#' the unconditional Bayes \eqn{R^2}. Defaults to \code{TRUE}
#' @param conditional Logical argument indicating whether to calculate
#' Bayes \eqn{R^2} conditional on the random effects. Defaults to
#' \code{FALSE}
#' @param effects A one-sided formula specifying the structure of
#' the random effects for the conditional Bayes \eqn{R^2}. Defaults
#' to \code{NULL}, which conditions on the full random effects structure.
#' @param save Logical argument indicating whether the results should be
#' saved to the local model file. Defaults to \code{FALSE} and is
#' currently only supported for \code{brmsfit} objects.
#' @param path The base directory for where the updated model object
#' should be saved. Defaults to the current working directory and is
#' ignored if save is \code{FALSE}.
#' @param overwrite Logical argument indicating whether to overwrite any
#' existing \code{Bayesian_R2} result saved in model$criteria. Defaults to
#' \code{FALSE}, which for \code{brmsfit} objects will simply return
#' the original model object supplied.
#' @param ... Additional arguments to be passed to either
#' \code{brms::bayes_R2} or \code{rstanarm::bayes_R2} depending on
#' the model class.
#'
#' @return The function returns the original model object with the
#' results stored in the \code{model$criteria$Bayesian_R2}.
#'
#' @details For additional details on the underlying calls see
#' \code{?brms::bayes_R2} or \code{?rstanarm::bayes_R2}. A technical
#' discussion of the Bayes \eqn{R^2} statistic can be found in
#' \insertCite{Gelman2019}{PoliSciTools}
#'
#' @references
#'     \insertAllCited{}
#'
#' @importFrom brms bayes_R2
#' @importFrom rstanarm bayes_R2
#' @export Bayesian_R2
#' @export
Bayesian_R2 <- function(model,
                        marginal = TRUE,
                        conditional = FALSE,
                        effects = NULL,
                        save = FALSE,
                        path = getwd(),
                        overwrite = FALSE,
                        ...
) {

  ### Initiate a list to store the results in
  out <- list()

  ## Check the class of the model object
  .mod_type <- class(model)

  ## Check that at least one of the arguments is specified
  .check_args <- ifelse(isTRUE(marginal) || isTRUE(conditional), 1, 0)

  ## Otherwise, return an error
  if (.check_args < 1) {
    stop("At least one of the 'conditional' or 'marginal' arguments must be TRUE")
  }

  ## If model is a brmsfit object call `brms::bayes_R2`
  else if (!is.na(match(.mod_type, "brmsfit"))) {

    ### Check if Bayesian_R2 is already stored to avoid redundant computation
    R2 <- .check_criteria(model, criterion = "Bayesian_R2", overwrite = overwrite)

    ### If Bayesian_R2 is already stored and overwrite is false, just return the existing object
    if (!is.null(R2)) {
      model$criteria$Bayesian_R2 <- R2
      return(model)
    }

    ## Otherwise, evaluate the arguments to the function
    else {

      ### Get the file path
      file <- model$file

      ### Marginal Bayes R2 for brmsfit objects
      if (isTRUE(marginal)) {

        ### Return the Marginal Bayes R2 to the temporary list object
        out$marginal_R2 <- brms::bayes_R2(model, re_formula = NA, ...)

        ### Assign a names attribute
        attr(out$marginal_R2, "names") <- "Marginal Bayes R2"

      }

      ### Conditional Bayes R2 for brmsfit objects
      if (isTRUE(conditional)) {

        ### Add Conditional Bayes R2 to the temporary list object
        out$conditional_R2 <- brms::bayes_R2(model, re_formula = effects, ...)

        ### Assign a names attribute
        attr(out$conditional_R2, "names") <- "Conditional Bayes R2"

      }

      ## Add the result to the original model object
      model$criteria$Bayesian_R2 <- out

      ## If save is true, save the model with the added criteria to the original file
      if (isTRUE(save)) {

        # Write using high file compression
        saveRDS(
          model,
          file = paste(path, file, sep = "/"),
          compress = TRUE
        )

        # Print the location the model was saved to
        print(paste("Saving the updated model object to ", file, sep = ""))
      }

      ## Return the updated model
      return(model)
    }
  }

  ## If model is a stanreg object call `rstanarm::bayes_R2`
  else if (!is.na(match(.mod_type, "stanreg"))) {

    ### Marginal Bayes R2 for brmsfit objects
    if (isTRUE(marginal)) {

      ### Return the Marginal Bayes R2 to the temporary list object
      out$marginal_R2 <- rstanarm::bayes_R2(model, re.form = NA, ...)

      ### Assign a names attribute
      attr(out$marginal_R2, "names") <- "Marginal Bayes R2"

    }

    ### Conditional Bayes R2 for brmsfit objects
    if (isTRUE(conditional)) {

      ### Add Conditional Bayes R2 to the temporary list object
      out$conditional_R2 <- rstanarm::bayes_R2(model, re.form = effects, ...)

      ### Assign a names attribute
      attr(out$conditional_R2, "names") <- "Conditional Bayes R2"

    }

    ## Add the result to the original model object
    model$criteria$Bayesian_R2 <- out

    ## Return the updated model
    return(model)
  }

  ## Return and error if model is not of class brmsfit or stanreg
  else {
    stop("This function only supports model objects of class brmsfit or stanreg")
  }
}

## Function to check whether model criterion is already stored
.check_criteria <- function (x, criterion, overwrite) {

  ## Check the class of the model object
  stopifnot(class(x) == "brmsfit")

  ## Retrive the names of the already stored criteria
  criteria <- names(x$criteria)

  ## Position of the stored criteria
  criteria_pos <- match(criterion, criteria)

  ## If overwrite is FALSE and criteria is found, return existing criteria
  if (!is.na(criteria_pos) & overwrite == FALSE) {
    return(x$criteria[[criteria_pos]])
  }
}