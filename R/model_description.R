#' Get the long description of a model
#'
#' @return Long description of a model as a character vector of length one.
#' @param model character vector of length one having the short name of the
#' model.
#' @export
model_description <- function(model) {
  if (model == "Berquist") {
    "Berquist-Sherman Incremental Severity"
  } else if (model == "CapeCod") {
    "Cape Cod"
  } else if (model == "Hoerl") {
    "Generalized Hoerl Curve Model with Trend"
  } else if (model == "Wright") {
    "Generalized Hoerl Curve with Individual Accident Year Levels"
  } else if (model == "CapeCod") {
    "Chain Ladder Model"
  }
}
