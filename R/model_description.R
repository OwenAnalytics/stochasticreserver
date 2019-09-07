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
