#' Round the numeric columns in a data.frame
#'
#' @param df a data.frame
#' @param digits the number of digits to round to
#'
#' @return a data.frame
#'
#' @examples
#' x <- matrix(runif(10), nrow=10, ncol=10)
#' x <- as.data.frame(x)
#' y <- roundDf(x, 2)
#''
#' @export
round_df <- function(df, digits=5) {
  nums <- vapply(df, is.numeric, FUN.VALUE=logical(1))
  df[,nums] <- round(df[,nums], digits=digits)
  return(df)
}
