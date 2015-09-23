#' Tcells2: an R package to estimate cell dynamics with mixed models.
#' 
#' Collection of functions for analysis of Tcell dynamics
#'
#' @section Functions:
#' \describe{
#'    \item{\code{\link{cond_ev}}}{Conditional expected value at time t+1 given value at time t}
#'    \item{\code{\link{cond_var}}}{Conditional variance at time t+1 given value at time t}
#'    \item{\code{\link{marg_ev}}}{Marginal expected value over time}
#'    \item{\code{\link{marg_var}}}{Marginal variance over time}
#'    \item{\code{\link{mean_cond_var}}}{Mean conditional variance over time}
#'    \item{\code{\link{mean_cond_var}}}{Mean conditional variance over time}
#'    \item{\code{\link{var_equation}}}{Solves for \eqn{V} in \eqn{ V = A + LVL'}
#'    where $L$ is lower triangular using the diagonalization of $L$}
#'    
#' }
#' @export
Tcells2 <- function() vignette("Tcells2")
