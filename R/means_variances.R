
#' Display value of an expression
#' 
#' Display value of an expression for debugging
#' 
#' @param x an expression
#' @export
disp <- function (x) 
{
  cat("\n::: ", deparse(substitute(x)), " :::\n")
  print(x)
}
#' Equilibrium expected value -- single cell type.
#' 
#' @param lambda autonomous rate of cell formation per unit of time
#' @param beta proportion of cells that replicate per unit of time
#' @param delta proportion of cells that die per unit of time
#' @return the equilibrium expected cell count
#' @examples
#' ceq1(10, .01, .02) 
#' @export
ceq1 <- function( lambda, beta, delta){
  theta <- delta - beta
  alpha <- delta + beta
  ret <- lambda / theta
  attr(ret, 'sd') <-  sqrt(lambda * (theta + alpha) / ( theta^2 * (2 - theta)))
  ret
}
#' Simulate single cell process
#' 
#' @param N number of observations
#' @param lambda autonomous rate of cell formation per unit of time
#' @param beta proportion of cells that replicate per unit of time
#' @param delta proportion of cells that die per unit of time
#' @param burnin starting the process at the equilibrium value, discard this many observations before generating the observation to be retured
#' @return N consecutive observations from the process
#' @examples
#' x <- sim1(10000, 10, .01, .02)
#' mean(x)
#' sd(x)
#' plot(x,pch='.')
#' @export
sim1 <- function(N, lambda, beta, delta, burnin = 1000){
  start <- ceq1(lambda, beta, delta)
  if(start <= 0) stop('series not stationary')
  x <- rep(NA,N+burnin)
  x[1] <- round(start)
  for ( i in 2:(N+burnin)) {
    n <- rpois(1,lambda)
    b <- rpois(1,beta*x[i-1])
    d <- rpois(1,delta*x[i-1])
    x[i] <- x[i-1] + n + b - d
  }
  attr(x,'pars') <- list(lambda=lambda,beta=beta,delta=delta,ex=start)
  x
}
#' Conditional expectation of next step of process
#' 
#' @param Y vector of cell counts at time t
#' @param lambda autonomous rate of cell formation for first cell type per unit of time
#' @param birth vector for each type of cell with proportions of cells that replicate per unit of time
#' @param death vector for each type of cell with proportions of cells that die per unit of time
#' @param tran vector of length equal to length(Y) -1 with proportions of cells that transform to the next cell type per unit of time
#' @param parms option list with parameters lambda, birth, death and tran
#' @param verbose (default FALSE) print additional output for debugging
#' @return conditional expectation at time t+1
#' @examples
#' cond_ev(c(100,100,100), 10, .01, .03, c(.01,.02))
#' parms <- list(lambda = 10, birth = rep(.01,3),
#'               death = .03, tran = c(.01,.02))
#' cond_ev(c(100,100,100), parms = parms)
#' 
#' Nsteps <- 2000
#' mat <- matrix(NA, 3, Nsteps)
#' mat[,1] <- c(100,100,100)
#' for ( i in 2:Nsteps) mat[,i] <- cond_ev(mat[,i-1],parms = parms)
#' # Equilibrium profile with these parameters:
#' head(t(mat),10)
#' tail(t(mat),10)
#' @export
cond_ev <- function(Y, lambda = parms$lambda, 
                    birth = parms$birth, 
                    death = parms$death,
                    tran = parms$tran,
                    parms, #' optional way of providing parameters
                    verbose = FALSE){
  # Returns: E( Y_{t+1} | Y_t)
  N <- length(Y)
  # birth, death and trans can be entered as single value
  birth <- rep(birth, length.out = N)
  death <- rep(death, length.out = N)
  if (N > 1) tran <- rep(tran, length.out = N-1)
  # net death rate
  theta <- numeric(N)
  if(N == 1) theta[1] <- death[1] - birth[1] 
  if(N > 1) {
    theta[1] <- death[1] - birth[1] + tran[1] 
    if (N > 2) {
      for (i in 2:(N-1)) theta[i] <- 
          death[i] + tran[i] - birth[i]
    }
    theta[N] <- death[N] - birth[N]
  }
  if(verbose) {
    disp(lambda)
    disp(birth)
    disp(death)
    disp(tran)
    disp(theta)
  }
  if(N == 1) L <- 1 - theta[1] else {
    L <- diag( 1 - theta)
    L[row(L)==(col(L)+1)] <- tran
  }
  if(verbose) disp(L)
  if(N == 1) return( lambda + L * Y)
  c(lambda,rep(0,N-1)) + L %*% Y
}
#' Mariginaal (equilibrium) expectation of process
#' 
#' @param lambda autonomous rate of cell formation for first cell type per unit of time
#' @param birth vector for each type of cell with proportions of cells that replicate per unit of time
#' @param death vector for each type of cell with proportions of cells that die per unit of time
#' @param tran vector of length equal to length(Y) -1 with proportions of cells that transform to the next cell type per unit of time
#' @param parms option list with parameters lambda, birth, death and tran
#' @param verbose (default FALSE) print additional output for debugging
#' @return marginal, i.e. long run, expectation of process
#' @examples
#' parms <- list(lambda = 10, birth = rep(.01,3),
#'               death = .03, tran = c(.01,.02))
#' marg_ev(parms=parms)
#' # comparing with asymptotic conditional expectation:
#' Nsteps <- 2000
#' mat <- matrix(NA, 3, Nsteps)
#' mat[,1] <- c(100,100,100)
#' for ( i in 2:Nsteps) mat[,i] <- cond_ev(mat[,i-1],parms = parms)
#' # Equilibrium profile with these parameters:
#' head(t(mat),10)
#' tail(t(mat),10)
#' @export
marg_ev <- function(lambda = parms$lambda, 
                    birth = parms$birth, 
                    death = parms$death,
                    tran = parms$tran,
                    parms, # optional way of providing parameters
                    verbose = FALSE) {
  # Returns: E( Y ) at equilibrium
  N <- length(birth)  #' must have correct length
  # death and trans can be entered as single value
  death <- rep(death, length.out = N)
  if (N > 1) tran <- rep(tran, length.out = N-1)
  # net death rate
  theta <- numeric(N)
  if(N == 1) theta[1] <- death[1] - birth[1]
  if(N > 1) {
    theta[1] <- death[1] - birth[1] + tran[1]
    if (N > 2) {
      for (i in 2:(N-1)) theta[i] <- 
          death[i] + tran[i] - birth[i]
    }
    theta[N] <- death[N] - birth[N]
  }
  ret <- numeric(N)
  ret[1] <- lambda / theta[1]
  if ( N > 1) for ( i in 2:N) ret[i] <- 
    ret[i-1] * tran[i-1]/theta[i]
  ret
}
#' Conditional variance
#' 
#' @param Y vector of cell counts at time t
#' @param lambda autonomous rate of cell formation for first cell type per unit of time
#' @param birth vector for each type of cell with proportions of cells that replicate per unit of time
#' @param death vector for each type of cell with proportions of cells that die per unit of time
#' @param tran vector of length equal to length(Y) -1 with proportions of cells that transform to the next cell type per unit of time
#' @param parms option list with parameters lambda, birth, death and tran
#' @param verbose (default FALSE) print additional output for debugging
#' @return conditional expectation at time t+1
#' @examples
#' parms <- list(lambda = 10, birth = rep(.01,3),
#'               death = .03, tran = c(.01,.02))
#' cond_ev(c(100,100,100), 10, .01, .03, c(.01,.02))
#' cond_var(c(100,100,100), parms = parms, verbose = T)  
#' @export
cond_var <- function(Y, 
                     lambda = parms$lambda, 
                     birth = parms$birth, 
                     death = parms$death,
                     tran = parms$tran,
                     parms, # optional way of providing parameters
                     verbose = FALSE) {
  # Returns: Var( Y_{t+1} | Y_t) = T D T'
  # where D is diagonal matrix of variances of new cell births, deaths and transitions
  N <- length(Y)
  # birth, death and trans can be entered as single value
  birth <- rep(birth, length.out = N)
  death <- rep(death, length.out = N)
  if (N > 1) tran <- rep(tran, length.out = N-1)
  # D
  D <- numeric(3*N)
  D <- c(lambda, birth[1]*Y[1], death[1]*Y[1])
  if(N > 1) {
    for (i in 2:N) {
      D[3*i-2] <- tran[i-1]*Y[i-1]
      D[3*i-1] <- birth[i]*Y[i]
      D[3*i] <- death[i]*Y[i]
    }
  }
  if(verbose) disp(D)
  # L matrix
  L <- matrix(0,N,3*N+1) # will drop last column
  for (i in 1:N) L[i,(3*i-2):(3*i+1)] <- c(1,1,-1,-1)
  L <- L[,-(3*N+1)] # drop last column
  if(verbose) disp(L)
  L %*% ( D * t(L))
}
#' Eigenvector of a lower triangular matrix
#' 
#' @param A a lower triangular matrix
#' @return the eigenvectors of A in the order corresponding the 
#' eigenvalues as they apper on the diagonal of A
#' @examples
#' A <- matrix(rnorm(25), 5,5)
#' A[col(A) > row(A)] <- 0
#' A
#' evec_lt(A) %*% diag(diag(A))- A %*% evec_lt(A)   # this should be machine 0
#' @export
evec_lt <- function(A){
  #' eigenvectors of a lower triangular matrix
  #' return E where A = E %*% diag(A) %*% E^{-1}
  #' Note that diag(A) is the diagonal matrix of eigenvectors (not ordered)
  N <- nrow(A)
  ret <- diag(N)
  for ( i in (N-1):1) {
    inds <- (i+1):N
    ret[inds,i] <- solve(A[inds,inds]-A[i,i]*diag(N-i),-A[inds,i])
  }
  ret
}
#' Solution to variance equation using diagonalization of lower triangular matrix
#' 
#' @param A,L matrices for which to find solution to \eqn{V = A +LVL'}
#' @return \eqn{v} where \eqn{V = A +LVL'} if a solution exists
#' @examples
#' A <- diag(3) + .1
#' A
#' L <- cbind( c(.8,.1,.1), c(0, .7, .1), c(0,0,.6))
#' L
#' V <- var_equation(A,L)
#' V
#' V - A - L%*% V %*% t(L)  #' should be machine 0
#' @export
var_equation <- function(A,L){
  #' Returns V such that V = A + LVL' if a solution exists ??? (check conditions)
  #' likely to return nonsense otherwise
  #' This does for matrices what the formula for a geometric series does for scalars
  G <- evec_lt(L)  #' this might not work although it probably should for reasonable
  #' parameters for T-cells  
  D <- diag(L)
  Ginv <- solve(G)
  #' disp( L - G %*% (D* Ginv)) #' ok
  B <- Ginv %*% A %*% t(Ginv)
  K <- B / (1 -outer(D,D))
  G %*% K %*% t(G)
}
#' Solution to variance equation using vec
#' 
#' Solves for \eqn{V} in the equation \eqn{V = A +LVL'} using the fact that
#' \eqn{vec(V) = vec(A) + (L \kronecker L) vec(V)}. This may work when 
#' \code{\link{var_equation}} does not because \eqn{L} is not diagonalizable.
#' 
#' @param A,L matrices for which to find solution to \eqn{V = A +LVL'}
#' @return \eqn{v} where \eqn{V = A +LVL'} if a solution exists
#' @examples
#' A <- diag(3) + .1
#' A
#' L <- cbind( c(.8,.1,.1), c(0, .7, .1), c(0,0,.6))
#' L
#' V <- var_equation_vec(A,L)
#' V
#' V - t(V)
#' V - A - L%*% V %*% t(L)  # should be machine 0
#' V - var_equation(A,L)
#' # L that can't be diagonalized
#' L <- cbind( c(.9,.1,0), c(0, .9, .1), c(0,0,.9))
#' var_equation(A,L) 
#' var_equation_vec(A,L)
#' eigen(var_equation_vec(A,L))
#' @export
var_equation_vec <- function(A,L){
  #' Returns V such that V = A + LVL' if a solution exists ??? (check conditions)
  #' likely to return nonsense otherwise
  N <- dim(A)[1]
  K <- N^2
  v <- solve( diag(K) - kronecker(L,L), c(A))
  matrix(v,N,N)
}
#' 
#' Mean conditional variance
#' 
#' @inheritParams cond_var
#' @return \eqn{E ( Var( Y_{t+1} | Y_t)) = T E(D) T'} where 
#' \eqn{D} is the diagonal matrix of variances of new cell births, 
#' deaths and transitions and \eqn{T} is the coding matrix of 0s and \eqn{\pm 1}s
#' to linearly transform the birth, death and transition components into the 
#' cell counts
#' @examples
#' parms <- list(lambda = 10, birth = rep(.01,3),
#'               death = .03, tran = c(.01,.02))
#' mean_cond_var(parms = parms) 
#' @export
mean_cond_var <- function( 
  lambda = parms$lambda, 
  birth = parms$birth, 
  death = parms$death,
  tran = parms$tran,
  parms, #' optional way of providing parameters
  verbose = FALSE) {
  # Returns: E ( Var( Y_{t+1} | Y_t)) = T E(D) T'
  # where D is diagonal matrix of variances of new cell births, deaths and transitions
  EY <- marg_ev(lambda=lambda,birth=birth,death=death,tran=tran) 
  cond_var(EY,lambda=lambda,birth=birth,death=death,tran=tran) 
}
#' Marginal variance
#' 
#' @inheritParams cond_var
#' @return marginal long-term variance of the process 
#' @examples
#' parms <- list(lambda = 10, birth = rep(.01,3),
#'               death = .03, tran = c(.01,.02))
#' marg_var(parms = parms)
#' @export
marg_var <- function(
  lambda = parms$lambda, 
  birth = parms$birth, 
  death = parms$death,
  tran = parms$tran,
  parms, #' optional way of providing parameters
  verbose = FALSE) {
  #' Returns: Var (Y)
  N <- length(birth)  #' must have correct length
  #' death and trans can be entered as single value
  death <- rep(death, length.out = N)
  if (N > 1) tran <- rep(tran, length.out = N-1)
  #' net death rate
  theta <- numeric(N)
  if(N == 1) theta[1] <- death[1] - birth[1]
  if(N > 1) {
    theta[1] <- death[1] - birth[1] + tran[1]
    if (N > 2) {
      for (i in 2:(N-1)) theta[i] <- 
          death[i] + tran[i] - birth[i]
    }
    theta[N] <- death[N] - birth[N]
  }
  P <- mean_cond_var(lambda=lambda,birth=birth,death=death,tran=tran)
  if(N == 1) return(P/(2*theta[1] - theta[1]^2))
  L <- diag(1 - theta)
  L[row(L)==(col(L)+1)] <- tran
  var_equation_vec(P, L)
}
#' Simulation step
#' 
#' @inheritParams cond_ev
#' @return randomly generated value at time t+1
#' @examples
#' parms <- list(lambda = 10, birth = rep(.01,3),
#'               death = .03, tran = c(.01,.02))
#' sim_step(c(100,100,100), parms = parms)
#' @export
sim_step <- function(Y, 
                     lambda = parms$lambda, 
                     birth = parms$birth, 
                     death = parms$death,
                     tran = parms$tran,
                     parms, #' optional way of providing parameters
                     verbose = FALSE) {
  N <- length(Y)  #' must have correct length
  #' death and trans can be entered as single value
  birth <- rep(birth, length.out = N)
  death <- rep(death, length.out = N)
  if (N > 1) tran <- rep(tran, length.out = N-1)
  D <- numeric(3*N)
  D <- c(lambda, birth[1]*Y[1], death[1]*Y[1])
  if(N > 1) {
    for (i in 2:N) {
      D[3*i-2] <- tran[i-1]*Y[i-1]
      D[3*i-1] <- birth[i]*Y[i]
      D[3*i] <- death[i]*Y[i]
    }
  }
  if(verbose) disp(D)
  #' L matrix
  L <- matrix(0,N,3*N+1) #' will drop last column
  for (i in 1:N) L[i,(3*i-2):(3*i+1)] <- c(1,1,-1,-1)
  L <- L[,-(3*N+1)] #' drop last column
  if(verbose) disp(L)
  Diffs <- rpois(length(D), lambda = D)
  Y + L %*% Diffs
}
#' Simulation sequence
#' 
#' @param N number of steps to simulate
#' @inheritParams cond_ev
#' @return sequence of randomly generated values at times t+1, ..., t+N
#' @examples
#' \dontrun{
#' parms <- list(lambda=10, birth = c(.01,.02,.01), death = .03, tran = c(.02, .04))
#' matplot(sim_seq(N = 50, parms = parms), type = 'l', lty = c(1,2,3), lwd = 2, col = c('black','red','green'))
#' legend("right", inset=.05, legend=c("1", "2", "3"), lty = c(1,2,3), lwd = 2, col = c('black','red','green') , horiz=TRUE)
#' matplot(sim_seq(N = 1000, parms = parms), type = 'l', lty = c(1,2,3), lwd = 2, col = c('black','red','green'))
#' legend("right", inset=.05, legend=c("1", "2", "3"), lty = c(1,2,3), lwd = 2, col = c('black','red','green') , horiz=TRUE)
#' }
#' @export
sim_seq <- function(N, Y = marg_ev(lambda,birth,death,tran), 
                    lambda = parms$lambda, 
                    birth = parms$birth, 
                    death = parms$death,
                    tran = parms$tran,
                    parms, #' optional way of providing parameters
                    verbose = FALSE) {
  ret <- matrix(0,nrow = N, ncol = length(Y))
  ret[1,] <- Y
  if (N > 1) for ( i in 2:N) ret[i,] <- sim_step(ret[i-1,],
                                                 lambda = lambda, 
                                                 birth = birth,
                                                 death = death,
                                                 tran = tran)
  ret
}


