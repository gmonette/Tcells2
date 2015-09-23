#' ---
#' title: "Conditional and Marginal Expectation and Variance"
#' date: "`r format(Sys.time(), '%H:%M %d %B %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      toc_depth: 4
#' ---
#' Generated:
#' 

MIGHT ABANDON IDEA OF VIGNETTE IN FAVOUR OF JUST COPYING STUFF TO 
HTML in blackwell



{{format(Sys.time(), '%H:%M %d %B %Y')}}
#' 
#' Derivation of function in Tcells package in file means_variances.R
#' 
#' <!--
#' * 2015 09 14 [GM] - initial version
#' * TO DO: report both Exp conditional variance and variance of conditional expectation separately
#' * Test variance algorithms more throroughly
#' -->
#' 
#+ opts,echo=FALSE
library(knitr)
opts_knit$set(width=120)
options(width=120)
opts_chunk$set(tidy=FALSE,comment='| ',fig.height=8,fig.width=10)
#'
#' The following are formulas and algorithms in R for the conditional and marginal expectations and variance of a T-cell process in which 
#' there is a hierarchy of cell types in which cells can transform from one type to the next type.
#' 
#' I am coding the algorithms in 'STAN', a probabilistic programming language for Bayesian inference which will allow an analysis in which
#' the process at the level of the population is treated as a hidden Markov process.  The resulting estimates, if the data is sufficient
#' for this approach to work, may have a number of advantages over the HLM based only on modelling observed concentrations.  
#' 
#' 1. The analysis will take into account that concentrations, when used as covariates, are measured with error with consequent biasing of estimates of
#'    regression coefficients,
#' 2. the HLM is implicitly based on the incorrect assumption that the conntrations are Markov thus ignoring potential information in concentrations
#'    measured prior to the previous occasion,
#' 3. when a previous concentration is missing, the regression of a concentration on its precursors is lost although there is potential information
#'    from prior concentrations,
#' 4. the analysis will combine both a 'static' model using equilibrium expectations and variances for the initial counts and the 'dynamic' model
#'    based on conditional expectations and variance for changes in population counts,
#' 5. the model is based on a linear approximation to what is more like an exponential process: the linear approximation may be valid if one models
#'    relatively short time intervals between measurement occasions.  The structural hidden Markov model for population counts can incorporate 
#'    additional occasions to 'bridge' situations in which there is long delay between measurement occasions.
#' 6. + others we'll think of.
#' 
#' The downside is that I'm not sure the method will work on this amount of data. That will take simulations to explore.
#' 
#' ## Expectation and variance for single cell type
#' 
#' \[\begin{aligned}
#'    {X_{t + 1}} =  & {X_t} + {N_t} + {B_t} - {D_t} \\ 
#'    {N_t} \sim & \operatorname{Poisson} (\lambda ) \\ 
#'    {B_t} \sim & \operatorname{Poisson} (\beta {X_t}) \\ 
#'    {D_t} \sim & \operatorname{Poisson} (\delta {X_t}) \\ 
#'    \operatorname{E} \left( {{X_{t + 1}}|{X_t}} \right) =  & \lambda  + \left( {1 + \beta  - \delta } \right){X_t} \\ 
#'    \operatorname{E} \left( {{X_{t + 1}}} \right) =  & \lambda  + \left( {1 + \beta  - \delta } \right)\operatorname{E} \left( {{X_t}} \right) \\ 
#'    \operatorname{E} (X) =  & \lambda  + \left( {1 + \beta  - \delta } \right)\operatorname{E} \left( X \right) \\ 
#'    \operatorname{E} (X) =  & \frac{\lambda }{{\delta  - \beta }} = \frac{\lambda }{\theta };\theta  = \delta  - \beta {\text{ is net death rate}} \\ 
#'    \operatorname{Var} \left( {{X_{t + 1}}|{X_t}} \right) =  & \lambda  + \left( {\beta  + \delta } \right){X_t} = \lambda  + \alpha {X_t};\alpha  = \beta  + \delta {\text{ is the 'activity rate'}} \\ 
#'    \operatorname{Var} \left( {{X_{t + 1}}} \right) =  & \operatorname{E} \left( {\lambda  + \alpha {X_t}} \right) + \operatorname{Var} \left( {\lambda  + \left( {1 - \theta } \right){X_t}} \right) \\ 
#'    =  & \lambda  + \alpha \operatorname{E} \left( {{X_t}} \right) + {\left( {1 - \theta } \right)^2}\operatorname{Var} \left( {{X_t}} \right) \\ 
#'    \operatorname{Var} \left( X \right) =  & \lambda  + \alpha \operatorname{E} \left( X \right) + {\left( {1 - \theta } \right)^2}\operatorname{Var} \left( X \right) \\ 
#'    \operatorname{Var} \left( X \right) =  & \frac{{\lambda  + \alpha \frac{\lambda }{\theta }}}{{1 - {{\left( {1 - \theta } \right)}^2}}} = \frac{{\frac{\lambda }{\theta }\left( {\theta  + \alpha } \right)}}{{2\theta  - {\theta ^2}}} = \frac{{\lambda \left( {\theta  + \alpha } \right)}}{{{\theta ^2}\left( {2 - \theta } \right)}} \\ 
#'    \end{aligned} \]
#'  
#'  

# Utility functions -----------------------------------------------------------------------------------

disp <- function (x) 
{
  cat("\n::: ", deparse(substitute(x)), " :::\n")
  print(x)
}

# Functions for a single cell type --------------------------------------------------------------------
ceq1 <- function( lambda, beta, delta){
  theta <- delta - beta
  alpha <- delta + beta
  ret <- lambda / theta
  attr(ret, 'sd') <-  sqrt(lambda * (theta + alpha) / ( theta^2 * (2 - theta)))
  ret
}

ceq1(10, .01, .02)

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

x <- sim1(10000, 10, .01, .02)
mean(x)
sd(x)
plot(x,pch='.')

ceq1(10,.01,.02)

# Simulate and test ---------------------------------------------------------------------------------
lam <- 10
beta <- .01
delta <- .02
ceq1(lam,beta,delta)
x <- sim1(10000, lam, beta, delta)
mean(x)
sd(x)

#'
#' ## Three cell types
#'
#' The following formulas and algorithms are general in that they can be applied to any number of cell types.
#' Three is a minimal number in order to represent the full model with 
#' 
#' 1. a 'Naive' type that are can be created autonomously at a rate that is independent of the current population size, 
#' 2. one or more 'transitory' types that can be transformed into and that can transform into the next type, and
#' 3. a 'terminal' type that does not transform.
#'
#' ### Definition of process
#'
#' \[\begin{aligned}
#'    {Y_{t + 1}} =  & \left[ {\begin{array}{*{20}{c}}
#'      {{X_{1,t + 1}}} \\ 
#'      {{X_{2,t + 1}}} \\ 
#'      {{X_{3,t + 1}}} 
#'      \end{array}} \right] = \left[ {\begin{array}{*{20}{c}}
#'        {{X_{1t}} + {N_t} + {B_{1t}} - {D_{1t}} - {T_{1t}}} \\ 
#'        {{X_{2t}} + {T_{1t}} + {B_{2t}} - {D_{2t}} - {T_{2t}}} \\ 
#'        {{X_{3t}} + {T_{2t}} + {B_{3t}} - {D_{3t}}} 
#'        \end{array}} \right] \\ 
#'    =  & {Y_t} + \left[ {\begin{array}{*{20}{c}}
#'      {{N_t} + {B_{1t}} - {D_{1t}} - {T_{1t}}} \\ 
#'      {{T_{1t}} + {B_{2t}} - {D_{2t}} - {T_{2t}}} \\ 
#'      {{T_{2t}} + {B_{3t}} - {D_{3t}}} 
#'      \end{array}} \right] \\ 
#'    {\text{Given }}{Y_t}: &  \\ 
#'    {N_t} \sim  & \operatorname{Poisson} \left( \lambda  \right) \\ 
#'    {B_{it}} \sim  & \operatorname{Poisson} \left( {{\beta _i}{X_{it}}} \right) \\ 
#'    {D_{it}} \sim  & \operatorname{Poisson} \left( {{\delta _i}{X_{it}}} \right) \\ 
#'    {T_{it}} \sim  & \operatorname{Poisson} \left( {{\tau _i}{X_{it}}} \right) \\ 
#'    \end{aligned} \]
#' 
#' In the following we define:
#' \[\begin{aligned}
#' {\theta _i} =  & {\delta _i} + {\tau _i} - {\beta _i} \\ 
#' {\alpha _i} =  & {\delta _i} + {\tau _i} + {\beta _i} \\ 
#' \end{aligned} \]
#' for $i = 1,...,I-1$, and
#' \[\begin{aligned}
#' {\theta _I} =  & {\delta _I} - {\beta _I} \\ 
#' {\alpha _I} =  & {\delta _I} + {\beta _I} \\ 
#' \end{aligned} \]
#' 
#' Note that $\theta_i$ and $\alpha_i$ represent the 'net conditional death rate' and a 'conditional activity rate' respectively.
#' 
#' ### Expected value of process
#' 
#' \[\begin{aligned}
#'    \operatorname{E} \left( {{Y_{t + 1}}|{Y_t}} \right) =  & \operatorname{E} \left( {\left. {\begin{array}{*{20}{c}}
#'      {{X_{1,t + 1}}} \\ 
#'      {{X_{2,t + 1}}} \\ 
#'      {{X_{3,t + 1}}} 
#'      \end{array}} \right|} \right.\left. {\begin{array}{*{20}{c}}
#'        {{X_{1t}}} \\ 
#'        {{X_{2t}}} \\ 
#'        {{X_{3t}}} 
#'        \end{array}} \right) = \operatorname{E} \left( {\left. {\begin{array}{*{20}{c}}
#'          {{X_{1t}} + {N_t} + {B_{1t}} - {D_{1t}} - {T_{1t}}} \\ 
#'          {{X_{2t}} + {T_{1t}} + {B_{2t}} - {D_{2t}} - {T_{2t}}} \\ 
#'          {{X_{3t}} + {T_{2t}} + {B_{3t}} - {D_{3t}}} 
#'          \end{array}} \right|} \right.\left. {\begin{array}{*{20}{c}}
#'            {{X_{1t}}} \\ 
#'            {{X_{2t}}} \\ 
#'            {{X_{3t}}} 
#'            \end{array}} \right) \\ 
#'    =  & \left[ {\begin{array}{*{20}{c}}
#'      \lambda  \\ 
#'      0 \\ 
#'      0 
#'      \end{array}} \right] + \left[ {\begin{array}{*{20}{c}}
#'        {1 - {\theta _1}}&0&0 \\ 
#'        {{\tau _1}}&{1 - {\theta _2}}&0 \\ 
#'        0&{{\tau _2}}&{1 - {\theta _3}} 
#'        \end{array}} \right]\left[ {\begin{array}{*{20}{c}}
#'          {{X_{1t}}} \\ 
#'          {{X_{2t}}} \\ 
#'          {{X_{3t}}} 
#'          \end{array}} \right] \\ 
#'    =  & \left[ {\begin{array}{*{20}{c}}
#'      \lambda  \\ 
#'      0 \\ 
#'      0 
#'      \end{array}} \right] + \left[ {\begin{array}{*{20}{c}}
#'        {1 - {\theta _1}}&0&0 \\ 
#'        {{\tau _1}}&{1 - {\theta _2}}&0 \\ 
#'        0&{{\tau _2}}&{1 - {\theta _3}} 
#'        \end{array}} \right]{Y_t} \\ 
#'    \operatorname{E} \left( {{Y_{t + 1}}} \right) =  & \left[ {\begin{array}{*{20}{c}}
#'      \lambda  \\ 
#'      0 \\ 
#'      0 
#'      \end{array}} \right] + \left[ {\begin{array}{*{20}{c}}
#'        {1 - {\theta _1}}&0&0 \\ 
#'        {{\tau _1}}&{1 - {\theta _2}}&0 \\ 
#'        0&{{\tau _2}}&{1 - {\theta _3}} 
#'        \end{array}} \right]\operatorname{E} \left( {{Y_t}} \right) \\ 
#'    \operatorname{E} (Y) =  & {\left\{ {I - \left[ {\begin{array}{*{20}{c}}
#'      {1 - {\theta _1}}&0&0 \\ 
#'      {{\tau _1}}&{1 - {\theta _2}}&0 \\ 
#'      0&{{\tau _2}}&{1 - {\theta _3}} 
#'      \end{array}} \right]} \right\}^{ - 1}}\left[ {\begin{array}{*{20}{c}}
#'        \lambda  \\ 
#'        0 \\ 
#'        0 
#'        \end{array}} \right] \\ 
#'    =  & {\left[ {\begin{array}{*{20}{c}}
#'      {{\theta _1}}&0&0 \\ 
#'      { - {\tau _1}}&{{\theta _2}}&0 \\ 
#'      0&{ - {\tau _2}}&{{\theta _3}} 
#'      \end{array}} \right]^{ - 1}}\left[ {\begin{array}{*{20}{c}}
#'        \lambda  \\ 
#'        0 \\ 
#'        0 
#'        \end{array}} \right] \\ 
#'         =  & \lambda \left[ {\begin{array}{*{20}{c}}
#'      {\theta _1^{ - 1}} \\ 
#'      {{\tau _1}\theta _1^{ - 1}\theta _2^{ - 1}} \\ 
#'      {{\tau _1}{\tau _2}\theta _1^{ - 1}\theta _2^{ - 1}\theta _3^{ - 1}} 
#'      \end{array}} \right]
#'    \end{aligned} \]
#' 
#' ### Variance of process
#' 
#' #### Mean conditional variance
#' 
#' \[\begin{aligned}
#'    \operatorname{Var} \left( {{Y_{t + 1}}|{Y_t}} \right) =  & \operatorname{Var} \left( {\left. {\begin{array}{*{20}{c}}
#'      {{X_{1,t + 1}}} \\ 
#'      {{X_{2,t + 1}}} \\ 
#'      {{X_{3,t + 1}}} 
#'      \end{array}} \right|} \right.\left. {\begin{array}{*{20}{c}}
#'        {{X_{1t}}} \\ 
#'        {{X_{2t}}} \\ 
#'        {{X_{3t}}} 
#'        \end{array}} \right)\\
#'         =  & \operatorname{Var} \left( {\left[ {\begin{array}{*{20}{c}}
#'          {{X_{1t}}} \\ 
#'          {{X_{2t}}} \\ 
#'          {{X_{3t}}} 
#'          \end{array}} \right] + \left[ {\begin{array}{*{20}{c}}
#'            {{N_t} + {B_{1t}} - {D_{1t}} - {T_{1t}}} \\ 
#'            {{T_{1t}} + {B_{2t}} - {D_{2t}} - {T_{2t}}} \\ 
#'            {{T_{2t}} + {B_{3t}} - {D_{3t}}} 
#'            \end{array}} \right]} \right| \left. {\begin{array}{*{20}{c}}
#'              {{X_{1t}}} \\ 
#'              {{X_{2t}}} \\ 
#'              {{X_{3t}}} 
#'              \end{array}} \right) \\ 
#'    =  & \operatorname{Var} \left( {\left[ {\begin{array}{*{20}{c}}
#'      1&1&{ - 1}&{ - 1}&0&0&0&0&0 \\ 
#'      0&0&0&1&1&{ - 1}&{ - 1}&0&0 \\ 
#'      0&0&0&0&0&0&1&1&{ - 1} 
#'      \end{array}} \right]} \right.\left[ {\begin{array}{*{20}{c}}
#'        {{N_t}} \\ 
#'        {{B_{1t}}} \\ 
#'        {{D_{1t}}} \\ 
#'        {{T_{1t}}} \\ 
#'        {{B_{2t}}} \\ 
#'        {{D_{2t}}} \\ 
#'        {{T_{2t}}} \\ 
#'        {{B_{3t}}} \\ 
#'        {{D_{3t}}} 
#'        \end{array}} \right]\left. {\left| {\begin{array}{*{20}{c}}
#'          {{X_{1t}}} \\ 
#'          {{X_{2t}}} \\ 
#'          {{X_{3t}}} 
#'          \end{array}} \right.} \right) \\ 
#'    =  & C\left( {\left[ {\begin{array}{*{20}{c}}
#'      1&1&{ - 1}&{ - 1}&0&0&0&0&0 \\ 
#'      0&0&0&1&1&{ - 1}&{ - 1}&0&0 \\ 
#'      0&0&0&0&0&0&1&1&{ - 1} 
#'      \end{array}} \right],diag\left( {\lambda ,{\beta _1}{X_{1t}},{\delta _1}{X_{1t}},{\tau _1}{X_{1t}},{\beta _2}{X_{2t}},{\delta _2}{X_{2t}},{\tau _2}{X_{2t}},{\beta _3}{X_{3t}},{\delta _3}{X_{3t}}} \right)} \right) \\ 
#'    \operatorname{E} \left( {\operatorname{Var} \left( {{Y_{t + 1}}|{Y_t}} \right)} \right) =  & C\left( {\left[ {\begin{array}{*{20}{c}}
#'      1&1&{ - 1}&{ - 1}&0&0&0&0&0 \\ 
#'      0&0&0&1&1&{ - 1}&{ - 1}&0&0 \\ 
#'      0&0&0&0&0&0&1&1&{ - 1} 
#'      \end{array}} \right],\operatorname{E} \left( {diag\left( {\lambda ,{\beta _1}{X_{1t}},{\delta _1}{X_{1t}},{\tau _1}{X_{1t}},{\beta _2}{X_{2t}},{\delta _2}{X_{2t}},{\tau _2}{X_{2t}},{\beta _3}{X_{3t}},{\delta _3}{X_{3t}}} \right)} \right)} \right) \\ 
#'    & {\text{where }}C\left( {A,B} \right) = ABA' \\ 
#'    \end{aligned} \]
#'    
#' 
#' #### Variance of conditional mean
#' 
#'  \[\begin{aligned}
#'    \operatorname{Var} \left( {\operatorname{E} \left( {{Y_{t + 1}}|{Y_t}} \right)} \right) =  & \operatorname{Var} \left( {\left[ {\begin{array}{*{20}{c}}
#'      {1 - {\theta _1}}&0&0 \\ 
#'      {{\tau _1}}&{1 - {\theta _2}}&0 \\ 
#'      0&{{\tau _2}}&{1 - {\theta _3}} 
#'      \end{array}} \right]{Y_t}} \right) = \left[ {\begin{array}{*{20}{c}}
#'        {1 - {\theta _1}}&0&0 \\ 
#'        {{\tau _1}}&{1 - {\theta _2}}&0 \\ 
#'        0&{{\tau _2}}&{1 - {\theta _3}} 
#'        \end{array}} \right]\operatorname{Var} \left( {{Y_t}} \right){\left[ {\begin{array}{*{20}{c}}
#'          {1 - {\theta _1}}&0&0 \\ 
#'          {{\tau _1}}&{1 - {\theta _2}}&0 \\ 
#'          0&{{\tau _2}}&{1 - {\theta _3}} 
#'          \end{array}} \right]^\prime } \\ 
#'  \operatorname{Var} \left( {{Y_{t + 1}}} \right) =  & P + Q\operatorname{Var} \left( {{Y_t}} \right)Q' \\ 
#'  \end{aligned} \]
#'    
#'  Thus the variance of the process, $\operatorname{Var} \left( {Y} \right)$, if it exists, is the solution to the equation:
#' 
#'  \[\operatorname{Var} \left( Y \right) = P + L\operatorname{Var} \left( Y \right)L'\]
#' 
#'  where $P$ is the mean conditional variance and $Q$ is a lower triangular matrix. In one dimension, 
#'  the solution is analogous to the well known solution for the sum of a geometric sequence.  
#'   
#'  To find a solution to the matrix version, we exploit the fact that $L$ is lower triangular with, 
#'  if the process is stationary, distinct eigenvalues. Thus, it is easily diagnalized, as shown below, with
#'  $L = GDG^{-1}$ where $D$ is diagonal.   
#'
#' #### Solution of *V = A + LVL'*
#'  
#'  I think that the following works (check details) if  $||L|| < 1$.  Consider:
#'  
#'  \[\begin{aligned} 
#'  V = & A + LVL' \\
#'  = & A + L \left( A + LVL' \right)L'\\
#'  = & A + LAL' + L^2\left( A + LVL' \right)L'^2\\
#'  = & A + LAL' + L^2AL'^2 + ... + L^nAL'^n + L^{n+1}VL'^{n+1}\\
#'  = & A + LAL' + L^2AL'^2 + ... + L^nAL'^n + ...\\
#'  \end{aligned}\]
#'  
#'  If $L = GDG^{-1}$, then
#'  
#'  $$V = A + GDG^{-1}AG'^{-1}DG'+ GD^2G^{-1}AG'^{-1}D^2G'+ GD^3G^{-1}AG'^{-1}D^3G'  +...$$
#'  
#'  Letting $B = G^{-1}AG'^{-1}$:
#'  
#'  $$V = G(B + DBD + D^2BD^2 + D^3BD^3 + D^4BD^4  +... )G'$$
#'  
#'  Each term of $K = B + DBD + D^2BD^2 + D^3BD^3 + D^4BD^4  +...$ has the form:
#'  
#'  \[\begin{aligned}  
#'  K_{jk} = & B_{jk} + D_{jj}B_{jk}D_{kk} + D_{jj}^2B_{jk}D_{kk}^2 + D^3_{jj}B_{jk}D^3_{kk} + D_{jj}^4B_{jk}D_{kk}^4  +...\\ 
#'  = & B_{jk} + D_{jj}D_{kk}B_{jk} + (D_{jj}D_{kk})^2B_{jk} + (D_{jj}D_{kk})^3B_{jk} + (D_{jj}D_{kk})^4B_{jk}  +...\\
#'  = & B_{jk} + \phi_{jk}B_{jk} + \phi_{jk}^2B_{jk} + \phi_{jk}^3B_{jk} + \phi_{jk}^4B_{jk}  +...\\
#'  \end{aligned}\]  
#'  
#'  Thus, if the geometric series converges:
#'  
#'  $$K_{jk} = \frac{B_{jk}}{1-D_{jj}D_{kk}}$$
#'  
#'  and
#'  
#'  $$V = GKG'$$
#'    
#'  Ultimately, this algorithm needs to be programmed in STAN since it needs to be used for the 
#'  evaluation of the 'equilibrium' likelihood in which likelihood is evaluated on the assumption that
#'  the starting profile reflects the process in equilibrium.
#'  
#'  To check the derivation, we program the algorithm in R and will check its accuracy through simulation.
#'  
# Conditional expectation ---------

cond_ev <- function(Y, lambda = parms$lambda, 
                    birth = parms$birth, 
                    death = parms$death,
                    tran = parms$tran,
                    parms, # optional way of providing parameters
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


cond_ev(c(100,100,100), 10, .01, .03, c(.01,.02))
parms <- list(lambda = 10, birth = rep(.01,3),
              death = .03, tran = c(.01,.02))
cond_ev(c(100,100,100), parms = parms)

Nsteps <- 2000
mat <- matrix(NA, 3, Nsteps)
mat[,1] <- c(100,100,100)
for ( i in 2:Nsteps) mat[,i] <- cond_ev(mat[,i-1],parms = parms)
# Equilibrium profile with these parameters:
head(t(mat),10)

tail(t(mat),10)
# matplot(t(mat), type = 'l')

# Marginal (equilibrium) expectation -------------------

marg_ev <- function(lambda = parms$lambda, 
                    birth = parms$birth, 
                    death = parms$death,
                    tran = parms$tran,
                    parms, # optional way of providing parameters
                    verbose = FALSE) {
  # Returns: E( Y ) at equilibrium
  N <- length(birth)  # must have correct length
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

marg_ev(parms=parms)

# By simulation

mat[,1]
mat[,10]
mat[,100]
mat[,500]
mat[,1000]

# Conditional variance ----------

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

cond_var(c(100,100,100), parms = parms, verbose = T)  

# Marginal (equilbrium) variance --------------------------------------------------------------------

# Eigenvectors of lower triangular matrix -----------------------------------------------------------
evec_lt <- function(A){
  # eigenvectors of a lower triangular matrix
  # return E where A = E %*% diag(A) %*% E^{-1}
  # Note that diag(A) is the diagonal matrix of eigenvectors (not ordered)
  N <- nrow(A)
  ret <- diag(N)
  for ( i in (N-1):1) {
    inds <- (i+1):N
    ret[inds,i] <- solve(A[inds,inds]-A[i,i]*diag(N-i),-A[inds,i])
  }
  ret
}
# Check evec_lt
A <- matrix(rnorm(25), 5,5)
A[col(A) > row(A)] <- 0
A
evec_lt(A) %*% diag(diag(A))- A %*% evec_lt(A)   # this should be machine 0

# Solution to variance equation ------------------------------------------------------------------

var_equation <- function(A,L){
  # Returns V such that V = A + LVL' if a solution exists ??? (check conditions)
  # likely to return nonsense otherwise
  # This does for matrices what the formula for a geometric series does for scalars
  G <- evec_lt(L)  # this might not work although it probably should for reasonable
  # parameters for T-cells  
  D <- diag(L)
  Ginv <- solve(G)
  # disp( L - G %*% (D* Ginv)) # ok
  B <- Ginv %*% A %*% t(Ginv)
  K <- B / (1 -outer(D,D))
  G %*% K %*% t(G)
}

A <- diag(3) + .1
A
L <- cbind( c(.8,.1,.1), c(0, .7, .1), c(0,0,.6))
L
V <- var_equation(A,L)
V
V - A - L%*% V %*% t(L)  # should be machine 0

# Variance of conditional mean --------------------------------------------------------------------

# Mean conditional variance

mean_cond_var <- function( 
  lambda = parms$lambda, 
  birth = parms$birth, 
  death = parms$death,
  tran = parms$tran,
  parms, # optional way of providing parameters
  verbose = FALSE) {
  # Returns: E ( Var( Y_{t+1} | Y_t)) = T E(D) T'
  # where D is diagonal matrix of variances of new cell births, deaths and transitions
  EY <- marg_ev(lambda=lambda,birth=birth,death=death,tran=tran) 
  cond_var(EY,lambda=lambda,birth=birth,death=death,tran=tran) 
}

mean_cond_var(parms = parms) 


# Marginal variance --------------------------------------------------------------

marg_var <- function(
  lambda = parms$lambda, 
  birth = parms$birth, 
  death = parms$death,
  tran = parms$tran,
  parms, # optional way of providing parameters
  verbose = FALSE) {
  # Returns: Var (Y)
  N <- length(birth)  # must have correct length
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
  P <- mean_cond_var(lambda=lambda,birth=birth,death=death,tran=tran)
  if(N == 1) return(P/(2*theta[1] - theta[1]^2))
  L <- diag(1 - theta)
  L[row(L)==(col(L)+1)] <- tran
  var_equation(P, L)
}

marg_var(parms = parms)

# Simulation tests --------------------------------------------------------------------------------------

sim_step <- function(Y, 
                     lambda = parms$lambda, 
                     birth = parms$birth, 
                     death = parms$death,
                     tran = parms$tran,
                     parms, # optional way of providing parameters
                     verbose = FALSE) {
  N <- length(Y)  # must have correct length
  # death and trans can be entered as single value
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
  # L matrix
  L <- matrix(0,N,3*N+1) # will drop last column
  for (i in 1:N) L[i,(3*i-2):(3*i+1)] <- c(1,1,-1,-1)
  L <- L[,-(3*N+1)] # drop last column
  if(verbose) disp(L)
  Diffs <- rpois(length(D), lambda = D)
  Y + L %*% Diffs
}

system.time(
  sim_step(c(100,100,100), parms = parms)
)
system.time(
  sim_step(c(100,100,100), parms = parms)
)

sim_seq <- function(N, Y = marg_ev(lambda,birth,death,tran), 
                    lambda = parms$lambda, 
                    birth = parms$birth, 
                    death = parms$death,
                    tran = parms$tran,
                    parms, # optional way of providing parameters
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

parms <- list(lambda=10, birth = c(.01,.02,.01), death = .03, tran = c(.02, .04))
matplot(sim_seq(N = 50, parms = parms), type = 'l', lty = c(1,2,3), lwd = 2, col = c('black','red','green'))
legend("right", inset=.05, legend=c("1", "2", "3"), lty = c(1,2,3), lwd = 2, col = c('black','red','green') , horiz=TRUE)
matplot(sim_seq(N = 1000, parms = parms), type = 'l', lty = c(1,2,3), lwd = 2, col = c('black','red','green'))
legend("right", inset=.05, legend=c("1", "2", "3"), lty = c(1,2,3), lwd = 2, col = c('black','red','green') , horiz=TRUE)


