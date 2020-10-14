library(bbmle)
library(devtools)
library(Deriv)
library(numDeriv)

## wherever the code is actually stored on your computer

## loglik_list <- list(
##     dpois = expression(x * log(lambda) - lambda - lfactorial(x)),
##     dnorm = expression(-log(2*pi*sd) - (x-mean)^2/(2*sd^2))
## )

loglik_list <- list(
    dpois = list(expr=expression(x * log(lambda) - lambda - lfactorial(x)),
                 params=c("lambda")),
    dnorm = list(expression(-log(2*pi*sd) - (x-mean)^2/(2*sd^2)),
                 params=c("mean","sd"))
)

## bbmle_home <- "~/R/pkgs/bbmle"
form <- y ~ dpois(lambda = b1 * latitude^2)
set.seed(101)
dd <- data.frame(
  y = rpois(100, lambda = 2),
  latitude = rnorm(100)
)

load_all(bbmle_home)
## calc_mle2_function
debug(bbmle:::calc_mle2_function)
debug(calc_mle2_function)

calc_mle2_function(form, start = list(b0 = 1), data = dd)

## form -> objective function (i.e. a function that computes -sum(dpois(y, lambda=b0*latitude^2))
##  where the variables in the formula are *evaluated* in an environment that
##  includes the current parameter
##  values and information stored in a 'data' variable


## data frame lat, long
y ~ dpois(exp(log_lambda), ...,
  parameters = list(log_lambda = ~ poly(lat, long, 2))
)

mkfun <- function(formula, data) {
  RHS <- formula[[3]]
  response <- formula[[2]]
  ddistn <- as.character(RHS[[1]]) ## get the name of distribution variable
  arglist <- as.list(RHS[-1]) ## delete function name
  arglist1 <- c(
    list(x = response),
    arglist, ##
    list(log = TRUE)
  )
  fn <- function(pars) {
    pars_and_data <- c(as.list(pars), data)
    r <- with(
      pars_and_data,
      -sum(do.call(ddistn, arglist1))
    )
    return(r)
  }
  gr <- function(pars) {
    pars_and_data <- c(as.list(pars), data)
    if (!ddistn %in% names(loglik_list)) {
      stop("I can't evaluate the derivative for ", sQuote(ddistn))
    }
    ## eventually we need to calculate partial derivatives of the log-likelihood
    ## with respect to all of its parameters
    LL <- loglik_list[[ddistn]]$expr
    mnames <- loglik_list[[ddistn]]$params
    ## ???
    ## setdiff(all.vars(LL), "x")  ## response var should be the only non-parameter
    d0 <- Deriv(LL, mnames)
    ## evaluate all of the arguments to the log-likelihood
    arglist_eval <- lapply(arglist, eval, pars_and_data)
    ## evaluate response variable and assign its value to 'x'
    arglist_eval$x <- eval(response, pars_and_data)
    ## derivative of log-lik wrt PDF parameters
    d1 <- eval(d0, arglist_eval)
    ## compute the deriv of log_lik with respect to its parameters
    parnames <- setdiff(all.vars(RHS), names(data))
    dlist <- list()
    glist <- list()
    for (p in parnames) {
        dlist[[p]] <- eval(Deriv(arglist$lambda, p), pars_and_data)
        glist[[p]] <- -sum(d1 * dlist[[p]])
    }
    return(unlist(glist))
    ## d(loglik_pois/d(lambda))* d(lambda)/d(b0)
  }
  return(list(fn = fn, gr = gr))
}

## debug(mkfun)
ff <- mkfun(form, data = dd)
## debug(ff)
ff$fn(c(b1 = 2))
debug(ff$gr)
ff$gr(c(b1 = 2))


## use finite differences and/or
##  numDeriv::grad() to check this calculation.

## TEST

## PDF param <- multiple model parameters
form2 <- y ~ dpois(lambda = b0 + b1 * latitude^2)
## ff <- mkfun(form2)
## what happens if we forget data? error ...
ff2 <- mkfun(form2, data=dd)
ff2$fn(c(b0=1,b1=2))
undebug(ff$gr)
ff2$gr(c(b0=1,b1=2))

## testing
ff$gr(c(b1=2)) == ff2$gr(c(b0=0, b1=2))[["b1"]]

-sum(dpois(dd$y, 1 + 2*dd$latitude^2, log=TRUE))


debug(ff$gr)
ff$gr(c(b0=1,b1=2))

ff(c(b0

form3 <- y ~ dnorm(mean = b0 * latitude^2, sd = sd)
form4 <- y ~ dnorm(mean = b0 + b1 * latitude^2, sd = sd)



###### JUNK

-sum(dpois(dd$y, lambda = 1 * dd$latitude^2, log = TRUE))

## deriv of -neg log likelihood of
-sum(dpois(y, lambda = b0 * latitude^2, log = TRUE))

## -sum(deriv(loglik(lambda,x))/lambda *
##     deriv(lambda/parameters)


##
deriv(lambda_nll, "lambda")
library(Deriv)

eval(d1, list(x = 1, lambda = 2))
