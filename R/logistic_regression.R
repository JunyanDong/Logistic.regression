#'Logistic regression

#'

#'Find the estimated parameters of the logistic regression model.

#'

#'@param X input predictors

#'@param y input binary responses taking value 0 and 1

#'@param maxi the maximum of iteration

#'@param threshold the tolerance of precision

#'@return the estimated parameters(beta)

#'

#'@examples

#'X <- cbind(rbinom(30,1,0.4),rnorm(30,3,1))

#'y <- rbinom(30, 1, 0.45)

#'logistic_regression(X,y)

#'

#'@export
#'
logistic_regression = function(X,y,maxi = 100,threshold = 1e-10){
  n = nrow(X)
  X = cbind(rep(1,n),X)
  expression_Xb = function(X,beta){ #expression of probabilities w.r.t X and beta
    b = as.vector(beta)
    return(exp(X%*%b) / (1+ exp(X%*%b)))
  }
  p = ncol(X)
  beta = rep(0,p)
  m = 10000
  i = 0
  while(m > threshold){
    p = as.vector(expression_Xb(X,beta))#estimated probabilities
    W =  diag(p*(1-p)) #weight matrix
    beta_change = solve(t(X)%*%W%*%X) %*% t(X)%*%(y - p) #change of beta
    beta = beta + beta_change#update value
    m = sum(beta_change^2)# change of beta in the iteration
    i = i + 1
    if(i> maxi) { # make sure not to stuck in loop
      stop("This doesn't converge.")
      }
   }
  return(beta)
}
