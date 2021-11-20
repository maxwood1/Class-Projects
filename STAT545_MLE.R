## Max Woodbury
## STAT 545: Introduction to Computational Statistics
## Methods for computing MLE 
## With proposed method to find the global optimum using the Bisection method



## Implementing Newton-Raphson Method to find MLE
## for Cauchy(mu,1) with unknown mu
## and sample size n

#Cauchy log-likelihood
llCauchy <- function(mu, x){
  sum(dt(x-mu, df=1, log=TRUE))
}

#Cauchy gradient (first derivative of log-likelihood)
g.llCauchy <- function(mu, x){
  2 * sum((x-mu)/(1+(x-mu)^2))
}

#Cauchy hessian (second derivative of log-likelihood)
h.llCauchy <- function(mu, x){
  2 * sum( (-(1+(x-mu)^2) + 2*(x-mu)^2) /(1+(x-mu)^2)^2)
}

#implement Newton-Raphson method to find MLE (with unknown mu)
#given Cauchy values with sample size n
newraph <- function(x, stop=1e-10) {
  #start at midpoint
  x <- sort(x)
  n <- length(x)
  if(length(x) %% 2 == 0) mu <- x[n/2]
  else mu <- x[(n+1)/2]
  newmu <- mu - g.llCauchy(mu,x) / h.llCauchy(mu,x)
  while(abs(newmu-mu) > stop) {
    mu <- newmu
    newmu <- mu - g.llCauchy(mu,x) / h.llCauchy(mu,x)
  }
  newmu
}

#test with random Cauchy sample (any mu) from size n
n1 <- 100
mu <- 10
randCauchy <- rt(n=n1, df=1) + mu
newraph(randCauchy)



## Implement bisection method to find global optimum

#return single iteration of bisection method
mybisect <- function(fn, a, b, eps = 1.0E-10, ...){
  fa <- fn(a, ...) 
  fb <- fn(b, ...) 
  if(fa*fb > 0) {
    warning("Invalid arguments: fn(a) * fn(b) > 0")
    return(NULL)
  }
  while(abs(b-a) > eps) {
    m <- (a+b)/2 
    fm <- fn(m, ...)
    if(fa*fm <= 0){
      b <- m
      fb <- fm
    } else {
      a <- m
      fa <- fm
    }
  }
  m <- (a+b)/2
  fm <- fn(m, ...) 
  m
}

#plot log-likelihood function
plot.ll <- function(theta, x){
  y <- numeric(length(theta))
  for(i in 1:length(theta)) {
    y[i] <- llCauchy(mu=theta[i], x=x)
  }
  plot(theta, y, type="l")
}

#iterate bisection method to and return global optimum
globaloptim <- function(x, graph=T) {
  x <- sort(x)
  
  # find three consective x-values to create search intervals for bisection
  f <- numeric(length(x))
  for(i in 1:length(x)) {
    f[i] <- llCauchy(mu=x[i], x=x)
  }
  s <- which.max(f) + c(-1, 0, 1)
  
  # try interval x[s[1:2]]
  rs12 <- mybisect(g.llCauchy, a=x[s[1]], b=x[s[2]], eps = 1.0E-10, x)
  # try interval x[s[2:3]]
  rs23 <- mybisect(g.llCauchy, a=x[s[2]], b=x[s[3]], eps = 1.0E-10, x)
  
  if (graph) {
    #plot the log likelihood
    theta <- seq(min(x), max(x), len=500)
    y <- numeric(length(theta))
    for(i in 1:length(theta)) {
      y[i] <- llCauchy(mu=theta[i], x=x)
    }
    plot.ll (theta, x=x)
  
    #plot global optimum
    if(!is.null(rs12)) {
      abline(v=rs12, col=2)
      abline(h= llCauchy(mu=rs12, x=x), col=2)
    }
    if(!is.null(rs23)) {
      abline(v=rs23, col=3)
      abline(h= llCauchy(mu=rs23, x=x), col=3)
    }
  }
  if(!is.null(rs12)) return(rs12)
  if(!is.null(rs23)) return(rs23)
}

x <- rt(n=5, df=1) + 8
globaloptim(x)



## Implementing Newton-Raphson Method to find MLE
## for Cauchy(mu,sigma^2) with unknown mu and sigma
## and sample size n

#Cauchy gradient unknown mu and sigma
g.llCauchy1 <- function(theta, x){
  mu <- theta[1]
  sigma <- theta[2]
  gmu <- 2 * sum((x-mu)/(sigma^2+(x-mu)^2))
  gsigma <- -length(x)/sigma + 2 * sum((x-mu)^2/(sigma^3+sigma*(x-mu)^2))
  return(as.matrix(c(gmu,gsigma)))
}

#Cauchy hessian unknown mu and sigma
h.llCauchy1 <- function(theta, x){
  mu <- theta[1]
  sigma <- theta[2]
  H <- matrix(c(0,0,0,0),nrow=2)
  H[1,1] <- 2 * sum(((x-mu)^2-sigma^2) / ((x-mu)^2+sigma^2)^2)
  H[1,2] <- -4 * sum(sigma*(x-mu)/(sigma^2+(x-mu)^2)^2)
  H[2,1] <- -4 * sum(sigma*(x-mu)/(mu^2-2*x*mu+x^2+sigma^2)^2)
  H[2,2] <- length(x)/sigma^2 - 2 * sum((x-mu)^2*(3*sigma^2+(x-mu)^2)/(sigma^3+(x-mu)^2*sigma)^2)
  return(H)
}

#implement Newton-Raphson method to find MLE (unknown mu and sigma)
newraph1 <- function(x, stop=1e-10) {
  #start at midpoint
  x <- sort(x)
  mu <- median(x)
  sigma <- sqrt(var(x))
  theta <- as.matrix(c(mu,sigma))
  newtheta <- theta - solve(h.llCauchy1(theta,x)) %*% g.llCauchy1(theta,x)
  while((abs(newtheta[1]-theta[1]) > stop && abs(newtheta[2]-theta[2]) > stop) && abs(det(h.llCauchy1(theta,x))) > 1e-10) {
    theta <- newtheta
    newtheta <- theta - solve(h.llCauchy1(theta,x)) %*% g.llCauchy1(theta,x)
  }
  newtheta
}

#test with random Cauchy sample from size n=5
mu0 <- 8
randCauchy1 <- rt(n=5, df=1) + mu0
newraph1(randCauchy1)

