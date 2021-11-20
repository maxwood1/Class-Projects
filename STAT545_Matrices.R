## Max Woodbury
## STAT 545: Introduction to Computational Statistics
## This code computes various matrix decompositions


## Calculating Moore-Penrose Inverse using SVD

mpinv <- function(X, e=1.0E-08){
  decomp <- svd(X)
  decomp$d[abs(decomp$d) < e] <- 0
  k1 <- which(decomp$d != 0)
  V1 <- decomp$v[,k1, drop=F]
  U1 <- decomp$u[,k1, drop=F]
  # Moore-Penrose inverse
  mp.X <- V1 %*% diag( 1/decomp$d[k1], nrow=length(k1)) %*% t(U1)
  mp.X
}

# Example 1.1
X1 <- matrix(c(1,2,3,4,5,6),nrow=2)
X1
mpinv(X1)

# Example 1.2
X2 <- matrix(c(8,4,2,6,3,5),nrow=3)
X2
mpinv(X2)


## Compute Cholesky decomposition for 2 x 2 matrices

cholesky <- function(X) {
  a <- sqrt(X[1,1])
  b <- X[1,2]/a
  c <- sqrt(X[2,2]-b^2)
  C <- matrix(c(a,b,0,c),nrow=2)
  C
}

# Example 2.1
X3 <- matrix(c(2,-1,-1,2),nrow=2)
X3
C1 <- cholesky(X3)
C1
# Check that you get original matrix back
C1 %*% t(C1)

# Example 2.2
X4 <- matrix(c(2,1,1,4),nrow=2)
X4
C2 <- cholesky(X4)
C2
# Check that you get original matrix back
C2 %*% t(C2)


## Define sweep operator and use it to calculate the inverse and determinant of a matrix

sweepy <- function(X) {
  n <- nrow(X)
  det <- 1
  #sweep on all diagonal elements
  for(k in 1:n) {
    temp <- X
    X[k,k]=-1/temp[k,k]
    for(i in 1:n) {
      if(i!=k) X[i,k]=temp[i,k]/temp[k,k]
    }
    for(j in 1:n) {
      if(j!=k) X[k,j]=temp[k,j]/temp[k,k]
    }
    for(i in 1:n) {
      for(j in 1:n) {
        if(i!=k & j!=k) X[i,j]=temp[i,j]-(temp[i,k]*temp[k,j])/temp[k,k]
      }
    }
    det = det*temp[k,k]
  }
  #return inverse and determinant
  return(list(-X,det))
}

# Example 3.1
X5 <- matrix(c(2,-1,-1,2), nrow=2)
sweepy(X5)
#Check with pre-defined R functions
solve(X5)
det(X5)

# Example 3.2
X6 <- matrix(c(1,5,3,4,6,6,7,8,1), nrow=3)
sweepy(X6)
#Check with pre-defined R functions
solve(X6)
det(X6)
