#' Function for creating chebyshev nodes in 1 dimension
#'
#' @param m a double - the number of nodes to create
#' @return a vector of m nodes on \[-1,1]
#' @export

chebnodes <- function(m, lb = -1, ub = 1){
  k = seq(1, m, 1)
  z = -cos(pi*(2*k - 1)/(2*m))
  z = (z+1)*(ub - lb)/2 + lb
  return(z)
}

#' Function for creating a 2d state space of nodes
#'
#' @param m a double - the number of nodes to create in each dimension
#' @param ub a vector of 2 upper bounds in each dimension
#' @param lb a vector of 2 lower bounds in each dimension
#' @return a dataframe with all combinations of nodes to pass for function evaluation using mapply or similar
#' @export

chebstatespace <- function(m, lb, ub){
  z = chebnodes(m)
  x = (z+1)*(ub[1] - lb[1])/2 + lb[1]
  y = (z+1)*(ub[2] - lb[2])/2 + lb[2]
  ss = tidyr::crossing("x" = x, "y" = y)
  return(ss)
}

#' Function for evaluating chebyshev polynomials at nodes
#'
#' @param z a vector of nodes length i
#' @param degree the degree j of the cheybyshev polynomial
#' @return a matrix with element (i,j) = T_j(z_i)
#' @export

chebeval <- function(z, degree = 5){
  Tvec = vector(mode = "list", length = length(z))
  j = 1
  for (x in z) {
    tv = c(1, x, rep(NA, degree - 1))
    i = 3
    while(i <= degree+1){
      tv[i] <- 2*x*tv[i-1] - tv[i-2]
      i = i+1
    }
    Tvec[[j]] = tv
    j = j+1
  }
  Tvec = do.call(rbind, Tvec)
  return(Tvec)
}

#' Function finding chebyshev coefficients
#'
#' @param fy a square matrix with element (i,j) corresponding to f(x_j, y_i)
#' @param degree the degree of the cheybyshev polynomial
#' @return a matrix of coefficients a_{i,j}
#' @export


chebcoefs <- function(fy, degree = 5){
  z = chebnodes(dim(fy)[1])
  Tmat = chebeval(z, degree)
  coefmat = matrix(NA, nrow = degree+1, ncol = degree+1)
  for (i in c(1:(degree+1))) {
    for (j in c(1:(degree+1))) {
      n = sum(sapply(c(1:length(z)), function(k) sum(fy[,k]*Tmat[k,i]*Tmat[,j])))
      d = sum(Tmat[,i]^2)*sum(Tmat[,j]^2)
      coefmat[i, j] = n/d
    }
  }
  return(coefmat)
}

#' Function for interpolating using chebyshev coefficients
#'
#' @param x0 x coordinate for interpolation
#' @param y0 y coordinate for interpolation
#' @param coefmat matrix of coefficients from chebcoefs
#' @param ub a vector of 2 upper bounds in each dimension
#' @param lb a vector of 2 lower bounds in each dimension
#' @return approximation for f(x,y)
#' @export

chebpred <- function(x0, y0, coefmat, lb, ub){
  xp = 2*(x0 - lb[1])/(ub[1] - lb[1]) - 1
  yp = 2*(y0 - lb[2])/(ub[2] - lb[2]) - 1
  degree = dim(coefmat)[1] - 1
  Tix = as.vector(chebeval(xp, degree))
  Tij = as.vector(chebeval(yp, degree))
  pxy = sum(sapply(c(1:(degree)), function(i) sum(coefmat[i,]*Tix[i]*Tij)))
  return(pxy)
}
