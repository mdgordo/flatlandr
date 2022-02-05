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

#' Function for bilinear and simplical interpolation in 2D
#'
#' @param x0 coordinate for interpolation
#' @param y0 coordinate for interpolation
#' @param x vector of x coordinates with known function values
#' @param y vector of y coordinates with known function values
#' @param vfx function values at each (x,y). Note x and y must form a crossing - must have function values at each possible combination
#' @param method string either "bilinear" or "simplical"
#' @return interpolated value at x0 and y0
#' @export

simplr <- function(x0, y0, x, y, vfx, method = "bilinear", extrapolation.warning = TRUE){
  ## find closest points
  xvals = sort(unique(x))
  xi = findInterval(x0, xvals, all.inside = TRUE)
  x1 = xvals[xi]
  x2 = xvals[xi+1]
  yvals = sort(unique(y))
  yi = findInterval(y0, yvals, all.inside = TRUE)
  y1 = yvals[yi]
  y2 = yvals[yi+1]
  if ((x0>max(x) | x0<min(x) | y0>max(y) | y0<min(y)) & extrapolation.warning) {
    print("Warning: Point outside boundary, applying constant extrapolation")}
  if (x0>max(x)) x0 = max(x) else if (x0<min(x)) x0 = min(x)
  if (y0>max(y)) y0 = max(y) else if (y0<min(y)) y0 = min(y)

  ## function vals at corners
  v1 = vfx[which(x==x1 & y==y1)]
  v2 = vfx[which(x==x1 & y==y2)]
  v3 = vfx[which(x==x2 & y==y1)]
  v4 = vfx[which(x==x2 & y==y2)]

  if (method=="bilinear") {
    ## change of variables
    xp = -1 + 2*(x0-x1)/(x2-x1)
    yp = -1 + 2*(y0-y1)/(y2-y1)
    vp = .25*(v1*(1-xp)*(1-yp) + v2*(1-xp)*(1+yp) + v3*(1+xp)*(1-yp) + v4*(1+xp)*(1+yp))
  } else {
    ### change of variables
    xp = (x0-x1)/(x2-x1)
    yp = (y0-y1)/(y2-y1)
    if (xp+yp <= 1) {
      vp = v1*(1-xp-yp) + v2*yp + v3*xp
    } else {
      vp = v2*(1-xp) + v3*(1-yp) + v4*(xp+yp-1)
    }
  }
  return(vp)
}

#' Function for simplical interpolation in 3
#'
#' @param x0 coordinates of point for interpolation
#' @param x matrix with x(i, j) the coordinate of point i in dimension j - must be a crossing
#' @param vfx function values at each x
#' @return interpolated value at x0
#' @export

complexr <- function(x0, x, vfx, extrapolation.warning = TRUE){
  closestpts = vector(mode = "list", length =3)
  ## find closest points in x
  for (i in c(1:3)) {
    xvals = sort(unique(x[,i]))
    xi = findInterval(x0[i], xvals, all.inside = TRUE)
    closestpts[[i]] = c(xvals[xi], xvals[xi+1])
    if ((x0[i]>max(x[,i]) | x0[i]<min(x[,i])) & extrapolation.warning) {
      print("Warning: Point outside boundary, applying constant extrapolation")}
    if (x0[i]>max(x[,i])) x0[i]=max(x[,i]) else if (x0[i]<min(x[,i])) x0[i]=min(x[,i])
  }
  ## function vals at vertices
  d = tidyr::crossing(x = closestpts[[1]], y = closestpts[[2]], z = closestpts[[3]])
  d$v = 0
  for (i in c(1:nrow(d))) {
    d$v[i] = vfx[x[,1]==d$x[i] & x[,2]==d$y[i] & x[,3]==d$z[i]]
  }

  ### change of variables
  xp = (x0[1]-closestpts[[1]][1])/(closestpts[[1]][2] - closestpts[[1]][1])
  yp = (x0[2]-closestpts[[2]][1])/(closestpts[[2]][2] - closestpts[[2]][1])
  zp = (x0[3]-closestpts[[3]][1])/(closestpts[[3]][2] - closestpts[[3]][1])
  xcv = c(xp,yp,zp)

  ### algorithm
  p = order(xcv)
  s = vector(mode = "list", length = 4)
  fs = vector(mode = "list", length = 4)
  s[[1]] = rep(1, 3)
  fs[[1]] = vfx[x[,1]==s[[1]][1] & x[,2]==s[[1]][2] & x[,3]==s[[1]][3]]
  for (i in c(2:4)) {
    svec = rep(0, 3)
    svec[p[i-1]] = 1
    s[[i]] = s[[i-1]] - svec
    f = vfx[x[,1]==s[[i]][1] & x[,2]==s[[i]][2] & x[,3]==s[[i]][3]]
    f2 = vfx[x[,1]==s[[i-1]][1] & x[,2]==s[[i-1]][2] & x[,3]==s[[i-1]][3]]
    fs[[i]] = fs[[i-1]] + (1 - xcv[p[i-1]])*(f - f2)
  }
  return(fs[[4]])
}
