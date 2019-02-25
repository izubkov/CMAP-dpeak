#' Euclidean distance
distEuclidean <- function(x, centers)
{
  if(ncol(x)!=ncol(centers))
    stop(sQuote("x")," and ",sQuote("centers"),
         " must have the same number of columns")
  z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
  for(k in 1:nrow(centers)){
    z[,k] <- sqrt( colSums((t(x) - centers[k,])^2) )
  }
  z
}

#' kmeans++ implementation
#'
#' @param x vector of data
#' @param k how many clusters
#' @param start random/peak start
#'
#' @return centers to initialise kmeans
#' @export
#'
#' @examples
kmeanspp <- function(x, k = 2, start = "rand")
{
  centers <- matrix(0, nrow=k, ncol=ncol(x))
  if(start == "max") {
    h <- hist(x, breaks = length(x), plot = F)
    mx <- h$breaks[h$counts == max(h$counts)]
    centers[1,] <- sample(mx, 1)
  } else if(start == "rand") {
    centers[1,] <- x[sample(1:nrow(x), 1), , drop=FALSE]
  }
  d <- distEuclidean(x, centers[1L,,drop=FALSE])^2
  for(l in 2:k){
    centers[l,] <- x[sample(1:nrow(x), 1, prob=d), , drop=FALSE]
    d <- pmin(d, distEuclidean(x, centers[l,,drop=FALSE])^2)
  }
  centers
}
