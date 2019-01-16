plot_clusters <- function(a, b, FI) {
  p.FI <- hist(FI, breaks = length(FI))
  p.a <- hist(a, breaks = p.FI$breaks)
  p.b <- hist(b, breaks = p.FI$breaks)
  plot( p.FI, col=rgb(1,1,1,1/4), xlim=c(0,max(FI)))
  plot( p.a,  col=rgb(1,0,0,1/4), xlim=c(0,max(FI)), add=T)
  plot( p.b,  col=rgb(0,0,1,1/4), xlim=c(0,max(FI)), add=T)
}

FI <- c(2469,1977,2410,2827,1826,2812,2213,2520,2328,2410,2481,2487,2574,2018,3027,2523,2304,2283,2180,2029,2625,2209,2509,2265,4723,2342,2481,2708,2887,2364,2967,2785,2597,2146,2304,66,3237,2913,2821,2650,2542,2339,2973,2725,2212,2442,2434,2983,2145,2730,2109,2230,2549,2317,2807,2447,2483,4723)

FI.s <- FI

a <- vector(mode = "integer") #, length = length(FI.s))
b <- vector(mode = "integer") #, length = length(FI.s))

# init
h <- hist(FI.s, breaks = length(FI.s), plot = F)
# select most frequent maximum
br <- h$breaks[2] - h$breaks[1]
mid <- h$mids[which(h$counts == max(h$counts))]
mx <- max(Filter(function(x) mid <= x && x <= mid + br, FI.s))
a[1] <- mx
idx <- which(FI.s == mx)[1]
# remove it from array
FI.s <- FI.s[-idx]
# latest point from current cluster
last <- a[1]

is.b <- T
b.skipping <- T

while(length(FI.s) > 1) {
  # next step
  p <- sample(FI.s, 1, prob = sapply(FI.s, function(x, c, xs) { (c - x)^2 / sum( (xs - c)^2 ) }, c = last, xs = FI.s) )
  p
  idx <- which(FI.s == p)[1]
  FI.s <- FI.s[-idx]
  if(is.b) {
    if(b.skipping) {
      # add to a
      p <- sample(FI.s, 1, prob = sapply(FI.s, function(x, c, xs) { (c - x)^2 / sum( (xs - c)^2 ) }, c = p, xs = FI.s) )
      idx <- which(FI.s == p)[1]
      FI.s <- FI.s[-idx]
      a[length(a) + 1] <- p
      #is.b <- F
      last <- p

      b.skipping <- F
    } else {
      b[length(b) + 1] <- p
      is.b <- F
      b.skipping <- T
      last <- p
    }
  } else {
    a[length(a) + 1] <- p
    is.b <- T
    last <- p
  }
  a; b
}

FI.s

if(is.a) {
  a[length(a) + 1] <- FI.s
} else {
  b[length(b) + 1] <- FI.s
}

plot_clusters(a, b, FI)
length(a); length(b)
