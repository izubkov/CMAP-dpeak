FI <- c(2469,1977,2410,2827,1826,2812,2213,2520,2328,2410,2481,2487,2574,2018,3027,2523,2304,2283,2180,2029,2625,2209,2509,2265,4723,2342,2481,2708,2887,2364,2967,2785,2597,2146,2304,66,3237,2913,2821,2650,2542,2339,2973,2725,2212,2442,2434,2983,2145,2730,2109,2230,2549,2317,2807,2447,2483,4723)

FI.s <- FI

a <- vector(mode = "integer", length = length(FI.s))
b <- vector(mode = "integer", length = length(FI.s))

# init
a[1] <- max(FI.s)
# select only one maximum
max.idx <- which(FI.s == max(FI.s))[1]
# remove it from array
FI.s <- FI.s[-max.idx]
# latest point from current cluster
last <- a[1]

while(length(FI.s) > 0) {
  # next step
  p <- sample(FI.s, 1,
              prob = sapply(FI.s, function(x, c, xs) {
                (c - x)^2 / sum( (xs - c)^2 )
                }, c = last, xs = FI.s)
              )
  idx <- which(FI.s == p)[1]
  b[1] <- p
  FI.s <- FI.s[-idx]
  last <- b[1]
}
