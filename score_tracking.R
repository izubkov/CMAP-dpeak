library(plot3Drgl)

f <- function(x,y) log(x^2 * y^2)
x <- y <- seq(-4,4,len=20)
z <- outer(x, y, f)
persp3Drgl(z=z)


f <- function(x,y)
  (x^2) * (y^2)
x <- y <- seq(0,1,len=20)
z <- outer(x, y, f)
persp3Drgl(z=z)

plot(function(x) exp(-x), xlim = c(0, 1))
