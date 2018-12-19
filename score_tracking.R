#library(plot3Drgl)
#f <- function(x,y)
#  (x^2) * (y^2)
#x <- y <- seq(0,1,len=20)
#z <- outer(x, y, f)
#persp3Drgl(z=z)

t.bench <- 162
exp_mul <- function(t) exp(-t/(3*t.bench))
plot(function(t) exp_mul(t), xlim = c(0, t.bench), ylim = c(0.7, 1.0),
     xlab = "t", ylab = "exponent multiplier")
abline(h = exp_mul(0.9*t.bench), col = "red")
abline(v = log(1/0.99)*3*t.bench, col = "green")
