#library(plot3Drgl)
#f <- function(x,y)
#  (x^2) * (y^2)
#x <- y <- seq(0,1,len=20)
#z <- outer(x, y, f)
#persp3Drgl(z=z)

t.bench.dpk <- 163
t.bench.lit <- 161
t.bench <- mean(c(t.bench.dpk, t.bench.lit))
COR.bench.dpk <- 0.701
COR.bench.lit <- 0.819
AUC.bench.dpk <- mean(c(0.904, 0.926))
AUC.bench.lit <- mean(c(0.920, 0.921))
exp_mul <- function(t) exp(-t/(3*t.bench))
par(mfrow = c(2, 2))
plot(function(t) exp_mul(t), xlim = c(0, t.bench), ylim = c(0.7, 1.0),
     xlab = "t", ylab = "exponent multiplier")
abline(h = exp_mul(0.9*t.bench), col = "red")
abline(v = log(1/0.99)*3*t.bench, col = "green")
plot(function(x) x^2, xlim = c(0, 1),
     xlab = "COR", ylab = "COR multiplier")
abline(v = COR.bench.dpk, col = "darkred")
abline(v = COR.bench.lit, col = "red")
plot(function(x) x^2, xlim = c(0, 1),
     xlab = "AUC", ylab = "AUC multiplier")
abline(v = AUC.bench.dpk, col = "darkgreen")
abline(v = AUC.bench.lit, col = "green")

score <- function(COR.dpk, COR.lit, AUC.dpk, AUC.lit, T.dpk, T.lit) {
  dpk <- 1e6 * COR.dpk^2 * AUC.dpk^2 * exp(-T.dpk / 3 / t.bench)
  lit <- 1e6 * COR.lit^2 * AUC.lit^2 * exp(-T.lit / 3 / t.bench)
  mean(c(dpk, lit))
}
# competitor_pack_2: 421009.415034549
# calculated: 351131.7

score(COR.bench.dpk, COR.bench.lit, AUC.bench.dpk, AUC.bench.lit, t.bench.dpk, t.bench.lit)
score(COR.bench.dpk, COR.bench.lit, AUC.bench.dpk, AUC.bench.lit, t.bench.dpk, t.bench.lit)

# COR dkp +10% (+0.07) => +8%
# COR lit +10% (+0.82) => +12%
# AUC dkp +10% (+0.09) => +8.7%
# AUC lit +1% (+0.009) => +13%
# time x163            => +58%
# time +10% of bench   => +16%
