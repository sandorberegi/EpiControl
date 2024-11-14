library(fastbeta)

set.seed(2L)
n <- 200L
d <- 50L
p <- 0.1
prob <- plogis(rlogis(d + n, location = qlogis(p), scale = 0.1))
delay <- diff(pgamma(0L:(d + 1L), 12, 0.4))

h <- function (x, a = 1, b = 1, c = 0) a * exp(-b * (x - c)^2)
ans <- floor(h(seq(-60, 60, length.out = d + n), a = 1000, b = 0.001))

x0 <- rbinom(d + n, ans, prob)
x <- tabulate(rep.int(1L:(d + n), x0) +
                sample(0L:d, size = sum(x0), replace = TRUE, prob = delay),
              d + n)[-(1L:d)]

str(D0 <- deconvolve(x, 0.1, delay, complete = FALSE))
str(D1 <- deconvolve(x, 0.1, delay, complete =  TRUE))

matplot(-(d - 1L):n,
        cbind(x0, c(rep.int(NA, d), x), prob * D0[["value"]], p * ans),
        type = c("p", "p", "p", "l"),
        col = c(1L, 1L, 2L, 4L), pch = c(16L, 1L, 16L, NA),
        lty = c(0L, 0L, 0L, 1L), lwd = c(NA, NA, NA, 3L),
        xlab = "time", ylab = "count")
legend("topleft", NULL,
       c("actual", "actual+delay", "actual+delay+deconvolution", "p*h"),
       col = c(1L, 1L, 2L, 4L), pch = c(16L, 1L, 16L, NA),
       lty = c(0L, 0L, 0L, 1L), lwd = c(NA, NA, NA, 3L),
       bty = "n")

plot(0L:D1[["iter"]], D1[["chisq"]], xlab = "iterations", ylab = quote(chi^2))
abline(h = 1, lty = 2L)
