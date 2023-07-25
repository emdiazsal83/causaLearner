

x <- rnorm(100)

dens <- density(x)


hist(x, prob=T)
lines(dens)

library(hdrcde)

n <- 1000

x <- rnorm(n)
y <- x*rnorm(n)

plot(x,y)


cdensXY <- cde(x,y)
cdensYX <- cde(y,x)

hist(y)

plot(cdensXY)
plot(cdensYX)

plot(cdensXY, plot.fn="hdr")
plot(cdensYX, plot.fn="hdr")
