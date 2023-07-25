
remove(list=ls())
setwd("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R/pkg_entropy")
dyn.load("./kdpee/kdpeer.so")

# void kdpeeR(int *n, int *d, double *x, double *mins, double *maxs, double *zcut, double *out)

n <- 10
d <- 2
X <- matrix(rnorm(n*d), n, 2)


n <- as.integer(nrow(X))
d <- as.integer(ncol(X))
x <- as.double(as.numeric(X))
  
mins <- as.double(apply(X, 2, min))
maxs <- as.double(apply(X, 2, max))
zcut <- as.double(1.96)
# void kdpeer(double *x, int *n, int *d, double *mins, double *maxs, double *zcut,  double *out)
res <- .C("kdpeer", x, n, d, mins, maxs, zcut, as.double(numeric(1)))
res[[length(res)]]


