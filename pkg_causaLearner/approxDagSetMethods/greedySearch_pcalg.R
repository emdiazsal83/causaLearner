
library(pcalg)
## Load predefined data
data(gmInt)
class(gmInt)
names(gmInt)
plot(gmInt$g)
class(gmInt$x)
dim(gmInt$x)

## Define the score (BIC)
score <- new("GaussL0penIntScore", gmInt$x, gmInt$targets, gmInt$target.index) 

gds <- function (score, labels = score$getNodes(), targets = score$getTargets(), 
          fixedGaps = NULL, phase = c("forward", "backward", "turning"), 
          iterate = length(phase) > 1, turning = TRUE, maxDegree = integer(0), 
          verbose = FALSE, ...) {
  if (!turning) {
    phase <- c("forward", "backward")
    iterate <- FALSE
    warning(paste("The argument 'turning' is deprecated; please use 'phase' instead", 
                  "(cf. ?ges)", sep = " "))
  }
  phase <- match.arg(phase, several.ok = TRUE)
  caus.inf(algorithm="GDS", score = score, labels = labels, targets = targets, 
           fixedGaps = fixedGaps, phase = phase, iterate = iterate, 
           maxDegree = maxDegree, verbose = verbose, ...)
}

# modify to give back dag and not just essential graph
caus.inf <- function (algorithm = c("GIES", "GDS", "SiMy"), score, labels = score$getNodes(), 
          targets = score$getTargets(), ...) {
  pars <- list(...)
  # pars <- list(fixedGaps=fixedGaps, phase=phase, iterate=iterate, maxDegree=maxDegree, verbose=verbose)
  algorithm <- match.arg(algorithm)
  if (is.numeric(score)) {
    p <- score
    if (is.list(labels) && is(targets, "Score")) {
      score <- targets
      targets <- labels
      labels <- as.character(1:p)
      warning(paste("You are using a DEPRECATED calling convention for", 
                    "gies(), gds() or simy(); please refer to the documentation", 
                    "of these functions to adapt to the new calling conventions."))
    }
    else if (is(labels, "Score")) {
      score <- labels
      labels <- as.character(1:p)
      warning(paste("You are using a DEPRECATED calling convention for", 
                    "ges(); please refer to the documentation", "to adapt to the new calling convention."))
    }
  }
  else if (is.numeric(labels) && length(labels) == 1) {
    labels <- as.character(1:labels)
    warning(paste("You are using a DEPRECATED calling convention for", 
                  "gies(), ges(), gds() or simy(); please refer to the documentation", 
                  "of these functions to adapt to the new calling conventions."))
  }
  if (!is(score, "Score")) {
    stop("'score' must be of a class inherited from the class 'Score'.")
  }
  if (!is.character(labels)) {
    stop("'labels' must be a character vector.")
  }
  if (!is.list(targets) || !all(sapply(targets, is.numeric))) {
    stop("'targets' must be a list of integer vectors.")
  }
  essgraph <- new("EssGraph", nodes = labels, targets = targets, score = score)
  
  pars <- c(list(algorithm=algorithm), pars)
  #if (essgraph$caus.inf(algorithm, ...)) {
  if (do.call(essgraph$caus.inf, pars)) {
    if (algorithm == "GIES") {
      list(essgraph = essgraph, repr = essgraph$repr())
    }
    else {
      list(essgraph = dag2essgraph(essgraph$repr(), targets = targets), 
           repr = essgraph$repr())
    }
  }
  else stop("invalid 'algorithm' or \"EssGraph\" object")
}


## Estimate the essential graph
gds.fit <- gds(score) 


## Plot the estimated essential graph and the true DAG
if (require(Rgraphviz)) {
  par(mfrow=c(1,2))
  plot(gds.fit$essgraph, main = "Estimated ess. graph")
  plot(gmInt$g, main = "True DAG")
}
