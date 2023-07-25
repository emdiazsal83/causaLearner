#test simulation functions
remove(list=ls())
source("/home/soulivanh/Documents/proyectos/indepReg/Mooij/R/pkg_causaLearner/genData/func_getData_v2.R")

library(ggplot2)
library(reshape)

##################################################################################################*
# Random Functions, random dags
##################################################################################################*

# simulate 100 pairs, 200 data points each, random additive functions, x->y dag, use applySimfxH
q <- 100
set.seed(4)
p <- 2
(ps <- rep(p, q))
(ns <- rep(200, q))
nodes <- list(dist="runif", pars=list(min=-3, max=3), a=1, b=1)
nodes <- rep(list(nodes),p)
names(nodes) <- as.character(seq(p))
nodess <- lapply(ps, function(p) nodes)
dagMat <- matrix(c(0,0,1,0), 2, 2)
rownames(dagMat) <- colnames(dagMat) <- names(nodes)

set.seed(5)
dataTestList <- simRandSEMs(q, ps, ns, nodess, sigma=3, dagMat=dagMat, markovEquiv=T)
dataTestList$dags
plotPairsList(dataTestList)

# simulate 100 pairs, 200 data points each, random non-additive functions, random dag, applySimfxH

q <- 100
set.seed(4)
p <- 2
(ps <- rep(p, q))
(ns <- rep(200, q))
nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
nodes <- rep(list(nodes),p)
names(nodes) <- as.character(seq(p))
nodess <- lapply(ps, function(p) nodes)

set.seed(5)
dataTestList <- simRandSEMs(q, ps, ns, nodess, sigma=5, sigmaErr=3, geU=function(y, nois, scale, constant) y)
plotPairsList(dataTestList)

# simulate 100 data sets, 4 vars each, 200 data points each, random non-additive functions, random dag, applySimfxH

q <- 100
set.seed(4)
p <- 4
(ps <- rep(p, q))
(ns <- rep(200, q))
nodes <- list(dist="runif", pars=list(min=-1, max=1), a=1, b=1)
nodes <- rep(list(nodes),p)
names(nodes) <- as.character(seq(p))
nodess <- lapply(ps, function(p) nodes)

set.seed(5)
dataTestList <- simRandSEMs(q, ps, ns, nodess, sigma=5, sigmaErr=3, geU=function(y, nois, scale, constant) y)

                     

# plot data
data <- mapply(function(xs,nm) cbind(xs, nm), xs=dataTestList$xs, nm=dataTestList$names, SIMPLIFY=FALSE)

data <- melt(dataTestList$xs)

colnames(data) <- c("numData","xY","val","dataset")
data <- cast(data, dataset+numData~xY, value="val")
colnames(data) <- c("dataset","numData","x","y","w","z")

p <- ggplot(data)
p <- p + geom_point(aes(x, y), alpha=0.2, colour="steelblue4")
p <- p + facet_wrap(dataset~., scales="free")
p <- p + theme(strip.background = element_blank(), strip.text = element_blank(),
               axis.text=element_blank(), axis.ticks=element_blank())
print(p)


library(GGally)
i <- 4
ggpairs(data[which(data$dataset==i),c("x","y","w","z")], aes(alpha = 0.4))


 
library(dplyr)
library(tidyr)

gatherpairs <- function(data, ..., 
                        xkey = '.xkey', xvalue = '.xvalue',
                        ykey = '.ykey', yvalue = '.yvalue',
                        na.rm = FALSE, convert = FALSE, factor_key = FALSE) {
  vars <- quos(...)
  xkey <- enquo(xkey)
  xvalue <- enquo(xvalue)
  ykey <- enquo(ykey)
  yvalue <- enquo(yvalue)
  
  data %>% {
    cbind(gather(., key = !!xkey, value = !!xvalue, !!!vars,
                 na.rm = na.rm, convert = convert, factor_key = factor_key),
          select(., !!!vars)) 
  } %>% gather(., key = !!ykey, value = !!yvalue, !!!vars,
               na.rm = na.rm, convert = convert, factor_key = factor_key)
}

data %>% 
  gatherpairs(w, x, y, z) %>% {
    ggplot(., aes(x = .xvalue, y = .yvalue)) +
      geom_point() + 
      geom_smooth(method = 'lm') +
      facet_wrap(dataset~ .xkey ~ .ykey, ncol = length(unique(.$.ykey)), scales = 'free', labeller = label_both) +
      scale_color_brewer(type = 'qual') + theme(strip.background = element_blank(), strip.text = element_blank(),
                                                axis.text=element_blank(), axis.ticks=element_blank())
  }

##################################################################################################*
# User provided functions, user provided dags
##################################################################################################*

