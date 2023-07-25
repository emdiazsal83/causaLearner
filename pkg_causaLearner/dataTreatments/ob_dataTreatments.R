# construct data treatments

# perm = FALSE o seed
norm <- list(perm=1234, scalingFunc="norml", maxPoints=20000)

stdr <- list(perm=1234, scalingFunc="stdrize", maxPoints=1000)

gaussianiz <- list(perm=1234, scalingFunc="gaussianize", maxPoints=1000)

none <- list(perm=1234, scalingFunc="iden", maxPoints=1000)
