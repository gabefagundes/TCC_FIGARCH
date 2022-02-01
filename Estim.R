#install.packages('StableEstim')

library(StableEstim)
library(stabledist)




xx = rstable(1000, 1.5, -1, 1, 0)
curve(dstable(x, 1.5, -1, 1, 0), from = -15, to = 10)


est1 = Estim('ML', xx)

est2 = Estim('GMM', xx)

est3 = Estim('Cgmm', xx)

est4 = Estim('Kout', xx)


est
