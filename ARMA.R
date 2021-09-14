{
    library(stabledist)
    library(tidyverse)
}


MA_proc = function(n, theta){
    
    #################
    n = 100
    theta = c(0.2, 0.5)
    ###################
    q = length(theta)
    theta = c(1, theta)
    
    e = rnorm(n + q)
    
    y = stats::filter(e, theta, 'convolution', sides = 1)
    
    y = y[(q+1):(n+q)]
    
}


AR_proc = function(p, n){
    
    
    
}


x = rep(1, 50)

y = stats::filter(x, -0.5, 'r', sides = 1)

