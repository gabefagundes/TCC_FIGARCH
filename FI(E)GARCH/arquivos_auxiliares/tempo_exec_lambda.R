
{
    library(tidyverse)
}

truncs = c(100, 1000, 10000)
lambda1_list = list()
lambda2_list = list()

tempos = tibble(
    delta1 = numeric(3),
    delta2 = numeric(3),
)


for (i in seq_along(truncs)){
    

    # Algoritmo com multiplicação de polinômios #
    
    message(paste("Algoritmo Polinomial com trunc =", truncs[i]))
    
    t00 = Sys.time()
    
    lambda1_list[[i]] = lambda(truncs[i], d = 0.8, 0.1, 0.1)
    
    t01 = Sys.time()
    
    tempos[i, 1] = as.numeric(t01 - t00)
    
    
    # Algoritmo recursivo
    
    message(paste("Algoritmo Recursivo com trunc =", truncs[i]))
    
    t10 = Sys.time()
    
    lambda2_list[[i]] = lambda2(truncs[i], d = 0.8, 0.1, 0.1)
    
    t11 = Sys.time()
    
    tempos[i, 2] = as.numeric(t11 - t10)
    
}


