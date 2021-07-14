
# Bibliotecas -------------------------------------------------------------

{
    library(stabledist)
    library(tidyverse)
}


# Auxiliar ----------------------------------------------------------------



f_alpha = function(x, alpha){
    
    m = length(x)
    out = 0
    for (i in 1:m){
        out = out + alpha[i]*x[i]
    }
    return(out)
}


# ARCH --------------------------------------------------------------------

arch_proc = function(mu, alpha0, alpha, y0, n){
    # Gera um processo ARCH(m) com inovações alpha-estaveis
    # mean = média do processo
    # alpha0 = constante de sig2
    # alpha: vetor de tamanho m >= 1
    # y0: vetor de tamanho m correspondente aos 'm' primeiros y's observados
    # n = tamanho do processo
    m = length(y0)
    e = rnorm(n)
    
    sig2 = numeric(n)
    
    y = numeric(n)
    y[1:m] = y0
    
    r0 = y0 - mu
    r = numeric(n)
    
    r[1:m] = r0
    
    for (i in (m+1):n){
        
        aux = (r[(i-m):(i-1)])^2
        
        sig2[i] = alpha0 + f_alpha(aux, alpha)
        
        r[i] = sqrt(sig2[i])*e[i]
        
        y[i] = mu + r[i]
    }
    out = list('y' = y, 'r' = r, 'sigma2' = sig2)
    
    return(out)
}
