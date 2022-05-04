{
    #library(SymTS)
    library(tidyverse)
    source('FI(E)GARCH/utils/poly_functions.R')
}




forecast_fiegarch = function(h = 1, x, d, omega, a, b, theta, gamma){
    
    
    n = length(x)
    sigma = numeric(n)
    g_z = numeric(2*n - 1)
    z = numeric(n)
    
    
    lambdas = lambda(n-1, d, a, b)
    
    # 'funcionou' porque gamma era pequeno .
    #omega = -3.084
    #gamma = 0.2
    
    sigma[1] = exp(0.5*omega)
    z[1] = x[1]/sigma[1]
    g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
    
    
    for(t  in 2:n){
        
        aux = (t-1):(t-1-(n-1)) + n
        
        (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
        x[t]
        (z[t] = x[t]/sigma[t])
        (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
    }
    
    
    plot.ts(z)
    hist(z)
}
