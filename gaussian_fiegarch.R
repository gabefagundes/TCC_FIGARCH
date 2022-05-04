{
    #library(SymTS)
    library(tidyverse)
    source('FI(E)GARCH/utils/poly_functions.R')
}



pi_func = function(k, d){
    d = -d
    
    for(i in 0:k){
        if(i == 0){
            out = 1
        } else{
            out = out*(i - 1 - d)/i
        }
        #print(out)
    }
    return(out)
}


# Função g(Z)

g_f_gauss = function(Z,
                  theta,
                  gamma){
    
    E_X = sqrt(2/pi)
    
    out = numeric(length(Z))
    
    for (t in seq_along(Z)){
        
        out[t] = theta*Z[t] + gamma*(abs(Z[t]) - E_X)
    }
    
    print("g(Z) Calculada.")
    return(out)
}



# Lambdas 
# usa funcoes do script poly_functions.R

lambda = function(m, d, a1, b1){
    
    p = c(1, a1)
    q = polyinv(c(1, b1), m = m)
    pq = poly.prod(a = p, b = q)
    
    s = sapply(0:m, pi_func, d)
    
    out = poly.prod(pq,s)[1:(m+1)]
    
    return(out)
    
}



# Geração Processo FIEGARCH -----------------------------------------------

fiegarch_gauss = function(N = 500,
                       d,
                       omega,
                       a,
                       b,
                       theta,
                       gamma,
                       m){
    print("Gerando variaveis iniciais.")
    
    sigma = numeric(N+m)
    Z = rnorm(N+m+1)
    x = numeric(N+m)
    
    
    print('Calculando g(Z)')
    g_z = g_f_gauss(Z,
                 theta = theta,
                 gamma = gamma)
    
    print('Gerando lambdas')
    lambdas = lambda(m, d, a, b)
    
    
    print("Gerando Serie")
    for (t in (m+1):(N+m)){
        #print(t)
        aux = (t):(t-m)
        
        sigma[t] = (omega + sum(lambdas*g_z[aux])) %>% 
            exp() %>% 
            sqrt()
        
        x[t] = sigma[t]*Z[t + 1]
        
    }
    print("Series e Sigmas Gerados")
    y = x[(m+1):(m + N)]
    sigma = sigma[(m+1):(m + N)]
    
    
    out = list(serie = y, sigma = sigma, res = y/sigma)
    return(out)
    
}




#### testes

{
    N = 5000
    d = .3
    theta = -0.0215
    gamma = 0.37
    omega = -5.79
    a = 0.1409
    b = -0.1611
    
    m = 50000
}




{
    t0 = Sys.time()
    teste = fiegarch_gauss(N, d, omega, a, b, theta, gamma, m)
    #teste1 = lambda(trunc = 50000, d = 0.3, 0.1, 0.68)
    t1 = Sys.time()
    print(round(t1-t0, 2))
    #print(tail(teste1))
}



{
    par(mfrow = c(3,1))
    plot.ts(teste$serie)
    plot.ts(teste$sigma)
    plot.ts(teste$res)
}


y = teste$serie
