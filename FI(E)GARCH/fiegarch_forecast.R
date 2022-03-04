{
    library(stabledist)
    library(tidyverse)
    source("FI(E)GARCH/stable_fiegarch.R")
}

fiegarch_forecast = function(x, sigma, Z, pars = numeric(7), h){
    #### x = serie de retornos
    #### sigma = volatilidade estimada
    #### Z = residuos do modelo ajustado
    #### params = vetor de parametros estimado 
    #### h = número de passos a frente para prever
    
    
    forecast = numeric(h)
    
    if(length(pars) != 7){
        stop("Vetor de parâmetros suprido tem comprimento diferente de 7.")
    }
    names(pars) = c('d', 'alpha', 'theta', 'gamma', 'omega', 'a1', 'b1')
    
    
    # nome dos parâmetros 
    N = length(x)
    d = pars['d']
    alpha = pars['alpha']
    theta = pars['theta']
    gamma = pars['gamma']
    omega = pars['omega']
    a1 = pars['a1'] 
    b1 = pars['b1']
    
    # variância de g(Z_t)
    sigma_g = sqrt(theta^2 + gamma^2 - gamma^2*(mean(abs(Z)))^2 + 2*theta*gamma*mean(Z*abs(Z)))
    
    # lambdas 
    # começando loop em 1 e terminando em n + h - 1
    lambdas = lambda(N + h - 1, d = d, alpha = a1, beta = b1 )
    
    
    for (t in (N+1):(N+h)){
        if(t == N+1){
            sigma_circ = sqrt(exp(omega + )) 
        }
    }
}
    