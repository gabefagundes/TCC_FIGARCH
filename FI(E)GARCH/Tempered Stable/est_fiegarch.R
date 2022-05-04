# script para estimar parametros de um fiegarch (p,d,q)

quasi_llf_fiegarch_1d1 = function(pars = numeric(6), x, lamb2 = NULL ){
    
    
    
    #x = y[1:2000]
    n = length(x)
    sigma = numeric(n)
    g_z = numeric(2*n - 1)
    z = numeric(n)
    
    
    if(length(pars) != 6){
        stop("Vetor de parâmetros suprido tem comprimento diferente de 6.")
    }
    names(pars) = c('d', 'theta', 'gamma', 'omega', 'a', 'b')
    

    
    d = pars['d'] #= 0.43
    theta = pars['theta'] #= 0.11
    gamma = pars['gamma'] #= -0.73
    omega = pars['omega'] #= -3
    a = pars['a'] #
    b = pars['b'] #= 0.54
    # 
    # if(d > 0.5){
    #     d = 0.49
    # }
    
    
    if( d <= -0.5){
        return(NA)
    }
    
    if (d >= 0.5){
        return(NA)
    }
    
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
    
    
    #plot.ts(z)
    #hist(z)
    plot.ts(sigma)
    
    log_func = -n/2*log(2/pi) - 1/2*sum(log(sigma[t]^2) + (x[t]^2)/(sigma[t]^2), na.rm = T)
    print(log_func)
    
    if(log_func == Inf | log_func == -Inf){
        return(NA)
    }
    
    return(log_func)
}


y = teste$serie



otimizacao = optim(c(0.45, -0.1, 0.3, -6, 0.5, -0.6), quasi_llf_fiegarch_1d1, x = y, method = 'Nelder-Mead', hessian = F,
                   control = list(fnscale = -1, reltol = 1e-5))



otimizacao$par

otimizacao2 = optim(c(0.45, 0.1, 0.2, -5, 0.5, -0.6), quasi_llf_fiegarch_1d1, x = y, method = 'BFGS', hessian = F,
                   control = list(fnscale = -1))

#otimizacao_boa = otimizacao2
otimizacao2$par

otimizacao3 = optim(c(0.45, 0.1, 0.2 , -6, 0.5,-0.6), quasi_llf_fiegarch_1d1, x = y, method = 'L-BFGS-B', hessian = F,
                    control = list(fnscale = -1))


otimizacao3$par
