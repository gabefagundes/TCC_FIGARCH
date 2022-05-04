source('FI(E)GARCH/Tempered Stable/tempered_stable_fiegarch_simulation.R', encoding = 'UTF-8')


{
    N = 5000
    d = 0.45
    theta = -0.22
    gamma = 0.37
    omega = -5.9
    a = 0.5
    b = -0.6
    m = 10000
    alpha = 1.5
}


lamb = lambda(m, d, a, b)

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
    
    
    if( d <= 0.35){
        return(NA)
    }

    if (d >= 0.45){
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
    #plot.ts(sigma)
    
    log_func = -n/2*log(2/pi) - 1/2*sum(log(sigma[t]^2) + (x[t]^2)/(sigma[t]^2), na.rm = T)
    print(log_func)
    
    if(log_func == Inf | log_func == -Inf){
        return(NA)
    }
    
    return(log_func)
}


re = 100
out = data.frame(
    d = numeric(),
    theta = numeric(),
    gamma = numeric(),
    omega = numeric(),
    a1 = numeric(),
    b1 = numeric()
)


arquivo = paste0('FI(E)GARCH/dados/estimacaoMC_N',N, '_d_', d, '_omega_', omega, '_a_', paste0(a, '_'),
                 'b_', paste0(b, '_', collapse = ''), 'theta_', theta, '_gamma_', gamma,
                 '_m_', m, '_alpha_', alpha, '.csv')

arquivo_series = paste0('FI(E)GARCH/dados/serie_N',N, '_d_', d, '_omega_', omega, '_a_', paste0(a, '_'),
                        'b_', paste0(b, '_', collapse = ''), 'theta_', theta, '_gamma_', gamma,
                        '_m_', m, '_alpha_', alpha, '.csv')


#write_csv2(out, arquivo)


for (i in 1:re){
    
    print(i)
    
    data = fiegarch_TS(N, d, omega, a, b, theta, gamma, m, alpha, lamb = lamb)
    
    est = try(optim(c(0.41, -0.1, 0.2, -5.5, 0.5, -0.6), quasi_llf_fiegarch_1d1, x = data$serie, method = 'BFGS', hessian = F,
                    control = list(fnscale = -1)))
    if(class(est) == "try-error") next 
    
    { #estimacao 
        auxx = read_delim(arquivo, delim = ';', locale = locale(decimal_mark = ',')) %>% 
            mutate_all(as.numeric)
        out[1,] = est$par
        print(est$par)
        final = bind_rows(auxx, out[1,])
        write_csv2(final,
                   arquivo)
    }
    
    { #serie gerada 
    aux_serie = read_delim(arquivo_series, delim = ';', locale = locale(decimal_mark = ','))
    
    aux_serie_2 = aux_serie %>% 
        cbind(data$serie)
    
    colnames(aux_serie_2) = paste0('x', 1:ncol(aux_serie_2))
    
    write_csv2(aux_serie_2,
               arquivo_series
               #'teste.csv'
               )
    }
    
}

