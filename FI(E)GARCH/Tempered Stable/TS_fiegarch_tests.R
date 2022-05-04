# Script para testar processos FIEGARCH com tempered stable dists
# Iniciado em 15/03

source('FI(E)GARCH/Tempered Stable/tempered_stable_fiegarch_simulation.R', encoding = 'UTF-8')

# Setando parametros
{
    N = 5000
    d = 0.3
    theta = 0.3
    gamma = -1.2
    omega = -4.5
    a = 0
    b = c(0.4, -0.21)
    m = 10000
    alpha = 0.7
}


lamb = lambda(m, d, a, b)


# Gerando Valores
{
    t0 = Sys.time()
    teste2 = fiegarch_TS(N, d, omega, a, b, theta, gamma, m, alpha, lamb = lamb)
    #teste1 = lambda(trunc = 50000, d = 0.3, 0.1, 0.68)
    t1 = Sys.time()
    print(round(t1-t0, 2))
    #print(tail(teste1))
}



# Resultados
{
    y = teste2$serie
    par(mfrow = c(3,1))
    plot.ts(teste2$serie)
    plot.ts(teste2$sigma)
    plot.ts(teste2$sigma^2 %>% exp)
    #plot.ts(teste$res)
}



Box.test(teste$serie, lag = 20, type = 'Box-Pierce')


saveRDS(teste,
        paste0('FI(E)GARCH/dados/TS_N',N, '_d_', d, '_omega_', omega, '_a_', paste0(a, '_'), 
               'b_', paste0(b, '_', collapse = ''), 'theta_', theta, '_gamma_', gamma, 
               '_m_', m, '_alpha_', alpha, '.rds')
)





teste5 = readRDS(paste0('FI(E)GARCH/dados/TS_N',N, '_d_', d, '_omega_', omega, '_a_', paste0(a, '_'), 
                       'b_', paste0(b, '_', collapse = ''), 'theta_', theta, '_gamma_', gamma, 
                       '_m_', m, '_alpha_', alpha, '.rds'))



Box.test(teste$serie, lag = 20, type = 'Box-Pierce')



acf(log(teste$serie^2), type = 'covariance')







# Lista de parametros que podem ser bons para mostrar ---------------------

{
    N = 2000
    d = 0.15
    theta = 0.5
    gamma = 0.3
    omega = -2
    a = 0.4
    b = c(0.55)
    m = 5000
    alpha = 1.5
}
