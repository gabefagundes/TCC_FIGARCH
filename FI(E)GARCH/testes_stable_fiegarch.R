source('FI(E)GARCH/stable_fiegarch.R', encoding = 'utf-8')


# testes ------------------------------------------------------------------

{
    N = 5000
    d = 0.4312
    alpha = 1.8
    theta = 0.1095
    gamma = -0.7376
    omega = -10
    a1 = 0
    b1 = 0.54
    m = 50000
    
}

{
    t0 = Sys.time()
    teste = fiegarch_1d1(N, d, a0, alpha, beta, theta, gamma, trunc)
    #teste1 = lambda(trunc = 50000, d = 0.3, 0.1, 0.68)
    t1 = Sys.time()
    print(round(t1-t0, 2))
    #print(tail(teste1))
}

{
    par(mfrow = c(3,1))
    plot.ts(teste$serie)
    plot.ts(teste$sigma)
    plot.ts(teste$residuos)
}


x = teste$serie

Box.test(teste$serie, lag = 20, type = 'Box-Pierce')


saveRDS(teste,
        paste0('FI(E)GARCH/dados/N',N, '_d_', d, '_omega_', a0, '_a1_', alpha, '_b1_', beta, '_theta_', theta, '_gamma_', gamma, '_m_', trunc, '_alpha_', alpha_stable, '.rds')
)
