library(stabledist)

# FUNÇÕES UTEIS -----------------------------------------------------------------
zeta0d0 = function(x, d, trunc = 1000){
    
    t = length(x)
    parcelas = numeric(trunc)
    
    for (j in 1:trunc){
        #print(j)
        parcelas[j] = gamma(j - d)/(gamma(j+1)*(gamma(-d))) * (x[t - j])^2
        
    }
    
    out = sum(parcelas)
    
    return(list('soma' = out, 'parcelas' = parcelas))
}

teste = zeta0d0(rnorm(1001), 0.1)

plot(teste$parcelas, type = 'h')

#####

figarch_sim_0d0 = function(N, d, x_init, trunc,  a0 = 0){
    
    z = rnorm(N+1+trunc)
    sigma = numeric(N+trunc)
    x = c(x_init, numeric(N))
    
    for(t in (trunc+1):(N+trunc)){
        print(t)
        
        aux = (t):(t-trunc)
        print(var(x[aux]))
        
        gammas = zeta(x[aux], d, trunc)$soma
        #print(gammas)
        
        sigma[t] = sqrt(a0 - gammas)
        
        
        x[t] = sigma[t]*z[t]
        
    }
    serie = x[(trunc+1): (N+trunc)]
    
    return(serie)
    
}


a = figarch_sim_0d0(1000, 0.5, rnorm(100) , trunc = 100, a0 = 0.5)


plot.ts(a)



# FIGARCH (1, d, 0) -------------------------------------------------------

zeta_1d0 = function(x,  d, trunc, alpha1){
    
    t = length(x)
    parcelas = numeric(trunc)
    
    for (j in 1:trunc){
        #print(j)
        parcelas[j] = gamma(j - d - 1)/(gamma(j)*(gamma(-d))) * (alpha1 - (j-1-d)/j)* (x[t - j])^2
        
    }
    
    out = sum(parcelas)
    
    return(list('soma' = out, 'parcelas' = parcelas))
}


teste = zeta_1d0(rt(51, df = 2), 0.4, 50, 0.5)

plot.ts(teste$parcelas)

#simulacao

figarch_sim_1d0 = function(N, d, x_init, trunc, a1, a0 = 0){
    
    z = rnorm(N+1+trunc)
    sigma = numeric(N+trunc)
    x = c(x_init, numeric(N))
    
    for(t in (trunc+1):(N+trunc)){
        print(t)
        
        aux = (t):(t-trunc)
        #print(var(x[aux]))
        
        gammas = zeta_1d0(x[aux], d, trunc, a1)$soma
        print(gammas)
        
        sigma[t] = sqrt(a0 + gammas)
        
        
        x[t] = sigma[t]*z[t]
        
    }
    serie = x[(trunc+1): (N+trunc)]
    
    return(serie)
    
}

a = figarch_sim_1d0(1000, 0.2, rnorm(100, sd = 5) , trunc = 100,  a1 = 0.4, a0 = 1)


plot.ts(a)



# FIGARCH(1,d,1) ----------------------------------------------------------
## ACHAR FORMA DE GERAR LISTA COM RECURSÃO

pi_k = function(d, k, alpha1, beta1){
    
    phi1 = alpha1 + beta1
    
    if (k == 1){
        return(gamma(k+1-d)/(gamma(k+2)*gamma(-d))*(phi1 - (k + 1 - d)/(k+2)))
    } 
    
    recurs = pi_k(d, k - 1, alpha1, beta1)
    #print(k)
    #print(recurs)
    
    out = recurs[k-1] * (phi1 - (k + 1 -d)/(k+2))*((k-d)/(phi1*(k+1) - (k - d)))
    
    vetor = c(recurs, out)
    
    return(vetor)
}

#plot.ts(pi_k(0.1, 50, 0.1, 0.1))


pi_k2 = function(d, k, phi1){
    
    out = gamma(k+1-d)/(gamma(k+2)*gamma(-d))*(phi1 - (k + 1 - d)/(k+2))
    
    return(out)
}



figarch_1d1 = function(N, # Número de observações que queremos gerar
                       d, # Parâmetro fracionário
                       trunc, #Onde truncaremos a série infinita sum_{k \geq 0} pi_k X_[t-2-k]^2
                       a1, #\alpha_1
                       b1, #\beta_1
                       a0, # \omega
                       x_init = rnorm(trunc + 2, sd = 5) # Valores iniciais para X, precisamos de trunc valores para que a série infinita seja calculada
                                                 # E dois valores a mais, pois a série começa em [t-2] (e o R não entende índices negativos)
                       ){
    
    # Sobre o trunc + 2 no argumento acima:
    # Para gerar cada \sigma^2_t da equação 2.11 do ProjetoTCC(090321), precisamos de \sigma^2_{t-1}, \pi_0, ... \pi_trunc,
    # e de X^2_{t-2} até X^2_{t-2-trunc}.
    # Suponha agora trunc = 10.
    # Como o R não suporta índices negativos, o primeiro valor 't' do nosso loop será {10 + 1 = 11}, para que possamos buscar os 10 valores anteriores.
    # Porém, como precisaremos de X^2_{t-2-trunc} = X^2_{11-2-10} = X^2_{-1}. Porém, como mencionado, o R não suporta índices negativos.
    # Por isso, são gerados 2 valores iniciais a mais, de forma que X_{t+2} = \sigma_t * Z_t
    
    z = rnorm(N+1+trunc) # Inicialmente, inovações N(0,1)
    sigma = numeric(N+trunc) # 50 valores a mais apenas por facilidade de leitura do código, serão removidos no fim da rotina.
    sigma[trunc] = 1 # Valor inicial para sigma
    x = c(x_init, numeric(N)) # Vetor de tamanho trunc + 2 + N; Função retornará apenas os N valores
    pi_ = pi_k2(d, 0:trunc, a1 + b1) # Função \pi, definida pela equação (2.9) do projetoTCC
    
    for(t in (trunc+1):(N+trunc)){
        
        
        aux = (t):(t-trunc) # No loop, sempre precisamos das quantidades X_[t - 2] até X_[t- trunc - 2] 
        
        
        gammas = sum(pi_*x[aux]^2)
       
        
        sigma[t] = sqrt(a0 + b1*sigma[t-1]^2 + (d + a1)*x[t+2-1]^2 + gammas)
        
        
        
        x[t+2] = sigma[t]*z[t]
        
    }
    serie = x[(trunc+1+2): (N+trunc+2)]
    dp = sigma[(trunc+1):(N+trunc)]
    
    return(list('y' = serie, 'se' = dp))
    
}

{
    N = 5000
    d = 0.5
    trunc = 100
    b1 = 0.3
    a1 = -b1
    phi1 = a1 + b1
    
    a0 = 0.1
    x_init = rnorm(trunc+2)
    
    cond1 = phi1 >= b1 - d & phi1 <= (2-d)/3
    cond2 = d*(phi1 - (1-d)/2) <= b1*(d+a1)
    cond3 = a1 + b1 < 1
    
    cond1 & cond2 & cond3
}




ret = figarch_1d1(N, d, trunc, a1, b1, a0)


plot.ts(ret$y)
plot.ts(ret$se)

acf1 = acf(ret$y, lwd = 1, lag.max = 100)
Box.test(ret$y, lag = 20, type = 'Ljung-Box')

pacf(ret$y, lwd = 5, lag.max = 200)
