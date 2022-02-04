#### stable_fiegarch #####



# libraries ---------------------------------------------------------------

{
    library(stabledist)
    library(tidyverse)
    medias_aux = read_csv2('FI(E)GARCH/arquivos_auxiliares/simulacoes_abs_alpha_estavel_medias.csv') %>% 
        select(alpha, media)
}


# funções auxiliares ------------------------------------------------------
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


E_X = medias_aux %>% 
    filter(alpha == 1.8) %>% 
    select(media) %>% 
    pull()

g_func_stable = function(z, 
                  theta, 
                  gamma,
                  alpha_stable = 2){
    
    # função feita para v.a.'s normais 
    
    if(z < 0){
        gamma = -gamma
    }
    
    out = (theta + gamma)*z - gamma*E_X
}


polyinv <- function(poli = NULL, m = 1000){
    phi = c(poli[-1], rep(0, m-length(poli)+1))
    psi = numeric(m+1)
    psi[1] = 1 # psi_0
    for(i in 2:(m+1)) psi[i] = -sum(rev(phi[1:(i-1)])*psi[1:(i-1)])
    return(psi)
}


poly.prod <- function(a = 1,b = 1){
    
    k1 = length(a)-1
    k2 = length(b)-1
    k = k1+k2
    p = numeric(k+1)
    
    ax = rep(0,k+1); ax[1:(k1+1)] = a
    bx = rep(0,k+1); bx[1:(k2+1)] = b
    
    for(i in 0:k){
        p[i+1] = sum(ax[1:(i+1)]*bx[(i+1):1])
    }
    
    return(p)
}


lambda = function(trunc, d, alpha, beta){
    
    p = c(1, alpha)
    q = polyinv(c(1, beta), m = trunc)
    pq = poly.prod(a = p, b = q)
    
    s = sapply(0:trunc, pi_func, d)
    
    out = poly.prod(pq,s)[1:(trunc+1)]
    
    return(out)
    
}


lambda2 = function(trunc, d, alpha, beta){
    ####### para debug #######
    #trunc = 2
    #d = 0.5
    #alpha = 0.1
    #beta = 0.1
    #######
    
    
    # Proposition 2 -> Lopes and Prass 
    out = numeric(trunc+1) #lambda_{d, k} = out[k+1]
    alpha = c(-1, alpha)
    beta = c(-1, beta, numeric(trunc-1))
    #beta = c(-1, beta)
    out[1] = 1
    
    for(k in 2:(trunc+1)){
        #print(paste0("k= ", k))
        parcelas_soma1 = numeric(length(1:(k-1)))
        for(i in 1:(k-1)){
            #print(paste0("i= ", i))
            parcelas_soma2 = numeric(length(1:(k - i + 1)))
            for (j in 1:(k - i + 1)){
                #print(paste0("j= ", j))
                parcelas_soma2[j] = beta[j]*pi_func(k-i-j+1, -d)
                #print(paste0('pi_func= ', pi_func(k-i-j+1, -d), ', beta_j = ', beta[j]))
            }
            soma_2 = sum(parcelas_soma2)
            parcelas_soma1[i] = out[i]*soma_2
            #print(paste0('multiplicado soma2 = ', soma_2, ' por lambda_',i, ' = ', out[i], ' resultando em ', parcelas_soma1[i]))
        }
        out[k] = sum(parcelas_soma1)
        if (k == 2){
            out[k] = out[k] - alpha[k]
        } else {
            out[k] = out[k]
        }
        #print(out[k])
    }
    return(out)
}


# função principal --------------------------------------------------------

fiegarch_1d1 = function(N = 500,
                        d, 
                        a0,
                        alpha,
                        beta,
                        theta, 
                        gamma,
                        trunc,
                        alpha_stable =  2){
    
    sigma = numeric(N+trunc)
    z = rstable(N + trunc + 1, #Z_0 = z[1] -> Z_t = z[t + 1]
                alpha = alpha_stable,
                beta = 0,
                gamma = 1, 
                delta = 0
                ) 
    x = numeric(N + trunc)
    print("Gerou variaveis iniciais, Começando a gerar valores da série.")
    
    
    for (t in (trunc+1):(N+trunc)){
        print(t)
        aux = (t):(t-trunc)
        
        sigma[t] = (a0 + sum(lambda(trunc, d, alpha, beta)*sapply(z[aux], g_func_stable, theta = theta, gamma = gamma, alpha_stable = alpha_stable))) %>% 
            exp() %>% 
            sqrt()
        
        x[t] = sigma[t]*z[t + 1]
        
    }
    print("Series e Sigmas Gerados")
    y = x[(trunc+1):(trunc + N)]
    sigma = sigma[(trunc+1):(trunc + N)]
    
    
    out = list(serie = y, sigma = sigma)
    return(out)
    
}


# testes ------------------------------------------------------------------


{
    N = 1000
    d = 0.7
    a0 = -9
    alpha = 0.1
    beta = 0.68
    theta = -0.21661
    gamma = 0.27
    trunc = 100
    alpha_stable = 1.8
}

teste = fiegarch_1d1(N, d, a0, alpha, beta, theta, gamma, trunc)

{
    par(mfrow = c(1,2))
    plot.ts(teste$serie)
    plot.ts(teste$sigma^2)
}


Box.test(teste$serie, lag = 20, type = 'Box-Pierce')
