#### stable_fiegarch #####



# libraries ---------------------------------------------------------------

{
    library(stabledist)
    library(tidyverse)
    #medias_aux = read_csv2('FI(E)GARCH/arquivos_auxiliares/simulacoes_abs_alpha_estavel_medias.csv') %>% 
    #    select(alpha, media)
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


g_func_stable = function(Z, 
                         theta, 
                         gamma,
                         alpha = 2,
                         method = c('sim', 'ana')
){

    
    # função feita para v.a.'s normais 
    if(method == 'sim'){
        E_X = medias_aux %>% 
            filter(alpha == alpha) %>% 
            select(media) %>% 
            pull()
    } else if (method == 'ana'){
        E_X = 2 * gamma(1 - 1/alpha)/pi
    }

    i = 0
    out = numeric(length(Z))
    
    for(z in Z){
        i = i+1
        if(z < 0){
            gamma = -gamma
        }
        out[i] = (theta + gamma)*z - gamma*E_X
    }
    print("g(Z) calculada")
    return(out)

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


lambda = function(m, d, a1, b1){
    
    p = c(1, a1)
    q = polyinv(c(1, b1), m = m)
    pq = poly.prod(a = p, b = q)
    
    s = sapply(0:m, pi_func, d)
    
    out = poly.prod(pq,s)[1:(m+1)]
    
    return(out)
    
}


lambda2 = function(m, d, a1, b1){
    ####### para debug #######
    #m = 2
    #d = 0.5
    #a1 = 0.1
    #b1 = 0.1
    #######
    
    
    # Proposition 2 -> Lopes and Prass 
    out = numeric(m+1) #lambda_{d, k} = out[k+1]
    a1 = c(-1, a1)
    b1 = c(-1, b1, numeric(m-1))
    #b1 = c(-1, b1)
    out[1] = 1
    
    for(k in 2:(m+1)){
        #print(paste0("k= ", k))
        parcelas_soma1 = numeric(length(1:(k-1)))
        for(i in 1:(k-1)){
            #print(paste0("i= ", i))
            parcelas_soma2 = numeric(length(1:(k - i + 1)))
            for (j in 1:(k - i + 1)){
                #print(paste0("j= ", j))
                parcelas_soma2[j] = b1[j]*pi_func(k-i-j+1, -d)
                #print(paste0('pi_func= ', pi_func(k-i-j+1, -d), ', b1_j = ', b1[j]))
            }
            soma_2 = sum(parcelas_soma2)
            parcelas_soma1[i] = out[i]*soma_2
            #print(paste0('multiplicado soma2 = ', soma_2, ' por lambda_',i, ' = ', out[i], ' resultando em ', parcelas_soma1[i]))
        }
        out[k] = sum(parcelas_soma1)
        if (k == 2){
            out[k] = out[k] - a1[k]
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
                        omega,
                        a1,
                        b1,
                        theta, 
                        gamma,
                        m,
                        alpha =  2){
    
    sigma = numeric(N+m)
    Z = rstable(N + m + 1, #Z_0 = z[1] -> Z_t = z[t + 1]
                alpha = alpha,
                beta = 0,
                gamma = 1, 
                delta = 0
    ) 
    x = numeric(N + m)
    print("Gerou variaveis iniciais, Começando a gerar valores da série.")

    g_z = g_func_stable(Z,
                        theta = theta, 
                        gamma = gamma, 
                        alpha = alpha, 
                        method = 'ana')
    
    lambdas = lambda(m, d, a1, b1)
    
    
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
    
    
    out = list(serie = y, sigma = sigma, residuos = y/sigma)
    return(out)
    
}



