# Script para gerar amostras de processos FIEGARCH com inovações TS #
# Script iniciado em 15/03 # 



# Bibliotecas -------------------------------------------------------------

{
    
    library(tidyverse)
    source('FI(E)GARCH/Tempered Stable/rCTS.R')
    source('FI(E)GARCH/utils/poly_functions.R')
}
library(readr)
TS_absmean_alpha1_5 <- read_delim("FI(E)GARCH/dados/TS_absmean_alpha1.5.csv", 
                                  ";", escape_double = FALSE, col_types = cols(X1 = col_skip()), 
                                  locale = locale(decimal_mark = ","), 
                                  trim_ws = TRUE)


TS_absmean_alpha0_7 = read_csv2("FI(E)GARCH/dados/TS_absmean_alpha0.7.csv")


# medias_aux = read_csv2('FI(E)GARCH/arquivos_auxiliares/medias_abs_tempered_stable2.csv') %>% 
#     rename('alpha' = X1,
#            value = ell_1) %>% 
#     mutate(alpha = str_remove(alpha, 'alpha_') %>% as.numeric)




media_alpha_1.5_ell_1 = mean(TS_absmean_alpha1_5$x)
media_alpha_0.7_ell_1 = mean(TS_absmean_alpha0_7$xx)



# Funcoes -----------------------------------------------------------------


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

g_f_ts = function(Z,
                  theta,
                  gamma,
                  alpha_ts){
    
    # E_X = medias_aux %>% 
    #     filter(alpha == alpha_ts) %>% 
    #     select(value) %>% 
    #     pull()
    
    E_X = media_alpha_1.5_ell_1
    #E_X = media_alpha_0.7_ell_1
    
    out = numeric(length(Z))
    
    for (t in seq_along(Z)){
        
        out[t] = theta*Z[t] + gamma*(abs(Z[t]) - E_X)
    }
    
    print("g(Z) Calculada.")
    return(out)
}



# Lambdas 
# usa funcoes do script poly_functions.R

lambda = function(m, d, a, b){
    
    p = c(1, a)
    q = polyinv(c(1, b), m = m)
    pq = poly.prod(a = p, b = q)
    
    s = sapply(0:m, pi_func, d)
    
    out = poly.prod(pq,s)[1:(m+1)]
    
    return(out)
    
}



# Geração Processo FIEGARCH -----------------------------------------------

fiegarch_TS = function(N = 500,
                       d,
                       omega,
                       a,
                       b,
                       theta,
                       gamma,
                       m,
                       alpha_ts,
                       lamb = NULL){
    print("Gerando variaveis iniciais.")
    
    
    sigma = numeric(N+m)
    Z = rCTS(N + m + 1,
            alpha = alpha_ts,
            ell_plus =  1,
            ell_minus = 1,
            c = 1^(alpha_ts-2)/(2*gamma(2-alpha_ts))
            )
    x = numeric(N+m)
    
    
    print('Calculando g(Z)')
    g_z = g_f_ts(Z,
                 theta = theta,
                 gamma = gamma,
                 alpha = alpha_ts)
    
    print('Gerando lambdas')
    if (is.null(lamb)){
        lambdas = lambda(m, d, a, b)
    } else {
        lambdas = lamb 
    }
    
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
