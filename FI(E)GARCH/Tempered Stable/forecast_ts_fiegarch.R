# Script para gerar forecasts para FIEGARCH(p,d,q)
# Iniciado em 22/03/2022

source('FI(E)GARCH/Tempered Stable/tempered_stable_fiegarch_simulation.R', encoding = 'UTF-8')

# Bibliotecas -------------------------------------------------------------

library(tidyverse)


# Modelo 1 ----------------------------------------------------------------



serie_m1 = read_delim(
    "FI(E)GARCH/dados/analise_final/series_N5000_d_0.4_omega_-5.9_a_0.5_b_-0.6_theta_-0.22_gamma_0.37_m_10000_alpha_1.5.csv", 
    ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
    trim_ws = TRUE
)


ruido_m1 = read_delim(
    "FI(E)GARCH/dados/analise_final/ruido_N5000_d_0.4_omega_-5.9_a_0.5_b_-0.6_theta_-0.22_gamma_0.37_m_10000_alpha_1.5.csv", 
    ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
    trim_ws = TRUE
)


sigma_m1 = read_delim(
    "FI(E)GARCH/dados/analise_final/sigma_N5000_d_0.4_omega_-5.9_a_0.5_b_-0.6_theta_-0.22_gamma_0.37_m_10000_alpha_1.5.csv", 
    ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
    trim_ws = TRUE
)


ests_m1 = read_delim(
    "FI(E)GARCH/dados/analise_final/MC_N5000_d_0.4_omega_-5.9_a_0.5_b_-0.6_theta_-0.22_gamma_0.37_m_10000_alpha_1.5.csv", 
    ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
    trim_ws = TRUE
)


ruido_treino = ruido_m1[1:4950,]
ruido_teste = ruido_m1[4951:5000,]
sigma_treino = sigma_m1[1:4950,]
sigma_teste = sigma_m1[4951:5000, ]


ruido = ruido_treino[,2] %>% pull
pars = ests_m1[2,] %>% as.matrix



forecast_fiegarch_1d1(h, ruido, volat, pars){
    
    n = length(ruido)
    
    pred = numeric(h)
    
    ruido_alinhado = ruido*(abs(ruido) < 50) 
    
    
    
    
    
    d = pars[1] #= 0.43
    theta = pars[2] #= 0.11
    gamma = pars[3] #= -0.73
    omega = pars[4] #= -3
    a = pars[5] #
    b = pars[6] 
    
    lambdas = lambda(n-2+h, d, a, b)
    
    mu_z = mean(abs(ruido_alinhado))
    
    ruido_pred = c(ruido, pred)
    sigma_pred = pred
    
    for(h_ in 1:h){
        
        aux_lambda = h_:(h_ + n - 1)
        aux_z = n:1
        
        g = theta*ruido[aux_z] + gamma*(abs(ruido_alinhado[aux_z]) - sqrt(2/pi))#mu_z)
        
        sigma_pred[t] = sqrt(exp(omega + sum(lambdas[aux_lambda] * g)))
        print(sigma_pred[t])
        
        
    }
    
}

