{
    library(tidyverse)
    source('FI(E)GARCH/Tempered Stable/tempered_stable_fiegarch_simulation.R', encoding = 'UTF-8')
    
}
# script para estimar valores de cts --------------------------------------


m1_ests <- read_delim(
    "FI(E)GARCH/dados/MC_N5000_d_0.45_omega_-5.9_a_0.5_b_-0.6_theta_-0.22_gamma_0.37_m_10000_alpha_1.5.csv", 
    ";", 
    escape_double = FALSE,
    locale = locale(decimal_mark = ","), 
    trim_ws = TRUE
)


m1_series <- read_delim(
    "FI(E)GARCH/dados/serie_N5000_d_0.45_omega_-5.9_a_0.5_b_-0.6_theta_-0.22_gamma_0.37_m_10000_alpha_1.5.csv", 
    ";", 
    escape_double = FALSE,
    locale = locale(decimal_mark = ","), 
    trim_ws = TRUE
)

m1_ests = m1_ests %>% 
    mutate(
        fl_exclui = ifelse(d < 0 | d > 0.7, 1, 0)
    )



m1_ests_final = m1_ests %>% 
    filter(!fl_exclui) %>% 
    select(-fl_exclui)
m1_series_final = m1_series[,!m1_ests$fl_exclui]


# write_csv2(
#     m1_ests_final,
#     "FI(E)GARCH/dados/analise_final/MC_N5000_d_0.45_omega_-5.9_a_0.5_b_-0.6_theta_-0.22_gamma_0.37_m_10000_alpha_1.5.csv"
# )
# 
# write_csv2(
#     m1_series_final,
#     "FI(E)GARCH/dados/analise_final/series_N5000_d_0.45_omega_-5.9_a_0.5_b_-0.6_theta_-0.22_gamma_0.37_m_10000_alpha_1.5.csv"
# )

hist(m1_ests_final$d)
hist(m1_ests_final$theta)
hist(m1_ests_final$gamma)
hist(m1_ests_final$omega)
hist(m1_ests_final$a1)
hist(m1_ests_final$b1)

arquivo_ruido = "FI(E)GARCH/dados/analise_final/ruido_N5000_d_0.45_omega_-5.9_a_0.5_b_-0.6_theta_-0.22_gamma_0.37_m_10000_alpha_1.5.csv"
arquivo_volatilidade = "FI(E)GARCH/dados/analise_final/sigma_N5000_d_0.45_omega_-5.9_a_0.5_b_-0.6_theta_-0.22_gamma_0.37_m_10000_alpha_1.5.csv"


# criacao ruidos ----------------------------------------------------------

for(i in 7:ncol(m1_series_final)){
    
    pars = m1_ests_final[i,] %>% as.matrix
    serie = m1_series_final[,i] %>% pull
    
    n = length(serie)
    d = pars[1] 
    theta = pars[2] 
    gamma = pars[3] 
    omega = pars[4] 
    a1 = pars[5] 
    b1 = pars[6]
    
    
    lambdas = lambda(n-1, d, a1, b1)
    
    sigma = numeric(n)
    g_z = numeric(2*n - 1)
    z = numeric(n)
    
    sigma[1] = exp(0.5*omega)
    z[1] = serie[1]/sigma[1]
    g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
    
    
    for(t  in 2:n){
        
        aux = (t-1):(t-1-(n-1)) + n
        
        (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
        #x[t]
        (z[t] = serie[t]/sigma[t])
        (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
    }
    
    ruido = read_delim(arquivo_ruido, delim = ';', locale = locale(decimal_mark = ','))
    volat = read_delim(arquivo_volatilidade, delim = ';', locale = locale(decimal_mark = ','))
    
    
    ruido2 = ruido %>% 
        cbind(z)
    volat2 = volat %>% 
        cbind(sigma)
    
    colnames(ruido2) =  colnames(volat2) = paste0('x', 1:ncol(ruido2))
    
   
    write_csv2(ruido2, arquivo_ruido)
    
    write_csv2(volat2, arquivo_volatilidade)
}




# Modelo 2 ----------------------------------------------------------------

m2_ests <- read_delim(
    "FI(E)GARCH/dados/analise_final/MC_M2.csv", 
    ";", 
    escape_double = FALSE,
    locale = locale(decimal_mark = ","), 
    trim_ws = TRUE
) 


set.seed(1234)
fl_exclui = sample(c(rep(0, 269),rep(1, 500)), size = nrow(m2_ests), replace = F)


m2_ests_final = m2_ests %>% 
    #mutate(
    #    d = d + rnorm(nrow(.), mean = 0, sd = 0.001),
    #    theta = theta + rnorm(nrow(.), mean = 0, sd = 0.001),
    #    gamma = gamma + rnorm(nrow(.), mean = 0, sd = 0.001),
    #    omega = omega + rnorm(nrow(.), mean = 0, sd = 0.001),
    #    b1 = b1 + rnorm(nrow(.), mean = 0, sd = 0.001),
    #    b2 = b2 + rnorm(nrow(.), mean = 0, sd = 0.001)
    #) %>% 
    mutate(fl_exclui = fl_exclui) %>% 
    filter(
        fl_exclui == 1
    ) %>% 
    select(
        -fl_exclui
    )


summary(m2_ests_final)

m2_series <- read_delim(
    "FI(E)GARCH/dados/analise_final/series_M2.csv", 
    ";", 
    escape_double = FALSE,
    locale = locale(decimal_mark = ","), 
    trim_ws = TRUE
) %>% 
    select(-x1)





m2_series_final = m2_series[, as.logical(fl_exclui)]

arquivo_ruido = "FI(E)GARCH/dados/analise_final/ruido_M2.csv"
arquivo_volatilidade = "FI(E)GARCH/dados/analise_final/sigma_M2.csv"


for(i in 7:ncol(m2_series)){
    
    print(i)
    
    pars = m2_ests[i,] %>% as.matrix
    serie = m2_series[,i] %>% pull
    
    n = length(serie)
    d = pars[1] 
    theta = pars[2] 
    gamma = pars[3] 
    omega = pars[4] 
    b1 = pars[5] 
    b2 = pars[6]
    
    
    lambdas = lambda(n-1, d, 0, c(b1, b2))
    
    sigma = numeric(n)
    g_z = numeric(2*n - 1)
    z = numeric(n)
    
    sigma[1] = exp(0.5*omega)
    z[1] = serie[1]/sigma[1]
    g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
    
    
    for(t  in 2:n){
        
        aux = (t-1):(t-1-(n-1)) + n
        #print(serie[t])
        (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
        #x[t]
        (z[t] = serie[t]/sigma[t])
        (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
        #print(z[t])
        #print(sigma[t])
        #readline()
    }
    
    #ruido = read_delim(arquivo_ruido, delim = ';', locale = locale(decimal_mark = ','))
    #volat = read_delim(arquivo_volatilidade, delim = ';', locale = locale(decimal_mark = ','))
    
    #ruido2 = as.data.frame(z)
    ruido2 = ruido2 %>% 
        cbind(z)
    #volat2 = as.data.frame(sigma)
    volat2 = volat2 %>% 
        cbind(sigma)
    
    colnames(ruido2) =  colnames(volat2) = paste0('x', 1:ncol(ruido2))
    
    
    write_csv2(ruido2, arquivo_ruido)
    
    write_csv2(volat2, arquivo_volatilidade)
}



fl_na = numeric(769)


for(i in 1:769){
    
    xx_ruido = ruido2[,i] 
    
    fl_na[i] = ifelse(any(is.na(xx_ruido)), 1, 0)
    
    
}


serie_final_m2 = m2_series[,!fl_na]
ruido_final_m2 = ruido2[,!fl_na]
sigma_final_m2 = volat2[,!fl_na]
ests_final_m2 = m2_ests[!fl_na,]

write_csv2(serie_final_m2, 'FI(E)GARCH/dados/analise_final/serie_M2.csv')
write_csv2(ruido_final_m2, 'FI(E)GARCH/dados/analise_final/ruido_M2.csv')
write_csv2(sigma_final_m2, 'FI(E)GARCH/dados/analise_final/sigma_M2.csv')
write_csv2(ests_final_m2 , 'FI(E)GARCH/dados/analise_final/MC_M2.csv')

