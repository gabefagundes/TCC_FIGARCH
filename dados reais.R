library(tidyverse)
library(optimization)
source('FI(E)GARCH/Tempered Stable/tempered_stable_fiegarch_simulation.R', encoding = 'UTF-8')
sp500_treated <- read_csv("FI(E)GARCH/dados/dados_reais/sp500_treated.csv", 
                          col_names = FALSE)
source('funcoes_llh.R')
# dados full --------------------------------------------------------------

{
    sp500_full = read_delim("FI(E)GARCH/dados/dados_reais/SP50015min0830to15_ComMissingData.csv", delim = ' ') %>% 
        select(
            dtbase = Data,
            hr = Hora,
            close = Close
        ) %>% 
        mutate(
            dtbase = as.Date(dtbase, format = '%d/%m/%Y')
        )
    
    sp500_full$close = sp500_treated$X1
    
    sp500_diario = sp500_full %>% 
        group_by(dtbase) %>% 
        filter(
            hr == max(hr)
        ) %>% 
        ungroup() %>% 
        mutate(
            returns = c(NA, diff(close, lag = 1)),
            log_returns = c(NA, diff(log(close)))
        ) %>% 
        select(
            -hr
        )
    
}

ggplot(sp500_diario)+
    geom_line(aes(x = dtbase, y = close))+
    theme_bw()+
    xlab('Year') +
    ylab('Closing Price')


ggplot(sp500_diario)+
    geom_line(aes(x = dtbase, y = log_returns*100))+
    theme_bw()+
    xlab('Year') +
    ylab('Log-Returns')
# 
# ggplot(sp500_diario)+
#     geom_line(aes(x = dtbase, y = log_returns))
# 
# 




# 0d1 --------------------------------------------------------------------

estim_0d1 = 0
class(estim_0d1) = 'try-error'
while(class(estim_0d1) == 'try-error'){
    pars = c(runif(1,0,0.5), runif(3,-3,3), runif(1,-0.5,0.5))
    #pars = c(0.3200371, -1.2002767,  1.6866940,  1.6297912,  0.1225334)
    estim_0d1 = try(optim(pars, quasi_llf_fiegarch_0d1, method = 'BFGS', hessian = T,
                          control = list(fnscale = -1)))
}



quasi_llf_fiegarch_0d1(estim_0d1$par)

# saveRDS(estim_0d1, 'FI(E)GARCH/dados/analise_realdata/estimacao_0d1_bfgs.RDS')
# write_csv2(as.data.frame(estim_0d1$par), 'FI(E)GARCH/dados/analise_realdata/estimacao_0d1_bfgs.csv')


#estim_0d1 = readRDS('FI(E)GARCH/dados/analise_realdata/estimacao_0d1_bfgs.RDS')

n = length(sp500_diario$returns[2:1204])


aic0d1 = -2*estim_0d1$value + 2*5
bic0d1 = -2*estim_0d1$value + 5*log(n)
hqc0d1 = -2*estim_0d1$value + 2*5*log(log(n))


# estima ruido

z_0d1 = ruido_0d1(estim_0d1$par)
sigma_0d1 = ruido_0d1(estim_0d1$par, T)


cts_0d1 = 0
class(cts_0d1) = 'try-error'
while(class(cts_0d1) == 'try-error'){
    pars_cts = c(runif(1, 0, 2), rexp(1, 1))
    print(pars_cts)
    #pars = c(0.3200371, -1.2002767,  1.6866940,  1.6297912,  0.1225334)
    cts_0d1 = try(optim(pars_cts, mle_std_cts, Z= z_0d1, method = 'BFGS', hessian = T,
                          control = list(fnscale = -1)))
}


# 0,d,2 -------------------------------------------------------------------
estim_0d2 = 0
class(estim_0d2) = 'try-error'
while(class(estim_0d2) == 'try-error'){
    pars = c(runif(1,0,0.5), runif(3,-3,3), runif(2,-0.5,0.5))
    # 0.2580414 1.2681108 2.4570429 1.9619513 0.3196015 0.0594135
    estim_0d2 = try(optim(pars, quasi_llf_fiegarch_0d2, method = 'BFGS', hessian = T,
                          control = list(fnscale = -1)))
}

quasi_llf_fiegarch_0d2(estim_0d2$par)

# saveRDS(estim_0d2, 'FI(E)GARCH/dados/analise_realdata/estimacao_0d2_bfgs.RDS')
# write_csv2(as.data.frame(estim_0d2$par), 'FI(E)GARCH/dados/analise_realdata/estimacao_0d2_bfgs.csv')

p = 6

aic0d2 = -2*estim_0d2$value + 2*p
bic0d2 = -2*estim_0d2$value + p*log(n)
hqc0d2 = -2*estim_0d2$value + 2*p*log(log(n))

# 1d1 -----------------------------------------------------------------


# saveRDS(estim_1d1, 'FI(E)GARCH/dados/analise_realdata/estimacao_1d1_bfgs.RDS')
# write_csv2(as.data.frame(estim_1d1$par), 'FI(E)GARCH/dados/analise_realdata/estimacao_1d1_bfgs.csv')

pars = c(runif(1,0,0.5), runif(3,-3,3), runif(2,-0.5,0.5))
#pars = c(0.35328756, 1.31824174, 2.06578728, 0.96917391, 0.04526795, 0.06138524)

estim_1d1 = 0
class(estim_1d1) = 'try-error'
while(class(estim_1d1) == 'try-error'){
    pars = c(runif(1,0,0.5), runif(3,-3,3), runif(2,-0.5,0.5))
    estim_1d1 = try(optim(pars, quasi_llf_fiegarch_1d1, method = 'BFGS', hessian = T,
                          control = list(fnscale = -1)))
}


#estim_1d1 = read_csv2('FI(E)GARCH/dados/analise_realdata/estimacao_1d1_bfgs.RDS')

#quasi_llf_fiegarch_1d1(estim_1d1$`estim_1d1$par`)


p = 6

aic1d1 = -2*estim_1d1$value + 2*p
bic1d1 = -2*estim_1d1$value + p*log(n)
hqc1d1 = -2*estim_1d1$value + 2*p*log(log(n))


# 1d2 ---------------------------------------------------------------------


#pars = c(runif(1,0,0.5), runif(3,-3,3), runif(2,-0.5,0.5))
#0.04167871  2.89115032  2.57083022  0.61617029  0.06389144 -0.05147172 -0.11415316

estim_1d2 = 0
class(estim_1d2) = 'try-error'
while(class(estim_1d2) == 'try-error'){
    pars = c(runif(1,0,0.5), runif(3,-3,3), runif(3,-0.5,0.5))
    print(pars)
    estim_1d2 = try(optim(pars, quasi_llf_fiegarch_1d2, method = 'BFGS', hessian = T,
                          control = list(fnscale = -1)))
}


# saveRDS(estim_1d2, 'FI(E)GARCH/dados/analise_realdata/estimacao_1d2_bfgs.RDS')
# write_csv2(as.data.frame(estim_1d2$par), 'FI(E)GARCH/dados/analise_realdata/estimacao_1d2_bfgs.csv')

p = 7

aic1d2 = -2*estim_1d2$value + 2*p
bic1d2 = -2*estim_1d2$value + p*log(n)
hqc1d2 = -2*estim_1d2$value + 2*p*log(log(n))


# 2d1 ---------------------------------------------------------------------


estim_2d1 = 0
class(estim_2d1) = 'try-error'
while(class(estim_2d1) == 'try-error'){
    pars = c(runif(1,0,0.5), runif(3,-3,3), runif(3,-0.5,0.5))
    #0.2184608  0.3120820  2.3140227  2.1145509  0.1325746 -0.1641791 -0.1041725
    print(pars)
    estim_2d1 = try(optim(pars, quasi_llf_fiegarch_2d1, method = 'BFGS', hessian = T,
                          control = list(fnscale = -1)))
}

# 
# saveRDS(estim_2d1, 'FI(E)GARCH/dados/analise_realdata/estimacao_2d1_bfgs.RDS')
# write_csv2(as.data.frame(estim_2d1$par), 'FI(E)GARCH/dados/analise_realdata/estimacao_2d1_bfgs.csv')


p = 7

aic2d1 = -2*estim_2d1$value + 2*p
bic2d1 = -2*estim_2d1$value + p*log(n)
hqc2d1 = -2*estim_2d1$value + 2*p*log(log(n))

# 2d2 ------

estim_2d2 = 0
class(estim_2d2) = 'try-error'
while(class(estim_2d2) == 'try-error'){
    pars = c(runif(1,0,0.5), runif(3,-3,3), runif(4,-0.5,0.5))
    #0.2184608  0.3120820  2.3140227  2.1145509  0.1325746 -0.1641791 -0.1041725
    print(pars)
    estim_2d2 = try(optim(pars, quasi_llf_fiegarch_2d2, method = 'BFGS', hessian = T,
                          control = list(fnscale = -1)))
}

# 
# saveRDS(estim_2d2, 'FI(E)GARCH/dados/analise_realdata/estimacao_2d2_bfgs.RDS')
# write_csv2(as.data.frame(estim_2d2$par), 'FI(E)GARCH/dados/analise_realdata/estimacao_2d2_bfgs.csv')

p = 8

aic2d2 = -2*estim_2d2$value + 2*p
bic2d2 = -2*estim_2d2$value + p*log(n)
hqc2d2 = -2*estim_2d2$value + 2*p*log(log(n))




# forecast ----------------------------------------------------------------

estim_0d1 = readRDS('FI(E)GARCH/dados/analise_realdata/estimacao_0d1_bfgs.RDS')

data_teste = sp500_diario %>% 
    filter(
        dtbase >= '2009-07-24'
    )


serie_teste = data_teste$log_returns
z_0d1 = ruido_0d1(estim_0d1$par)


forecast_fiegarch_0d1= function (h, ruido, pars){
    
    n = length(ruido)
    
    pred = numeric(h)
    
    ruido_alinhado = ruido*(abs(ruido) < 50) 
    
    
    
    
    
    d = pars[1] #= 0.43
    theta = pars[2] #= 0.11
    gamma = pars[3] #= -0.73
    omega = pars[4] #= -3
    #a = pars[5] #
    b = pars[5] 
    
    lambdas = lambda(n-2+h, d, 0, b)
    
    mu_z = mean(abs(ruido_alinhado))
    
    ruido_pred = c(ruido, pred)
    sigma_pred_2 = pred
    
    for(h_ in 1:h){
        
        aux_lambda = h_:(h_ + n - 1)
        aux_z = n:1
        
        g = theta*ruido[aux_z] + gamma*(abs(ruido_alinhado[aux_z]) - mu_z)#sqrt(2/pi))
        
        sigma_pred_2[h_] = exp(omega + sum(lambdas[aux_lambda] * g))
        print(sigma_pred_2[h_])
        
        
    }
    
    return(sigma_pred_2)
}



pred_sp500 = forecast_fiegarch_0d1(h = length(serie_teste), ruido = z_0d1, pars = estim_0d1$par)


data_teste = data_teste[1:5,] %>% 
    mutate(
        log_returns2 = (log_returns*100)^2,
        pred_sigma2 = pred_sp500[1:5],
        bias = log_returns2 - pred_sigma2
    )



mse = var(data_teste$pred_sigma2) + mean(data_teste$bias)^2 

z_0d1 %>% shapiro.test

ks.test(z_0d1, dnorm)


