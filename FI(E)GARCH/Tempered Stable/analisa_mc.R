{
    library(tidyverse)
}



# modelo1 -----------------------------------------------------------------


est_mc_m1 = read_delim(
    "FI(E)GARCH/dados/analise_final/MC_N5000_d_0.4_omega_-5.9_a_0.5_b_-0.6_theta_-0.22_gamma_0.37_m_10000_alpha_1.5.csv"
    , delim = ';'
    , locale = locale(decimal_mark = ',')
)

calc_mae = function(x, pars){
    
    out = x %>% 
        mutate(
            d = abs(d - pars[1]),
            theta = abs(theta - pars[2]),
            gamma = abs(gamma - pars[3]),
            omega = abs(omega - pars[4]),
            a1 = abs(a1 - pars[5]),
            b1 = abs(b1 - pars[6]),
            alpha = abs(alpha - pars[7]),
            ell = abs(ell - pars[8])
        )
    
    return(colMeans(out))
    
}

tabela_m1 = tibble(
    par = names(est_mc_m1),
    media = colMeans(est_mc_m1),
    sd = apply(est_mc_m1, 2, sd),
    bias = colMeans(est_mc_m1) - c(0.4, -0.22, 0.37, -5.9, 0.5, -0.6, 1.5, 1),
    mse = sd^2 + bias^2
)




# modelo2 -----------------------------------------------------------------


est_mc_m2 = read_csv2('FI(E)GARCH/dados/analise_final/MC_M2_CTS.csv')[1:500, ]



tabela_m2 = tibble(
    par = names(est_mc_m2),
    media = colMeans(est_mc_m2),
    sd = apply(est_mc_m2, 2, sd),
    bias = colMeans(est_mc_m2) - c(0.25, 0.3 , -1.2, -4.5, 0.4, -0.21, 0.7, 1),
    mse = sd^2 + bias^2
)

