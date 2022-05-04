library(tidyverse)
#source('FI(E)GARCH/Tempered Stable/rCTS.R')


#alpha = 1.5
#ell = 1

#Z = rCTS(1000, alpha, ell, ell,  c = (gamma(2-alpha)*(ell^(alpha-2) + ell^(alpha-2)))^(-1))

cf_std_cts3 = function(t, alpha, ell){
    
    c = (gamma(2-alpha)*(ell^(alpha-2)+ell^(alpha-2)))^(-1)
    
    
    p1 = gamma(-alpha)*c*ell^alpha*((1 - (1i*t)/ell)^alpha - 1 + (1i*t*alpha)/ell)
    p2 = gamma(-alpha)*c*ell^alpha*((1 - (1i*t)/ell)^alpha - 1 + (1i*t*alpha)/ell)
    
    
    out = exp(p1+p2)
    
    return(out)
}

#tau = c(1.5, 1)


mle_std_cts = function(tau = numeric(2), Z, N = 2^15){
    
    alpha = tau[1]
    ell = tau[2]
    
    zmin = min(Z) - 0.1
    zmax = max(Z) + 0.1
    
    h = (zmax - zmin)/N
    
    k = 1:N
    n = 1:N
    
    x = (k - 1 - N/2)*h
    
    s = 1/(h*N)
    
    
    aux_p = s*(-1)^(k - 1 - (N/2))
    
    
    
    cf_calc_cts = (-1)^(n-1)*cf_std_cts3(2*pi*s*(n - 1 - (N/2)), alpha = 1.5, ell = 1)
    
    fft_cf_cts = fft(cf_calc_cts)
    
    
    final_cts = Re(aux_p*fft_cf_cts) #parte real para tirar problemas de precisão
    
    #plot(x, Re(final_cts), type = 'l' )
    
    
    lin_inter = approx(x, final_cts, xout = Z)

    
    
    out = sum(log(lin_inter$y), na.rm = T)
    
    #print(out)
    return(out)
}



#optim(c(1.5, 1), mle_std_cts, Z = Z, method = 'Nelder-Mead', hessian = F,
#      control = list(fnscale = -1))

arquivo_ruido = "FI(E)GARCH/dados/analise_final/ruido_N5000_d_0.4_omega_-5.9_a_0.5_b_-0.6_theta_-0.22_gamma_0.37_m_10000_alpha_1.5.csv"
ruido = read_delim(arquivo_ruido, delim = ';', locale = locale(decimal_mark = ',')) %>% 
    mutate_all(as.numeric)

estimativas = data.frame(
    alpha = numeric(ncol(ruido)),
    ell = numeric(ncol(ruido))
)


for(i in 1:ncol(ruido)){
    
    print(i)
    
    Z_aux = ruido[,i] %>% pull
    
    alpha_init = rnorm(1, mean = 1.37, sd = 0.1)
    ell_init = rnorm(1, mean = 1.08, sd = 0.17)
    
    par_init = c(alpha_init, ell_init)
    
    est_final = try(optim(par_init, mle_std_cts, Z = Z_aux, method = 'BFGS', hessian = F,
                          control = list(fnscale = -1)))
    
    if(class(est_final) == 'try-error'){
        estimativas[i, ] = par_init 
        print('erro')
        next 
    }
    
    estimativas[i,] = est_final$par
    
    
}


arquivo_mc = "FI(E)GARCH/dados/analise_final/MC_N5000_d_0.4_omega_-5.9_a_0.5_b_-0.6_theta_-0.22_gamma_0.37_m_10000_alpha_1.5.csv"

est_mc = read_delim(arquivo_mc, delim = ';', locale = locale(decimal_mark = ','))



teste_bind = est_mc[, 1:6] %>% 
    cbind(estimativas)



write_csv2(teste_bind, arquivo_mc)




# M2 ----------------------------------------------------------------------


arquivo_ruido = "FI(E)GARCH/dados/analise_final/ruido_M2.csv"
ruido = read_delim(arquivo_ruido, delim = ';', locale = locale(decimal_mark = ',')) %>% 
    mutate_all(as.numeric)

estimativas = data.frame(
    alpha = numeric(ncol(ruido)),
    ell = numeric(ncol(ruido))
)


for(i in 1:ncol(ruido)){
    
    print(i)
    
    Z_aux = ruido[,i] %>% pull
    
    alpha_init = rnorm(1, mean = 0.792, sd = 0.07)
    ell_init = rnorm(1, mean = 1.23, sd = 0.14)
    
    par_init = c(alpha_init, ell_init)
    
    est_final = try(optim(par_init, mle_std_cts, Z = Z_aux, method = 'BFGS', hessian = F,
                          control = list(fnscale = -1)))
    
    if(class(est_final) == 'try-error'){
        estimativas[i, ] = par_init 
        print('erro')
        next 
    }
    
    estimativas[i,] = est_final$par
    
    
}

ests_m2 = read_csv2('FI(E)GARCH/dados/analise_final/MC_M2.csv')


ests_m2_cts = ests_m2 %>% 
    bind_cols(estimativas)


write_csv2(ests_m2_cts, 'FI(E)GARCH/dados/analise_final/MC_M2_CTS.csv')
