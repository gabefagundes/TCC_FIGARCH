


# algoritmo mittkikk ------------------------------------------------------

cf_gauss = function(t, mu, sd){
    
    return(exp(1i*mu*t -(t*sd)^2/2))
}


cf_std_cts = function(s, alpha, ell){
    
    
    p1 = ell^2/(alpha^2 - alpha)
    p2 = (1 - (1i*s)/ell)^alpha - 1
    
    expoente = p1*p2
    
    out = exp(expoente)
    
    return(out)
} 

cf_std_cts2 = function(t, alpha, ell){
    
    c = (gamma(2-alpha)*(ell^(alpha-2)+ell^(alpha-2)))^(-1)
    
    p1 = 1i*t*gamma(1-alpha)*(c*ell^(alpha-1) - c*ell^(alpha-1))
    p2 = gamma(-alpha)*c*ell^alpha*((1 - (1i*t)/ell)^alpha - 1)
    p3 = gamma(-alpha)*c*ell^alpha*((1 - (1i*t)/ell)^alpha - 1)
 
    
    out = exp(p1+p2+p3)
    
    return(out)
}

cf_std_cts3 = function(t, alpha, ell){
    
    c = (gamma(2-alpha)*(ell^(alpha-2)+ell^(alpha-2)))^(-1)
    
    
    p1 = gamma(-alpha)*c*ell^alpha*((1 - (1i*t)/ell)^alpha - 1 + (1i*t*alpha)/ell)
    p2 = gamma(-alpha)*c*ell^alpha*((1 - (1i*t)/ell)^alpha - 1 + (1i*t*alpha)/ell)
    
    
    out = exp(p1+p2)
    
    return(out)
}

cf_std_cts3(0, 1.5, 1)

mu = 0
sd = 1 

N = 2^15

h = 8/N

k = 1:N

n = 1:N

x = sd*(k - 1 - (N/2))*h + mu


s = 1/(h*N)


aux_p = s*(-1)^(k - 1 - (N/2))

# normal
cf_calc = (-1)^(n-1)*cf_gauss(2*pi*s*(n - 1 - (N/2)))

fft_cf = fft(cf_calc) 


final = aux_p*fft_cf
plot(x, Re(final), type = 'l')


# tempered stable

cf_calc_cts = (-1)^(n-1)*cf_std_cts3(2*pi*s*(n - 1 - (N/2)), alpha = 1.5, ell = 1)

fft_cf_cts = fft(cf_calc_cts)


final_cts = aux_p*fft_cf_cts

plot(x, Re(final_cts), type = 'l' )

