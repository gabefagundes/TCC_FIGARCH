library(magrittr)



#pi_func = function(k, d){
#    
#    # A PROXIMA LINHA FAZ COM QUE O SINAL DA FUNÇÃO FIQUE MUDANDO -- ESTÁ ERRADO! 
#    d = -d # \delta_{-d, k} = \pi_{d, k} 
#    
#    
#    if(k == 0){
#        return(1)
#    }
#    
#    out = pi_func(k - 1, d) * (k - 1 - d)/k
#    print(out)
#    return(out)
#}


pi_func = function(k, d){
    d = -d
    
    for(i in 0:k){
        if(i == 0){
            out = 1
        } else{
            out = out*(i - 1 - d)/i
        }
        print(out)
    }
    return(out)
}

g_func = function(z, theta, gamma){
    
    # função feita para v.a.'s normais 
    
    
    if(z < 0){
        gamma = -gamma
    }
    
    out = (theta + gamma)*z - gamma*sqrt(2/pi)
}


fiegarch_0d0 = function(N = 500,
                        d, 
                        a0,
                        theta, 
                        gamma,
                        trunc){
    
    sigma = numeric(N+trunc)
    z = rnorm(N + trunc + 1) #Z_0 = z[1] -> Z_t = z[t + 1]
    x = numeric(N + trunc)
    
    for (t in (trunc+1):(N+trunc)){
        
        aux = (t):(t-trunc)
        
        sigma[t] = (a0 + sum( sapply(0:trunc, pi_func, d = d)* 
                                  sapply(z[aux], g_func, theta = theta, gamma = gamma))
        ) %>% 
            exp() %>% 
            sqrt()
        
        x[t] = sigma[t]*z[t + 1]
        
    }
    
    y = x[(trunc+1):(trunc + N)]
    sigma = sigma[(trunc+1):(trunc + N)]
    
    
    out = list(serie = y, sigma = sigma)
    return(out)
    
}

# testes

{
    N = 2000
    d = 0.35
    a0 = 5
    theta = -0.25
    gamma = 0.24
    trunc = 20
}

teste = fiegarch_0d0(N, d, a0, theta, gamma, trunc)


par(mfrow = c(1,2))
plot.ts(teste$serie)
plot.ts(teste$sigma)



Box.test(teste$serie, lag = 20, type = 'Ljung-Box')

plot.ts(lapply(0:100, delta_func, d= -0.15))


# FIEGARCH (1,d,1) --------------------------------------------------------
polyinv <- function(poli = NULL, m = 1000){
    phi = c(poli[-1], rep(0, m-length(poli)+1))
    psi = numeric(m+1)
    psi[1] = 1 # psi_0
    for(i in 2:(m+1)) psi[i] = -sum(rev(phi[1:(i-1)])*psi[1:(i-1)])
    return(psi)
}



#-----------------------------------------------------------------------------
#
#  Funçõees para calcular o produto de dois polinômios 
#
#   a(z) = 1 + a1*z + a2*z^2 + ... + ak*z^k1
# 
#   b(z) = 1 + b1*z + b2*z^2 + ... + bk*z^k2
#
#   a(z)b(z) = p(z) = 1 + p1*z + ... + pk*zk, onde k = k1+k2
#
#-----------------------------------------------------------------------------
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



lambda = function(alpha, beta, trunc, d){
    
    p = c(1, alpha)
    q = polyinv(c(1, beta), m = trunc)
    pq = poly.prod(a = p, b = q)
    
    s = sapply(0:trunc, pi_func, d)
    
    out = poly.prod(pq,s)[1:(trunc+1)]
    
}

lambda2 = function(trunc, d, alpha, beta){
    # Proposition 2 -> Lopes and Prass 
    out = numeric(trunc+1) #lambda_{d, k} = out[k+1]
    alpha = c(-1, alpha)
    beta = c(-1, beta)
    out[1] = 1
    
    for(k in 2:(trunc+1)){
        for(i in 1:(k-1)){
            for (j in 1:(k - i + 1)){
                out[k] = out[k] + beta[j]*pi_func(-d, k-i-j+1)
            }
            out[k] = out[k]*out[i]
        }
        if (k == 2){
            out[k] = out[k] + alpha[k]
        } else {
            out[k] = 0
        }
    }
}


fiegarch_1d1 = function(N = 500,
                        d, 
                        a0,
                        alpha,
                        beta,
                        theta, 
                        gamma,
                        trunc){
    
    sigma = numeric(N+trunc)
    z = rnorm(N + trunc + 1) #Z_0 = z[1] -> Z_t = z[t + 1]
    x = numeric(N + trunc)
    
    for (t in (trunc+1):(N+trunc)){
        
        aux = (t):(t-trunc)
        
        sigma[t] = (a0 + sum( lambda(alpha, beta, trunc, d)*sapply(z[aux], g_func, theta = theta, gamma = gamma))
        ) %>% 
            exp() %>% 
            sqrt()
        
        x[t] = sigma[t]*z[t + 1]
        
    }
    
    y = x[(trunc+1):(trunc + N)]
    sigma = sigma[(trunc+1):(trunc + N)]
    
    
    out = list(serie = y, sigma = sigma)
    return(out)
    
}


# testes ------------------------------------------------------------------


{
    N = 400
    d = 0.8
    a0 = -7
    alpha = 0
    beta = 0.68
    theta = -0.21661
    gamma = 0.27
    trunc = 100
}

teste = fiegarch_1d1(N, d, a0, alpha, beta, theta, gamma, trunc)

{
    par(mfrow = c(1,2))
    plot.ts(teste$serie)
    plot.ts(teste$sigma^2)
}


Box.test(teste$serie, lag = 20, type = 'Ljung-Box')

