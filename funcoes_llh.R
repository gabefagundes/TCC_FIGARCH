
# 1,d,q -------------------------------------------------------------------


quasi_llf_fiegarch_1d1 = function(pars = numeric(6)){
    
    
    x =  sp500_diario$log_returns[-1][1:1204]
    
    n = length(x)
    sigma = numeric(n)
    g_z = numeric(2*n - 1)
    z = numeric(n)
    
    
    if(length(pars) != 6){
        stop("Vetor de parâmetros suprido tem comprimento diferente de 6.")
    }
    names(pars) = c('d', 'theta', 'gamma', 'omega', 'a', 'b')
    
    
    
    d = pars['d'] #= 0.43
    theta = pars['theta'] #= 0.11
    gamma = pars['gamma'] #= -0.73
    omega = pars['omega'] #= -3
    a = pars['a'] #
    b = pars['b'] #= 0.54
    # 
    # if(d > 0.5){
    #     d = 0.49
    # }
    
    
    if( d <= 0){
        return(NA)
    }
    
    if (d >= 0.5){
        return(NA)
    }
    

    lambdas = lambda(n-1, d, a, b)
    plot.ts(lambdas)
    
    

    sigma[1] = exp(0.5*omega)
    z[1] = x[1]/sigma[1]
    g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
    
    
    for(t  in 2:n){
        
        aux = (t-1):(t-1-(n-1)) + n
        
        (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
        x[t]
        (z[t] = x[t]/sigma[t])
        (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
    }
    
    if( is.na(sum(sigma)) ){
        return(NA)
    }
    
    #plot.ts(z)
    #hist(z)
    plot.ts(sigma )
    
    log_func = -n/2*log(2/pi) - 1/2*sum(log(sigma[t]^2) + (x[t]^2)/(sigma[t]^2))
    print(log_func)
    
    if(log_func == Inf | log_func == -Inf | is.na(log_func)){
        return(NA)
    }
    
    return(log_func)
}


ruido_1d1 = function(pars){
  
  x =  sp500_diario$log_returns[-1][1:1204]
  
  n = length(x)
  sigma = numeric(n)
  g_z = numeric(2*n - 1)
  z = numeric(n)
  
  
  if(length(pars) != 6){
    stop("Vetor de parâmetros suprido tem comprimento diferente de 6.")
  }
  names(pars) = c('d', 'theta', 'gamma', 'omega', 'a', 'b')
  
  
  
  d = pars['d'] #= 0.43
  theta = pars['theta'] #= 0.11
  gamma = pars['gamma'] #= -0.73
  omega = pars['omega'] #= -3
  a = pars['a'] #
  b = pars['b'] #= 0.54
  # 
  # if(d > 0.5){
  #     d = 0.49
  # }
  
  
  if( d <= 0){
    return(NA)
  }
  
  if (d >= 0.5){
    return(NA)
  }
  
  
  lambdas = lambda(n-1, d, a, b)
  plot.ts(lambdas)
  
  
  
  sigma[1] = exp(0.5*omega)
  z[1] = x[1]/sigma[1]
  g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
  
  
  for(t  in 2:n){
    
    aux = (t-1):(t-1-(n-1)) + n
    
    (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
    x[t]
    (z[t] = x[t]/sigma[t])
    (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
  }
  
  return(z)
}

quasi_llf_fiegarch_1d2 = function(pars = numeric(7)){
    
    
    x =  sp500_diario$log_returns[-1][1:1204]
    
    n = length(x)
    sigma = numeric(n)
    g_z = numeric(2*n - 1)
    z = numeric(n)
    
    
    if(length(pars) != 7){
        stop("Vetor de parâmetros suprido tem comprimento diferente de 6.")
    }
    names(pars) = c('d', 'theta', 'gamma', 'omega', 'a', 'b1', 'b2')
    
    
    
    d = pars['d'] #= 0.43
    theta = pars['theta'] #= 0.11
    gamma = pars['gamma'] #= -0.73
    omega = pars['omega'] #= -3
    a = pars['a'] #
    b1 = pars['b1'] #= 0.54
    b2 = pars['b2']
    b = c(b1,b2)
    # 
    # if(d > 0.5){
    #     d = 0.49
    # }
    
    
    if( d <= 0){
        return(NA)
    }
    
    if (d >= 0.5){
        return(NA)
    }
    
    
    lambdas = lambda(n-1, d, a, b)
    plot.ts(lambdas)
    
    
    
    sigma[1] = exp(0.5*omega)
    z[1] = x[1]/sigma[1]
    g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
    
    
    for(t  in 2:n){
        
        aux = (t-1):(t-1-(n-1)) + n
        
        (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
        x[t]
        (z[t] = x[t]/sigma[t])
        (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
    }
    
    if( is.na(sum(sigma)) ){
        return(NA)
    }
    
    #plot.ts(z)
    #hist(z)
    plot.ts(sigma)
    
    log_func = -n/2*log(2/pi) - 1/2*sum(log(sigma[t]^2) + (x[t]^2)/(sigma[t]^2))
    print(log_func)
    
    if(log_func == Inf | log_func == -Inf | is.na(log_func)){
        return(NA)
    }
    
    return(log_func)
}


ruido_1d2 = function(pars){
  
  x =  sp500_diario$log_returns[-1][1:1204]
  
  n = length(x)
  sigma = numeric(n)
  g_z = numeric(2*n - 1)
  z = numeric(n)
  
  
  if(length(pars) != 7){
    stop("Vetor de parâmetros suprido tem comprimento diferente de 6.")
  }
  names(pars) = c('d', 'theta', 'gamma', 'omega', 'a', 'b1', 'b2')
  
  
  
  d = pars['d'] #= 0.43
  theta = pars['theta'] #= 0.11
  gamma = pars['gamma'] #= -0.73
  omega = pars['omega'] #= -3
  a = pars['a'] #
  b1 = pars['b1'] #= 0.54
  b2 = pars['b2']
  b = c(b1,b2)
  # 
  # if(d > 0.5){
  #     d = 0.49
  # }
  
  
  if( d <= 0){
    return(NA)
  }
  
  if (d >= 0.5){
    return(NA)
  }
  
  
  lambdas = lambda(n-1, d, a, b)
  plot.ts(lambdas)
  
  
  
  sigma[1] = exp(0.5*omega)
  z[1] = x[1]/sigma[1]
  g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
  
  
  for(t  in 2:n){
    
    aux = (t-1):(t-1-(n-1)) + n
    
    (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
    x[t]
    (z[t] = x[t]/sigma[t])
    (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
  }
  
  return(z)
}

# 0,d,q -------------------------------------------------------------------



quasi_llf_fiegarch_0d1 = function(pars = numeric(5)){
    
    
    x =  sp500_diario$log_returns[-1][1:1204]
    
    n = length(x)
    sigma = numeric(n)
    g_z = numeric(2*n - 1)
    z = numeric(n)
    
    
    if(length(pars) != 5){
        stop("Vetor de parâmetros suprido tem comprimento diferente de 6.")
    }
    names(pars) = c('d', 'theta', 'gamma', 'omega', 'b1')
    
    
    
    d = pars['d'] 
    theta = pars['theta'] 
    gamma = pars['gamma'] 
    omega = pars['omega'] 
    #a = pars['a'] 
    b1 = pars['b1'] 
 
    
    
    if( d <= 0){
        return(NA)
    }
    
    if (d >= 0.5){
        return(NA)
    }
    
    
    lambdas = lambda(n-1, d, 0, b1)
    #plot.ts(lambdas)
    
    
    
    sigma[1] = exp(0.5*omega)
    z[1] = x[1]/sigma[1]
    g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
    
    
    for(t  in 2:n){
        
        aux = (t-1):(t-1-(n-1)) + n
        
        (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
        x[t]
        (z[t] = x[t]/sigma[t])
        (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
    }
    
    if( is.na(sum(sigma)) ){
        return(NA)
    }
    
    #plot.ts(z)
    #hist(z)
    plot.ts(sigma)
    
    log_func = -n/2*log(2/pi) - 1/2*sum(log(sigma[t]^2) + (x[t]^2)/(sigma[t]^2))
    print(log_func)
    
    if(log_func == Inf | log_func == -Inf | is.na(log_func)){
        return(NA)
    }
    
    return(log_func)
}

ruido_0d1 = function(pars, volat = F){
  
  x =  sp500_diario$log_returns[-1][1:1204]
  
  n = length(x)
  sigma = numeric(n)
  g_z = numeric(2*n - 1)
  z = numeric(n)
  
  
  if(length(pars) != 5){
    stop("Vetor de parâmetros suprido tem comprimento diferente de 6.")
  }
  names(pars) = c('d', 'theta', 'gamma', 'omega', 'b1')
  
  
  
  d = pars['d'] 
  theta = pars['theta'] 
  gamma = pars['gamma'] 
  omega = pars['omega'] 
  #a = pars['a'] 
  b1 = pars['b1'] 
  
  
  
  if( d <= 0){
    return(NA)
  }
  
  if (d >= 0.5){
    return(NA)
  }
  
  
  lambdas = lambda(n-1, d, 0, b1)
  #plot.ts(lambdas)
  
  
  
  sigma[1] = exp(0.5*omega)
  z[1] = x[1]/sigma[1]
  g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
  
  
  for(t  in 2:n){
    
    aux = (t-1):(t-1-(n-1)) + n
    
    (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
    x[t]
    (z[t] = x[t]/sigma[t])
    (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
  }
  
  if (volat){
    return(sigma)
  }
  
  return(z)
}

quasi_llf_fiegarch_0d2 = function(pars = numeric(6)){
    
    
    x =  sp500_diario$log_returns[-1][1:1204]
    
    n = length(x)
    sigma = numeric(n)
    g_z = numeric(2*n - 1)
    z = numeric(n)
    
    
    if(length(pars) != 6){
        stop("Vetor de parâmetros suprido tem comprimento diferente de 6.")
    }
    names(pars) = c('d', 'theta', 'gamma', 'omega', 'b1', 'b2')
    
    
    
    d = pars['d'] 
    theta = pars['theta'] 
    gamma = pars['gamma'] 
    omega = pars['omega'] 
    #a = pars['a'] 
    b1 = pars['b1'] 
    b2 = pars['b2']
    
    b = c(b1,b2)
    
    if( d <= 0){
        return(NA)
    }
    
    if (d >= 0.5){
        return(NA)
    }
    
    
    lambdas = lambda(n-1, d, 0, b)
    #plot.ts(lambdas)
    
    
    
    sigma[1] = exp(0.5*omega)
    z[1] = x[1]/sigma[1]
    g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
    
    
    for(t  in 2:n){
        
        aux = (t-1):(t-1-(n-1)) + n
        
        (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
        x[t]
        (z[t] = x[t]/sigma[t])
        (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
    }
    
    if( is.na(sum(sigma)) ){
        return(NA)
    }
    
    #plot.ts(z)
    #hist(z)
    plot.ts(sigma)
    
    log_func = -n/2*log(2/pi) - 1/2*sum(log(sigma[t]^2) + (x[t]^2)/(sigma[t]^2))
    print(log_func)
    
    if(log_func == Inf | log_func == -Inf | is.na(log_func)){
        return(NA)
    }
    
    return(log_func)
}

ruido_0d2 = function(pars){
  
  x =  sp500_diario$log_returns[-1][1:1204]
  
  n = length(x)
  sigma = numeric(n)
  g_z = numeric(2*n - 1)
  z = numeric(n)
  
  
  if(length(pars) != 6){
    stop("Vetor de parâmetros suprido tem comprimento diferente de 6.")
  }
  names(pars) = c('d', 'theta', 'gamma', 'omega', 'b1', 'b2')
  
  
  
  d = pars['d'] 
  theta = pars['theta'] 
  gamma = pars['gamma'] 
  omega = pars['omega'] 
  #a = pars['a'] 
  b1 = pars['b1'] 
  b2 = pars['b2']
  
  b = c(b1,b2)
  
  if( d <= 0){
    return(NA)
  }
  
  if (d >= 0.5){
    return(NA)
  }
  
  
  lambdas = lambda(n-1, d, 0, b)
  #plot.ts(lambdas)
  
  
  
  sigma[1] = exp(0.5*omega)
  z[1] = x[1]/sigma[1]
  g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
  
  
  for(t  in 2:n){
    
    aux = (t-1):(t-1-(n-1)) + n
    
    (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
    x[t]
    (z[t] = x[t]/sigma[t])
    (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
  }
  
  
}



# 2,d,q -------------------------------------------------------------------

quasi_llf_fiegarch_2d1 = function(pars = numeric(7)){
    
    
    x =  sp500_diario$log_returns[-1][1:1204]
    
    n = length(x)
    sigma = numeric(n)
    g_z = numeric(2*n - 1)
    z = numeric(n)
    
    
    if(length(pars) != 7){
        stop("Vetor de parâmetros suprido tem comprimento diferente de 6.")
    }
    names(pars) = c('d', 'theta', 'gamma', 'omega', 'a1', 'a2', 'b1')
    
    
    
    d = pars['d']
    theta = pars['theta'] 
    gamma = pars['gamma'] 
    omega = pars['omega'] 
    a1 = pars['a1']
    a2 = pars['a2']
    b1 = pars['b1'] 
    a = c(a1,a2)
    b = b1
  
    
    
    if( d <= 0){
        return(NA)
    }
    
    if (d >= 0.5){
        return(NA)
    }
    
    
    lambdas = lambda(n-1, d, a, b)
    plot.ts(lambdas)
    
    
    
    sigma[1] = exp(0.5*omega)
    z[1] = x[1]/sigma[1]
    g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
    
    
    for(t  in 2:n){
        
        aux = (t-1):(t-1-(n-1)) + n
        
        (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
        x[t]
        (z[t] = x[t]/sigma[t])
        (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
    }
    
    if( is.na(sum(sigma)) ){
        return(NA)
    }
    
    #plot.ts(z)
    #hist(z)
    plot.ts(sigma)
    
    log_func = -n/2*log(2/pi) - 1/2*sum(log(sigma[t]^2) + (x[t]^2)/(sigma[t]^2))
    print(log_func)
    
    if(log_func == Inf | log_func == -Inf | is.na(log_func)){
        return(NA)
    }
    
    return(log_func)
}

ruido_2d1 = function(pars){
  
  
  x =  sp500_diario$log_returns[-1][1:1204]
  
  n = length(x)
  sigma = numeric(n)
  g_z = numeric(2*n - 1)
  z = numeric(n)
  
  
  if(length(pars) != 7){
    stop("Vetor de parâmetros suprido tem comprimento diferente de 6.")
  }
  names(pars) = c('d', 'theta', 'gamma', 'omega', 'a1', 'a2', 'b1')
  
  
  
  d = pars['d']
  theta = pars['theta'] 
  gamma = pars['gamma'] 
  omega = pars['omega'] 
  a1 = pars['a1']
  a2 = pars['a2']
  b1 = pars['b1'] 
  a = c(a1,a2)
  b = b1
  
  
  
  if( d <= 0){
    return(NA)
  }
  
  if (d >= 0.5){
    return(NA)
  }
  
  
  lambdas = lambda(n-1, d, a, b)
  plot.ts(lambdas)
  
  
  
  sigma[1] = exp(0.5*omega)
  z[1] = x[1]/sigma[1]
  g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
  
  
  for(t  in 2:n){
    
    aux = (t-1):(t-1-(n-1)) + n
    
    (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
    x[t]
    (z[t] = x[t]/sigma[t])
    (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
  }
  
  
  return(z)
}


quasi_llf_fiegarch_2d2 = function(pars = numeric(8)){
    
    
    x =  sp500_diario$log_returns[-1][1:1204]
    
    n = length(x)
    sigma = numeric(n)
    g_z = numeric(2*n - 1)
    z = numeric(n)
    
    
    if(length(pars) != 8){
        stop("Vetor de parâmetros suprido tem comprimento diferente de 8.")
    }
    names(pars) = c('d', 'theta', 'gamma', 'omega', 'a1', 'a2', 'b1', 'b2')
    
    
    
    d = pars['d']
    theta = pars['theta'] 
    gamma = pars['gamma'] 
    omega = pars['omega'] 
    a1 = pars['a1']
    a2 = pars['a2']
    b1 = pars['b1']
    b2 = pars['b2']
    a = c(b1,b2)
    b = c(b1,b2)
    
    
    
    if( d <= 0){
        return(NA)
    }
    
    if (d >= 0.5){
        return(NA)
    }
    
    
    lambdas = lambda(n-1, d, a, b)
    plot.ts(lambdas)
    
    
    
    sigma[1] = exp(0.5*omega)
    z[1] = x[1]/sigma[1]
    g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
    
    
    for(t  in 2:n){
        
        aux = (t-1):(t-1-(n-1)) + n
        
        (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
        x[t]
        (z[t] = x[t]/sigma[t])
        (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
    }
    
    if( is.na(sum(sigma)) ){
        return(NA)
    }
    
    #plot.ts(z)
    #hist(z)
    plot.ts(sigma)
    
    log_func = -n/2*log(2/pi) - 1/2*sum(log(sigma[t]^2) + (x[t]^2)/(sigma[t]^2))
    print(log_func)
    
    if(log_func == Inf | log_func == -Inf | is.na(log_func)){
        return(NA)
    }
    
    return(log_func)
}


ruido_2d2 = function(pars){
  
  x =  sp500_diario$log_returns[-1][1:1204]
  
  n = length(x)
  sigma = numeric(n)
  g_z = numeric(2*n - 1)
  z = numeric(n)
  
  
  if(length(pars) != 8){
    stop("Vetor de parâmetros suprido tem comprimento diferente de 8.")
  }
  names(pars) = c('d', 'theta', 'gamma', 'omega', 'a1', 'a2', 'b1', 'b2')
  
  
  
  d = pars['d']
  theta = pars['theta'] 
  gamma = pars['gamma'] 
  omega = pars['omega'] 
  a1 = pars['a1']
  a2 = pars['a2']
  b1 = pars['b1']
  b2 = pars['b2']
  a = c(b1,b2)
  b = c(b1,b2)
  
  
  
  if( d <= 0){
    return(NA)
  }
  
  if (d >= 0.5){
    return(NA)
  }
  
  
  lambdas = lambda(n-1, d, a, b)
  plot.ts(lambdas)
  
  
  
  sigma[1] = exp(0.5*omega)
  z[1] = x[1]/sigma[1]
  g_z[1+n] = theta*z[1] + gamma*(abs(z[1]) - sqrt(2/pi)) 
  
  
  for(t  in 2:n){
    
    aux = (t-1):(t-1-(n-1)) + n
    
    (sigma[t] = exp(omega/2 + 1/2*sum(lambdas*g_z[aux])))
    x[t]
    (z[t] = x[t]/sigma[t])
    (g_z[t+n] = theta*z[t] + gamma*(abs(z[t]) - sqrt(2/pi))) 
  }
  
  
  return(z)
}


# estima CTS ----

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