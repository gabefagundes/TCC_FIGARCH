
# bibliotecas -------------------------------------------------------------

{
    library(stabledist)
    library(tidyverse)
    source("FI(E)GARCH/stable_fiegarch.R")
}



# funcoes_auxiliares ------------------------------------------------------


f_loglikelihood = function(x, pars = numeric(7)){
    # Função de log-verossimilhança de um processo FIEGARCH com inovações alpha-estáveis
    # @pars : vetor de parâmetro com 7 entradas a serem otimizadas: d, alpha, theta, gamma, omega, a1, b1
    # @x : amostra observada
    #browser()
    {
        x = teste$serie[1:2000]*1
        #x = rstable(1000, alpha = 1.5, beta = 0)
        #x = rnorm(1000, sd = 0.1)
        sd_x = sd(x)
        
        if(length(pars) != 7){
            stop("Vetor de parâmetros suprido tem comprimento diferente de 7.")
        }
        names(pars) = c('d', 'alpha', 'theta', 'gamma', 'omega', 'a1', 'b1')
        
        # parâmetros iniciais
        
        N = length(x)
        d = pars['d'] = 0.43
        alpha = pars['alpha'] = 1.8
        theta = pars['theta'] = 0.11
        gamma = pars['gamma'] = -0.73
        omega = pars['omega'] = -3
        a1 = pars['a1'] 
        b1 = pars['b1'] = 0.54
    }
    
    # Condições para valores amostrais abaixo de zero
    {
        m = N
        
        sigma = numeric(2*N+1)
        x = c(numeric(N + 1), x)
        Z = numeric((2*N+1))
        gZ = numeric(2*N+1)
        lambdas = lambda2(m-1, d, a1, b1)
    }
    
    for(t in (m+2):(2*m+1)){
        
        print(t-m-1)
        
        aux = (t-1):(t-m)
        
        sigma[t] = sqrt(exp(omega + sum(lambdas*gZ[aux])))
        Z[t] = x[t]/sigma[t]
        #(gZ[t] = g_func_stable(Z = Z[t], theta = theta, gamma = gamma, alpha_stable = alpha, method = 'ana'))
        #gZ[t] = theta*Z[t] + gamma*(abs(Z[t]) - 2*gamma(1 - 1/alpha)/pi)
        gZ[t] = theta*Z[t] + gamma*(abs(Z[t]) - sqrt(2/pi)) 
        
        
        {
            print(paste('x_t = ', round(x[t],2)))
            print(paste('sigma_t = ', round(sigma[t],2)))
            print(paste('Z[t] = ', round(Z[t],2)))
            print(paste('g(Z[t]) = ', round(gZ[t],2)))
            #Sys.sleep(3)
        }
        
        #readline()
        
    }
    
    {
        index = (m + 2):(2*m+1)
        
        x = x[index]
        sigma = sigma[index]
        Z = Z[index]
        gZ = gZ[index]
        
        par(mfrow = c(2,2))
        plot.ts(x)
        plot.ts(sigma)
        plot.ts(Z)
        plot.ts(gZ)
    }
    
    k = 1:m
    t = 1:2
    denom_sigma = 1/(pi*alpha*sigma)
    quociente_gamas = gamma(alpha/k)/gamma(k)*sin(k*pi/2) 
    parcelas = log(denom_sigma*sum(quociente_gamas*(-Z)^(k-1)))
    out = sum(parcelas)
}



debug(f_loglikelihood)
undebug(f_loglikelihood)

f_loglikelihood(x = 1)














# testes optim ------------------------------------------------------------


