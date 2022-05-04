# Script que gera variaveis aleatorias temperadas classicas  # 
# Desenvolvimento começou em 16/03/2022 # 


rCTS = function(n, alpha, ell_plus, ell_minus, c, M = 10000){
    
    out = numeric(n)
    
    for(i in 1:n){
        v = sample(c(-ell_minus, ell_plus), size = M, replace = T)
        
        
        u = runif(M)
        
        e = rexp(M)
        
        e_prime = rexp(M)
        
        gama = cumsum(e_prime)
        
        
        S = 0
        
        aux2 = - gamma(1-alpha)*c*(ell_plus^(alpha-1) - ell_minus^(alpha-1))
        
        
        for (j in 1:M){
            
            aux1 = min((alpha*gama[j]/(2*c))^(-1/alpha), e[j]*u[j]^(1/alpha)*abs(v[j])^(-1))*sign(v[j])
            #print(aux1)
            
            
            S = S + (aux1 + aux2)
            #print(S)
        }
        
        out[i] = S
    }
    return(out)
}



# testes

# alpha = 1.5
# ell_plus = 1
# ell_minus = 1
# c = (gamma(2-alpha)*(ell_minus^(alpha-2) + ell_plus^(alpha-2)))^(-1)
# 
# 
# 
# y = rCTS(1, alpha, ell_plus, ell_minus, c, M = 10000)
# 
# 
# var(y)
# hist(y)



# Geração de banco auxiliar -----------------------------------------------


#alphas = seq(0.05, 1.95, by = 0.05)
#alphas = 0.7


#yy = matrix(numeric(1000*length(alphas)), ncol = length(alphas))
xx = numeric(1000)

# for (j in 1:19){
#     print(alphas[j])
#     for(i in 1:1000){
#         print(i)
#         yy[i,j] = mean(abs(rCTS(50000, alpha = alphas[j],ell_plus = 1, ell_minus = 1, c = (gamma(2-alphas[j])*(ell_minus^(alphas[j]-2) + ell_plus^(alphas[j]-2)))^(-1))))
#     }
# }

# alpha = 0.7
# ell =  1
# 
# for(i in 1:1000){
#     print(i)
#     t0 = Sys.time()
#     xx[i] = mean(abs(rCTS(5000, alpha = alpha ,ell_plus = ell, ell_minus = ell , c = (gamma(2-alpha)*(ell^(alpha-2) + ell_plus^(alpha-2)))^(-1), M = 5000)))
#     t1 = Sys.time()
#     print(t1 - t0)
# }
# 
# write_csv2(as.data.frame(xx), 'FI(E)GARCH/dados/TS_absmean_alpha0.7.csv')
# yy