b_j = function(d, j){
    #expressão (5)
    out = gamma(j+d)/(gamma(d)*gamma(j+1))
    return(out)
}


c_j = function(d, j, c_ar, phi, theta){
    # expressão (8)
    # Função que calcula c_j de acordo com suas três parcelas.
    # p1 = b_j(d)
    # p2 = \sum_{i=1}^{p}(\phi_i c_{j-i}(d))
    # p3 = \sum_{k=1}^{q}(\theta_i b_{j-k)(d)})
    p1  = b_j(d, j)
    p2 = sum(phi*c_ar)
    p3 = 0
    for (k in seq_along(theta)){
        p3 = p3 + theta[k]*b_j(d, j-k)
    }
    
    out = p1 + p2 + p3
    
    return(out)
    
    
}


arfima_proc = function(phi, 
                    d, 
                    theta, 
                    n = 500, 
                    iter = 50, 
                    innov = rnorm(n+max(length(phi),length(theta))+iter),
                    init = NULL
                    ){
    # ARGUMENTOS
    # phi: vetor autorregressivo de tamanho p
    # theta: vetor de média móvel de tamanho q
    # n : tamanho final da série temporal a ser gerada
    # iter: número de iterações para a expressão (7)
    # innov: por padrão, vetor de v.a.'s gaussianas padrão de tamanho (n + max(p,q) + iter)
    # vetor de valores iniciais para c_j
    
    
    ### Teste lógico para valor inicial
    # Se não for dado um valor inicial pelo usuário, gera um vetor de tamanho p com c_j = 0.5 para j = 1, 2, ... p
    if (is.null(init)){
        c_ar = rep(0.5, length(phi))
    } else if (length(init) == 1){
        c_ar = rep(init, length(phi))
    } else if (length(init) != length(phi)){
        stop("Valores iniciais incompatíveis")
    } 
    
    p = length(phi) # tamanho do vetor AR
    if(phi == 0) p = 0               # 
    q = length(theta)                # Linhas para que se phi_1 = 0 e/ou theta_1 = 0, p e ou q vão ser 0
    if(theta == 0) q = 0             #
    
    y = numeric(n) # Gera um vetor y de tamanho n que será usado para retornar a série gerada                  
    C = c(c_ar, numeric(iter)) # Vetor de tamanho 'iter' com os valores de c_j a serem utilizados na expressão (7)
    
    
    for(i in 1:n){ #gera-se uma observação do processo para cada tempo de 1 até n 
        
        for (j in (p+1):(iter+p)){ #Esse vetor começa em p+1, pois precisamos dos p valores anteriores caso p != 0
            
            C[j] = c_j(d, j, C[(j-p):(j-1)], phi, theta)
            
        }
            y[i] = sum(C*innov[i:(iter+i)]) # Expressão (7)
    }
    
    
    return(y)
}
    

# teste 

set.seed(123)
par(mfrow = c(1,2))
innov = rnorm(550)
a1 = arfima_proc(0, 0.3, 0, init = 1, innov = innov)
plot.ts(a1)
acf(a1)


set.seed(123)
a2 = arfima_proc(0, 0.3, 0, init = 1, innov = innov)
plot.ts(a2)
acf(a2)


a3 = arfima_proc(0, 0.47, 0)
plot.ts(a3)
acf(a3)


par(mfrow = c(2,5))

set.seed(123)
innov = rnorm(550)
d = 0.45

par(mfrow = c(1,1))
for(i in 1:10) {
    c_i = 0.1*i
    a = arfima_proc(0, d, 0, innov = innov, init = c_i)
    #plot.ts(a, main = paste("Valor inicial = ", c_i, ", d =", d))
    acf(a, main = paste("Valor inicial = ", c_i, ", d =", d), lag.max = 100 )
    #plot.ts(a)
    #title(paste("Valor inicial = ", c_i))
}


# teste stable

library(stabledist)
