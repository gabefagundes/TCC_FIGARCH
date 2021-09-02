b_j = function(d, j){
    out = gamma(j+d)/(gamma(d)*gamma(j+1))
    return(out)
}


c_j = function(d, j, c_ar, phi, theta){
    p1  = b_j(d, j)
    p2 = sum(phi*c_ar)
    p3 = 0
    for (k in seq_along(theta)){
        p3 = p3 + theta[k]*b_j(d, j-k)
    }
    
    out = p1 + p2 + p3
    
    return(out)
    
    
}

b_j(0.3, 10)
c_j(0.3, 10, 0.1, 0, 0)



arfima_proc = function(phi, 
                    d, 
                    theta, 
                    n = 500, 
                    iter = 50, 
                    innov = rnorm(n+max(length(phi),length(theta))+iter),
                    init = NULL
                    ){
    if (is.null(init)){
        c_ar = rep(0.5, length(phi))
    } else if (length(init) == 1){
        c_ar = rep(init, length(phi))
    } else if (length(init) != length(phi)){
        stop("Valores iniciais incompatíveis")
    } 
    
    p = length(phi)
    q = length(theta)
    y = numeric(n + max(p,q))
    C = c(c_ar, numeric(iter))
    
    
    for(i in 1:n){
        
        for (j in (p+1):(iter+q)){
            
            C[j] = c_j(d, j, C[(j-p):(j-1)], phi, theta)
            
        }
            y[i] = sum(C*innov[i:(iter+i)])
    }
    
    return(y)
}
    

# teste 
par(mfrow = c(1,2))
a1 = arfima_proc(0, 0.3, 0)
plot.ts(a1)
acf(a1)



a2 = arfima_proc(0, 0.4, 0)
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
