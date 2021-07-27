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

d = 0.3
#======
phi = c(0.5)
c_ar = c(1)
theta = c(0.5)


arfima_2 = function(phi, d, theta, n = 500, iter = 50, innov = rnorm(n+max(p,q)+iter), c_ar = rep(1, length(phi))){
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
}
    