
# Bibliotecas -------------------------------------------------------------

{
    library(stabledist)
    library(tidyverse)
}



# Aux ---------------------------------------------------------------------

comb = function(d,k){
    gamma(1+d)/(gamma(1+k)*gamma(1+d-k))
}


f_aux1 = function(d,p){
    
    if(p == 0){
        return(1)
    }
    
    out = ((p-1) - d)*f_aux1(d, p - 1)
    
    return(out)
}


delta_d = function(y, t, d, iter = 5) {
    
    out = 0
    
    ####
    d = 0.3
    t = 1000
    iter = 150
    ####
    
    for (i in 0:(iter-1)){
        out = out + f_aux1(d,i)*y[t-i]/factorial(i)
    }
    
    return(out)
}
  


# teste -------------------------------------------------------------------


delta_d(y, 5, 0.3)
    
