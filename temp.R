
# Bibliotecas -------------------------------------------------------------

{
    library(stabledist)
    library(tidyverse)
}




# funções -----------------------------------------------------------------



# debug -------------------------------------------------------------------

mu = 2
alpha0 = 3
alpha = c(0.5, 0.3, 0.4)
y0 = c(3, 2, 1)
n = 100


# Teste 

proc = arch_proc(mu, alpha0, alpha, y0, n)


plot.ts(proc$y)
plot.ts(proc$sigma2)
