
# Bibliotecas -------------------------------------------------------------

{
    library(stabledist)
    library(tidyverse)
}




# Teste 

proc = arch_proc(alpha = c(0.3, 1.5), y0 = c(0.2, 1.5), alpha0 = 0.5)


plot.ts(proc$y)
plot.ts(proc$sigma2)
