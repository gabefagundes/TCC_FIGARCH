
# Tempered Stable Testes --------------------------------------------------

{
    source('FI(E)GARCH/Tempered Stable/rCTS.R')
    library(tidyverse)
}




# Valor para padronizar alpha estável 
alpha = 1.9
ell = 1
c = ell^(alpha-2)/(2*gamma(2-alpha))


y = rCTS(10000, alpha = alpha, c = c, ell = ell, mu = 0)


var_real = ell^(2-alpha)*c*gamma(2 - alpha)*2


#hist(y)


mean(y)
var(y)




# Teste para valor absoluto da tempered stable ----------------------------

N = 50000
alphas_aux = seq(0.1, 1.9, 0.1)
ells_aux = 1
EX_hat = matrix(numeric(length(alphas_aux)*length(ells_aux)), ncol = length(ells_aux))
rownames(EX_hat) = paste0('alpha_',alphas_aux)
colnames(EX_hat) = paste0('ell_',ells_aux)

abs_values = tibble(
    alphas = rep(alphas_aux, each = length(ells_aux)),
    ells = rep(ells_aux, length(alphas_aux)),
    E_X = 0
)

t0 = Sys.time()


for(i in seq_along(alphas_aux)){
    print(paste0('i=',i))
    for(j in seq_along(ells_aux)){
        print(paste0('j=', j))
        alpha = alphas_aux[i]
        ell = ells_aux[j]
        c = ell^(alpha-2)/(2*gamma(2-alpha))
        amostra = rCTS(N, 
                  alpha = alphas_aux[i],
                  ell = ells_aux[j], 
                  c = c,
                  mu = 0
        ) 
        EX_hat[i,j] = mean(abs(amostra), na.rm = T)
    }
}

t1 = Sys.time()

print(t1 - t0)

write.csv2(EX_hat, 'FI(E)GARCH/dados/medias_abs_tempered_stable2.csv')

#write.csv2(EX_hat3, 'FI(E)GARCH/dados/medias_abs_tempered_stable_pivot.csv')
#EX_hat3 = read_csv2('FI(E)GARCH/dados/medias_abs_tempered_stable_pivot.csv')


ggplot(EX_hat3) +
    geom_line(aes(as.numeric(ell), E_X, color = as.factor(alpha), group = as.factor(alpha)))







# Explore ell -------------------------------------------------------------

cols = RColorBrewer::brewer.pal(9, 'Spectral')

x = seq(-10, 10, 0.1)
curve(dCTS(x, alpha = 1.5, ell = 1, c = ell^(alpha-2)/gamma(2-alpha), mu = 0), from = -10, to = 10 )
lambdas_graf = 2:10

for (i in 1:9){
    lines(x, dCTS(x, alpha = 1.5, ell = lambdas_graf[i], c = lambdas_graf[i]^(alpha-2)/gamma(2-alpha), mu = 0), col = cols[i])
}

