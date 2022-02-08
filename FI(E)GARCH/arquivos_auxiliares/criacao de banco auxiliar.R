library(tidyverse)


alphas = seq(1.1, 2, by = 0.05)

y = matrix(numeric(1000*19), ncol = 19)


for (j in 1:19){
    print(alphas[j])
    for(i in 1:1000){
        print(i)
        y[i,j] = mean(abs(rstable(100000, alpha = alphas[j], 0, 1, 0)))
    }
}
    




colunas = paste0('alpha_', alphas)

df = as.data.frame(y) 


colnames(df) = colunas


write.csv2(df, 'FI(E)GARCH/arquivos_auxiliares/simulacoes_alpha_estavel_completo.csv')


df2 = df %>%
    summarise_all(mean) %>% 
    pivot_longer(cols = alpha_1.1:alpha_2,
                 names_to = 'alpha',
                 names_prefix = 'alpha_',
                 values_to = 'media') %>% 
    mutate(alpha = as.numeric(alpha))


write.csv2(df2, 'FI(E)GARCH/arquivos_auxiliares/simulacoes_abs_alpha_estavel_medias.csv')
