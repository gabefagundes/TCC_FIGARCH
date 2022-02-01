{
    library(stabledist)
    library(StableEstim)
}


##################################
#### Estimação alpha-Estáveis ####
##################################

# Primeira proposta de resolver exemplo 4.1 do Documento de Estimação Clássica e Bayesiana

#Gerar uma amostra i.i.d, de uma S(alpha, beta, gamma, delta; 0) (parametrização 0) com
# alpha = 1.7, beta =0, gamma = 1, delta = 0 (alpha estavel simétrica)

set.seed(1234)
x_1.7 = rstable(1000, alpha = 1.7, beta = 0, pm = 0)
teste = rstable(1000, alpha = 1.7, beta = 0, pm = 1)

# Teste para afirmação na página 7 de Nolan(2020), falando que a pm = 1 é igual a pm = 0 para 
# alpha estáveis simétricas

print('parametrização 0 vs parametrização 1')
plot(ecdf(x_1.7), xlim = c(-5, 5))
lines(ecdf(teste), col = 'red', lty = 2)

print('teste KS')
ks.test(x_1.7, teste) # Afirmação se sustenta com KS-test

# Formas diferentes de estimar:

# McCulloch(1986):

estim_mcculloch = McCullochParametersEstim(x_1.7)

# MLE:

estim_MLE = Estim('ML', 
                  data = x_1.7, 
                  ComputeCov = T)

# parâmetros estimados: 

estim_MLE@par

# matriz de variancia e covariancia dos parametros

estim_MLE@vcov

# intervalo de confiança para os parametros estimados

estim_MLE@confint




# Exemplo 4.2; alpha < 0  -------------------------------------------------------------

x_0.7 = rstable(1000, alpha = 0.7, beta = 0, pm = 0)

# McCulloch(1986):

estim_mcculloch_2 = McCullochParametersEstim(x_0.7)

# MLE:

estim_MLE_2 = Estim('ML', 
                  data = x_0.7, 
                  ComputeCov = T)




# Estimação de Monte Carlo ------------------------------------------------
B = 500
alpha = 1.7
df_estim = data.frame(
    alpha_mcc = c(),
    alpha_mle = c()
)



for(i in 1:B){
    print(i)
    x = rstable(1000, alpha, beta = 0, pm = 0)
    
    df_estim[i, 1] = c(McCullochParametersEstim(x)[1])
}




# Medida Espectral --------------------------------------------------------
library(stabledist)
library(pracma)
library(nnls)
library(StableEstim)


metodo = 'ML'
alpha = 1.4
m = 3 
pontos = 12
n = 1000

gamma = rep(1/3, 3)
S = numeric(2)

grid = matrix(numeric(2*m), ncol = 2)

for(j in 1:m){
    grid[j, ] = c(cos(2*pi*(j-1)/m), sin(2*pi*(j-1)/m))
    S = S + gamma[j] * grid[j,]
}


N <-pontos
s<-matrix(numeric(N*2), ncol = N)
for (j in 1:(N/2)) {
    s[,j]<- c(cos(2*pi*(j-1)/N),sin(2*pi*(j-1)/N))
    s[,j+N/2]<- - Conj(s[,j])}


stv<-matrix(rstable(n*m,alpha,1,1,0,pm=1),ncol = m)
X = matrix(numeric(n*2),ncol = 2)
for (l in 1:n){
    for(j in 1:m){X[l,]<- X[l,] + gamma[j]^(1/alpha)*stv[l,j]*grid[j,]
    }
    X[l,] <- X[l,] - S*tan(pi*alpha/2)
}



Y = X%*%s


if(metodo == 'mc'){
    KLS = apply(Y, 2, McCullochParametersEstim)
} else if (metodo == 'ML') {
    for (i in 1:pontos){
        print(i)
        KLS[,i] = Estim("ML", data = Y[,i])@par
    }
}



# Pooled estimate of alpha
aest <- numeric(1)
for(j in 1:N){
    aest<-aest+as.numeric(KLS[1, j])/N
}


# Construction of the matrix Phi, given in page 1116 in Nolan, Panorska and McCulloch (2001)
phi <- function(a,u){
    if(a==1){
        output<-abs(u)*(1+2i/pi*sign(u)*log(abs(u)))
    } else{
        output<-abs(u)^a*(1-1i*sign(u)*tan(pi*a/2))}
        return(output)
}

Phi <- matrix(complex(N*N/4),ncol=N/2)

for (i in 1:(N/2)) {
    for(j in 1:(N/2)){
        Phi[i,j]<-phi(aest,dot(s[,i],s[,j]))
        }
    }


# Estimated exponent functions
I<-numeric(N)
for (j in 1:N) {
    beta<- as.numeric(KLS[2,j])
    sigma<-as.numeric(KLS[3,j])
    delta<-as.numeric(KLS[4,j])
    if(aest==1){
        I[j]<-sigma*(1-1i*delta)
    } else{
            I[j]<-sigma^aest*(1-1i*beta*tan(pi*aest/2))
            }
    }
# Construction of vector c, given in eq. (3) in Nolan, Panorska and McCulloch (2001)
c <- numeric(N)
for (j in 1:(N/2)) {
    c[j] <- Re(I[j])
    c[j+N/2]<-Im(I[j])
    }
# Construction of matrix A, given in eq. (3) in Nolan, Panorska and McCulloch (2001)
A <- matrix(numeric(N*N), ncol=N)
for (j in 1:(N/2)) {
    for(l in 1:(N/2)){
        A[l,j] <- Re(Phi[l,j])
        A[l+N/2,j] <- Im(Phi[l,j])
        A[l,j+N/2] <- Re(Phi[l,j])
        A[l+N/2,j+N/2] <- -Im(Phi[l,j])
        } 
    }
KLSm <- numeric(N)
for (j in 1:N) {
    KLSm[j]<- as.numeric(KLS[4,j])
    }
    

# Mean of location estimates
delta = mean(KLSm)
# List of the Spectral measures' estimated weights
Phi_ = nnls(A,c)
Phi_$x
# Estimated value of alpha
alpha_hat = aest




# Estimação Final da medida espectral

apply(Phi_$x * t(s), 2, sum)

