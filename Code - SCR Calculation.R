

folder = getwd()
paste(folder, sep = '/', "chun")
######################################################
#PARTIE I) Installing packages and parameters
#####################################################
#install.packages("orthopolynom")
library(orthopolynom)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("gridExtra")
library(gridExtra)

#############################################################

maturite = 2 #la maturité du contrat
sigma_hi = 0.15 #la volatilité historique 
sigma_rn = 0.25 #la volatilité dqns l'environement risque neutre

mu_es = 0.02 # l'esperence de rendrement

S0 = 1 #valeur initial de S 

r = 0.01  #le taux sans risque
r_g = 0.005 #le taux technique garantie
VM_0 = 110   #la valeur initiale du fond
x = 0.3 #la proportion investit dans l'actif risqué
PM_0 = 100
R = PM_0/VM_0


N0 = 1000 #Nombre de simulation pour NAV0

#function simuler St with t = 1
S_1 = function(Nsim = N0, mu = mu_es){
  return( S0 * exp(( mu - ( sigma_hi^2 ) / 2) * 1 + sigma_hi * rnorm(Nsim,0,1)))
}

#function simuler St with t = 2

S_2 = function(S_init,Nsim = N0, mu = r ){
  return(S_init * exp(( mu - ( sigma_rn^2 ) / 2) * 1 + sigma_rn * rnorm(Nsim,0,1)))
}


#simuler S2
S1 = S_1()
S2 = S_2(S1)


VM = function(St , VM0 = VM_0, taux = r, t = maturite){
  return(VM0 * ((1 - x) * exp(taux * t) + x * St))}



#la fonction permettant de calculer la valeur de la partie garantie BEt
BE = function(VMt){
  return(R*VMt)
}

#la fonction permettant de calculer la valeur de l'option en t

#valeur initiale 
K0 = x * PM_0 
Option = function(St , rg = r_g){
  K = R * VM_0 * (exp( rg * maturite) - (1 - x)*exp(r * maturite))
  (K - K0*St)*( K > K0 * St)    
}


actu = function(taux=r, t1 = 0, t2 = maturite){
  return(exp(-taux*(t2 - t1)))
} 

################################################################
#PARTIE 2) CALCULATION
################################################################

#P2.I) CACUL DE NAV0

VM_T = VM(S2)
BE_T = BE(VM_T)
O_T = Option(S2)
NAV0 = actu() * mean(VM_T - BE_T - O_T)
NAV0



######################################################################

NAV1 = function(S_debut, taux_nav = r, rg_nav = r_g, Nsim_nav = M ){
  S2 = S_2(S_init = S_debut, Nsim = Nsim_nav)
  VM_nav = VM(S2, taux = taux_nav)
  BE_nav = BE(VM_nav)
  O_nav = Option(S2, rg = rg_nav)
  actu_nav = actu(taux = taux_nav, t1 = 1 , t2 = maturite)
  return( actu_nav * mean (VM_nav - BE_nav - O_nav))
}

#P2.II) Methode de calcul NAV1

N = 25000 #nombre de simulation primaire
M = 1000 #nombre de simulation secondaire


#1.Methode SdS

#1.1.Simulate Simuler N trajectoires du sous-jacent à la date t = 1

set.seed(2)
S1 = S_1(Nsim = N)

#1.2.Pour chaque simulation k, estimer par Monte-Carlo NAVk1

NAV1s <- sapply(S1, function(x) NAV1(S_debut = x, Nsim_nav = M))

SCR_SdS  = NAV0 - quantile(NAV1s, p = 0.005)/(1 + r) 
SCR_SdS

#2.Methode LMSC


N_lmsc = 10000
M_lmsc = 1000

#2.1  Simuler N trajectoires du sous-jacent à la date t = 1

set.seed(2)
S1 = S_1(Nsim = N_lmsc)

#2.2 Pour chaque simulation k, estimer par Monte-Carlo NAVk

NAV1s <- sapply(S1, function(x) NAV1(S_debut = x, Nsim = M_lmsc))



#2.3 Utiliser une base orthonormée et estimer par moindres carrées la régression linéaire

d = 7

#2.3.1 modele canonique

modele_canonique <- lm(NAV1s ~ poly (S1 ,d,raw=F))
summary(modele_canonique)

NAV1_canonique = modele_canonique$fitted.values
SCR_poly_canonique = NAV0 - quantile(NAV1_canonique,p=0.005)/(1+r)
SCR_poly_canonique


#scale X pour les polynomials qui requissent les variables dans la plage -1 et 1

x_scale = scaleX(S1, u=-1, v=1)

#2.3.2 modele Jaccobi
jacobi = jacobi.g.polynomials(n = d, p = 3, q = 2)

x_jacobi = as.matrix(as.data.frame(polynomial.values(jacobi, x_scale))[, 2: (d+1)])
modele_jacobi = lm(NAV1s ~ x_jacobi)
NAV1_jacobi = modele_jacobi$fitted.values
SCR_poly_jacobi = NAV0 - quantile(NAV1_jacobi, p = 0.005)/(1+r)
SCR_poly_jacobi

#2.3.3 modele Chebysev
chebyshev = chebyshev.t.polynomials(n = d)
x_chebyshev = as.matrix(as.data.frame(polynomial.values(chebyshev, S1)))

modele_chebyshev = lm(NAV1s ~ x_chebyshev)
NAV1_chebyshev = modele_chebyshev$fitted.values
SCR_poly_chebyshev = NAV0 - quantile(NAV1_chebyshev, p = 0.005)/(1+r)
SCR_poly_chebyshev


#2.3.4 Gegenbauer polynomials

gegenbauer = gegenbauer.polynomials(n = d, alpha = 3)

x_gegenbauer = as.matrix(as.data.frame(polynomial.values(gegenbauer, S1 )))
modele_gegenbauer  = lm(NAV1s ~ x_gegenbauer)
NAV1_gegenbauer = modele_gegenbauer$fitted.values
SCR_poly_gegenbauer = NAV0 - quantile(NAV1_gegenbauer, p = 0.005)/(1+r)
SCR_poly_gegenbauer


#2.3.5 Hermite polynomials

hermite = hermite.h.polynomials(n = d)
x_hermite = as.matrix(as.data.frame(polynomial.values(hermite, S1)))
modele_hermite = lm(NAV1s ~ x_hermite)
NAV1_hermite = modele_hermite$fitted.values
SCR_poly_hermite = NAV0 - quantile(NAV1_hermite, p = 0.005)/(1+r)
SCR_poly_hermite


#2.3.6 Laguerre polynomials

laguerre = laguerre.polynomials(n = d)
x_laguerre = as.matrix(as.data.frame(polynomial.values(laguerre, S1)))

modele_laguerre = lm(NAV1s ~ x_laguerre)
NAV1_laguerre = modele_laguerre$fitted.values
SCR_poly_laguerre = NAV0 - quantile(NAV1_laguerre, p = 0.005)/(1+r)
SCR_poly_laguerre


#III) Modele tuning

#1. Tuning SdS

#1.1 Nombre de simulation primaire

SCR_sensi_N = c()

M = 1000 #fix M pour chercher N optimal
N_sen <- seq(20000, 30000, 5)

for (i in N_sen ){
  S1 = S_1(Nsim = i)
  NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
  SCR = NAV0 - quantile(NAV1s, p = 0.005)
  SCR_sensi_N = c(SCR_sensi_N, SCR)
}

png(file = "Sensibilité en fonction du nombre de simulations primaires N.png")
plot(SCR_sensi_N~N_sen, type = 'o', pch=20, col= "red", xlab="N", ylab="SCR", main="Sensibilité en fonction du nombre de simulations primaires N")
dev.off()

write.csv(SCR_sensi_N, "D:/M2/Technique/SCR_sensi_N", row.names=FALSE)


#lower bound N
SCR_sensi_N = c()
for (i in seq(1000,20000, 10)){
  S1 = S_1(Nsim = i)
  NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
  SCR = NAV0 - quantile(NAV1s, p = 0.005)
  SCR_sensi_N = c(SCR_sensi_N, SCR)
}

plot(SCR_sensi_N~seq(1000,20000, 10), type = 'o', pch=20, col= "red", xlab="N", ylab="SCR", main="Sensibilité en fonction du nombre de simulations primaires N")
write.csv(SCR_sensi_N, paste(folder, sep = '/', "SCR_sensi_N"), row.names=FALSE)




plot(sensi_N~N_sen, type = 'o', pch=20, col= "red", xlab="N", ylab="SCR", main="Sensibilité en fonction du nombre de simulations primaires N")

#--------------------------------------#

#1.2 Nombre de simulation secondaire

SCR_sensi_M = c()

N = 25000 #fix N pour chercher M optimal
M_sen <- seq(1, 2000, 1)

for (i in M_sen ){
  S1 = S_t(Nsim = N)
  NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = i))
  SCR = NAV0 - quantile(NAV1s, p = 0.005)
  SCR_sensi_M = c(SCR_sensi_M, SCR)
}


plot(SCR_sensi_M ~ M_sen, type = 'o', pch=20, col= "blue", xlab="M", ylab="SCR", main="Tuning en fonction du nombre de simulations secondaires M")

write.csv(SCR_sensi_M, paste(folder, sep = '/', "SCR_sensi_M"), row.names=FALSE)

plot(SCR_sensi_M ~ M_sen, type = 'o', pch=20, col= "blue", xlab="M", ylab="SCR", main="Tuning en fonction du nombre de simulations secondaires M")


#Choisir N = 25000 avec M = 250


#-----------------------------------------------------------#
# tuning LMSC

lsmc_calcul = function(d){

  #canonique
  modele_canonique <- lm(NAV1s ~ poly (S1 ,d,raw=F))
  NAV1_canonique = modele_canonique$fitted.values
  SCR_canonique = NAV0 - quantile(NAV1_canonique,p=0.005)/(1+r)
  
  #jacobi
  x_scale = scaleX(S1, u=-1, v=1)
  jacobi = jacobi.g.polynomials(n = d, p = 3, q = 2)
  x_jacobi = as.matrix(as.data.frame(polynomial.values(jacobi, x_scale))[, 2: (d+1)])
  modele_jacobi = lm(NAV1s ~ x_jacobi)
  NAV1_jacobi = modele_jacobi$fitted.values
  SCR_jacobi = NAV0 - quantile(NAV1_jacobi,p=0.005)/(1+r)
  
  #Hermite polynomials
  
  hermite = hermite.h.polynomials(n = d)
  x_hermite = as.matrix(as.data.frame(polynomial.values(hermite, S1)))
  modele_hermite = lm(NAV1s ~ x_hermite)
  NAV1_hermite = modele_hermite$fitted.values
  SCR_hermite = NAV0 - quantile(NAV1_hermite,p=0.005)/(1+r)
  
  return(c(SCR_canonique, SCR_jacobi, SCR_hermite))
}

#sensi N
M = 1000
d = 7

df = data.frame()

for ( i in 1000:10000){
  S1 = S_1(Nsim = i)
  NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
  #LMSC
  df = rbind(df, lsmc_calcul(d))

}

colnames(df) = c('Canonique', 'Jacobi', 'Hermite')
write.csv(df, paste(folder, sep = '/', "df"), row.names=FALSE)


par(mfrow=c(3,1))
plot(df$Canonique ~ seq(1000,10000,1), xlab="N", ylab="SCR", main="Tuning nombre de simulation primaire N - Cano")
plot(df$Jacobi ~ seq(1000,10000,1), xlab="N", ylab="SCR", main="Tuning nombre de simulation primaire N - Jaco")
plot(df$Hermit ~ seq(1000,10000,1), xlab="N", ylab="SCR", main="Tuning nombre de simulation primaire N - Hermit")


  

#sensi M

N = 4000
d = 7
df_M = data.frame()

for ( i in 1:1000){
  S1 = S_1(Nsim = N)
  NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = i))
  #LMSC
  df_M = rbind(df_M, lsmc_calcul(d))
}

write.csv(df_M, paste(folder, sep = '/', "df_M"), row.names=FALSE)

df_M = read.csv("D:/M2/Technique/df_M")
colnames(df_M) = c('Canonique', 'Jacobi', 'Hermite')

plot(df_M$Canonique ~ seq(1,1000,1), xlab="M", ylab="SCR", main="Tuning nombre de simulation seconddaire M - Cano")




#sensi d

N = 10000
M = 1000
df_d = data.frame()

for ( i in 1:20){
  S1 = S_1(Nsim = N)
  NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = i))
  #LMSC
  df_d = rbind(df_d, lsmc_calcul(d = i))
}

colnames(df_d) = c('Canonique', 'Jacobi', 'Hermite')
write.csv(df_d, paste(folder, sep = '/', "df_d"), row.names=FALSE)

plot(df_d$Canonique ~ seq(1,20,1), xlab="d", ylab="SCR", main="Tuning nombre de simulation polynomials- Cano")




####################################################################
#Confidence Interval and Time of Execution
####################################################################
table1 = data.frame()
colnames(table1) = c("Methode", "M", "N", "d","time", "SCR_inf", "SCR_sup")
#Interval de confiance SdS

M = 250
N = 25000

t1 = Sys.time()
S1 = S_1(Nsim = N)
NAV1s <- sapply(S1, function(x) NAV1(S_debut = x, Nsim_nav = M))
SCR_SdS   = NAV0 - quantile(NAV1s, p = 0.005)/(1 + r) 

t2 = Sys.time()

timeSdS = t2 - t1



hist_SdS <- sapply(1:1000, function(i) {
  S1 = S_1(Nsim = N)
  NAV1s <- sapply(S1, function(x) NAV1(S_debut = x, Nsim_nav = M))
  SCR   = NAV0 - quantile(NAV1s, p = 0.005)/(1 + r) 
  return(SCR)
})

hist(hist_SdS,probability = T,main = "Distribution de l'estimateur du SCR avec la méthode SDS", ylab="Densité")
table1 = rbind(table1, c("SdS", M, N, d = 0, timeSdS, t.test(hist_SdS)$conf.int ))

table1

#Interval confidence

N = 4000
M = 200
d = 7

#Canonique

t1 = Sys.time()
S1 = S_1(Nsim = N)
NAV1s <- sapply(S1, function(x) NAV1(S_debut = x, Nsim_nav = M))
SCR_SdS   = NAV0 - quantile(NAV1s, p = 0.005)/(1 + r) 
modele_canonique <- lm(NAV1s ~ poly (S1 ,d,raw=F))
NAV1_canonique = modele_canonique$fitted.values
NAV0 - quantile(NAV1_canonique,p=0.005)/(1+r)

t2 = Sys.time()

time = t2 - t1

SCR_cano = sapply(1:1000, function(i) {
  S1 = S_1(Nsim = N)
  NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
  modele_canonique <- lm(NAV1s ~ poly (S1 ,d,raw=F))
  NAV1_canonique = modele_canonique$fitted.values
  
  return(NAV0 - quantile(NAV1_canonique,p=0.005)/(1+r))
})

hist(SCR_cano , probability = T, main = "Méthode LSMC avec la base canonique", ylab="Densité", xlab="SCR")

table1 = rbind(table1, c("Canonique", M, N, d , time, t.test(SCR_cano)$conf.int ))
table1

#Jacobi

t1 = Sys.time()
S1 = S_1(Nsim = N)
NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
x_scale = scaleX(S1, u=-1, v=1)
jacobi = jacobi.g.polynomials(n = d, p = 3, q = 2)
x_jacobi = as.matrix(as.data.frame(polynomial.values(jacobi, x_scale))[, 2: (d+1)])
modele_jacobi = lm(NAV1s ~ x_jacobi)
NAV1_jacobi = modele_jacobi$fitted.values
NAV0 - quantile(NAV1_jacobi,p=0.005)/(1+r)
t2 = Sys.time()

time = t2 - t1

SCR_jacobi = sapply(1:1000, function(i){
  S1 = S_1(Nsim = N)
  NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
  x_scale = scaleX(S1, u=-1, v=1)
  jacobi = jacobi.g.polynomials(n = d, p = 3, q = 2)
  x_jacobi = as.matrix(as.data.frame(polynomial.values(jacobi, x_scale))[, 2: (d+1)])
  modele_jacobi = lm(NAV1s ~ x_jacobi)
  NAV1_jacobi = modele_jacobi$fitted.values
  
  return(NAV0 - quantile(NAV1_jacobi,p=0.005)/(1+r))
})

hist(SCR_jacobi , probability = T, main = "Méthode LSMC Jacobi", ylab="Densité", xlab="SCR")


table1 = rbind(table1, c("Jacobi", M, N, d , time, t.test(SCR_jacobi)$conf.int ))
table1


#Chebyshev

t1 = Sys.time()
S1 = S_1(Nsim = N)
NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
S1 = S_1(Nsim = N)
NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
chebyshev = chebyshev.t.polynomials(n = d)
x_chebyshev = as.matrix(as.data.frame(polynomial.values(chebyshev, S1)))

modele_chebyshev = lm(NAV1s ~ x_chebyshev)
NAV0 - quantile(modele_chebyshev$fitted.values,p=0.005)/(1+r)
t2 = Sys.time()

time = t2 - t1

SCR_chebyshev = sapply(1:1000, function(i){
  S1 = S_1(Nsim = N)
  NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
  chebyshev = chebyshev.t.polynomials(n = d)
  x_chebyshev = as.matrix(as.data.frame(polynomial.values(chebyshev, S1)))
  
  modele_chebyshev = lm(NAV1s ~ x_chebyshev)
  
  return(NAV0 - quantile(modele_chebyshev$fitted.values,p=0.005)/(1+r))
  
})
hist(SCR_chebyshev , probability = T, main = "Méthode LSMC Chebyshev", ylab="Densité", xlab="SCR")

table1 = rbind(table1, c("Chebyshev", M, N, d , time, t.test(SCR_chebyshev)$conf.int ))
table1



#Gegenbauer

t1 = Sys.time()
S1 = S_1(Nsim = N)
S1 = S_1(Nsim = N)
NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
gegenbauer = gegenbauer.polynomials(n = d, alpha = 3)

x_gegenbauer = as.matrix(as.data.frame(polynomial.values(gegenbauer, S1 )))
modele_gegenbauer  = lm(NAV1s ~ x_gegenbauer)
NAV0 - quantile(modele_gegenbauer$fitted.values,p=0.005)/(1+r)
t2 = Sys.time()

time = t2 - t1

SCR_gegenbauer = sapply(1:1000, function(i){
  S1 = S_1(Nsim = N)
  NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
  gegenbauer = gegenbauer.polynomials(n = d, alpha = 3)
  
  x_gegenbauer = as.matrix(as.data.frame(polynomial.values(gegenbauer, S1 )))
  modele_gegenbauer  = lm(NAV1s ~ x_gegenbauer)
  
  return(NAV0 - quantile(modele_gegenbauer$fitted.values,p=0.005)/(1+r))
})

hist(SCR_gegenbauer , probability = T, main = "Méthode LSMC Gegenbauer", ylab="Densité", xlab="SCR")

table1 = rbind(table1, c("gegenbauer", M, N, d , time, t.test(SCR_gegenbauer)$conf.int ))
table1

#Hermite

t1 = Sys.time()
S1 = S_1(Nsim = N)
NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
hermite = hermite.h.polynomials(n = d)
x_hermite = as.matrix(as.data.frame(polynomial.values(hermite, S1)))
modele_hermite = lm(NAV1s ~ x_hermite)

NAV0 - quantile(modele_hermite$fitted.values,p=0.005)/(1+r)
t2 = Sys.time()

time = t2 - t1


SCR_hermite = sapply(1:1000, function(i){
  S1 = S_1(Nsim = N)
  NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
  hermite = hermite.h.polynomials(n = d)
  x_hermite = as.matrix(as.data.frame(polynomial.values(hermite, S1)))
  modele_hermite = lm(NAV1s ~ x_hermite)
  return(NAV0 - quantile(modele_hermite$fitted.values,p=0.005)/(1+r))
  
})

hist(SCR_hermite , probability = T, main = "Méthode LSMC Hermite", ylab="Densité", xlab="SCR")


table1 = rbind(table1, c("Hermite", M, N, d , time, t.test(SCR_hermite)$conf.int ))
table1



#Laguerre

t1 = Sys.time()

S1 = S_1(Nsim = N)
NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
laguerre = laguerre.polynomials(n = d)
x_laguerre = as.matrix(as.data.frame(polynomial.values(laguerre, S1)))

modele_laguerre = lm(NAV1s ~ x_laguerre)
NAV0 - quantile(modele_laguerre$fitted.values,p=0.005)/(1+r)

t2 = Sys.time()

time = t2 - t1

SCR_laguerre = sapply(1:1000, function(i){
  S1 = S_1(Nsim = N)
  NAV1s  = sapply(S1,function(x) NAV1(S_debut = x, Nsim_nav = M))
  laguerre = laguerre.polynomials(n = d)
  x_laguerre = as.matrix(as.data.frame(polynomial.values(laguerre, S1)))
  
  modele_laguerre = lm(NAV1s ~ x_laguerre)
  
  return(NAV0 - quantile(modele_laguerre$fitted.values,p=0.005)/(1+r))
  
})

hist(SCR_laguerre , probability = T, main = "Méthode LSMC Laguerre", ylab="Densité", xlab="SCR")


table1 = rbind(table1, c("Laguerre", M, N, d , time, t.test(SCR_laguerre)$conf.int ))
table1


write.csv(table1, "D:/M2/Technique/confi et temps")




###################################################################
#III) Sensibilite
##################################################################

#1.sensibilite taux sans riques

N = 25000
M = 500
d = 7

#1.1
sensi_taux = function(r_){
  S1 <- St(Nsim = N)
  NAV1s = sapply(S1, function(x) NAV1(x, Nsim_nav = M, taux_nav  = r_))
  
  #canonique
  modele_canonique <- lm(NAV1s ~ poly (S1 ,d,raw=F))
  NAV1_canonique = modele_canonique$fitted.values
  SCR_canonique = NAV0 - quantile(NAV1_canonique,p=0.005)/(1+r)
  
  #jacobi
  x_scale = scaleX(S1, u=-1, v=1)
  jacobi = jacobi.g.polynomials(n = d, p = 3, q = 2)
  x_jacobi = as.matrix(as.data.frame(polynomial.values(jacobi, x_scale))[, 2: (d+1)])
  modele_jacobi = lm(NAV1s ~ x_jacobi)
  NAV1_jacobi = modele_jacobi$fitted.values
  SCR_jacobi = NAV0 - quantile(NAV1_jacobi,p=0.005)/(1+r)
  
  #Hermite polynomials
  
  hermite = hermite.h.polynomials(n = d)
  x_hermite = as.matrix(as.data.frame(polynomial.values(hermite, S1)))
  modele_hermite = lm(NAV1s ~ x_hermite)
  NAV1_hermite = modele_hermite$fitted.values
  SCR_hermite = NAV0 - quantile(NAV1_hermite,p=0.005)/(1+r)
  
  return(c(SCR_canonique, SCR_jacobi, SCR_hermite))
}


SCR_taux_r = data.frame()
for(r_ in c(0, 0.015, 0.02, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075)){
  SCR_taux_r = rbind(SCR_taux_r, c(r_, sensi_taux(r_)))
}

colnames(SCR_taux_r) = c('taux sans risque', 'SCR_canonique', 'SCR_jacobi', 'SCR_hermite')
write.csv(SCR_taux_r, paste(folder, sep = '/', "SCR_taux_r"), row.names=FALSE)

plot1 = read.csv(paste(folder, sep = '/', "SCR_taux_r"))

plot(plot1$SCR_canonique ~ plot1$taux.sans.risque,type = 'o', xlab="taux sans risque", ylab="SCR", main="Sensibilite en fonction de taux sans risque")




#1.2code for illustrations with more iterations

N = 4000
M = 500
d = 5
SCR_taux_r = data.frame()
for(r_ in (1: 3000)/30000){
  SCR_taux_r = rbind(SCR_taux_r, c(r_, sensi_taux(r_)))
  print(r_)
}


colnames(SCR_taux_r) = c('taux sans risque', 'SCR_canonique', 'SCR_jacobi', 'SCR_hermite')
write.csv(SCR_taux_r, paste(folder, sep = '/', "SCR_taux_r_illu"), row.names=FALSE)

plot2 = read.csv(paste(folder, sep = '/', "SCR_taux_r_illu"))

plot(plot2$SCR_canonique ~ plot2$taux.sans.risque, xlab="taux sans risque", ylab="SCR", main="Sensibilite en fonction de taux sans risque")




#2)sensibilite taux garanti
N = 25000
M = 500
d = 7

sensi_taux_g = function(r_g){
  S1 <- St(Nsim = N)
  NAV1s = sapply(S1, function(x) NAV1(x, Nsim_nav = M, rg_nav  = r_g))
  
  #canonique
  modele_canonique <- lm(NAV1s ~ poly (S1 ,d,raw=F))
  NAV1_canonique = modele_canonique$fitted.values
  SCR_canonique = NAV0 - quantile(NAV1_canonique,p=0.005)/(1+r)
  
  #jacobi
  x_scale = scaleX(S1, u=-1, v=1)
  jacobi = jacobi.g.polynomials(n = d, p = 3, q = 2)
  x_jacobi = as.matrix(as.data.frame(polynomial.values(jacobi, x_scale))[, 2: (d+1)])
  modele_jacobi = lm(NAV1s ~ x_jacobi)
  NAV1_jacobi = modele_jacobi$fitted.values
  SCR_jacobi = NAV0 - quantile(NAV1_jacobi,p=0.005)/(1+r)
  
  #Hermite polynomials
  
  hermite = hermite.h.polynomials(n = d)
  x_hermite = as.matrix(as.data.frame(polynomial.values(hermite, S1)))
  modele_hermite = lm(NAV1s ~ x_hermite)
  NAV1_hermite = modele_hermite$fitted.values
  SCR_hermite = NAV0 - quantile(NAV1_hermite,p=0.005)/(1+r)
  
  return(c(SCR_canonique, SCR_jacobi, SCR_hermite))
}


SCR_taux_rg = data.frame()
for(rg in c(0, 0.003, 0.005,  0.01, 0.02, 0.025, 0.035, 0.045, 0.055, 0.065)){
  SCR_taux_rg = rbind(SCR_taux_rg, c(rg, sensi_taux_g(rg)))
}


colnames(SCR_taux_rg) = c('taux garantie', 'SCR_canonique', 'SCR_jacobi', 'SCR_hermite')


write.csv(SCR_taux_rg, paste(folder, sep = '/', "SCR_taux_rg"), row.names=FALSE)

plot3 = read.csv(paste(folder, sep = '/', "SCR_taux_rg"))


plot(plot3$SCR_canonique ~ plot3$taux.garantie,type = 'o', xlab="rg", ylab="SCR", main="Sensibilite en fonction de taux garanti")




#illu for rg 
N = 4000
M = 500
d = 5
SCR_taux_rg = data.frame()
for(r_ in (1: 3000)/30000){
  SCR_taux_rg = rbind(SCR_taux_rg, c(r_, sensi_taux_g(r_)))
  print(r_)
}


colnames(SCR_taux_rg) = c('taux sans risque', 'SCR_canonique', 'SCR_jacobi', 'SCR_hermite')
write.csv(SCR_taux_rg, paste(folder, sep = '/', "SCR_taux_rg_illu"), row.names=FALSE)

plot4 = read.csv(paste(folder, sep = '/', "SCR_taux_rg_illu"))
plot(plot4$SCR_canonique ~ plot4$taux.sans.risque, xlab="taux garanti", ylab="SCR", main="Sensibilite en fonction de taux garanti")




























