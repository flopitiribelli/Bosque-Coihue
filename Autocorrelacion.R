setwd("C:/Users/flor/Google Drive/Fuego/Datos")
library("R2WinBUGS")
library("plyr")
datos <- read.table("multicoihue.txt", header=T)
#a?o <- datos$A?o
#transecta <- datos$Transecta
#punto <- datos$Punto
#estado <- datos$Estado
#suma <- ddply(datos, .(a?o, transecta, punto, estado), numcolwise(sum))
#write.table (suma, "autocor_toques_vivo_seco.txt")
suma <- read.table("autocor_toques_vivo_seco.txt", header=T)

TSF <- 2015 - suma$A?o
W <- list()

for(i in 1:12) {
      subtrans <- subset(suma, suma$transecta == i)
      a <- subtrans$punto
      w <- numeric(length(subtrans$tq_s))
      for(j in 1:(length(subtrans$tq_s))) {            
            w1 <- subtrans$tq_s[which(subtrans$punto == (subtrans$punto[j] + 1))]
            w2 <- subtrans$tq_s[which(subtrans$punto == (subtrans$punto[j] - 1))]
            if (length(w1) > 0) {
                  w[j] <- sum(w1) 
            }
            if (length(w2) >  0) {
                  w[j] <- w[j] + sum(w2) 
            }
      } 
      
      
      W[[i]] <- cbind(a, w)
      
}

y_y <- do.call(rbind,W)

p <- cbind(suma$tq_s, 8)
y <- suma$tq_s
tsf <- TSF

autocov <- y_y[,2]

n <- length(TSF)  
N <- p[,2]
transecta <- numeric(n)
transecta[which(suma$transecta==1)] <- 1
transecta[which(suma$transecta==2)] <- 2
transecta[which(suma$transecta==3)] <- 3
transecta[which(suma$transecta==4)] <- 4
transecta[which(suma$transecta==5)] <- 5
transecta[which(suma$transecta==6)] <- 6
transecta[which(suma$transecta==7)] <- 7
transecta[which(suma$transecta==8)] <- 8
transecta[which(suma$transecta==9)] <- 9
transecta[which(suma$transecta==10)] <- 10
transecta[which(suma$transecta==11)] <- 11
transecta[which(suma$transecta==12)] <- 12

J <- max(transecta)



cat(file="autocorr.bug",
    "
    model { 
    for (i in 1:n){
    y[i] ~ dbin (p.bound[i], N[i]) 
    p.bound[i] <- max(0, min(1, p[i]))
    logit(p[i]) <- Xbeta[i]
    Xbeta[i] <- b0 + b1 * tsf[i] + b2[transecta[i]]*autocov[i]
    }
    
    for (j in 1:J) {
    b2 [ j ] ~ dnorm (mu.b2, tau.b2)
    } 
    
    b0 ~ dnorm (0, .0001)
    b1 ~ dnorm (0, .0001)
    mu.b2 ~ dnorm (0, .0001)
    sd.b2 ~ dunif (0, 100)
    tau.b2 <- 1/pow(sd.b2,2)
    }
    
    ")



library(jagsUI)


data <- list("tsf","autocov","n", "transecta", "J", "y","N")
inits <- function() list(b0 = runif(1,1,10),
                         b1 = rnorm(n=1,mean=0,sd=1),
                         b2 = rnorm(n=J,mean=0,sd=1), 
                         mu.b2 = rnorm(n=1,mean=0,sd=1),
                         sd.b2 = runif(1,1,10))

params <- c("b0","b1","b2", "mu.b2", "sd.b2")

mod <- jags(data, inits, params, model.file = "autocorr.bug", n.chains = 3,
            n.iter = 10000, n.burnin = 3000)#, n.thin = 4)


#mod <- bugs(data,inits,params, model.file="autocorr.txt", debug = T,
# n.chains=3, n.iter=30000, n.sims=3000) 
#DIC modeo lineal 1777.825 


print(mod)
mod$summary 
save(mod,file="ND_tqs_autocor.RData")

plot(jitter(TSF), jitter(suma$tq_s/8))
lines(TSF, plogis(mod$mean$b0 + mod$mean$b1*TSF), col="red", lwd=2)
curve(plogis(mod$mean$b0*(mod$mean$b1*x)), add=T, col="blue", lwd=2)


##### modelo exponencial #####


cat(file="autocorr.bug",
    "
    model { 
    for (i in 1:n){
    y[i] ~ dbin (p.bound[i], N[i]) 
    p.bound[i] <- max(0, min(1, p[i]))
    logit(p[i]) <- Xbeta[i]
    Xbeta[i] <- (b0* exp(b1 * tsf[i])) + b2[transecta[i]]*autocov[i]
    }
    
    for (j in 1:J) {
    b2 [ j ] ~ dnorm (mu.b2, tau.b2)
    } 
    
    b0 ~ dnorm (0, .0001)
    b1 ~ dnorm (0, .0001)
    mu.b2 ~ dnorm (0, .0001)
    sd.b2 ~ dunif (0, 100)
    tau.b2 <- 1/pow(sd.b2,2)
    }
    
    ")



library(jagsUI)


data <- list("tsf","autocov","n", "transecta", "J", "y","N")
inits <- function() list(b0 = runif(1,1,10),
                         b1 = rnorm(n=1,mean=0,sd=1),
                         b2 = rnorm(n=J,mean=0,sd=1), 
                         mu.b2 = rnorm(n=1,mean=0,sd=1),
                         sd.b2 = runif(1,1,10))

params <- c("b0","b1","b2", "mu.b2", "sd.b2")

mod <- jags(data, inits, params, model.file = "autocorr.bug", n.chains = 3,
            n.iter = 10000, n.burnin = 3000)#, n.thin = 4)

print(mod)
mod$summary 
save(mod,file="ND_tqs_exp_autocor.RData")