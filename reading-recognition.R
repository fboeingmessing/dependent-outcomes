##### set working directory to source file location:
setwd()


##### load packages:
library(tidyr)
library(lme4)
library(matrixsampling)
library(MASS)


##### source functions:
source("functions.R")


##### import data:
reading.full <- read.csv("reading-recognition.csv", header=TRUE, sep=";", dec=".")
reading.full <- reading.full[, 5:8]
reading.full$id <- 1:nrow(reading.full)
reading.full[, 1:4] <- reading.full[, 1:4]*10


##### settings:
n.start <- 10
nsim <- 1e5
burnin <- 100
counter <- 1e4
# set hyperparameters:
epsilon <- 0.001
p <- 4
nu <- p-1+2*epsilon
Lambda <- diag(rep(2*epsilon, p))


##### preliminaries:
niter <- nsim + burnin
p <- 4
n.full <- nrow(reading.full)
ns <- c((1:(n.full%/%10))*10, n.full)
set.seed(5814)
seeds <- sample(1:1e4, size=length(ns))
prior.samp <- array(NA, c(nsim, p, length(ns)), dimnames=list(c(), c(), ns))
post.samp <- array(NA, c(niter, p, length(ns)), dimnames=list(c(), c(), ns))


##### sample from posterior:
for (j in which(ns==n.start):length(ns)) {
  set.seed(seeds[j])
  n <- ns[j]
  reading <- reading.full[1:n, ]
  # compute missing data information:
  Y <- reading[, 1:4]
  M <- is.na(Y)
  comcases <- which(complete.cases(Y))
  incomcases <- which(!complete.cases(Y))
  Obs <- apply(Y, 1, function(x) unname(which(!is.na(x))))
  Mis <- apply(Y, 1, function(x) unname(which(is.na(x))))
  q <- sapply(Mis, length)
  # reshape data from wide to long format:
  reading <- gather(reading, key="time", value="read", read1:read4, factor_key=TRUE)
  reading <- reading[order(reading$id, reading$time), ]
  # prepare predictors:
  reading$const <- rep(1, nrow(reading))
  reading$time <- as.numeric(reading$time)# - 1
  reading$time <- reading$time - mean(reading$time)
  reading$time2 <- (reading$time)^2
  reading <- reading[, c("id", "read", "const", "time", "time2")]
  # create data matrices:
  id <- reading$id
  p <- length(unique(reading$time))
  y <- as.matrix(reading[, "read"])
  X <- unname(as.matrix(reading[, c("const", "time", "time2")]))
  k <- ncol(X)
  # create matrices for storing the draws from the full conditionals:
  Sigma <- array(NA, dim=c(p, p, niter))
  beta <- matrix(NA, nrow=niter, ncol=k)
  # obtain starting values using lme4:
  fit.lmer <- lmer(read ~ (1 | id) + time + time2, data=reading)
  Sigma[, , 1] <- summary(fit.lmer)$varcor$id[1, 1] + diag(rep(sigma(fit.lmer)^2, p))
  post.samp[1, , j] <- diag(Sigma[, , 1])
  beta[1, ] <- fixef(fit.lmer)
  ymism <- is.na(y)*X%*%as.matrix(beta[1, ])
  ymism[ymism==0] <- NA
  # gibbs sampler:
  for (m in 2:niter) {
    yimp <- replace(y, is.na(y), ymism[!is.na(ymism)])
    # draw Sigma:
    S <- matrix(0, nrow=p, ncol=p)
    e <- yimp-X%*%as.matrix(beta[m-1, ])
    for (i in 1:n) {
      S <- S + e[id==i, ]%*%t(e[id==i, ])
    }
    Sigma[, , m] <- rinvwishart(1, nu=nu+n, Omega=Lambda+S)[, , 1]
    post.samp[m, , j] <- diag(Sigma[, , m])
    # draw beta:
    Sigmainv <- solve(Sigma[, , m])
    Sigmaninv <- kronecker(diag(n), Sigmainv)
    Sigmabeta <- solve(t(X)%*%Sigmaninv%*%X)
    mubeta <- Sigmabeta%*%t(X)%*%Sigmaninv%*%yimp
    beta[m, ] <- mvrnorm(1, mu=mubeta, Sigma=Sigmabeta)
    # draw missing values:
    for (i in incomcases) {
      Obsi <- Obs[[i]]
      Misi <- Mis[[i]]
      Sigmareord <- Sigma[c(Misi, Obsi), c(Misi, Obsi), m]
      Sigmamm <- as.matrix(Sigmareord[1:q[i], 1:q[i]])
      Sigmamo <- matrix(Sigmareord[1:q[i], (q[i]+1):p], nrow=q[i])
      Sigmaooinv <- solve(as.matrix(Sigmareord[(q[i]+1):p, (q[i]+1):p]))
      yimpi <- as.matrix(yimp[id==i, ])
      Xi <- X[id==i, ]
      mumis <- Xi[Misi, ]%*%as.matrix(beta[m, ]) + Sigmamo%*%Sigmaooinv%*%
      (as.matrix(yimpi[Obsi, ])-Xi[Obsi, ]%*%as.matrix(beta[m, ]))
      Sigmamis <- Sigmamm - Sigmamo%*%Sigmaooinv%*%t(Sigmamo)
      ymism[id==i, ][Misi] <- mvrnorm(1, mu=mumis, Sigma=Sigmamis)
    }
    # print iteration number:
    if(m%%counter==0) {
      cat(paste("Iteration", m, "complete.\n"))
    }
  }
  cat(paste("\nSample size", n, "complete.\n\n"))
}
post.samp <- post.samp[-(1:burnin), , ]


##### sample from prior:
set.seed(9076)
for (j in 1:length(ns)) {
  count <- 1
  while (count<=nsim) {
    Sigma.samp <- try(rinvwishart(1, nu=nu, Omega=Lambda)[, , 1], silent=TRUE)
    if (class(Sigma.samp)=="try-error" | any(!is.finite(Sigma.samp))) {
      next
    } else {
      prior.samp[count, , j] <- diag(Sigma.samp)
      count <- count + 1
    }
  }
}


##### compute bayes factors and posterior probabilities:
PR <- array(NA, dim=c(2, 4, length(ns)),
            dimnames=list(c("Posterior", "Prior"),
                          c("H_1", "H_2", "H_3", "H_4"), ns))
WOE <- matrix(NA, nrow=length(ns), ncol=3,
              dimnames=list(c(), c("WOE_21", "WOE_23", "WOE_24")))
PMP <- matrix(NA, nrow=length(ns), ncol=4,
              dimnames=list(c(), c("H_1", "H_2", "H_3", "H_4")))
for (j in 1:length(ns)) {
  # posterior:
  sigma2 <- post.samp[, , j]
  simpleorder <- (sigma2[, 1]<sigma2[, 2])*
                 (sigma2[, 2]<sigma2[, 3])*
                 (sigma2[, 3]<sigma2[, 4])
  ratioincrease <- (sigma2[, 2]/sigma2[, 1]<sigma2[, 3]/sigma2[, 2])*
                   (sigma2[, 3]/sigma2[, 2]<sigma2[, 4]/sigma2[, 3])
  ratiodecrease <- (sigma2[, 4]/sigma2[, 3]<sigma2[, 3]/sigma2[, 2])*
                   (sigma2[, 3]/sigma2[, 2]<sigma2[, 2]/sigma2[, 1])
  PR[1, 1, j] <- sum(simpleorder*ratioincrease)/nsim
  PR[1, 2, j] <- sum(simpleorder*ratiodecrease)/nsim
  PR[1, 4, j] <- 1-sum(simpleorder)/nsim
  PR[1, 3, j] <- 1-sum(PR[1, c(1, 2, 4), j])
  # prior:
  sigma2 <- prior.samp[, , j]
  simpleorder <- (sigma2[, 1]<sigma2[, 2])*
                 (sigma2[, 2]<sigma2[, 3])*
                 (sigma2[, 3]<sigma2[, 4])
  ratioincrease <- (sigma2[, 2]/sigma2[, 1]<sigma2[, 3]/sigma2[, 2])*
                   (sigma2[, 3]/sigma2[, 2]<sigma2[, 4]/sigma2[, 3])
  ratiodecrease <- (sigma2[, 4]/sigma2[, 3]<sigma2[, 3]/sigma2[, 2])*
                   (sigma2[, 3]/sigma2[, 2]<sigma2[, 2]/sigma2[, 1])
  PR[2, 1, j] <- sum(simpleorder*ratioincrease)/nsim
  PR[2, 2, j] <- sum(simpleorder*ratiodecrease)/nsim
  PR[2, 4, j] <- 1-sum(simpleorder)/nsim
  PR[2, 3, j] <- 1-sum(PR[2, c(1, 2, 4), j])
  # bayes factors and posterior probabilities:
  WOE[j, ] <- bayes_factors(PR[, , j], log.BF=TRUE)[c(2, 10, 14)]
  PMP[j, ] <- posterior_probabilities(PR[, , j])
}


##### plot weights of evidence:
ticks <- (0:8)*50
xlimits <- c(0, n.full)
ylimits <- c(min(WOE[is.finite(WOE)]), max(WOE[is.finite(WOE)])+1.0)
par(mar=c(5.1, 5.1, 1.1, 1.1))
plot(c((1:(n.full%/%10))*10, n.full), WOE[, 1], type="l", lty=1, lwd=2,
     main="", xaxt="n", xlab="", ylab="Weight of evidence",
     xlim=xlimits, ylim=ylimits, cex.lab=1.75, cex.axis=1.5)
for (i in 2:ncol(WOE)) lines(c((1:(n.full%/%10))*10, n.full), WOE[, i],
                             type="l", lty=c(1, 3, 4)[i], lwd=c(2, 3, 2)[i])
axis(side=1, at=ticks, labels=FALSE)
text(ticks, par("usr")[3]-4*(ylimits[2]-ylimits[1])/100, srt=45, adj=1,
     labels=ticks, xpd=TRUE, cex=1.5)
title(xlab="n", line=4, cex.lab=1.75)
abline(h=0)
leg <- as.expression(c(bquote(WOE[21]), bquote(WOE[23]), bquote(WOE[24])))
legend("topleft", leg, bg="white", lty=c(1, 3, 4), lwd=c(2, 3, 2), cex=1.5)


##### plot posterior probabilities:
ticks <- (0:8)*50
xlimits <- c(0, n.full)
ylimits <- c(0, 1)
par(mar=c(5.1, 5.1, 1.1, 1.1))
plot(c((1:(n.full%/%10))*10, n.full), PMP[, 1], type="l", lty=1, lwd=2,
     main="", xaxt="n", xlab="", ylab="Posterior probability",
     xlim=xlimits, ylim=ylimits, cex.lab=1.75, cex.axis=1.5)
for (i in 2:ncol(PMP)) lines(c((1:(n.full%/%10))*10, n.full), PMP[, i],
                             type="l", lty=i, lwd=c(2, 2, 3, 2)[i])
axis(side=1, at=ticks, labels=FALSE)
text(ticks, par("usr")[3]-4*(ylimits[2]-ylimits[1])/100, srt=45, adj=1,
     labels=ticks, xpd=TRUE, cex=1.5)
title(xlab="n", line=4, cex.lab=1.75)
leg <- as.expression(c(bquote(H[1]), bquote(H[2]), bquote(H[3]), bquote(H[4])))
legend("topright", leg, bg="white", lty=1:ncol(PMP), lwd=c(2, 2, 3, 2), cex=1.5)
