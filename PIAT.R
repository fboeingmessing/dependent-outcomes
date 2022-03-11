##### set working directory to source file location:
setwd()


##### load packages:
library(matrixsampling)


##### source functions:
source("functions.R")


##### import data:
PIAT <- read.csv("PIAT.txt", header=FALSE, sep="\t")
PIAT <- PIAT[, 4:6]


##### sequential analysis:
set.seed(5781)
n <- nrow(PIAT)
ns <- c((1:(n%/%10))*10, n)
p <- ncol(PIAT)
k <- 1
epsilon <- 0.001
nu <- p-1+2*epsilon
Lambda <- diag(rep(2*epsilon, p))
nsim <- 1e5
PR <- array(NA, dim=c(2, 3, length(ns)),
            dimnames=list(c("Posterior", "Prior"),
                          c("H_1", "H_2", "H_3"), ns))
WOE <- matrix(NA, nrow=length(ns), ncol=2,
              dimnames=list(c(), c("WOE_12", "WOE_13")))
PMP <- matrix(NA, nrow=length(ns), ncol=3,
              dimnames=list(c(), c("H_1", "H_2", "H_3")))
for (i in ns) {
  Y <- PIAT[1:i, ]
  S <- (i-1)*cov(Y)
  post.samp <- rinvwishart(nsim, nu=nu+i-k, Omega=Lambda+S)
  post.samp <- t(sapply(1:nsim, function(x) diag(post.samp[ , , x])))
  prior.samp <- matrix(NA, nrow=nsim, ncol=p)
  count <- 1
  while (count<=nsim) {
    Sigma.samp <- try(rinvwishart(1, nu=nu, Omega=Lambda)[, , 1], silent=TRUE)
    if (class(Sigma.samp)=="try-error" | any(!is.finite(Sigma.samp))) {
      next
    } else {
      prior.samp[count, ] <- diag(Sigma.samp)
      count <- count + 1
    }
  }
  PR[1, 1, which(ns==i)] <- sum((post.samp[, 1]<post.samp[, 2])*
                                (post.samp[, 1]<post.samp[, 3]))/nsim
  PR[1, 2, which(ns==i)] <- sum((post.samp[, 2]<post.samp[, 1])*
                                (post.samp[, 3]<post.samp[, 1]))/nsim
  PR[1, 3, which(ns==i)] <- 1-sum(PR[1, c(1, 2), which(ns==i)])
  PR[2, 1, which(ns==i)] <- sum((prior.samp[, 1]<prior.samp[, 2])*
                                (prior.samp[, 1]<prior.samp[, 3]))/nsim
  PR[2, 2, which(ns==i)] <- sum((prior.samp[, 2]<prior.samp[, 1])*
                                (prior.samp[, 3]<prior.samp[, 1]))/nsim
  PR[2, 3, which(ns==i)] <- 1-sum(PR[2, c(1, 2), which(ns==i)])
}
for (i in 1:length(ns)) {
  WOE[i, ] <- bayes_factors(PR[, , i], log.BF=TRUE)[c(4, 7)]
  PMP[i, ] <- posterior_probabilities(PR[, , i])
}


##### plot weights of evidence:
ticks <- (0:10)*100
xlimits <- c(0, n)
ylimits <- c(min(WOE[is.finite(WOE)]), max(WOE[is.finite(WOE)]))
par(mar=c(5.1, 5.1, 1.1, 1.1))
plot(c((1:(n%/%10))*10, n), WOE[, 1], type="l", lty=2, lwd=2, main="", xaxt="n",
     xlab="", ylab="Weight of evidence", xlim=xlimits, ylim=ylimits,
     cex.lab=1.75, cex.axis=1.5)
for (i in 2:ncol(WOE)) lines(c((1:(n%/%10))*10, n), WOE[, i],
                             type="l", lty=c(2, 3)[i], lwd=c(2, 3)[i])
axis(side=1, at=ticks, labels=FALSE)
text(ticks, par("usr")[3]-4*(ylimits[2]-ylimits[1])/100, srt=45, adj=1,
     labels=ticks, xpd=TRUE, cex=1.5)
title(xlab="n", line=4, cex.lab=1.75)
abline(h=0)
leg <- as.expression(c(bquote(WOE[12]), bquote(WOE[13])))
legend("topleft", leg, bg="white", lty=c(2, 3), lwd=c(2, 3), cex=1.5)


##### plot posterior probabilities:
ticks <- (0:10)*100
xlimits <- c(0, n)
ylimits <- c(0, 1)
par(mar=c(5.1, 5.1, 1.1, 1.1))
plot(c((1:(n%/%10))*10, n), PMP[, 1], type="l", lty=1, lwd=2, main="", xaxt="n",
     xlab="", ylab="Posterior probability", xlim=xlimits, ylim=ylimits,
     cex.lab=1.75, cex.axis=1.5)
for (i in 2:ncol(PMP)) lines(c((1:(n%/%10))*10, n), PMP[, i],
                             type="l", lty=i, lwd=c(2, 2, 3)[i])
axis(side=1, at=ticks, labels=FALSE)
text(ticks, par("usr")[3]-4*(ylimits[2]-ylimits[1])/100, srt=45, adj=1,
     labels=ticks, xpd=TRUE, cex=1.5)
title(xlab="n", line=4, cex.lab=1.75)
leg <- as.expression(c(bquote(H[1]), bquote(H[2]), bquote(H[3])))
legend("right", leg, bg="white", lty=1:ncol(PMP), lwd=c(2, 2, 3), cex=1.5)
