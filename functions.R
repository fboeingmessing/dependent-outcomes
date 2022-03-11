########## posterior and prior probability that the inequality constraints hold:
inequality_probabilities <- function(S, n, k, hypotheses, complement=FALSE,
                                     Lambda="default", nu="default", nsim=1e5) {
  P <- ncol(S)
  if (all(Lambda=="default")) {
    Lambda <- diag(rep(1e-6, P))
  } else if (!isSymmetric(Lambda)) {
    warning("Lambda is not symmetric")
  } else if (any(eigen(Lambda)$values<=0)) {
    warning("Lambda is not positive definite")
  }
  if (nu=="default") {
    nu <- P
  } else if (!nu>P-1) {
    warning(paste("nu must be greater than", as.character(P-1)))
  }
  H <- length(hypotheses)
  Pr <- matrix(0, 2, H, dimnames=list(c("Posterior", "Prior"), hypotheses))
  compsymm <- length(unique(diag(Lambda)))==1 &&
    length(unique(Lambda[lower.tri(Lambda) | upper.tri(Lambda)]))==1
  post.samp <- matrix(NA, nrow=nsim, ncol=P)
  Psiinv <- solve(Lambda+S)
  for (i in 1:nsim) {
    Sigmainv <- rWishart(1, nu+n-k, Psiinv)[, , 1]
    post.samp[i, ] <- diag(solve(Sigmainv))
  }
  if (!compsymm) {
    prior.samp <- matrix(NA, nrow=nsim, ncol=P)
    Psiinv <- solve(Lambda)
    for (i in 1:nsim) {
      Sigmainv <- rWishart(1, nu, Psiinv)[, , 1]
      prior.samp[i, ] <- diag(solve(Sigmainv))
    }
  }
  for (i in 1:H) {
    hypothesis <- hypotheses[i]
    if (hypothesis=="order") {
      hypothesis <- paste(as.character(1:P), collapse="<")
    }
    if (hypothesis=="unconstrained" || !grepl("<", hypothesis)) {
      Pr[, i] <- c(1, 1)
      next
    }
    complement.h <- if (grepl("not", hypothesis)) {TRUE} else {FALSE}
    hypothesis <- gsub("[not ]", "", hypothesis)
    if (complement.h && grepl("r", hypothesis)) {
      hypothesis <- substr(hypothesis, start=2, stop=nchar(hypothesis)-1)
    }
    hypothesis <- unlist(strsplit(hypothesis, split="r"))
    hypothesis <- hypothesis[grepl("<", hypothesis)]
    for (j in 1:length(hypothesis)) {
      inequalities <- hypothesis[j]
      if (grepl("[(]", inequalities)) {
        brackets <- c(unlist(gregexpr(pattern="[(]", inequalities)),
                      unlist(gregexpr(pattern="[)]", inequalities)))
        brackets <- matrix(brackets, nrow=length(brackets)/2, ncol=2)
        commas <- unlist(apply(brackets, 1, function(x) unlist(gregexpr(pattern=",",
                         substr(inequalities, start=x[1], stop=x[2])))+(x[1]-1)))
        for (l in 1:length(commas)) {substring(inequalities, commas[l]) <- ";"}
      }
      inequalities <- unlist(strsplit(inequalities, split=","))
      inequalities <- as.list(inequalities[grepl("<", inequalities)])
      for (l in 1:length(inequalities)) {
        if (grepl("[(]", inequalities[[l]])) {
          combinations <- unlist(strsplit(inequalities[[l]], split="<"))
          combinations <- gsub("[()]", "", combinations)
          combinations <- strsplit(combinations, split=";")
          inequalities[[l]] <- apply(as.matrix(expand.grid(combinations)), 1,
                                     function(x) paste(x, collapse="<"))
        }
      }
      inequalities <- unlist(inequalities)
      ninequalities <- length(unlist(gregexpr(pattern="<", inequalities)))
      inequalities <- lapply(strsplit(unlist(inequalities), split="<"), as.numeric)
      variances <- unique(unlist(inequalities))
      ineqs <- rep(0, nsim)
      for (l in 1:length(inequalities)) {
        for (m in 1:(length(inequalities[[l]])-1)) {
          ineqs <- ineqs+(post.samp[, inequalities[[l]][m]]<
                            post.samp[, inequalities[[l]][m+1]])
        }
      }
      Pr[1, i] <- Pr[1, i]+sum(ineqs==ninequalities)/nsim
      if (!compsymm) {
        ineqs <- rep(0, nsim)
        for (l in 1:length(inequalities)) {
          for (m in 1:(length(inequalities[[l]])-1)) {
            ineqs <- ineqs+(prior.samp[, inequalities[[l]][m]]<
                              prior.samp[, inequalities[[l]][m+1]])
          }
        }
        Pr[2, i] <- Pr[2, i]+sum(ineqs==ninequalities)/nsim
      } else if (length(inequalities)==1) {
        Pr[2, i] <- Pr[2, i]+1/factorial(length(inequalities[[1]]))
      } else {
        intersections <- sapply(inequalities, paste, collapse=",")
        intersections <- combn(intersections, 2)
        intersections <- apply(intersections, 2,
            function(x) intersect(as.numeric(unlist(strsplit(x[1], split=","))),
                                  as.numeric(unlist(strsplit(x[2], split=",")))))
        if (length(intersections)==0) {
          Pr[2, i] <- Pr[2, i]+1/prod(factorial(sapply(inequalities, length)))
        } else {
          prior.samp <- matrix(NA, nrow=nsim, ncol=P)
          for (l in 1:length(variances)) {
            prior.samp[, variances[l]] <- runif(nsim)
          }
          ineqs <- rep(0, nsim)
          for (l in 1:length(inequalities)) {
            for (m in 1:(length(inequalities[[l]])-1)) {
              ineqs <- ineqs+(prior.samp[, inequalities[[l]][m]]<
                                prior.samp[, inequalities[[l]][m+1]])
            }
          }
          Pr[2, i] <- Pr[2, i]+sum(ineqs==ninequalities)/nsim
        }
      }
    }
    if (complement.h) {Pr[, i] <- 1-Pr[, i]}
  }
  if (complement) {
    Pr <- cbind(Pr, 1-rowSums(Pr))
    colnames(Pr)[ncol(Pr)] <- "complement"
  }
  return(Pr)
}


########## bayes factors:
bayes_factors <- function(Pr, log.BF=FALSE) {
  Btu <- Pr[1, ]/Pr[2, ]
  B <- unname(Btu%*%t(1/Btu))
  diag(B) <- 1
  if (log.BF) {B <- log(B)}
  return(B)
}


########## posterior probabilities:
posterior_probabilities <- function(Pr, prior.probabilities="default") {
  if (all(prior.probabilities=="default")) {
    prior.probabilities <- rep(1/ncol(Pr), ncol(Pr))
  } else if (!identical(sum(prior.probabilities), 1)) {
    warning("the prior probabilities do not sum to 1")
  }
  Btu <- Pr[1, ]/Pr[2, ]
  PP <- Btu*prior.probabilities/sum(Btu*prior.probabilities)
  names(PP) <- paste("H", as.character(1:length(PP)), sep="")
  return(PP)
}
