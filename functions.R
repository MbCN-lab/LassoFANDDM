#-----------------------------------------------------
#----- Title: Lasso FA NDDM (functions)
#----- Author: Inhan Kang, Woojong Yi, and Dr. Brandon Turner
#----- Model-based Cognitive Neuroscience Lab
#----- The Ohio State University
#-----------------------------------------------------


single.trial.ddm <- function(t, DE.list){
#-- a function to calculate the log likelihood value of a given set of single-trial ddm parameters
  rt <- DE.list$rt[t]
  resp <- DE.list$resp[t]
  resp12 <- 1;
  if(resp == "lower") resp12 <- 2;
  
  n.trials <- length(rt);
  neural <- DE.list$neural;
  
  # indicies
  k <- DE.list$k; 
  q <- DE.list$q;
  p <- k + q;
  
  
  # read the single-trial parameters and the boundary separation
  ter1 <- DE.list$ter[t]
  bias1 <- DE.list$bias[t]
  xi1 <- DE.list$xi[t]
  a <- DE.list$a
  
  #----- DDM likelihood
  if(ter1 > rt | ter1 < 0.0001 | bias1 > 0.99 | bias1 < 0.01 | xi1 > 15 | xi1 < -15){
    DE.lik <- -Inf
  }else{
    DE.lik <-  dwiener(rt, alpha = a, tau = ter1, beta = bias1, delta = xi1, resp = resp, give_log = TRUE) 
  }
  
  #----- FA likelihood (in fact, prior for the DDM parameters)
  Lam <- matrix(0, nrow = k+q, ncol = k)
  Lam[1:k, 1:k] <- diag(1, k);
  Lam[-(1:k), ] <- DE.list$Lam
  Phi <- matrix(0, k, k);
  diag(Phi) <- DE.list$Phi[1:k];
  Phi[upper.tri(Phi)] <- DE.list$Phi[-c(1:k)]
  Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)]
  Psi <- diag(c(rep(0, k), DE.list$Psi));
  mu <- DE.list$mu;
  sig.mat <- Lam %*% Phi %*% t(Lam) + Psi
  dat <- c(log(ter1), log(bias1/(1-bias1)), xi1, neural[t,])
  SEM.lik <- dmvnorm(dat, mean = mu, sigma = sig.mat, log = TRUE)
  
  return(list(DE.lik = DE.lik, SEM.lik = SEM.lik))
}



fa.nddm <- function(DE.list){ 
#-- a function to calculate the likelihood given all the data and parameters
#-- mostly used to track likelihood values

  rt <- DE.list$rt;
  resp <- DE.list$resp;
  resp12 <- rep(1, length(resp));
  resp12[resp == "lower"] <- 2;
  n.trials <- length(rt);
  neural <- DE.list$neural;
  
  # indicies
  k <- DE.list$k; 
  q <- DE.list$q;
  p <- k + q; 
  
  
  #----- DDM likelihood
  a <- DE.list$a;
  bias <- DE.list$bias;
  xi <- DE.list$xi;
  ter <- DE.list$ter;
  
  ddm.log.lik.vec <- numeric(n.trials)
  if(any(ter > rt) | any(bias > 0.99) | any(bias < 0.01)  | any(xi > 15) | any(xi < -15)){
    ddm.log.lik.vec <- -Inf
    ddm.log.lik <- -Inf
  }else{ 
    for(t in 1:n.trials){ ddm.log.lik.vec[t]  <- dwiener(rt[t], alpha = a, tau = ter[t], beta = bias[t], delta = xi[t], resp = resp[t], give_log = TRUE) }
    ddm.log.lik <- sum(ddm.log.lik.vec);
  }
  
  #----- FA likelihood (in fact, prior for the DDM parameters)
  Lam <- matrix(0, nrow = k+q, ncol = k)
  Lam[1:k, 1:k] <- diag(1, k);
  Lam[-(1:k), ] <- DE.list$Lam
  Phi <- matrix(0, k, k);
  diag(Phi) <- DE.list$Phi[1:k];
  Phi[upper.tri(Phi)] <- DE.list$Phi[-c(1:k)]
  Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)]
  Psi <- diag(c(rep(0, k), DE.list$Psi));
  mu <- DE.list$mu;
  sig.mat <- Lam %*% Phi %*% t(Lam) + Psi
  dat <- cbind(log(ter), log(bias/(1-bias)), xi, neural)
  
  sig.det <- det(sig.mat)   # to check the positive definiteness of the implied cov matrix
  if(sig.det <= 0){ 
    sem.log.lik <- -Inf;
  }else{
    sem.log.lik.vec <- dmvnorm(dat, mean = mu, sigma = sig.mat, log = TRUE)
    sem.log.lik <- sum(sem.log.lik.vec ) 
  }
  log.lik <- sum(ddm.log.lik + sem.log.lik); # the total log lik.
  return(list(log.lik = log.lik, ddm.log.lik = ddm.log.lik, sem.log.lik = sem.log.lik));
  
}

fa.nddm.prior <- function(DE.list){
#-- a function to calculate the log prior of the parameters
#-- some uniform priors are redundant (e.g., since some parameers have the linking function as their priors / parameters are sampled from their conditional posterior)
#-- but left them to monitor the log prior values of the current sample and check if they are within their valid ranges.
  
  #----- DDM part
  a <- DE.list$a;
  bias <- DE.list$bias;
  xi <- DE.list$xi;
  ter <- DE.list$ter;
 
  # indicies
  k <- DE.list$k; 
  q <- DE.list$q;
  p <- k + q;
  
  ddm.log.prior.vec <- c(log(dunif(a, 0, 10)), log(dunif(bias,0, 1)), log(dunif(xi, -15, 15)), log(dunif(ter, 0, max(rt))))
  ddm.log.prior <- sum(ddm.log.prior.vec)
  if(any(ter > rt)) ddm.log.prior <- -Inf;
  
  #----- SEM part
  #-- hyperparameters for priors
  hyp.alpha.psi <- 10;  hyp.beta.psi <- 4; 
  hyp.alpha.kap <- 1;   hyp.beta.kap <- 0.05; 

  hyp.rho <- 8;
  hyp.phi <- diag(5, 3);
  hyp.H <- 0.25;  
  hyp.Lam0 <- 0;
  
  Lam <- matrix(0, nrow = k+q, ncol = k)
  Lam[1:k, 1:k] <- diag(1, k);
  Lam[-(1:k), ] <- DE.list$Lam
  Lam.valid <- Lam[-(1:k), ]
  Psi <- DE.list$Psi
  Phi <- matrix(0, k, k);
  diag(Phi) <- DE.list$Phi[1:k];
  Phi[upper.tri(Phi)] <- DE.list$Phi[-c(1:k)]
  Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)]
  mu <- DE.list$mu;
  
  lam.log.prior <- matrix(0, q, k)
  for(j in 1:q){ lam.log.prior[j,] <- dnorm(Lam.valid[j,], mean = hyp.Lam0, sd = sqrt(hyp.H * Psi[j]), log = TRUE)  }  # not the lasso prior, but just for the monitoring purpose // not used for the updating, and can be removed.
  
  sem.log.prior.vec <- c(lam.log.prior,
                         log(dinvwishart(Phi, nu = hyp.rho, S = hyp.phi)),
                         log(dinvgamma(DE.list$Psi, shape = hyp.alpha.psi, scale = 1/hyp.beta.psi)),
                         log(dnorm(mu, 0, 100)))
  sem.log.prior <- sum(sem.log.prior.vec)
  if(is.na(sem.log.prior)) sem.log.prior <- -Inf
  if(any(DE.list$Phi[1:k] < 0) | any(DE.list$Psi < 0)) sem.log.prior <- -Inf; 
  
  return(list(log.prior = ddm.log.prior + sem.log.prior, ddm.log.prior = ddm.log.prior, sem.log.prior = sem.log.prior))  
}



DEMC.FA <- function(j, theta, theta.ddm.lik, theta.sem.lik, theta.ddm.prior, theta.sem.prior, idx.list, data.list, de.ac.ddm1, de.ac.ddm2, b = 0.005){
#-- a function to do Bayesian updating: a combination of the conditional posterior sampling and the DEMC sampling
  
  # import indices
  n.chains <- idx.list$n.chains; n.trials <- idx.list$n.trials; npar <- idx.list$npar;
  sem.idx <- idx.list$sem.idx; ddm.idx <- idx.list$ddm.idx; reg.idx <- idx.list$reg.idx;
  
  Lam.idx <- idx.list$Lam.idx; Phi.idx <- idx.list$Phi.idx; Psi.idx  <- idx.list$Psi.idx; mu.idx  <- idx.list$mu.idx;
  kap.sq.idx <- idx.list$kap.sq.idx; tau.sq.idx <- idx.list$tau.sq.idx
  a.idx <- idx.list$a.idx; tau.idx <- idx.list$tau.idx; ter.idx <- tau.idx;
  w.idx <- idx.list$w.idx; bias.idx <- w.idx; xi.idx <- idx.list$xi.idx; par.idx <- 1:npar;
  kap.sq.idx <- idx.list$kap.sq.idx; tau.sq.idx <- idx.list$tau.sq.idx
  
  
  k <- idx.list$k; q <- idx.list$q; p <- k+q;

  # import data
  rt <- data.list$rt; resp <- data.list$resp; neural <- data.list$neural;
  
  
  #------------------------------------------------------------------------
  #----- update the FA and regularization parameters with posterior distributions
  #------------------------------------------------------------------------
  
  N <- length(rt)
  dat <- cbind(theta[tau.idx, j], theta[w.idx, j], theta[xi.idx, j], neural)
  
  Lam <- matrix(0, nrow = k+q, ncol = k)
  Lam[1:k, 1:k] <- diag(1, k);
  Lam[-(1:k), ] <- theta[Lam.idx, j]
  Phi <- matrix(0, k, k);  Phi0 <- theta[Phi.idx, j]
  diag(Phi) <- Phi0[1:k];
  Phi[upper.tri(Phi)] <- Phi0[-c(1:k)]
  Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)]
  Psi <- diag(c(rep(0, k), theta[Psi.idx,j]));
  kap.sq <- theta[kap.sq.idx, j]
  tau.sq <- matrix(theta[tau.sq.idx, j], nrow = q, ncol = k)
  
  #----- update mu
  curr.mean <- cbind(apply(dat, 2, mean))   # sample mean of the current sample (ddm single params and neural data)
  curr.sig <- Lam %*% Phi %*% t(Lam) + Psi  # implied cov from the current sample
  
  inv.curr.sig <- tryCatch({
    inv.curr.sig <- solve(curr.sig)
  }, warning = function(w) {
    w = solve(curr.sig + diag(rnorm(1, 0, 0.01), q+k))
    return(w)
  }, error = function(e) {
    e = solve(curr.sig + diag(rnorm(1, 0, 0.01), q+k))
    return(e)
  })
  
  prior.mean <- cbind(rep(0, k+q));
  prior.sig <- diag(100, k+q);
  inv.prior.sig <- diag(1/diag(prior.sig))
  
  mu.var <- tryCatch({
    mu.var <- solve(inv.prior.sig + N*inv.curr.sig)
  }, warning = function(w) {
    w = solve(inv.prior.sig + N*inv.curr.sig + diag(rnorm(1, 0, 0.01), q+k))
    return(w)
  }, error = function(e) {
    e = solve(inv.prior.sig + N*inv.curr.sig + diag(rnorm(1, 0, 0.01), q+k))
    return(e)
  })
  
  mu.mean <- mu.var %*% (inv.prior.sig %*% prior.mean + N*inv.curr.sig %*% curr.mean)
  new.mu <- rmvnorm(1, mean = mu.mean, sigma = mu.var)
  old.mu <- theta[mu.idx, j] 
  theta[mu.idx, j] <- new.mu
  
  #----- calculate the factor scores
  sig.omega <- Phi
  old.mu.mat <- t(matrix(rep(old.mu[1:k], n.trials), k, n.trials))
  omega.sample <- dat[, 1:k] - old.mu.mat
  omega.sq <- t(omega.sample) %*% omega.sample
  
  
  #----- update psi and lambda 
  #----- hyperparameter values for the priors
  hyp.alpha.psi <- 10;  hyp.beta.psi <- 4;
  hyp.alpha.kap <- 1;   hyp.beta.kap <- 0.05; 
  hyp.rho <- 8
  hyp.phi <- diag(5, 3);
  

  Lam.valid <- Lam[-c(1:k), ]
  Psi.valid <- Psi[-c(1:k), -c(1:k)]
  Lam.valid.new <- matrix(0, q, k)
  Psi.valid.new <- numeric(q)
  hyp.Lam.post.mat <- matrix(0, 3, q)
  hyp.alpha.psi.post <- hyp.alpha.psi + N/2;
  
  for(i in 1:q){
    hyp.H <- diag(tau.sq[i,])
    inv.hyp.H.post <-  diag(1/diag(hyp.H)) + omega.sq 
    hyp.H.post <- tryCatch({
      hyp.H.post <- solve(inv.hyp.H.post)
    }, warning = function(w) {
      w = solve(inv.hyp.H.post + diag(rnorm(1, 0, 0.01), k))
      return(w)
    }, error = function(e) {
      e = solve(inv.hyp.H.post + diag(rnorm(1, 0, 0.01), k))
      return(e)
    })
    hyp.Lam.post <- hyp.H.post %*% (t(omega.sample) %*% neural[,i]) # note Lam0 = 0
    Lam.valid.new[i, ] <- mvrnorm(1, mu = hyp.Lam.post, Sigma = diag(Psi.valid)[i] * hyp.H.post)  # sample lambda from the condi.post
    hyp.beta.psi.post <- hyp.beta.psi + 1/2*( sum(neural[,i]^2) - t(hyp.Lam.post) %*% inv.hyp.H.post %*% hyp.Lam.post) # sample psi from the condi.post
    Psi.valid.new[i] <- 1/rgamma(1, shape = hyp.alpha.psi.post, rate = hyp.beta.psi.post)
  }  
  
  #-- Solutions for the sign-switching problem should be used here.
  #-- e.g.,
  if(Lam.valid.new[33, 1] < 0) Lam.valid.new[ , 1] <- -Lam.valid.new[ , 1]
  if(Lam.valid.new[59, 2] < 0) Lam.valid.new[ , 2] <- -Lam.valid.new[ , 2]
  if(Lam.valid.new[57, 3] < 0) Lam.valid.new[ , 3] <- -Lam.valid.new[ , 3]

  theta[Lam.idx, j] <- Lam.valid.new
  theta[Psi.idx, j] <- Psi.valid.new
  
  
  #----- update tau.sq
  mu.tau.sq <- tau.sq.new <- matrix(0, q, k)
  for(i in 1:q){ 
    mu.tau.sq[i, ] <- sqrt(kap.sq * Psi.valid.new[i] / (Lam.valid.new[i,]^2))
    tau.sq.new[i, ] <- 1/rinvgauss(k, mean = mu.tau.sq[i,], dispersion = 1/kap.sq) }
  theta[tau.sq.idx, j] <- tau.sq.new

  #----- update kappa
  hyp.alpha.kap.post <- hyp.alpha.kap + k * q
  hyp.beta.kap.post <- hyp.beta.kap + 1/2 * sum(tau.sq.new)
  kap.sq.new <- rgamma(1, shape = hyp.alpha.kap.post, rate = hyp.beta.kap.post)
  kap.sq.new[which(kap.sq.new < 1e-04)] <- 1e-04  # an arbitrary constraint for numerical stability (kap.sq > 0)
  theta[kap.sq.idx, j] <- kap.sq.new


  #----- update phi 
  Phi.new <- rinvwishart(nu = N + hyp.rho, S = hyp.phi + omega.sq)
  theta[Phi.idx, j] <- c(diag(Phi.new), Phi.new[upper.tri(Phi.new)])

  
  #----- re-evaluate the log lik & prior values 
  DE.list <- list(rt = rt, resp = resp, neural = neural, k = k, q = q, p = k+q,
                  Lam = theta[Lam.idx, j], Phi = theta[Phi.idx, j], Psi = theta[Psi.idx, j], mu = theta[mu.idx, j],
                  a = theta[a.idx, j], ter = exp(theta[tau.idx, j]), bias = exp(theta[w.idx, j])/(exp(theta[w.idx, j])+1), xi = theta[xi.idx, j],
                  kap.sq = theta[kap.sq.idx, j], tau.sq = theta[tau.sq.idx, j]);
  DE.log.lik.list <- fa.nddm(DE.list)
  DE.log.prior <- fa.nddm.prior(DE.list)
  theta.sem.lik[j] <- DE.log.lik.list$sem.log.lik
  theta.sem.prior[j] <- DE.log.prior$sem.log.prior;
  
  
  #------------------------------------------------------------------------
  #----- update the DDM parameters with posterior distributions
  #------------------------------------------------------------------------
  
  #----- update a
  gam <- runif(1, 0.5, 1)
  DE.proposal <- theta[,j]
  DE.proposal[a.idx] <- -999
  
  #-- make a DEMC proposal
  while(any(DE.proposal[a.idx] <= 0)){
    chain.numbers <- (1:n.chains)[-j] # chain numbers except the current one
    chain.idx <- sample(chain.numbers, 2);
    chain.sample <- theta[a.idx, chain.idx];
    DE.proposal[a.idx] <- theta[a.idx, j] + gam * (chain.sample[2]-chain.sample[1]) + runif(1, -b, b);
  }
  
  #-- Calculate lik and priors to evaluate the new proposal
  DE.list <- list(rt = rt, resp = resp, neural = neural, k = k, q = q, p = k+q,
                  Lam = DE.proposal[Lam.idx], Phi = DE.proposal[Phi.idx], Psi = DE.proposal[Psi.idx], mu = DE.proposal[mu.idx],
                  a = DE.proposal[a.idx], ter = exp(DE.proposal[tau.idx]), bias = exp(DE.proposal[w.idx])/(exp(DE.proposal[w.idx])+1), xi = DE.proposal[xi.idx],
                  kap.sq = DE.proposal[kap.sq.idx], tau.sq = DE.proposal[tau.sq.idx]);
  DE.log.lik.list <- fa.nddm(DE.list)
  DE.log.prior <- fa.nddm.prior(DE.list)
  ll.proposal <- DE.log.lik.list$ddm.log.lik + DE.log.prior$ddm.log.prior;
  ll.theta <- theta.ddm.lik[j] + theta.ddm.prior[j]  
  
  #-- Decide if accept or reject the new proposal
  log.r <- ll.proposal - ll.theta
  if(is.na(log.r)) log.r <- -10000
  
  if( log.r > log(runif(1)) ){
    theta[a.idx, j] <- DE.proposal[a.idx];
    de.ac.ddm1[j] <- de.ac.ddm1[j]+1;
    theta.ddm.lik[j] <- DE.log.lik.list$ddm.log.lik;
    theta.ddm.prior[j] <- DE.log.prior$ddm.log.prior;
  }
  
  
  #----- update single trial parameters
  update.idx <- c(tau.idx, w.idx, xi.idx);  # indicies for all single-trial parameters -- parameters to be updated
  single.trial.theta <- theta[update.idx, j]
  gam <- runif(1, 0.5, 1)  
  DE.proposal <- theta[, j]
  single.accepted <- 0
  
  for(t in 1:n.trials){
    single.update.idx <- c(tau.idx[t], w.idx[t], xi.idx[t]);
    single.ter <- -999
    single.bias <- -999
    
    while.idx <- 0
    while(single.ter < min(rt[t]/5, 0.1) | single.ter > rt[t] | single.bias > 0.99 | single.bias < 0.01){
      if(while.idx > 20){
        DE.proposal[single.update.idx] <- theta[single.update.idx, j] 
        sfCat(t,":", "Bad behavior in the while loop. Keep the current sample. \n")
        break;
      }
      
      chain.numbers <- (1:n.chains)[-j]     # chain numbers except the current one
      chain.idx <- sample(chain.numbers, 2);
      chain.sample <- theta[single.update.idx, chain.idx];
      DE.proposal[single.update.idx] <- theta[single.update.idx, j] + gam * (chain.sample[,2]-chain.sample[,1]) + runif(length(single.update.idx), -b, b);
    
      single.DE.par <- DE.proposal[single.update.idx]
      single.ter <- exp(single.DE.par[1])
      single.bias <- exp(single.DE.par[2])/(exp(single.DE.par[2])+1)
      while.idx <- while.idx + 1;
    }
  } # end of t-loop
  
  #-- Calculate log lik and prior values for single trials
  # current sample
  theta.list <- DE.list 
  theta.list$a <- theta[a.idx, j]  # take the boundary separation updated in the previous section
  theta.ddm.res <- matrix(unlist(lapply(1:n.trials, single.trial.ddm, DE.list = theta.list)), nrow = 2, ncol = n.trials)
  theta.ddm.val <- colSums(theta.ddm.res)
  
  # new proposal
  DE.list <- list(rt = rt, resp = resp, neural = neural, k = k,  q = q, p = k+q,
                  Lam = DE.proposal[Lam.idx], Phi = DE.proposal[Phi.idx], Psi = DE.proposal[Psi.idx], mu = DE.proposal[mu.idx],
                  a = DE.proposal[a.idx], ter = exp(DE.proposal[tau.idx]), bias = exp(DE.proposal[w.idx])/(exp(DE.proposal[w.idx])+1), xi = DE.proposal[xi.idx],
                  kap.sq = DE.proposal[kap.sq.idx], tau.sq = DE.proposal[tau.sq.idx]);
  DE.ddm.res <- matrix(unlist(lapply(1:n.trials, single.trial.ddm, DE.list = DE.list)), nrow = 2, ncol = n.trials)
  DE.ddm.val <- colSums(DE.ddm.res)
  
  #--- accept or reject?  
  log.r.vec <- DE.ddm.val - theta.ddm.val
  log.r.vec[which(is.na(log.r.vec))] <- -10000
  log.r.vec[which(log.r.vec == -Inf)] <- -10000
  
  
  single.trial.par <- DE.proposal[update.idx]
  
  change.idx <- which(log.r.vec > log(runif(n.trials)))   # indicate which samples should be accepted
  single.trial.theta[c(change.idx, change.idx + n.trials, change.idx + 2*n.trials)] <- single.trial.par[c(change.idx, change.idx + n.trials, change.idx + 2*n.trials)] 
  single.accepted <- length(change.idx)
  theta[update.idx, j]  <- single.trial.theta;
  de.ac.ddm2[j] <- de.ac.ddm2[j] + single.accepted/n.trials    # average acceptance of the single-trial DDM parameters
  
  
  #----- re-evaluate the log lik & prior after the updating
  DE.list <- list(rt = rt, resp = resp, neural = neural, k = k, q = q, p = k+q,
                  Lam = theta[Lam.idx, j], Phi = theta[Phi.idx, j], Psi = theta[Psi.idx, j], mu = theta[mu.idx, j],
                  a = theta[a.idx, j], ter = exp(theta[tau.idx, j]), bias = exp(theta[w.idx, j])/(exp(theta[w.idx, j])+1), xi = theta[xi.idx, j],
                  kap.sq = theta[kap.sq.idx, j], tau.sq = theta[tau.sq.idx, j]);
  DE.log.lik.list <- fa.nddm(DE.list)
  DE.log.prior <- fa.nddm.prior(DE.list)
  
  theta.ddm.lik[j] <- DE.log.lik.list$ddm.log.lik;
  theta.ddm.prior[j] <- DE.log.prior$ddm.log.prior;  
  theta.sem.lik[j] <- DE.log.lik.list$sem.log.lik;
  theta.sem.prior[j] <- DE.log.prior$sem.log.prior;
  
  return(list(theta[,j], theta.ddm.lik[j], theta.sem.lik[j], theta.ddm.prior[j], theta.sem.prior[j], de.ac.ddm1[j], de.ac.ddm2[j]))
} # End of function  




fa.nddm.migrate <- function(theta, theta.ddm.lik, theta.sem.lik, theta.ddm.prior, theta.sem.prior, idx.list, data.list){
  
  # import indices
  n.chains <- idx.list$n.chains; n.trials <- idx.list$n.trials; npar <- idx.list$npar;
  sem.idx <- idx.list$sem.idx; ddm.idx <- idx.list$ddm.idx; reg.idx <- idx.list$reg.idx;
  Lam.idx <- idx.list$Lam.idx; Phi.idx <- idx.list$Phi.idx; Psi.idx  <- idx.list$Psi.idx; mu.idx  <- idx.list$mu.idx;
  kap.sq.idx <- idx.list$kap.sq.idx; tau.sq.idx <- idx.list$tau.sq.idx
  a.idx <- idx.list$a.idx; tau.idx <- idx.list$tau.idx; ter.idx <- tau.idx;
  w.idx <- idx.list$w.idx; bias.idx <- w.idx; xi.idx <- idx.list$xi.idx; par.idx <- 1:npar;
  las.idx <- c(sem.idx, reg.idx)  
  k <- idx.list$k; q <- idx.list$q; p <- k+q;

  # import data
  rt <- data.list$rt; resp <- data.list$resp; neural <- data.list$neural;
  
  # tuning param
  b.ddm <- 0.005; b.sem <- 0.005;
  
 
  #----- start migration
  lnum1 <- sample(c(1:n.chains),1)										   # determine how many groups to work with
  lnum2 <- sort(sample(c(1:n.chains),lnum1,replace=F))	 # which groups we will work with
  thetaset <- matrix(NA, npar, lnum1)									   # initialize
  currentset.ddm <- propset.ddm <- propw.ddm <- currw.ddm <- propw.prior.ddm <- currw.prior.ddm <- numeric(lnum1)
  currentset.sem <- propset.sem <- propw.sem <- currw.sem <- propw.prior.sem <- currw.prior.sem <- numeric(lnum1)
  
  index <- numeric(lnum1)
  for(i in 1:lnum1){
    index[i] <- sample(1:n.chains,1,replace=F)	
    thetaset[,i] <- theta[ , lnum2[i]] + c(runif(length(sem.idx), -b.sem, b.sem), runif(length(ddm.idx),-b.ddm,b.ddm), rep(0, times = length(reg.idx)))				# create a set of these particles to swap
    DE.proposal <- thetaset[, i];
    
    Lam.DE <- DE.proposal[Lam.idx];
    Phi.DE <- DE.proposal[Phi.idx];
    Psi.DE <- DE.proposal[Psi.idx];
    mu.DE <- DE.proposal[mu.idx];
    
    a.DE <- DE.proposal[a.idx];
    ter.DE <- exp(DE.proposal[tau.idx]); 
    bias.DE <-  exp(DE.proposal[w.idx])/(exp(DE.proposal[w.idx])+1);
    xi.DE <- DE.proposal[xi.idx];
    kap.sq = DE.proposal[kap.sq.idx]
    tau.sq = DE.proposal[tau.sq.idx]
    
    DE.list <- list(rt = rt, resp = resp, neural = neural, k = k, q = q, p = k+q,
                        Lam = Lam.DE, Phi = Phi.DE, Psi = Psi.DE, 
                        mu = mu.DE, a = a.DE, ter = ter.DE, bias = bias.DE, xi = xi.DE,
                        kap.sq = kap.sq, tau.sq = tau.sq);
    
    DE.log.lik.list <- fa.nddm(DE.list)
    DE.log.prior <- fa.nddm.prior(DE.list)
    
    for(j in 1:2){
      if(j == 1){
        propset.ddm[i] <- DE.log.lik.list$ddm.log.lik; # likelihood of migrated theta
        currentset.ddm[i] <- theta.ddm.lik[lnum2[i]]   # likelihood of current theta
        propw.ddm[i] = propset.ddm[i]
        currw.ddm[i] = currentset.sem[i]
        propw.prior.ddm[i] <- DE.log.prior$ddm.log.prior; # not used
        currw.prior.ddm[i] <- theta.ddm.prior[lnum2[i]]   # not used
      } # end of j==1
      
      if(j == 2){
        propset.sem[i] <- DE.log.lik.list$sem.log.lik;
        currentset.sem[i] <- theta.sem.lik[lnum2[i]]
        propw.sem[i] = propset.sem[i]
        currw.sem[i] = currentset.sem[i]
        propw.prior.sem[i] <- DE.log.prior$sem.log.prior;
        currw.prior.sem[i] <- theta.sem.prior[lnum2[i]]
      }
    }
    
  } # end of i
  
  
  for(j in 1:2){  # j = 1, ddm parameters, j = 2, sem (FA) parameters
     
    if(j == 1){
      if(runif(1) < exp(propw.ddm[lnum1]+propw.prior.ddm[lnum1] - (currw.ddm[1]+currw.prior.ddm[1]))){
        theta[ddm.idx, lnum2[1]] <- thetaset[ddm.idx, lnum1]							# swap the first with the last (creating a circle)
        theta.ddm.lik[lnum2[1]] <- propset.ddm[lnum1]
        theta.ddm.prior[lnum2[1]] <- propw.prior.ddm[lnum1]
      }
      if(lnum1!=1){											
        for(i in 1:(lnum1-1)){		
          if(runif(1) < exp(propw.ddm[i]+propw.prior.ddm[i] - (currw.ddm[i+1]+currw.prior.ddm[i+1]) )){
            theta[ddm.idx, lnum2[i+1]] <- thetaset[ddm.idx, i]							# swap the first with the last (creating a circle)
            theta.ddm.lik[lnum2[i+1]] <- propset.ddm[i]
            theta.ddm.prior[lnum2[i+1]] <- propw.prior.ddm[i]
          }}}
    }
    if(j == 2){
      if(runif(1) < exp( (propw.sem[lnum1]+propw.prior.sem[lnum1]) - (currw.sem[1]+currw.prior.sem[1]))){
        theta[las.idx, lnum2[1]] <- thetaset[las.idx, lnum1]							# swap the first with the last (creating a circle)
        theta.sem.lik[lnum2[1]] <- propset.sem[lnum1]
        theta.sem.prior[lnum2[1]] <- propw.prior.sem[lnum1]
      }
      if(lnum1!=1){											
        for(i in 1:(lnum1-1)){		
          if(runif(1) < exp( (propw.sem[i]+propw.prior.sem[i]) - (currw.sem[i+1]+currw.prior.sem[i+1]))){
            theta[las.idx, lnum2[i+1]] <- thetaset[las.idx, i]							# swap the first with the last (creating a circle)
            theta.sem.lik[lnum2[i+1]] <- propset.sem[i]
            theta.sem.prior[lnum2[i+1]] <- propw.prior.sem[i]
          }}}
    }
  }
  
  list(theta = theta, theta.ddm.lik = theta.ddm.lik, theta.sem.lik = theta.sem.lik, theta.ddm.prior = theta.ddm.prior, theta.sem.prior = theta.sem.prior)
}


