#-----------------------------------------------------
#----- Title: Lasso FA NDDM (simulation or model fitting)
#----- Author: Inhan Kang, Woojong Yi, and Dr. Brandon Turner
#----- Model-based Cognitive Neuroscience Lab
#----- The Ohio State University
#-----------------------------------------------------
source("functions.R")
sfExportAll()

#-----------------------------------------------


#------ data set
#-- behavioral data
rt <- data$rt
n.trials <- length(data$rt)  # the number of trials
resp.old <- data$resp
resp <- rep("upper", n.trials); resp[resp.old == 0] <- "lower";
resp <- as.factor(resp);
dec.dat <- data.frame(rt = rt, resp = resp.old)

#-- neural data
neural <- data$neural



#----- define some indices
# note) notations here (k, q, p) are different from those in the article ( k -> q, q -> p, p -> q+p in the article)
k <- 3   # the number of behavioral factors
q <- ncol(data$neural)
p <- k+q


# indicies that indicate where the parameters are in "theta", a vector of all parameters aggregated.
Lam.idx <- 1:(q*k)
Phi.idx <- (max(Lam.idx)+1):(max(Lam.idx)+k+k*(k-1)/2);
Psi.idx <- (max(Phi.idx)+1):(max(Phi.idx) + q);
mu.idx <- (max(Psi.idx)+1):(max(Psi.idx) + q + k);
sem.idx <- c(Lam.idx, Phi.idx, Psi.idx, mu.idx);

a.idx <- max(sem.idx)+1;
tau.idx <- ter.idx <- (max(a.idx)+1):(max(a.idx)+n.trials);
w.idx <- bias.idx <- tau.idx + n.trials; 
xi.idx <- w.idx + n.trials;
ddm.idx <- c(a.idx, tau.idx, w.idx, xi.idx)

kap.sq.idx <- max(xi.idx)+1
tau.sq.idx <- (max(kap.sq.idx)+1):(max(kap.sq.idx) + k*q)
reg.idx <- c(kap.sq.idx, tau.sq.idx)

npar <- length(sem.idx) + length(ddm.idx) + length(reg.idx) # the number of parameters

# the list of indicies
idx.list <- list(n.chains = n.chains, n.trials = n.trials, npar = npar, sem.idx = sem.idx, ddm.idx = ddm.idx, reg.idx = reg.idx,
                 Lam.idx = Lam.idx, Phi.idx = Phi.idx, Psi.idx = Psi.idx, mu.idx = mu.idx, 
                 a.idx = a.idx, tau.idx = tau.idx, w.idx = w.idx, xi.idx = xi.idx, kap.sq.idx = kap.sq.idx, tau.sq.idx = tau.sq.idx,
                 k = k, q = q, p = p);

data.list <- list(rt = rt, resp = resp, neural = neural)

# the frequency of printing trace out
trace.out <- 10
lik.out <- 100


#----- storage
sample.array <- array(data = NA, dim = c(n.iter, npar, n.chains));
de.ac.ddm1 <- de.ac.ddm2  <- rep(0, n.chains);  # storage for #(accepted)



#----- initialize chains
theta <- matrix(NA, nrow = npar, ncol = n.chains);
for(l in 1:n.chains) {
  
  #----- DDM parameters
  a0 <- runif(1, 1, 3); 
  ter0 <- runif(n.trials, rt/3, rt*(2/3)); # should be smaller than RT
  tau0 <- log(ter0);
  w0 <- rnorm(n.trials, 0, 0.3)
  xi0 <- rnorm(n.trials, 1.2, 0.3)
  ddm0 <- c(a0, tau0, w0, xi0)
  
  #----- FA parameters
  lam0 <- runif(length(Lam.idx), 0.2, 0.8);
  phi0 <- c(rep(1, k) + rnorm(k, 0, 0.1), rep(0, k*(k-1)/2) + rnorm(k*(k-1)/2, 0, 0.05))
  psi0 <- runif(length(Psi.idx), 0.01, 0.3);
  
  mu0 <- c(mean(tau0), mean(w0), mean(xi0), apply(neural, 2, mean) + rnorm(q, 0, 0.01));
  fa0 <- c(lam0, phi0, psi0, mu0);
  
  #----- Regularization parameters
  kap.sq <- runif(1, 0, 1); #rgamma(1, shape = 10, scale = 1/4)
  tau.sq <- matrix(0, q, k);
  for(m in 1:q) {tau.sq[m, ] <- rexp(k, rate = kap.sq/2)}
  reg0 <- c(kap.sq, as.vector(t(tau.sq)))
  
  #----- a matrix of all parameters used to iterate chains
  theta[,l] <- c(fa0, ddm0, reg0); 
}


#-- check if initial values are okay by calculating their lik and prior
theta.ddm.lik <- theta.sem.lik <- theta.ddm.prior <- theta.sem.prior <- numeric(n.chains); # these are the storage for the log lik/prior values of the current sample

for(j in 1:n.chains){
  Lam.theta <- theta[Lam.idx, j]; 
  Phi.theta <- theta[Phi.idx, j]; 
  Psi.theta <- theta[Psi.idx, j]; 
  mu.theta <- theta[mu.idx, j];
  a.theta <- theta[a.idx, j]; 
  
  ter.theta <- exp(theta[tau.idx, j]); 
  bias.theta <-  exp(theta[w.idx, j])/(exp(theta[w.idx, j])+1); 
  xi.theta <- theta[xi.idx, j];
  
  kap.sq.theta <- theta[kap.sq.idx, j]
  tau.sq.theta <- theta[tau.sq.idx, j]
  
  theta.ddm.list <- list(rt = rt, resp = resp, neural = neural, k = k, q = q, p = k+q,
                         Lam = Lam.theta, Phi = Phi.theta, Psi = Psi.theta, mu = mu.theta, 
                         a = a.theta, ter = ter.theta, bias = bias.theta, xi = xi.theta,
                         kap.sq = kap.sq.theta, tau.sq = tau.sq.theta);
  
  # Calculate (initial) log likelihood and prior values 
  theta.log.lik.list <- fa.nddm(theta.ddm.list);
  theta.log.prior <- fa.nddm.prior(theta.ddm.list);
  
  theta.ddm.lik[j] <- theta.log.lik.list$ddm.log.lik;
  theta.sem.lik[j] <- theta.log.lik.list$sem.log.lik;
  theta.ddm.prior[j] <- theta.log.prior$ddm.log.prior;
  theta.sem.prior[j] <- theta.log.prior$sem.log.prior;
}





#------------ Run the Bayesian sampling procedure

sfExportAll(except = "sample.array")
start.time <- Sys.time();
iter.start.time <- Sys.time()
for(i in 1:n.iter){
  
  DEMC.result <- matrix(unlist(sfLapply(1:n.chains, DEMC.FA, theta = theta, theta.ddm.lik = theta.ddm.lik, theta.sem.lik = theta.sem.lik, 
                                        theta.ddm.prior = theta.ddm.prior, theta.sem.prior = theta.sem.prior,
                                        idx.list = idx.list, data.list = data.list, de.ac.ddm1 = de.ac.ddm1, de.ac.ddm2 = de.ac.ddm2)), nrow = npar+6, ncol = n.chains)
  
  # update the values
  theta <- DEMC.result[-c(npar+c(1:7)), ]
  theta.ddm.lik <- DEMC.result[npar+1, ]
  theta.sem.lik <- DEMC.result[npar+2, ]
  theta.ddm.prior <- DEMC.result[npar+3, ]
  theta.sem.prior <- DEMC.result[npar+4, ]
  de.ac.ddm1 <- DEMC.result[npar+5, ]
  de.ac.ddm2 <- DEMC.result[npar+6, ]
  
  sample.array[i , , ] <- theta;
  
   
   #----- migration
   if(i < migrate.duration){
     if(runif(1) < migrate.prob){
      out <- fa.nddm.migrate(theta, theta.ddm.lik, theta.sem.lik, theta.ddm.prior, theta.sem.prior, idx.list, data.list)
      theta <- out$theta
      theta.ddm.lik <- out$theta.ddm.lik
      theta.sem.lik <- out$theta.sem.lik
      theta.ddm.prior <- out$theta.ddm.prior
      theta.sem.prior <- out$theta.sem.prior
    }}
   
  
  #----- print out the elapsed time and acceptance ratios every -trace.out- iteration
  if(i %in% c(1, seq(trace.out, n.iter, by = trace.out))){ 
      iter.end.time <- Sys.time();
      iter.elapsed <- iter.end.time - iter.start.time;
      sfCat("     Iteration ", i, " (sampling) finished. Elapsed: ", round(iter.elapsed, 2), " ", units(iter.elapsed), ".\n", sep="") 
      sfCat("       DDM1 acceptance ratio: ", round(de.ac.ddm1/(i), 3), "\n")   # a
      sfCat("       DDM2 acceptance ratio: ", round(de.ac.ddm2/(i), 3), "\n")   # single-trial parameters
      # note) FA parameters are updated by their conditional posteriors and so their acceptance ratios are 1
    iter.start.time <- Sys.time()
  }
  
  #----- print out current log likelihood and prior values 
  if(i %% lik.out == 0 ){
    if(i == n.burnin) sfCat("     Burn-in finished.\n", sep="") 
    sfCat("       DDM log likelihood: ", round(theta.ddm.lik, 3), "\n")
    sfCat("       SEM log likelihood: ", round(theta.sem.lik, 3), "\n")
    sfCat("       DDM log prior: ", round(theta.ddm.prior, 3), "\n")
    sfCat("       SEM log prior: ", round(theta.sem.prior, 3), "\n")
    #save.image(save.name)
  }    
} # END OF i-loop


end.time <- Sys.time();
(elapsed <- end.time - start.time);






