#-----------------------------------------------------
#----- Title: Lasso FA NDDM (main)
#----- Author: Inhan Kang, Woojong Yi, and Dr. Brandon Turner
#----- Model-based Cognitive Neuroscience Lab
#----- The Ohio State University
#-----------------------------------------------------


rm(list = ls())
setwd("-----path-of-the-files-----")


save.name <- '-----name-of-the-result-file-----'   # file name '*.rdata': see line 86
load('-----load *.rdata -----')   
# the dataset should include a list variable 'data' which includes
# 1) rt (data$rt): a vector containing all RT observations across trials (length = number of observations).
# 2) resp (data$resp): a vector containing all responses (e.g., 1: correct, 0: incorrect) observations across trials (length = number of observations).
# 3) neural (data$neural): a matrix (nrow = number of observations, ncol = number of neural features) containing neural data.
# 4) etc.

# load("toy_data.rdata")
# names(data)


#----- require packages
req.pack <- function(pkg){ # require packages
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("MASS", "snowfall", "rjags", "RWiener", "mvtnorm", "truncnorm", "rtdists", "statmod", "LaplacesDemon")
req.pack(packages)
load.module("wiener")  # Wabersich & Vandekerckhove (2014), Behav Res., for wiener process functions => see the article 
set.seed(1)


#----- some simulation settings  (for Bayesian sampling)
# the number of burn-in period
n.burnin <- 2000;
# the number of sampling after burn-in
n.sampling <- 18000;
# the number of total iterations
n.iter <- n.burnin + n.sampling;
# the number of chains for DEMC
n.chains <- 12; 
# the index for the iterations after burn-in
take.idx <- (n.burnin + 1):n.iter;

# migration
migrate.duration <- round(n.burnin*.25)+1
migrate.prob <- 0.2



#----- snowfall setting & start multicore processing
cores <- n.chains;

sfInit(parallel = TRUE, cpus = cores);
sfLibrary(MASS);
sfLibrary(rjags);
sfLibrary(RWiener)
sfLibrary(mvtnorm)
sfLibrary(truncnorm)
sfLibrary(statmod)
sfLibrary(LaplacesDemon)
sfLibrary(snowfall)
load.module("wiener")

sfClusterSetupRNG()

#----- start simulation
start.time <- Sys.time();
cpu.start.time <- proc.time()[3];
source("sim.R") 
cpu.end.time <- proc.time()[3];
end.time <- Sys.time();
cpu.elapsed <- cpu.end.time - cpu.start.time;
elapsed <- end.time - start.time;
#----- end simulation
#----- stop multicore processing
sfStop()


#----- save the result
save.image(save.name)

