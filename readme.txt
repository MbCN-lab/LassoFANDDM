This Github includes R scripts and a toy data associated with the manuscript: Kang, I., Yi, W., & Turner, B. M. (Submitted). A Regularization Method for Linking Brain and Behavior.

- main.R: specifies settings for model-fitting (e.g., data file, number of iterations, number of chains, etc) and calls sim.R
- sim.R: imports functions.R and fits the Lasso FA NDDM with the settings specified in main.R
- functions.R: includes R functions to fit the model
- toy_data.rdata: a toy example data (based on Simulation 1 in the manuscript)

Users need to specify data file (*.rdata; see comments in main.R) and then fit the model in R, by 
	sources("~/main.R") 


Users only need to work on main.R and do not need to change anything in sim.R and functions.R except for the solutions to the sign-switching problem (Equation 9 in the manuscript). Users can find this in lines 284-286 in functions.R. The default is to monitor the factor loading of the first factor (corresponding to the non-decision time) on the first neural covariates (\lambda_{1,1}), that of the second factor (corresponding to the initial bias) on the second neural covariates (\lambda_{2,2}), and that of the third factor (corresponding to the drift rate) on the third neural covariate (\lambda_{3,3}). However, this should be set based on a priori knowledge or a preliminary analysis result (e.g., we chose ROIs 33, 59, and 57 for the non-decision time, initial bias, and drift rate, respectively, for our application with data under the accuracy-stressed condition; see Section 6: Application).