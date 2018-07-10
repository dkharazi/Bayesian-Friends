## Import libraries
library(knitr)
library(readr)
library(rjags)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

## Import data
friends <- read_csv("~/Desktop/friends.csv", col_types = c(.default = col_integer()))
friends.df <- data.frame(friends)

## Set the seed
set.seed(252)

## Initialize some variables
isFriends <- friends.df$isFriends
school <- friends.df$school
isSameGender <- friends.df$isSameGender
nSchools <- length(unique(school))
nSurvey <- length(isFriends)

## Initialize objects for JAGS
dataList <- list("isFriends" = isFriends,
                 "isSameGender" = isSameGender,
                 "school" = school,
                 "nSchools" = nSchools,
                 "nSurvey" = nSurvey)

## List of parameters to be monitored  
parameters <- c("alpha", 
                "beta",
                "mu_alpha",
                "mu_beta",
                "sigma2_alpha",
                "sigma2_beta")

## Set initial values
initsValues <- list("alpha" = rep(0, nSchools), 
                    "beta" = rep(0, nSchools),
                    "mu_alpha" = 0,
                    "mu_beta" = 0,
                    "sigma2_alpha" = 1,
                    "sigma2_beta" = 1)

## Number of iteration for "tuning" 
adaptSteps <- 5000 

## Number of iterations for "burn-in" 
burnInSteps <- 5000   

## Number of chains to run
nChains <- 2          

## Total number of iterations to save
numSavedSteps <- 5000           

## Thinning to keep every iteration
thinSteps <- 1                  

## Iterations per chain
ITER <- ceiling( (numSavedSteps * thinSteps) / nChains )

## Create, initialize, and adapt the model
jagsModel <- jags.model("friendsmodel.txt", 
                        data = dataList, 
                        inits = initsValues, 
                        n.chains = nChains, 
                        n.adapt = adaptSteps )

## Burn-in the algorithm
update(jagsModel, 
       n.iter = burnInSteps)

## Run algorithm to get interations for inference
codaSamples <- coda.samples(jagsModel, 
                            variable.names = parameters, 
                            n.iter = ITER, 
                            thin = thinSteps)
                            
## Make a dataframe with the posterior samples
mcmcChainDF <- data.frame(as.matrix(codaSamples, 
                                    iters = T, 
                                    chains = T ) )

## Create a vector with the variable names
varNames <- names(mcmcChainDF)[3:(dim(mcmcChainDF)[2])]

## Number of variables
nVars <- length(varNames)
mcmcChainDF$CHAIN <- as.factor(mcmcChainDF$CHAIN)

## Construct trace plots
p <- list()
for(k in 1:nVars) {
  plot_frame <- mcmcChainDF
  plot_frame$dep_var <- mcmcChainDF[,varNames[k]]
  p[[k]] <- ggplot(plot_frame, aes(x = ITER, y = dep_var)) +
    geom_line(aes(color = CHAIN)) + 
    labs(y = varNames[k])
}

## Trace plots
do.call(grid.arrange, c(p, list("ncol" = 1)))

## Setup data frame with posterior distributions for alphas
alphaPostDFreshape <- melt(mcmcChainDF, id.vars = "ITER", measure.vars = c("alpha.1.", "alpha.2.", "alpha.3.",
                                                                           "alpha.4.", "alpha.5.", "alpha.6.",
                                                                           "alpha.7.", "alpha.8.", "mu_alpha"))

## Plot posterior distributions for alphas
ggplot(alphaPostDFreshape, aes(x = variable, y = value)) +
  geom_boxplot() +
  ylab("Posterior") +
  xlab("")
  
## Setup data frame with posterior distributions for betas
betaPostDFreshape <- melt(mcmcChainDF, id.vars = "ITER", measure.vars = c("beta.1.", "beta.2.", "beta.3.", 
                                                                          "beta.4.", "beta.5.", "beta.6.",
                                                                          "beta.7.", "beta.8.", "mu_beta"))

## Plot posterior distributions for betas
ggplot(betaPostDFreshape, aes(x = variable, y = value)) +
  geom_boxplot() +
  ylab("Posterior") +
  xlab("")
  
## Predictive probabilit that they're the same gender
mcmcChainDF$alphaPred <- rnorm(numSavedSteps, mcmcChainDF$mu_alpha, sqrt(mcmcChainDF$sigma2_alpha))
mcmcChainDF$betaPred <- rnorm(numSavedSteps, mcmcChainDF$mu_beta, sqrt(mcmcChainDF$sigma2_beta))

mcmcChainDF$predSame <- exp(mcmcChainDF$alphaPred) / (1+exp(mcmcChainDF$alphaPred))
mean(mcmcChainDF$predSame)

## Predictive probability that they're different genders
mcmcChainDF$predDiff <- exp(mcmcChainDF$alphaPred+mcmcChainDF$betaPred ) / (1+exp(mcmcChainDF$alphaPred+mcmcChainDF$betaPred))
mean(mcmcChainDF$predDiff)