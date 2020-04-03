
library(gems)
hf <- generateHazardMatrix(3)
hf[[1,2]] <- "Exponential"
hf[[1,3]] <- "Exponential"
hf[[2,3]] <- "Exponential"

par <- generateParameterMatrix(hf)
par[[1, 2]] <- list( 0.02)
par[[1, 3]] <- list( 0.05)
par[[2, 3]] <- list( 0.1)


cohortSize <- n_i
 cohort <- simulateCohort(transitionFunctions = hf,
                             parameters = par,
                            cohortSize = cohortSize,
                             to = n_t)
head(cohort)
post <- transitionProbabilities(cohort, times = seq(0,n_t, 1))

matplot(post@probabilities, type="l")
