#########################################################################################################
##
## Illustration
##
## Supplementary to the manuschipt:
##
##        Comparing strategies for modeling competing risks in discrete event simulations: 
##        a simulation study and illustration in colorectal cancer
##
##        by K. Degeling,  H. Koffijberg, M.D. Franken, M. Koopman, M.J. IJzerman
##
##        DOI: https://doi.org/10.1177/0272989X18814770
## 
## This code illustrates how all functions defined in the provided R scripts can be used to implement
## individual patient time-to-event data in a discrete event simulation. The illustration includes both 
## examples on data analysis for the base-case analysis, as well as an example of performing a 
## probabilistic sensitivity analysis by bootstrapping the dataset to reflect parameter uncertainty in 
## the estimated time-to-event related proabilities and distribution parameters. The data used is an 
## extension of the exemplary data presented in the manuscript. This code also illustrates how the 
## results of the functions can be used to perform discrete event simulations.
##
#########################################################################################################

## Initialization ----

# Clear the workspace
rm(list=ls());

# Set the working directory to the folder which contains all R files
setwd("...SET DIRECTORY INCLUDING ALL PROVIDED R FILES...");

# Load all functions
source("function ESD.R");
source("function ESPD.R");
source("function UDR.R");
source("function MDR.R");

# Install the required Packages
install.packages("fitdistrplus");   # Package for fitting unimodal distributions
install.packages("mixtools");       # Package for fitting multimodal mixture distributions

# Load the exemplary data object 'IPD' 
load("data IPD.RData");


## Base-case analysis: ESD Modeling Approach ----

## Data analysis

# Call the 'ESD' function
ESD.in <- ESD(data=IPD, events=c("Progression","Death"), dist="weibull");

# Assess the data analysis results
ESD.in;

## Simulation

# Number of hypothetical patients
n.sim <- 1000;

# Initialize hypothetical patients
ESD.sim <- data.frame(id=c(1:n.sim));

# Draw times to each event
ESD.sim$time.to.death <- rweibull(n=n.sim, shape=ESD.in["shape","Death"], scale=ESD.in["scale","Death"]);
ESD.sim$time.to.progression <- rweibull(n=n.sim, shape=ESD.in["shape","Progression"], scale=ESD.in["scale","Progression"]);

# Select the event with the lowest drawn time
ESD.sim$time <- apply(cbind(ESD.sim$time.to.death, ESD.sim$time.to.progression), 1, min);
ESD.sim$event <- c("death","progressoin")[apply(cbind(ESD.sim$time.to.death, ESD.sim$time.to.progression), 1, which.min)];

# View results
View(ESD.sim);



## Base-case analysis: ESPD Modeling Approach ----

## Data analysis

# Call the 'ESPD' function
ESPD.in <- ESPD(data=IPD, events=c("Progression","Death"), dist="weibull");

# Assess the data analysis results
ESPD.in;

## Simulation

# Number of hypothetical patients
n.sim <- 1000;

# Initialize hypothetical patients
ESPD.sim <- data.frame(id=c(1:n.sim));

# Select the event to occur according to the event-specific probabilities
ESPD.sim$event <- ifelse(runif(nrow(ESPD.sim))<ESPD.in["prob","Progression"], "progression", "death");

# Draw a time to the corresponding event
ESPD.sim$time <- ifelse(ESPD.sim$event=="progression", rweibull(n=n.sim, shape=ESPD.in["shape","Progression"], scale=ESPD.in["scale","Progression"]),
                        rweibull(n=n.sim, shape=ESPD.in["shape","Death"], scale=ESPD.in["scale","Death"]));

# View results
View(ESPD.sim);




## Base-case analysis: UDR Modeling Approach ----

## Data analysis

# Call the 'UDR' function
UDR.in <- UDR(data=IPD, events=c("Progression","Death"), dist="weibull");

# Assess the data analysis results, i.e. the regression model and time-to-event distribution parameters
UDR.in;

## Simulation

# Number of hypothetical patients
n.sim <- 1000;

# Initialize hypothetical patients
UDR.sim <- data.frame(id=c(1:n.sim));

# Draw time to event
UDR.sim$time <- rweibull(n=n.sim, shape=UDR.in$tte.dist["shape"], scale=UDR.in$tte.dist["scale"]);

# Select the corresponding event
UDR.sim$event <- ifelse(runif(nrow(UDR.sim))<predict(UDR.in$regression.model, newdata=UDR.sim, type="response"),"progression","death");
  
# View results
View(UDR.sim);




## Base-case analysis: MDR Modeling Approach ----

## Data analysis

# Call the 'MDR' function
MDR.in <- MDR(data=IPD, events=c("Progression","Death"), maxit=500, maxit.survreg=500);

# Assess the data analysis results, i.e. the regression model and time-to-event mixture distribution parameters
# Disabled due to the large output for the mixture model
#MDR.in;

## Simulation

# Number of hypothetical patients
n.sim <- 1000;

# Initialize hypothetical patients
MDR.sim <- data.frame(id=c(1:n.sim));

# Draw time to event
MDR.sim$time <- rweibullmix(n=n.sim, lambda=MDR.in$tte.dist$lambda, shape=MDR.in$tte.dist$shape, scale=MDR.in$tte.dist$scale);

# Select the corresponding event
MDR.sim$event <- ifelse(runif(nrow(MDR.sim))<predict(MDR.in$regression.model, newdata=MDR.sim, type="response"),"progression","death");

# View results
View(MDR.sim);






## Probabilistic sensitivity analysis: ESPD Modeling Approach ----

# The ESPD modeling was chosen to illustrate the application of the approaches in a probabilistic
# sensitivity analysis, rather than the MDR approach, as the MDR approach requires substantial
# computer resources for estimating the mixture distribution. Using the MDR approach in this 
# illustration would, therefore, result in long waiting times.

## Data analysis

# Number of (probabilistic sensitivity analysis) runs to be performed
n.bootstraps <- 100;

# Set seed for reproducibility
set.seed(34512);

# Perform data analysis and put the estimated parameter values in a list for each run
ESPD.psa.in <- lapply(1:n.bootstraps, function(run) {
  
  # Draw a bootstrap sample from the original data
  bootstrap.sample <- IPD[sample(nrow(IPD), size=nrow(IPD), replace=T),];
  
  # Analyze the bootstrap sample
  ESPD(data=bootstrap.sample, events=c("Progression","Death"), dist="weibull");
  
}) # sapply


## Simulation

# Number of hypothetical patients
n.sim <- 1000;

# Perform the simulation for all (probabilistic sensitivity analysis) runs
ESPD.psa.sim <- sapply(ESPD.psa.in, function(parameters) {
  
  # Initialize hypothetical patients
  sim <- data.frame(id=c(1:n.sim));
  
  # Select the event to occur according to the event-specific probabilities
  sim$event <- ifelse(runif(nrow(sim))<parameters["prob","Progression"], "progression", "death");
  
  # Draw a time to the corresponding event
  sim$time <- ifelse(sim$event=="progression", rweibull(n=n.sim, shape=parameters["shape","Progression"], scale=parameters["scale","Progression"]),
                     rweibull(n=n.sim, shape=parameters["shape","Death"], scale=parameters["scale","Death"]));
  
  # Return the outcome of interest
  c(mean.tte=mean(sim$time));
  
}) # sapply

# Summarize and present the result of interest
summary(ESPD.psa.sim);