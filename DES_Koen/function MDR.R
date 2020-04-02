#########################################################################################################
##
## MDR Function
##
## Supplementary to the manuschipt:
##
##        Comparing strategies for modeling competing risks in discrete event simulations: 
##        a simulation study and illustration in colorectal cancer
##
##        by K. Degeling,  H. Koffijberg, M.D. Franken, M. Koopman, M.J. IJzerman
##
##        DOI: 
## 
## This code illustrates how the MDR approach, as defined in the paper, can be implemented in a general
## function in R Statistical Software. Depending on specific case study needs, modifications to the 
## code may be required.
##
## The functions uses the 'weibullRMM_SEM' function, which is part of the 'mixtools' package, to fit 
## a Weibull mixture distribution. The starting values provided to this function are estimated using 
## the 'fitdist' function that is part of the 'fitdistrplus' package.
##
## Function inputs
## - data             matrix or data.frame including two columns:
##                      > time:   containing the time observations for the individuals
##                      > event:  containing the event corresponding to 'time' observations
## - events           vector containing the events that should be included in the analysis
## - maxit            maximum number of iterations for the weibullRMM_S function (see documentation)
## - maxit.survreg    maximum number of iterations for the weibullRMM_S function (see documentation)
##
## Function outputs
## - the function returns the fitted logistic regression model and mixture distribution for in two list
##   items: 'regression.model' and 'tte.dist', respectively
##
#########################################################################################################


## Function ----
MDR <- function(data, events, maxit=200, maxit.survreg=200) {
  
  
  ## Initialization
  
  # General checks for the provided input
  if(!(is.data.frame(data) | is.matrix(data))) stop("The 'data' object should be a matrix or a data.frame!"); # The 'data' object is matrix or data.frame
  if(!("time" %in% colnames(data))) stop("There is no 'time' column in 'data'!");                             # A 'time' column should be present in the 'data' object
  if(!(is.numeric(data[,"time"]))) stop("The provided values in the 'time' column are not numerical!");       # Values in the 'time' column should be numerical
  if(any(is.na(data[,"time"]))) stop("NA values are present in the 'time' column!");                          # No missing values are allowed in the 'time' column
  if(any(data[,"time"]<=0)) stop("Negative values or zeros are present in the 'time' column!");               # No non-positive values are allowed in the 'time' column    
  if(!("event" %in% colnames(data))) stop("There is no 'event' column in 'data'!");                           # An 'event' column should be present in the 'data' object
  if(any(is.na(data[,"event"]))) stop("NA values are present in the 'event' column!");                        # No missing values are allowed in the 'event' column
  if(!all(events %in% data[,"event"])) stop("Not all events are present in the data!");                       # All requested events in 'events' should be present in the 'event' column

  # Load the 'fitdistrplus' package
  library(fitdistrplus);
  
  # Load the 'mixtools' package
  library(mixtools);
  
  # Convert 'data' object into data.frame if necessary  
  if(is.matrix(data)) {
    data <- as.data.frame(data);
  }
  
  # Make sure events are listed as factor for the analysis
  if(!is.factor(data$event)) {
    data$event <- as.factor(data$event);
  }
  
  # Make sure only the times and events of the included events are present
  data <- subset(data, event %in% events);
  
  
  ## Estimate the parameter values
  
  # Estimate the logistic regression model
  regression.model <- glm(event~time, family=binomial(link='logit'), data=data);
      
  # Estimate the starting values for the algorithm estimating the mixtools package
  start.values <- sapply(events, function(e) fitdist(data$time[data$event==e], "weibull")$estimate);
  
  # Estimate the distribution parameters
  tte.dist <- weibullRMM_SEM(x=data$time, shape=start.values["shape",], scale=start.values["scale",],
                             maxit=maxit, maxit.survreg=maxit.survreg);
    
  
  ## Return parameter values
  return(list(regression.model=regression.model, tte.dist=tte.dist));

}