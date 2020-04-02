#########################################################################################################
##
## ESPD Function
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
## This code illustrates how the ESPD approach, as defined in the paper, can be implemented in a general
## function in R Statistical Software. Depending on specific case study needs, modifications to the 
## code may be required. Separate code is available in the 'function bootESPD.R' file, including a 
## sophisticated error handling, which is useful when repeatedly applying to bootstrap samples in a
## probabilistic sensitivity analysis, for example.
##
## The functions uses the 'fitdist' function, which is part of the 'fitdistrplus' package, to fit 
## parametrical distributions.
##
## Function inputs
## - data     matrix or data.frame including two columns:
##            > time:   containing the time observations for the individuals
##            > event:  containing the event corresponding to 'time' observations
## - events   vector containing the events that should be included in the analysis
## - dist     the distribution to be fitted to the time-to-event data, as supported by the
##            'fitdist' function
##
## Function outputs
## - the function returns the parameter values according to the distribution provided in 'dist'
## - the parameter values are returned in a matrix for which the rows and columns are named
##   according the distribution parameter names and the event names provided in 'events'
##
#########################################################################################################


## Function ----
ESPD <- function(data, events, dist) {
  
  
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
  if(!is.character(dist)) stop("The 'dist' input should be a character!");                                    # The value provided to 'dist' should be a character
    
  # Load the 'fitdistrplus' package
  library(fitdistrplus);
  
  # Convert 'data' object into data.frame if necessary  
  if(is.matrix(data)) {
    data <- as.data.frame(data);
  }
  
  
  ## Estimate the parameter values
  
  # Loop over all events
  parameters <- sapply(events, function(e) {
    
    # Calculate the event incidence, i.e. probability
    prob <- sum(data$event==e)/sum(data$event %in% events);
      
    # Estimate the distribution parameters
    fit <- fitdist(data$time[data$event==e], dist)$estimate;
    
    # Return parameter values
    c(prob=prob, fit);
  
  }) # sapply
    
  
  ## Return the event-specific probabilities and distribution parameters
  return(parameters)

}