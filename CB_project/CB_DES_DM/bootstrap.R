##Bootstrap function

# Set seed for reproducibility
set.seed(1234);
bootstrap <- function(data,nIteration){
# Perform data analysis and put the estimated parameter values in a list for each run
   # Draw a bootstrap sample from the original data
    bootstrap.sample <- data[sample(nrow(data), size=(nIteration), replace=T),];
  # sapply
}




