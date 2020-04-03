
###Function to estimate time to event

TTE <- function(data){
  
  # Convert 'data' object into data.frame if necessary  
  if(!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  with(data,{
  TTE_data <- vector("numeric",nn) 
    for(j in 1:length(TTE_data)){
     # browser()
      if(dist[j]=="weibull"){
      TTE_data[j] <- ((1 / exp((lambda[j])))*log(2))^(1/exp(rho[j]))
      } else if(dist[j]=="gompertz") {
      TTE_data[j] <- log(((1 / exp((lambda[j])))*log(2))^(1/exp(rho[j])))
    }
    }
  return(TTE_data=TTE_data)
  
  })
}

TTE_CHF <- as.data.frame(TTE(data =param_CHF))
TTE_IHD <- as.data.frame(TTE(data =param_IHD))
TTE_death <- as.data.frame(TTE(data =param_death))

