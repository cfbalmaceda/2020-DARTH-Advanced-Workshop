

###UKPDS equation time to event

#Set working directory
setwd("/Users/cfbalmac/Desktop/Github/2020-DARTH-Advanced-Workshop/CB_project/data")
synthetic_data <- read.csv('synthetic_sample_data.csv', header = T, sep = ",")
ukpds_param <- read.csv('ukpds_param_eq.csv', header = T, sep = ";")
nn <- nrow(synthetic_data)

param2 <- ukpds_param$rho
ukpds_param <- dplyr::select(ukpds_param, -rho)
dist_cq <- as.character(ukpds_param$distribution)
ukpds_param <- dplyr::select(ukpds_param, -distribution)
cons <- ukpds_param$cons
cons <- as.matrix(cons)
ukpds_param <- dplyr::select(ukpds_param, -cons)
conseq <- as.character(ukpds_param$X)
ukpds_param <- dplyr::select(ukpds_param, -X)
nconsq <- nrow(ukpds_param)
ukpds_param <- as.matrix(ukpds_param)

#gen variable to estimate every lamda
age <- synthetic_data$age
afro <- rep(0,nn)
age_dx <- synthetic_data$age_dx
gender <- synthetic_data$gender
indian <- rep(0,nn)
atfib <- rep(0,nn)
bmi <- synthetic_data$bmi
bmi_cat1 <- ifelse(bmi<=18.5,1,0)
bmi_cat3 <- ifelse(bmi>=25,1,0)
egfr <- rgamma(nn,((110)^2)/((20)^2),(110)/((20)^2))
egfr <- egfr/10
egfrL60 <- ifelse(egfr<60,egfr,0)
egfrM60 <- ifelse(egfr>=60,egfr,0)
haem <- (synthetic_data$haem)/10
hba1c <- synthetic_data$hba1c
hdl <- synthetic_data$hdl
hr <- synthetic_data$hr
ldl <- synthetic_data$ldl
ldl35 <- ifelse(ldl>35,ldl,0)
mmalb <- rep(0,nn)
pvd <- rep(0,nn)
sbp <- synthetic_data$sbp
smoke <- synthetic_data$smoke
wbc <- rnorm(nn,6.6,1)
amp_event <- rep(0,nn)    #This variable should be a tracker or something like that
amp_hist <- rep(0,nn)
blind_hist <- rep(0,nn)
chf_hist <- rep(0,nn)
ihd_event <- rep(0,nn)    #This variable should be a tracker or something like that
ihd_hist <- rep(0,nn)
mi_event <- rep(0,nn)     #This variable should be a tracker or something like that
mi_hist <- rep(0,nn)
renal_event <- rep(0,nn)  #This variable should be a tracker or something like that
renal_hist <- rep(0,nn)
stroke_event <- rep(0,nn) #This variable should be a tracker or something like that
stroke_hist <- rep(0,nn)
ulcer_hist <- rep(0,nn)

##extra variables
years_dx <- synthetic_data$years_dx
height <- synthetic_data$height
weight <- synthetic_data$weight


#Create matrix for individual data
indiv_data <- c(age,afro,age_dx,gender,indian,atfib,bmi,bmi_cat1,bmi_cat3,egfr,egfrL60,egfrM60,haem,
              hba1c,hdl,hr,ldl,ldl35,mmalb,pvd,sbp,smoke,wbc,amp_event,amp_hist,blind_hist,chf_hist,
              ihd_event,ihd_hist,mi_event,mi_hist,renal_event,renal_hist,stroke_event,stroke_hist,
              ulcer_hist)
indiv_data <- matrix(indiv_data,nn,36)

###Its not working###
lamda <- matrix(c(rep(NA,nconsq),param2,dist_cq),nconsq,3)
for (k in 1:nconsq) {
  for(j in 1:nn ){  
    lamda[j,] <- cons[k,]+sum(indiv_data[j,] %*% ukpds_param[k,])
}
}
lamda2 <- data.frame(lamda)
lamda2$X1 <- is.numeric(lamda2$X1)

#####################################################################################
#Create parameters array (lambda, shape and distribution for every patients)

###I realize that i need an array where row should be for patients data, the columns for the parameters and the matrix for the consequences. I dont know if i really need this or im lost!


row.names <- "row"
column.names <- c("lamda","shape","distribution")
matrix.names <- c("CHF","IHD","MI1_male","MI1_female","MI2","STROKE1","STROKE2","BLINDNESS","ULCER",
                  "AMPUTATION1 w/ULCER","AMPUTATION w/o ULCER","AMPUTATION2","RENAL FAILURE",
                  "DEATH w/NO EVENTS","DEATH 1stYEAR 1stEVENT","DEATH over 1EVENT","DEATH 2+EVENTS")
lamda <- array(c(rep(NA,nn),param2,dist_cq),dim = c(nn,3,nconsq))
  for (k in 1:nconsq) {
    for(j in 1:nn) {
    lamda[j,1,k] <- cons[k,]+sum(indiv_data[i,k] %*% ukpds_param[k,])
    }
  }
lamda2 <- data.frame(lamda)
lamda2$X1 <- is.numeric(lamda2$X1)


#####################################################################################

#########Function for time to event#################    I dont think that its works. 
##Time to event for every consequence
#param1 correspond to lamda[,1] - It is the scale parameter
#param2 correspond to lamda[,2] - It is the shape parameter
#distribution correspond to lamda[,3] - It is the distribution used to estimate the value for every consequence
timeto_event <- function(param1, param2, distribution) {
  timeto_event <- matrix(c(conseq,rep(NA,nconsq)),nconsq,2)
  for(i in 1:length(nconsq)){
    if(distribution=="exponential"){
      timeto_event[i,] <- (1/exp(exp(param1[i,]))*ln(2))
    } else if(distribution=="weibull"){
      timeto_event[i,] <- ((1/exp(exp(param1[i,])))*ln(2))^(1/exp(param2[i,]))
    } else if(distribution=="logistic"){
      timeto_event[i,] <- (1/0.5/(1-exp(-param1[i,])/(1+exp(-param1[i,])))*ln(2))
    } else 
      timeto_event[i,] <- (ln(((param1[i,]/param2[i,])*ln(2))+1)/param1[i,]) #Gompertz distribution
  }
  return(timeto_event=timeto_event)
}

























