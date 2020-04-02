








setwd("/Users/cfbalmac/Desktop/Github/2020-DARTH-Advanced-Workshop/CB_project/data/");

synthetic_data <- read.csv('synthetic_sample_data.csv', header = T, sep = ",")
ukpds_param <- read.csv('ukpds_param_eq.csv', header = T, sep = ";")
nn <- nrow(synthetic_data)

param2 <- as.matrix(ukpds_param[c(1:2,14),"rho"])     #Shape parameters
dist_cq <-as.matrix(as.character(ukpds_param[c(1:2,14),"distribution"]))
cons <- as.matrix(ukpds_param[c(1:2,14),"cons"])
conseq <-  c("CHF","IHD","Death")  
nconsq <- nrow(cons)
ukpds_param <- dplyr::select(ukpds_param, -c(rho,cons,X,distribution))
ukpds_param_v2 <- ukpds_param[c(1:2,14),]
ukpds_param_v2 <- as.matrix(ukpds_param_v2)

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

###Parameter for consequences
colnames <- c("lamda","rho","dist")

param_CHF <- matrix(rep(c(NA,param2[1,],dist_cq[1,]),nn,each=nn),nn,3)
for(i in 1:nn){
  param_CHF[i,1] <- cons[1,]+sum(indiv_data[i,] %*% ukpds_param_v2[1,])
}
param_IHD <- matrix(rep(c(NA,param2[2,],dist_cq[2,]),nn,each=nn),nn,3)
for(i in 1:nn){
  param_IHD[i,1] <- cons[1,]+sum(indiv_data[i,] %*% ukpds_param_v2[2,])
}
param_death <- matrix(rep(c(NA,param2[3,],dist_cq[3,]),nn,each=nn),nn,3)
for(i in 1:nn){
  param_death[i,1] <- cons[1,]+sum(indiv_data[i,] %*% ukpds_param_v2[3,])
}
colnames(param_CHF) <- c("lamda","rho","dist")
colnames(param_IHD) <- c("lamda","rho","dist")
colnames(param_death) <- c("lamda","rho","dist")

param_CHF <- as.data.frame(param_CHF)

###Function to estimate time to event

TTE <- function(data,lamda,rho,dist){
  
  # Convert 'data' object into data.frame if necessary  
  if(is.matrix(data)) {
    data <- as.data.frame(data);
  }
  
  TTE_data <- data.frame(rep(NA,nn)) 
    for(j in 1:nrow(TTE_data)){
      if(dist[j,]=="weibull"){
      TTE_data[j,] <- ((1/exp(exp(lamda[j,])))*ln(2))^(1/exp(rho[j,]))
      } else if(dist[j,]=="gompertz") {
      TTE_data[j,] <- (ln(((lamda[j,]/rho[j,])*ln(2))+1)/lamda[j,])
    }
    }
  return(TTE_data=TTE_data)
}

TTE_CHF <- TTE(param_CHF,lamda= param_CHF$lamda,param_CHF$rho,param_CHF$dist)



