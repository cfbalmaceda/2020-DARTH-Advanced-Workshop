###DES code DM model


rm(list = ls())      # clear memory (removes all the variables from the workspace)

# 00 Input dataset

source("input.R")

# 01 Load packages
library(survival)

# 02 Load functions
source("TTE_function.R")
source("bootstrap.R")

# 03 Input model parameters

# Strategy names
v_names_str <- c("Base Case")  

# Number of strategies
n_str <- length(v_names_str)

# Markov model parameters
v_n  <- c("NoEvent","CHF", "IHD", "Dead")    # state names

n_states  <- length(v_n)                # number of states
n_t  <- 60                              # number of cycles
n_i <- 10000



# Number of (simulated patients) runs to be performed
n.bootstraps <- n_i;

# Set seed for reproducibility
set.seed(1234);

# Perform data analysis and put the estimated parameter values in a list for each run
#Draw the values for all state transitions. We ll check later if they actuyally transitioned
t_NoEvent_CHF <- bootstrap(TTE_CHF,n.bootstraps)
t_NoEvent_IHD <- bootstrap(TTE_IHD,n.bootstraps)
t_NoEvent_Death <- bootstrap(TTE_death,n.bootstraps)


##Time to a second event
#Excludes negative differences (We dont want negative time values)
diff_IHD_CHF <- ifelse(t_NoEvent_IHD-t_NoEvent_CHF<0,0,t_NoEvent_IHD-t_NoEvent_CHF)
diff_Death_CHF <- ifelse(t_NoEvent_Death-t_NoEvent_CHF<0,0,t_NoEvent_Death-t_NoEvent_CHF)
diff_CHF_IHD <- ifelse(t_NoEvent_CHF-t_NoEvent_IHD<0,0,t_NoEvent_CHF-t_NoEvent_IHD)
diff_Death_IHD <- ifelse(t_NoEvent_Death-t_NoEvent_IHD<0,0,t_NoEvent_Death-t_NoEvent_IHD)

t_2ndEvent <- ifelse(diff_Death_CHF>diff_IHD_CHF,diff_IHD_CHF,ifelse(diff_Death_IHD>diff_CHF_IHD,
                                                                     diff_CHF_IHD,0))
#Time to Death
#Time from first event to death (No second event)
t_CHF_Death <- ifelse(t_2ndEvent!=0,0,diff_Death_CHF)
t_IHD_Death <- ifelse(t_2ndEvent!=0,0,diff_Death_IHD)
#Time from second event to death
t_2E_Death <- ifelse(t_2ndEvent!=0,t_NoEvent_Death-t_2ndEvent,0)


t_NoEvent <- cbind(t_NoEvent_CHF,t_NoEvent_IHD,t_NoEvent_Death)
next.state <- apply(t_NoEvent,1, which.min)
next.time  <- apply(t_NoEvent,1, min)

# look at the first "jump" and assign the corresponding state
t_1 <- ifelse(t_NoEvent_Death > t_NoEvent_IHD & t_NoEvent_CHF, t_NoEvent_Death ,
              ifelse(t_NoEvent_CHF>t_NoEvent_IHD,t_NoEvent_CHF,t_NoEvent_IHD))
m_1 <- ifelse(t_NoEvent_Death > t_NoEvent_IHD & t_NoEvent_CHF, "Dead", 
              ifelse(t_NoEvent_CHF>t_NoEvent_IHD,"CHF","IHD"))

#for those one consequence (CHF or IHD) draw a time of Death

# look for a  second "jump". if the patient already went to the event state in the previous jump then their t2 is the sum of t1 and t_Event_Death ot t_2ndEvent
t_2 <-  ifelse(m_1=="CHF",ifelse(t_2ndEvent>t_CHF_Death,t_1+t_CHF_Death,t_1+t_2ndEvent),
               ifelse(m_1=="IHD",ifelse(t_2ndEvent>t_IHD_Death,t_1+t_IHD_Death,t_1+t_2ndEvent),t_1))# the time at the second jump for the 
m_2  <- ifelse(m_1=="CHF",ifelse(t_2ndEvent>t_CHF_Death,"Death","IHD"),
               ifelse(m_1=="IHD",ifelse(t_2ndEvent>t_IHD_Death,"Death","CHF"),m_1))
  
 
# look for a  third "jump". if the patient already went to the event state in the previous jump then their t3 is the sum of t2 and t_2ndEvent_Death
t_3 <-  ifelse(m_2== "CHF", t_2 + t_2E_Death, ifelse(m_2=="IHD",t_2+t_2E_Death,t_2)) # the time at the second jup for the 
m_3  <- ifelse(m_2=="CHF", "Dead", ifelse(m_2=="IHD","Dead",m_2))

time <- cbind( t_0 = 0, t_1,  t_2, t_3)
state <- cbind(m_0 = "No Event", m_1, m_2, m_3)
res <- list("t"= time, "st" = state)


