####################################################
# Trend test power and size simulation code
# Western EcoSystems Technology, Inc.
# LA Starcevich 2024

####################################################
setWindowTitle("MD Test Power and Size Simulation")
#rm(list=ls())

#################################
# Load packages
#################################
library(MASS)
library(doParallel)
library(glmmTMB)
require(insight)

######################################
# define paths
######################################
path <- getwd() 
setwd(path)
codePath <- file.path(path, 'code')

######################################
# Read R scripts
######################################
file.sources <- list.files(path=codePath,pattern="*.R")
setwd(codePath)
sapply(file.sources,source,.GlobalEnv)
setwd(path)

####################################################
# Specify power/size simulation inputs
####################################################

# Define each variance components as defined in Piepho and Ogutu (2002)
sig2b_sim<-c(0.001,0.1,0.001,0.001,0.001)  # year-to-year variation
sig2a_sim<-c(0.001,0.001,0.1,0.001,0.001)  # site-to-site variation   
sig2t_sim<-c(0.001,0.001,0.001,0.1,0.001)  # site-level slope variation
sig2c_sim<-c(0.001,0.001,0.001,0.001,0.1)  # site-by-year interaction variation


# function to compute covariance based on variation and assumed correlation
calcCov <- function(sig2a,sig2t,cor_at) {
  return(cor_at*sqrt(sig2a*sig2t))
}

# use mean corr(ai,ti) from MD analysis to compute covariance
sigat_sim <- calcCov(sig2a_sim,sig2t_sim,-0.3517)  

# Use mean values from Marine Debris data to inform intercept/status
trendResults_sim <- data.frame(DataSet='simulated',
                               outcome=c('y1','y2','y3','y4','y5'),
                               region=1,
                               beta0 = 1.585,              # mean across all metrics and regions
                               logStatus = 1.6,            # mean across all metrics and regions
                               year.to.year.var = sig2b_sim,
                               site.to.site.var = sig2a_sim,
                               site.trend.var = sig2t_sim,
                               site.slope.intercept.cov = sigat_sim,
                               site.by.year.var = sig2c_sim,
                               Tpower = 1.59,  
                               phi = 1.045   
)


######################################
# set cores for parallel processing
detectCores()
coreNum = 4
logit <- function(x) log(x/(1-x))

######################################
# General settings

N <- 200     # population size from which sample is drawn
its <- 500   # number of iterations in each simulation

######################################
# Calculate annual trend
######################################

CalcAnnualTrend <- function(net,mb) {
# Takes net trend as a proportion and number of years
# and returns the associated annual trend. 
  return(((1+net)^(1/(mb-1)))-1)
}

# Calculate the annual trend corresponding to a 20% net decline over 11 years
p_20n <-CalcAnnualTrend(-.2,11) # -0.02206723

######################################
# Enumerate all simulation scenarios
######################################

simScenarios_sim <- expand.grid(p = c(0,p_20n), 
                           ma = c(50,62),         # Total sample size
                           mb = 11,  		  # Number of years
                           revisit=1:5,           # Revisit codes
                           annual = c(NA,12))     # Number of sites in annual panel, if needed

# remove the scenarios with revisit = 3 or 6 and annual = 12,
#                           revisit = 2 and annual = NA
simScenarios_sim <- simScenarios_sim[-which(simScenarios_sim$revisit %in% c(1,3,5) & simScenarios_sim$annual == 12),]
simScenarios_sim <- simScenarios_sim[-which(simScenarios_sim$revisit %in% c(2,4)  & is.na(simScenarios_sim$annual)),]
simScenarios_sim <- simScenarios_sim[-which(simScenarios_sim$ma==50 & simScenarios_sim$revisit %in% c(2,4)),]
simScenarios_sim <- simScenarios_sim[-which(simScenarios_sim$ma==62 & simScenarios_sim$revisit %in% c(1,3,5)),]
simScenarios_sim$m_alt <- 1
simScenarios_sim$m_alt[which(simScenarios_sim$revisit %in% c(1,3,5))] <- NA
simScenarios_sim   


################################################################
# Run power sims for all metrics and scenarios
################################################################

ans_sim <- NULL

# i indexes each metric
for (i in 1:nrow(trendResults_sim)) {       
   simInputs_sim.i <- trendResults_sim[(trendResults_sim$region == trendResults_sim$region[i]) & 
                               (trendResults_sim$outcome == trendResults_sim$outcome[i]) &
                               (trendResults_sim$DataSet == trendResults_sim$DataSet[i]),]
   if(nrow(simInputs_sim.i)>1) print('error - more than one set of inputs flagged') 
# j indexes each simulation scenario
   for (j in 1:nrow(simScenarios_sim)) {  
print(c(i,j))
# Calculate the power to detect a given trend for metric i and power/size scenario j
     powerSim_sim.ij <-  CalcPower_Tweedie(
           its=its,
           coreNum=coreNum,
           alpha=0.15,
           sampsize=simScenarios_sim$ma[j], 
           yrs=simScenarios_sim$mb[j], 
           p=simScenarios_sim$p[j], 
           revisit=simScenarios_sim$revisit[j], 
           N=N,
           StartYear=2023,
           mu=simInputs_sim.i$logStatus,
           sig2b=simInputs_sim.i$year.to.year.var,
           sig2a=simInputs_sim.i$site.to.site.var,
           sig2t=simInputs_sim.i$site.trend.var,
           sigat=simInputs_sim.i$site.slope.intercept.cov,
           sig2c=simInputs_sim.i$site.by.year.var,
           annual=simScenarios_sim$annual[j], 
           m=c(1,4,12),
           m_alt=simScenarios_sim$m_alt[j],
           Tpower=simInputs_sim.i$Tpower,
           phi=simInputs_sim.i$phi) 
# Remove commenting when running simulations
     #saveRDS(powerSim_sim.ij, paste0(
     #        'Output/power_',trendResults_sim$outcome[i],
     #        '_Region',trendResults_sim$region[i],
     #        '_revisit',simScenarios_sim$revisit[j],
     #        '_mb',simScenarios_sim$mb[j],
     #        '_ma',simScenarios_sim$ma[j],
     #        '_m',simScenarios_sim$m[j],
     #        '_p',simScenarios_sim$p[j],'.rds'))

     #save.image("PowerSims.RData")

     ans_sim <- rbind(ans_sim, data.frame(outcome=trendResults_sim$outcome[i],
             revisit=simScenarios_sim$revisit[j],
             mb=simScenarios_sim$mb[j],
             ma=simScenarios_sim$ma[j],
             m=c(1,4,12),
             p=simScenarios_sim$p[j],
             power_z=powerSim_sim.ij[[1]][,4],
             power_LRT=powerSim_sim.ij[[2]][,4],
             power_BW=powerSim_sim.ij[[3]][,4]))

      gc()

     }

}

save.image("PowerSims.RData")
     

