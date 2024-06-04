#' @export
#' 
#' @name CalcPower_Tweedie
#'
#' @description Conduct a Monte Carlo power simulation for trend detection in a
#'   Tweedie random variable.
#'
#' @title Calculate Trend Analysis for Power Simulation 
#'
#' @param its Number of iterations (independent populations generated) used to evaluate each simulation scenario.
#'
#' @param coreNum  Number of cores to use in parallel processing.
#' 
#' @param alpha Type I error rate.
#'
#' @param sampsize  Number of unique sites in the entire sample.
#' 
#' @param yrs  Number of years in the monitoring period; may be a vector of years.
#' 
#' @param p  Annual trend (e.g. p = -0.01 is a 1% annual decrease).
#'
#' @param revisit A numerical value from 1-4 to set sampling revisit design.
#'  revisit=1 = [1-4] with a single alternating panel
#'  revisit=2 = [1-0,1-4] with a single alternating panel
#'  revisit=3 = [1-0] 
#'  revisit=4 = [1-0,1-2] with a single alternating panel
#'  revisit=5 = [1-2] with a single alternating panel
#'    
#' @param N Number of sites in the population. 
#' 
#' @param StartYear Calendar year for beginning of monitoring period.
#'
#' @param mu  Fixed intercept effect on link (log) scale.
#' 
#' @param sig2b Year to year variation.
#' 
#' @param sig2a Site to site variation.
#' 
#' @param sig2t Variation among site trends.
#' 
#' @param sigat Covariance between site intercept and slope
#' 
#' @param sig2c Interaction variation.
#' 
#' @param sig2e Overdispersion random effect (Lee and Nelder 2000).
#' 
#' @param annual The number of sites in the annual panel (NA if there is no annual panel). 
#'   
#' @param m The number of within-year replications in annual/main panels.
#'   
#' @param m_alt The number of within-year replications in alternating panels.
#' 
#' @param Tpower The Tweedie power parameter. 
#' 
#' @param phi The Tweedie dispersion parameter.
#' 
#' @details Conduct a Monte Carlo power simulation to detect trend in a
#'   Tweedie random variable specified by fixed and random effects, the Tweedie
#'   power parameter, and the Tweedie dispersion parameter. Power is calculated
#'   for specified temporal revisit designs, sample sizes of sites, and monitoring
#'   period lengths. Power is assessed for a two-sided trend test based on the 
#'   z-test, t-test with between-within df (Schluchter and Elashoff 1990), and
#'   the likelihood ratio test.   
#'
#' @return A list of power analysis results for three trend tests and simulation results by iteration.
#'
#' @author Leigh Ann Starcevich, WEST, Inc.
#' 

CalcPower_Tweedie <-function (its,coreNum,alpha,sampsize, yrs, p,
                                                 revisit, N,StartYear,mu=0,
                                                 sig2b=0,sig2a=0,sig2t=0,sigat=0,
                                                 sig2c=0,sig2e=0,annual=NA, 
                                                 m,m_alt,
                                                 Tpower=1,phi=NA)  {
# define logit
    logit <- function(x) log(x/(1-x))
# define clusters for parallel processing
    cl <- makeCluster(coreNum)
    registerDoParallel(cl)
# run power simulatio
    ans.list<-foreach(ii=icount(its),.export=c("Sim_Wrapper_Tweedie",
                                                "SimPopnChangeLogY_Tweedie",
                                                "GetRevisitSamp","FitTrendModel_Tweedie",
                                                "its","logit"),.combine=list,
                       .maxcombine=its, .inorder=TRUE,.multicombine=TRUE,
                       .packages=c('MASS','glmmTMB','mgcv','reshape','lmtest','insight')) %dopar%  {

                        Sim_Wrapper_Tweedie(alpha=alpha,sampsize=sampsize,yrs=yrs,
                                             p=p,revisit=revisit, N=N,StartYear=StartYear,
                                             mu=mu,sig2b=sig2b,sig2a=sig2a,sig2t=sig2t,
                                             sigat=sigat, sig2c=sig2c,sig2e=sig2e,
                                             annual=annual,m=m,m_alt=m_alt,
                                             Tpower=Tpower,phi=phi)
                         }

stopCluster(cl)

# Format output
    trendTest_z <- do.call(rbind, lapply(1:length(ans.list), function(u, x) x[[u]][[1]], x=ans.list)) # save trend test results
    trendTest_LRT <- do.call(rbind, lapply(1:length(ans.list), function(u, x) x[[u]][[2]], x=ans.list)) # save trend test results
    trendTest_BW <- do.call(rbind, lapply(1:length(ans.list), function(u, x) x[[u]][[3]], x=ans.list)) # save trend test results
    trendEst <- do.call(rbind, lapply(1:length(ans.list), function(u, x) x[[u]][[4]], x=ans.list)) # save trend test results
    trendSE <- do.call(rbind, lapply(1:length(ans.list), function(u, x) x[[u]][[5]], x=ans.list)) # save trend test results
    ans_z <- aggregate(trendTest_z[,4],list(sampsize=trendTest_z[,1], yrs=trendTest_z[,2],m=trendTest_z[,3]),mean)
    ans_LRT <- aggregate(trendTest_LRT[,4],list(sampsize=trendTest_LRT[,1], yrs=trendTest_LRT[,2],m=trendTest_z[,3]),mean)
    ans_BW <- aggregate(trendTest_BW[,4],list(sampsize=trendTest_BW[,1], yrs=trendTest_BW[,2],m=trendTest_z[,3]),mean)
    names(ans_z) <- names(ans_LRT) <- names(ans_BW) <- c('sampsize','yrs','m','power')
    return(list(ans_z,ans_LRT,ans_BW,trendTest_z,trendTest_LRT,trendTest_BW,trendEst,trendSE,ans.list))
}

