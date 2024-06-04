#' @export
#' 
#' @name Sim_Wrapper_Tweedie
#'
#' @description Wrapper function for power simulation to obtain the power to 
#' detect a trend in a Tweedie random variable with a two-sided trend test.
#' The trend model on the link scale takes the form of the trend model of 
#' Piepho and Ogutu (2002).
#'
#' @title  Wrapper function for power simulation.
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
#' @details Simulate a population specified by fixed and random effects, draw a sample given
#'   a specific sample size and temporal revisit design, fit a trend model, and conduct
#'   checks for model convergence with model reduction as required. 
#'
#' @return The results from power analyses of a two-sided trend test based on the z-test, 
#'   t-test with between-within df (Schluchter and Elashoff 1990), and the likelihood ratio test.  
#'   Additionally, trend estimates and standard errors for every iteration of the simulation are returned. 
#'
#' @author Leigh Ann Starcevich, WEST, Inc.
#' 
#' citEntry(entry="article",
#'   author   = "J. Lee and J.A. Nelder",
#'   title    = "Two ways of modelling overdispersion in non-normal data",
#'   journal  = "Journal of the Royal Statistical Society Series C: Applied Statistics",
#'   year     = 2000,
#'   volume   = 49,
#'   number   = 4,
#'   pages    = "591–598",
#'   url      = "https://doi.org/10.1111/1467-9876.00214",
#'   textVersion=
#'     paste("Lee and Nelder (2000),",
#'           "Two ways of modelling overdispersion in non-normal data,", 
#'           "JRSS C 49(4): 591-598",
#'           sep=" ")
#' )
#' citEntry(entry="article",
#'   author   = "H.P. Piepho and J.O. Ogutu",
#'   title    = "A simple mixed model for trend analysis in wildlife populations",
#'   journal  = "Journal of Agricultural, Biological, and Environmental Statistics",
#'   year     = 2002,
#'   volume   = 7,
#'   pages    = "350-360",
#'   url      = "https://doi.org/10.1198/108571102366",
#'   textVersion=
#'     paste("Piepho and Ogutu (2002),",
#'           "A simple mixed model for trend analysis in wildlife populations", 
#'           "JABES 7: 350-360",
#'           sep=" ")
#' )
#' citEntry(entry="article",
#'   author   = "M.D. Schluchter and J. T. Elashoff",
#'   title    = "Small-sample adjustments to tests with unbalanced repeated measures assuming several covariance structures",
#'   journal  = "Journal of Statistical Computation and Simulation",
#'   year     = 1990,
#'   volume   = 37,
#'   number   = "1-2",
#'   pages    = "69–87",
#'   url      = "https://doi.org/10.1198/108571102366",
#'   textVersion=
#'     paste("Schluchter and Elashoff (1990),",
#'           "Small-sample adjustments to tests with unbalanced repeated measures assuming several covariance structures", 
#'           "JSCS 37(1-2): 69–87",
#'           sep=" ")
#' )

Sim_Wrapper_Tweedie <- function (alpha=alpha,sampsize=sampsize, 
                                             yrs=yrs,p=p,revisit=revisit, 
                                             N=N,StartYear=StartYear,mu=mu,
                                             sig2b=sig2b,sig2a=sig2a,sig2t=sig2t,
                                             sigat=sigat, sig2c=sig2c,sig2e=sig2e,
                                             annual=annual,m=m,m_alt=m_alt,
                                             Tpower=Tpower,phi=phi)  {

    annual_orig <- annual
    nsamp <- length(sampsize)
    nyrs <- length(yrs) 
    ans_test_z <- ans_test_LRT <- ans_test_BW <- ans_trend <- ans_SEtrend <-data.frame(matrix(0,nsamp*nyrs*length(m),4))

    # Simulate population
    simpop <- SimPopnChangeLogY_Tweedie(N, w = 0:(max(yrs) - 1),StartYear,p,mu,sig2b,
                                sig2a,sig2t,sigat,sig2c,sig2e,
                                m=max(m),Tpower=Tpower,phi=phi)
    index <- 0
    for (i in 1:nsamp) {         # sample size loop
      annual <- annual_orig

      # Apply revisit design and draw sample
      if(length(dim(annual_orig))>1) annual <- unlist(annual[i,]) # get panel sizes, if specified
      # Draw sample for all years
      samp_ALL<- GetRevisitSamp(sampsize[i], max(yrs), simpop, revisit, annual, m=max(m), m_alt=m_alt)

      for (j in 1:nyrs) {     # years loop
        for (k in 1:length(m)) {     # revisits loop
          index <- index + 1
          samp<- samp_ALL[samp_ALL$WYear %in% 0:(yrs[j]-1) & samp_ALL$Month <= m[k],]  
          samp$Site <- as.factor(as.character(samp$Site))
          samp$Year <- as.factor(as.character(samp$Year))

          # Fit trend model
          fit <- NULL
          slope <- int <- TRUE      # Start with full model
          if(m[k]==1) int <- FALSE  # If have no Site-by-Year replication, remove random intercept effect

          fitModel <- FitTrendModel_Tweedie(fit,samp,slope,int)
          fit <- fitModel[[1]]
          slope <- fitModel[[2]]
          int <- fitModel[[3]]

          fit_red <- NULL
          tryCatch({
            fit_red <- update(fit, ~.-WYear)
          }, error = function(e)cat("ERROR :",conditionMessage(e), "\n"))

   # Check convergence of the reduced model.
   # If the reduced model did not converge, then simplify full model and refit reduced model
       if(is.null(fit_red) | is.na(summary(fit_red)$coef$cond[1, 2]) | is.na(AIC(fit_red))) {                                      # if any issues with red model, refit full model
          fitModel <- FitTrendModel_Tweedie(fit,samp,slope,int)
          fit <- fitModel[[1]]
          slope <- fitModel[[2]]
          int <- fitModel[[3]]

     # Assume we have a final model, and fit the reduced model
          tryCatch({
          fit_red <- update(fit, ~.-WYear)
       }, error = function(e)cat("ERROR :",conditionMessage(e), "\n"))
     }   # end checks for reduced model

  # If LRT p-value is NA, then reduce model again
    lrt_p <- lmtest::lrtest(fit, fit_red)$"Pr(>Chisq)"[2]
    if(is.na(lrt_p)) { # if LRT had error
       fitModel <- FitTrendModel_Tweedie(fit,samp,slope,int)
       fit <- fitModel[[1]]
       slope <- fitModel[[2]]
       int <- fitModel[[3]]

# Assume we have a final model, and fit the reduced model
      tryCatch({
        fit_red <- update(fit, ~.-WYear)
      }, error = function(e)cat("ERROR :",conditionMessage(e), "\n"))
    }   # end checks for LRT

# Assume converging full and reduced models have been obtained.
# Capture trend test results
        coefs<-as.data.frame(summary(fit)$coef$cond)
        ans_test_z[index,]<-c(sampsize[i],yrs[j],m[k],ifelse(coefs[2,4]>alpha,0,1)) # 1 if reject Ho
        ans_trend[index,] <- c(sampsize[i],yrs[j],m[k],coefs[2,1])
        ans_SEtrend[index,] <- c(sampsize[i],yrs[j],m[k],coefs[2,2])
        ans_test_LRT[index,]<-c(sampsize[i],yrs[j],m[k],ifelse(lrt_p > alpha,0,1)) # 1 if reject Ho
        ans_test_BW[index,]<-c(sampsize[i],yrs[j],m[k],ifelse(2*pt(abs(coefs[2,3]),insight::get_df(fit,type = "betwithin")[2],lower.tail=F)>alpha,0,1)) # 1 if reject Ho
} } }  # end i j k loops

# Delete unnecessary objects
rm(simpop,samp_All,samp,fit)
gc()

return(list(data.frame(ans_test_z),data.frame(ans_test_LRT),data.frame(ans_test_BW),data.frame(ans_trend),data.frame(ans_SEtrend)))
}
