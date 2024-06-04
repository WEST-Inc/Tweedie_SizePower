#' @export
#' 
#' @name SimPopnChangeLogY_Tweedie
#'
#' @description Simulate a population of random variables, Y, that follow the Tweedie distribution.
#'  Generate Y with known trend with estimates of status and variation from modeled data.
#'
#' @title  Simulate population of Tweedie random variables with known trend.
#' 
#' @param N Number of units in the population. 
#' 
#' @param w Vector of survey years for the power simulation, where the first
#' year is 0 and the last year is the number of years minus 1.
#' 
#' @param StartYear Calendar year for beginning of monitoring period.
#'
#' @param p Annual trend (e.g. p = -0.01 is a 1% annual decrease).
#' 
#' @param mu Fixed intercept effect on link (log) scale.
#' 
#' @param sig2b Year to year variation.
#' 
#' @param sig2a Site to site variation.
#' 
#' @param sig2t Variation among site trends.
#' 
#' @param sigat Covariance between site intercept and slope.
#' 
#' @param sig2c Interaction variation.
#' 
#' @param sig2e Overdispersion random effect (Lee and Nelder 2000).  
#' 
#' @param m The number of within-year replications in annual/main panels.
#'   
#' @param m_alt The number of within-year replications in alternating panels.
#' 
#' @param Tpower The Tweedie power parameter. 
#' 
#' @param phi The Tweedie dispersion parameter.
#' 
#' 
#' @details Generates a population in decline/growth measured in annual percent change
#' (p) over a number of years (w) for a population of size N.  The trend model of 
#' Piepho and Ogutu (2002) is used assuming an "inter-site" analysis which can be 
#' modeled as replicated (m>1, sig2c>0) or unreplicated (m=1, sig2c=0).
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
#'
#'
#' @return 
#'
#' @author Leigh Ann Starcevich, WEST, Inc.

SimPopnChangeLogY_Tweedie<-function(N,w,StartYear,p,mu=0,sig2b=0,sig2a=0,sig2t=0,
                            sigat=0,sig2c=0,sig2e=0,m=12,Tpower=NA,phi=NA) {
	
noyrs<-length(w)							     # number of years 
index<-1
p<-p+1	                                               # Annual change as a multiplier

# Generate random effects
bj<- rnorm(noyrs,0, sqrt(sig2b))		                 # random effect of jth year
at<- mvrnorm(N,c(0,0), matrix(c(sig2a, sigat, sigat, sig2t),2,2), tol=1)
ai<-at[,1]                                                 # random intercept of ith site
ti<-at[,2]                                                 # random trend effect of ith site
cij<- matrix(rnorm(N*noyrs,0, sqrt(sig2c)),N,noyrs)	     # site by yr interaction

# Create outcome matrix
yij <- data.frame(matrix(rep(0,6*N*noyrs*m),N*noyrs*m,6))

for (i in 1:N) {                                            # Site loop
  for (j in 1:noyrs){                                       # Year loop
       yij[index:(index+m-1),1]<-StartYear+ w[j]            # Year
       yij[index:(index+m-1),2]<-w[j]                       # WYear
       yij[index:(index+m-1),3]<-i                          # Site
       yij[index:(index+m-1),4]<-1:m                        # Visit
# Generate mean of Tweedie RV on logged scale
       yij[index:(index+m-1),5]<-(mu			      # Fixed intercept  
          +(w[j]*log(p))			                  # trend = log(p)
          +bj[j]	
          +ai[i]
          +ti[i]	
          +cij[i,j]
      )	
			

# Generate Tweedie Y with mu = exp(yij)
    yij[index:(index+m-1),6]<-rTweedie(mu=exp(yij[index:(index+m-1),5]),p=Tpower, phi=phi)		# original scale
    index <- index + m                                      # overall index of sites x years x seasons x weeks
    }
  }

yij<- data.frame(yij)
colnames(yij)<-c("Year","WYear","Site","Month","LogY","Y")

yij$Year <- as.factor(yij$Year)
yij$Site <- as.factor(yij$Site)

return(yij)
}
