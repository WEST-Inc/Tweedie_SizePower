#' @export
#' 
#' @name FitTrendModel_Tweedie
#'
#' @description Fit trend model with Tweedie GLMM using the 
#'   glmmTMB package (Brooks et al. 2017) and the form of the linked function 
#'   of fixed and random effects following that of Piepho and Ogutu (2002). 
#'   Check for model convergence and simplify the trend model in by removing
#'   random effects in the following order: random site-by-year interaction
#'   effect, random slope effect, random year effect. 
#'
#' @title  Wrapper function for power simulation.
#'
#' @param fit Trend model fit; NULL if no model has been successfully fit yet.
#'
#' @param samp The sample data.
#' 
#' @param slope  Logical indicator that random slope effect is currently in the trend model.
#' 
#' @param int Logical indicator that random site-by-year interaction effect is currently in the trend model.
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
#'    author = {Mollie E. Brooks and Kasper Kristensen and Koen J. {van Benthem} and Arni Magnusson and Casper W. Berg and Anders Nielsen and Hans J. Skaug and Martin Maechler and Benjamin M. Bolker},
#'    title = {{glmmTMB} Balances Speed and Flexibility Among Packages for Zero-inflated Generalized Linear Mixed Modeling},
#'    year = {2017},
#'    journal = {The R Journal},
#'    doi = {10.32614/RJ-2017-066},
#'    pages = {378--400},
#'    volume = {9},
#'    number = {2},
#'   url      = "https://doi.org/10.1198/108571102366",
#'   textVersion=
#'     paste("Brooks et al. (2017),",
#'           "glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling", 
#'           "The R Journal 9(2): 378-400",
#'           sep=" ")
#' )

FitTrendModel_Tweedie <- function(fit,samp,slope,int) {
  
# If model null and int = TRUE, fit full model. 
    if(is.null(fit) & int) {
      tryCatch({
        fit <- glmmTMB(Y ~ WYear + (1+WYear|Site) + (1|Year) + (1|Year:Site),
                 data = samp,
                 family = tweedie(link = 'log'),
                 control=glmmTMBControl(optCtrl = list(iter.max=500,eval.max=500),
                 profile=TRUE,collect=FALSE))
        }, error = function(e)cat("ERROR :",conditionMessage(e), "\n"))
       }

   # Check for perfectly correlated site int and slope REs
   # check for convergence - if not, remove random slope effect
   if(is.null(fit)&int) {
     tryCatch({
        fit <- glmmTMB(Y ~ WYear + (1|Site) + (1|Year) + (1|Year:Site), # drop slope RE
                 data = samp,
                 family = tweedie(link = 'log'),
                 control=glmmTMBControl(optCtrl = list(iter.max=500,eval.max=500),
                 profile=TRUE,collect=FALSE))
        }, error = function(e)cat("ERROR :",conditionMessage(e), "\n"))
     slope <- FALSE
  }

  if(is.null(fit)| !int) {   # if still NULL, drop random site:year interaction RE
     tryCatch({
        fit <- glmmTMB(Y ~ WYear + (1|Site) + (1|Year) , # drop int RE
                 data = samp,
                 family = tweedie(link = 'log'),
                 control=glmmTMBControl(optCtrl = list(iter.max=500,eval.max=500),
                 profile=TRUE,collect=FALSE))
        }, error = function(e)cat("ERROR :",conditionMessage(e), "\n"))
     slope <- FALSE
     int <- FALSE
  }
   if(is.null(fit)) {   # if still NULL, drop year RE
     tryCatch({
       fit <- glmmTMB(Y ~ WYear + (1|Site) , # drop int RE
                 data = samp,
                 family = tweedie(link = 'log'),
                 control=glmmTMBControl(optCtrl = list(iter.max=500,eval.max=500),
                 profile=TRUE,collect=FALSE))
        }, error = function(e)cat("ERROR :",conditionMessage(e), "\n"))
    }

  # check for convergence - if not, remove random slope effect
  if(slope) {   # if have slope, then also have int - check for convergence
     absCorr <- abs(attr(VarCorr(fit)$cond$Site,"correlation")[1,2])
     SEtrend <- summary(fit)$coef$cond[2, 2]
     if(is.na(AIC(fit))|is.na(SEtrend)|absCorr>0.95) {
      tryCatch({
         fit <- glmmTMB(Y ~ WYear + (1|Site) + (1|Year) + (1|Year:Site),  # drop slope RE
                  data = samp,
                  family = tweedie(link = 'log'),
                  control=glmmTMBControl(optCtrl = list(iter.max=500,eval.max=500),
                  profile=TRUE,collect=FALSE))
         }, error = function(e)cat("ERROR :",conditionMessage(e), "\n"))
        slope <- FALSE
     }
  }   # end if slope TRUE
  if((!slope)) {  # & int) {   # if no slope, could have int or not - check for convergence via SEtrend
     SEtrend <- summary(fit)$coef$cond[2, 2]
     if(is.na(AIC(fit))|is.na(SEtrend)) {
     tryCatch({
         fit <- glmmTMB(Y ~ WYear + (1|Site) + (1|Year),  # drop interaction RE
                   data = samp,
                   family = tweedie(link = 'log'),
                   control=glmmTMBControl(optCtrl = list(iter.max=500,eval.max=500),
                   profile=TRUE,collect=FALSE))
         }, error = function(e)cat("ERROR :",conditionMessage(e), "\n"))
      int <- FALSE
      }
  }   # end if slope FALSE & int TRUE 

  if(!slope & !int) {   # if no slope or int - check for convergence via SEtrend
    SEtrend <- summary(fit)$coef$cond[2, 2]
    if(is.na(AIC(fit))|is.na(SEtrend)) {
      tryCatch({
        fit <- glmmTMB(Y ~ WYear + (1|Site) ,  # drop year RE
                  data = samp,
                  family = tweedie(link = 'log'),
                  control=glmmTMBControl(optCtrl = list(iter.max=500,eval.max=500),
                  profile=TRUE,collect=FALSE))
        }, error = function(e)cat("ERROR :",conditionMessage(e), "\n"))
     }
  }   # end if !slope & !int

return(list(fit,slope,int))


}