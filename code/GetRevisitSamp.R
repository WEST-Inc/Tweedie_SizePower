#' @export
#' 
#' @name GetRevisitSamp
#'
#' @description Draw a sample with specified revisit design.
#'
#' @title Generate sample data for each revisit.
#'
#' @param sampsize Number of unique sites in the entire sample.
#'
#' @param yrs Number of years to sample.
#'
#' @param simpop Simulated population from which to draw sample.
#' 
#' @param revisit A numerical value from 1-4 to set sampling revisit design.
#'  revisit=1 = [1-4] with a single alternating panel
#'  revisit=2 = [1-0,1-4] with a single alternating panel
#'  revisit=3 = [1-0] 
#'  revisit=4 = [1-0,1-2] with a single alternating panel
#'  revisit=5 = [1-2] with a single alternating panel
#'    
#' @param annual The number of sites in the annual panel (NA if there is no annual panel). 
#'   
#' @param m The number of within-year replications in annual/main panels
#'   
#' @param m_alt The number of within-year replications in alternating panels.
#'   
#' @details Select a simple random sample of specified size for a given revisit design 
#'   and a number of years in the monitoring design. 
#'
#' @return A data frame object containing a sample of the specified size
#'  for the specified revisit design.
#'
#' @author Leigh Ann Starcevich, WEST, Inc.


GetRevisitSamp<-function(sampsize, yrs, simpop, revisit=1,annual=NA,m,m_alt) {
	
UniqueWYear<-unique(simpop$WYear)
N<-length(unique(simpop$Site))

# Select the total number of sites then apportion over the panels
Years<-sort(unique(simpop$WYear))[1:yrs]

# [1-4] One panel visited every 5 years
if(revisit==1) {   
	panelsperyr<-1
	panels<-1
	panelsizes<-rep(sampsize,1)
	panelcutoffs<-c(0,cumsum(panelsizes))
	sampsites<-sample(1:N,sum(panelsizes))
	panel1<-sampsites[1: panelcutoffs[2]]
	years1<-sort(Years[seq(1,yrs,5)])
	samp<-simpop[(is.element(simpop$Site,panel1))&(is.element(simpop$WYear,years1)),]
}

# [1-0, 1-4] One panel visited annually, one panel visited every 5 years
if(revisit==2) {   
	panelsperyr<-2
	panels<-2
      # if annual panel size not specified, allocate sites equally across panels
	if(is.na(annual)) {  # if annual panel size not specified, allocate evenly across panels
          panel2num<-floor(sampsize/panelsperyr)
	    panel1num<- sampsize-((panelsperyr-1)*floor(sampsize/panelsperyr)) 
      } else {
        panel1num <- annual
	  panel2num <- sampsize - annual[1]
      }
	flag<- (panel2num==panel1num)
	nosites<- (panel2num*(panels-1))+ panel1num
	sampsites<-sample(1:N,nosites)
	if(flag==TRUE) panelcutoffs<-seq(0,nosites,nosites/panels)
	if(flag==FALSE) panelcutoffs<-c(0,panel1num,seq(panel1num+panel2num,nosites, panel2num))

	panel1<-sampsites[1: panelcutoffs[2]]
	panel2<-sampsites[(panelcutoffs[2]+1): panelcutoffs[3]]

	years2<- UniqueWYear[seq(1,yrs,5)]

	samp<-simpop[is.element(simpop$Site,panel1)&(is.element(simpop$WYear, Years))&(simpop$Month <= m),]
      samp<-rbind(samp,
               simpop[(is.element(simpop$Site,panel2))&(is.element(simpop$WYear, years2))&(simpop$Month <= m_alt),])
}

# [1-0] One panel visited annually
if(revisit==3) {  
   nosites<- sampsize
   panels<-1
   panel1<-sample(1:N,nosites)
   samp<-simpop[is.element(simpop$Site,panel1) &(is.element(simpop$WYear, UniqueWYear[1:yrs])),] 
}

# [1-0, 1-2] One panel visited annually, one panel visited every 3 years
if(revisit==4) {   
	panelsperyr<-2
	panels<-2
      # if annual panel size not specified, allocate sites equally across panels
	if(is.na(annual)) {
          panel2num<-floor(sampsize/panelsperyr)
	    panel1num<- sampsize-((panelsperyr-1)*floor(sampsize/panelsperyr)) 
      } else {
        panel1num <- annual
	  panel2num <- sampsize - annual[1]
      }
	flag<- (panel2num==panel1num)
	nosites<- (panel2num*(panels-1))+ panel1num
	sampsites<-sample(1:N,nosites)
	if(flag==TRUE) panelcutoffs<-seq(0,nosites,nosites/panels)
	if(flag==FALSE) panelcutoffs<-c(0,panel1num,seq(panel1num+panel2num,nosites, panel2num))

	panel1<-sampsites[1: panelcutoffs[2]]
	panel2<-sampsites[(panelcutoffs[2]+1): panelcutoffs[3]]

	years2<- UniqueWYear[seq(1,yrs,3)]

	samp<-simpop[is.element(simpop$Site,panel1)&(is.element(simpop$WYear, Years))&(simpop$Month <= m),]
      samp<-rbind(samp,
                simpop[(is.element(simpop$Site,panel2))&(is.element(simpop$WYear, years2))&(simpop$Month <= m_alt),])
}

# [1-2] One panel visited every 3 years
if(revisit==5) {   
	panelsperyr<-1
	panels<-1
	panelsizes<-rep(sampsize,1)
	panelcutoffs<-c(0,cumsum(panelsizes))
	sampsites<-sample(1:N,sum(panelsizes))

	panel1<-sampsites[1: panelcutoffs[2]]
	
	years1<-sort(Years[seq(1,yrs,3)])

	samp<-simpop[(is.element(simpop$Site,panel1))&(is.element(simpop$WYear,years1)),]

}

return(samp)
}

