#' calculate_k_function, calculates Ripley's K-function with border correction 
#' with options to specify cell phenotype, plot, plot confidence envelope and plot confidence interval
#' of K function. Plotted K-functions utilise an inhomogenous Poisson distrubution with 
#' point intensity produced to reflect the original data as the null distribution. Unable to plot both condience interval and envelope
#' simultaneously 
#' 
#' Inputs - 
#' point_pattern: ppp object (marked) to determine K function 
#' determine_confidence: assign any integer value to indicate to plot the confidence interval for the K-function itself. Set to NULL as default, calculating this can be computationally expensive
#' phenotypes_of_interest: the cell phenotype/phenotypes to be considered in the K-function calculation if not the entire ppp object. Set to NULL as default
#' plot_envelope: assign any integer value to indicate to plot the envelope of the simulated inhomogenous poisson distribution
#' plot: assign any integer value to indicate to plot the K-function itself, default is to plot. Set to NULL if plot is not to be produced
#' Outputs - 
#' k_function: Array containing values required to plot K function and the associated inhomogenous Poisson distribution (null distribution)
#' @description Returns K function of cell distribution
#' 
#' @param point_pattern ppp object used as input data for K function 
#' @param determine_confidence Indicate to plot confidence interval of the K function calculated (adds significant run time)
#' @param phenotypes_of_interest String supplying the phenotype of interest, if unspecified the entire point pattern will be considered
#' @param plot_envelope Indicates whether the significance envelope of simulated data should be plotted (based on inhomogenous Poisson distribution)
#' @param plot Indicates whether to plot K-function
#' @importFrom spatstat.core lohboot Kinhom envelope

calculate_k_function <- function(point_pattern, phenotypes_of_interest = NULL, determine_confidence = NULL, plot_envelope = NULL, plot = 1){
  #check for phenotype of interest and replace point pattern object as required
  if(!is.null(phenotypes_of_interest)) {
    point_pattern <- find_subset(point_pattern, phenotypes_of_interest)
  } 
  
  #calculate K function
  D <- density(point_pattern)
  k_function <- Kinhom(point_pattern, correction = "border", D)
  
  #plot confidence interval if required
  if (!is.null(determine_confidence)) {
    confidence <- lohboot(point_pattern, "Kinhom", nsim = 50)
    plot(confidence)
  }
  
  #plot significance envelope if required
  if (!is.null(plot_envelope)) {
    envelope_k <- envelope(point_pattern, Kinhom, funargs = c(D, correction = "border"),  simulate = expression(rpoispp(D)), nsim= 49)
    plot(envelope_k, main = "K-function", xlab = "Radius")
  } 
  #plot K-function if required
  else if (!is.null(plot)){
    plot(k_function, main = "K-function", xlab = "Radius")
  }
  
  #return K function plot data
  return(k_function)
}
