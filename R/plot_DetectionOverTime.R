#' Diagnostics for the detection model with respect to Length 
#' 
#' Creates a plot of detectability by year for differing list lengths from an occupancy model output.

#' @param model a fitted sparta model of class \code{OccDet}.
#' @param spname optional name of the species (used for plotting)
#' @param min.yr optional first year of time series (used for plotting)
#' @param legend_labels optional names for legend labels. Should be a character vector with three elements if the model is fitted with categorical or continuous list length specifications, and four elements if the model is fitted with a mixed list length specification
#' @param legend_title optional name for legend title. Should be a character vector.
#' 
#' @details 
#' Takes a object of \code{OccDet}
#' 
#' Calculates the detection probability and produces a plot of detectability over time for the reference data type.
#'
#' @return This function returns plot showing the detection probability on the y axis and year on the x.
#'
#' @importFrom reshape2 melt 
#' @importFrom plyr ddply
#' @importFrom ggplot2 ggplot
#' @importFrom boot inv.logit
#' @export


plot_DetectionOverTime <- function(model, spname = NULL, min.yr = NULL, legend_labels = NULL, legend_title = NULL){

  sims_list <- model$BUGSoutput$sims.list
  
  # the base: alpha.p is common to all models: 
  # it's the logit probability of detection on a single species list
  pDet1 <- sims_list$alpha.p
  # pDet1 is an array of dims equal to (niter, nyr)
  
  if("beta1" %in% names(sims_list)){
    # we ran the Julian Date option
    # So let's scale the detection probabilities to end June (day 180)
    pDet1 <- apply(pDet1, 2, function(x) 
      x + 180 * sims_list$beta1[,1] +  180^2 * sims_list$beta2[,1]
      )
  }
  
  pDet5 <- NULL # Set pDet5 to NULL as it is only used when plotting mixed list lengths, and therefore would return an error when melt-ing pDet lists for continous and categorical list lengths
  
  # now calculate the equivalent values for different lists length specifications
  if("LL.p" %in% names(sims_list) && "dtype2.p"  %in% names(sims_list)){
    # the model was fitted with a mixed list length
    pDet2 <- pDet1 + sims_list$dtype2.p[,1] 
    pDet4 <- pDet1 + sims_list$dtype3.p[,1]
    pDet5 <- pDet1 + sims_list$LL.p[,1]*log(5) # data type 1 with list length 5
    } else if("LL.p" %in% names(sims_list)){
      # the model was fitted with continuous list length
      pDet2 <- pDet1 + sims_list$LL.p * log(2) # list length 2
      pDet4 <- pDet1 + sims_list$LL.p * log(4) # list length 4
      } else if("dtype2.p" %in% names(sims_list)){
        # the model was fitted with categorical list length
        pDet2 <- pDet1 + sims_list$dtype2.p[,1]
        pDet4 <- pDet1 + sims_list$dtype3.p[,1]
        } 
  # there is also an option to ignore list length, 
  # in which case the probability of detection is assumed to be constant across surveys
  # i.e. if the survey was systematic
  
  pDet <- melt(list(pDet1, pDet2, pDet4, pDet5))
  names(pDet) <- c("it", "year", "lgt_pDet", "ListLength")
  pDet$ListLength[pDet$ListLength>=3] <- pDet$ListLength[pDet$ListLength>=3] + 1 # the "third" category is for a list of length 4 for continuous LLs (>=4 for categorical) and the "fourth" is for a list length of 5 for mixed LLs
  
  pDet$pDet <- inv.logit(pDet$lgt_pDet)
  
  # now summarize these posterior distributions
  pDet_summary <-ddply(
        pDet, .(year, ListLength), summarise, 
        mean_pDet = mean(pDet),
        lower95CI = quantile(pDet, 0.025),
        upper95CI = quantile(pDet, 0.975))
  
  # if the user has supplied a year then switch the x axis to start at that minimum
  if(!is.null(min.yr)) pDet_summary$year <- pDet_summary$year + min.yr - 1
  
  # now plot the detection over time
  gp <- ggplot(data=pDet_summary, x=year, y=mean_pDet) +
    geom_line(aes(x=year, y=mean_pDet, col=factor(ListLength))) +
    geom_ribbon(aes(x=year, ymin=lower95CI, ymax=upper95CI, fill=factor(ListLength)), alpha=0.2) +
    ylab("Detection probability") +
    xlab("Year") +
    ggtitle(spname) +
    theme_bw()
  
    if(!is.null(legend_title)){
      if(!(is.character(legend_title))){
        stop('legend_title is not a character vector')
      } else{
        gp <- gp + labs(color=legend_title, fill = legend_title)
      }}
  
  if(!is.null(legend_title) && !is.null(legend_labels)){
    if(!(is.character(legend_labels))){
      stop('legend_labels is not a character vector')
    } else{
      gp <- gp + scale_fill_discrete(name = legend_title, 
                                   labels = legend_labels) +
        scale_color_discrete(name = legend_title, 
                             labels = legend_labels)
    }}
  
  gp 

}


