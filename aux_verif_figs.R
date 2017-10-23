##################################################################################
#
## Auxiliary functions for plotting of verification figures
#
# For informations on using ggplot2:
#  - http://t-redactyl.io/blog/2016/04/creating-plots-in-r-using-ggplot2-part-10-boxplots.html
#  - http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization
#
#
##################################################################################

## install needed packages:
if(FALSE) {
    install.packages("ggplot2")
}

## load libraries:
library(dplyr)
library("ggplot2")

    
## Auxiliary function for plotting of sharpness boxplots:

sharp_boxplots <- function(datavec, alpha, mods, permut = NULL, lts, lts_type, ncol = NULL, horizlab = TRUE) {
  ## - datavec: numeric vector of prediction interval widths
  ## - alpha the width of the prediction interval of interest
  ## - mods:  character vector of model names
  ## - permut: permutation of the models (needed for custom ordering of ggplot2 boxplots)
  ## - lts: numberic vector of lead times (in an arbitrary unit (i.e., hour, day, month))
  ## - lts_type: character vector that defines the unit of the lead times
  ## - ncol: number of columns (defaults to NULL)
  ## - horizlab: logical argument, defines direction of x labels (defaults to TRUE)
  
  tdata <- data.frame(y = datavec,
                      mod = factor(rep(rep(mods, length(lts)),
                                       each = length(datavec)/length(lts)/length(mods))),
                      lt =  factor(rep(paste("lead month", lts, sep = " "), each = length(datavec)/length(lts))))
  
  df <- data.frame(matrix(unlist(by(tdata[,1], tdata[,c("mod", "lt")], quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))), nrow = nlevels(tdata$lt) * nlevels(tdata$mod), byrow = T))                                                                                                                                            
  colnames(df) <- c("y05", "y25", "y50", "y75", "y95")
  
  df1 <- mutate(df, lt = rep(unique(tdata$lt), each = nlevels(tdata$mod)),
                mod = rep(unique(tdata$mod), nlevels(tdata$lt)))
  
  if(!is.null(permut)) {
      df1$mod <- factor(df1$mod, levels = mods[permut])
  }  

  p <- ggplot(df1, aes(x=mod)) +
    geom_boxplot(aes(ymin = y05, lower = y25, middle = y50, upper = y75,    max = y95),
                 stat = "identity") +
    facet_wrap(~lt, ncol = ncol) + theme_bw() + xlab("") +
    ylab(paste(alpha, "% prediction interval width", sep = " "))
  
  if(horizlab) {
    p
  } else {
    p + theme(axis.text.x = element_text(angle = 90, hjust = 0))
  }
  
}


    
## testing of the boxplot function
if(FALSE) {
    sharp_boxplots(datavec = rnorm(1500, rep(seq(0,0.5*30, length.out = 15), each = 100)),
                   alpha = 75,
                   mods = paste("model_", 1:5, sep = ""),
                   permut = NULL,
                   lts = c(1,3,7),
                   lts_type = "lead month",
                   ncol = 3,
                   horizlab = FALSE)
}





