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
# install.packages("ggplot2")

## Sharpness boxplot:

library("ggplot2")


## Create test data for function development

testdata <- data.frame(y = c(rnorm(100, 0), rnorm(100, 0.1),
                             rnorm(100, 0.2), rnorm(100, 0.3)),
                       x = factor(rep(c(1:2, 1:2), each = 100),
                                  labels = c("mod1", "mod2")),
                       z = factor(rep(1:2, each = 200),
                                  labels = c("lag_01", "lag_02")))


## start of actual function
ggplot(testdata, aes(x = x, y = y)) +
    geom_boxplot() + facet_grid(. ~ z)




