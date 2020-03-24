# install Dev tools
install.packages("devtools")

# load the software
library(devtools)

#Now create a new packages
create("procPharm")
create_package("C:/Users/leele/Desktop/procPharm")

# This is how we load in all the functions
load_all()

# This is how we add documentation to a function in R

#' Bind two factors
#'
#' Create a new factor from two existing factors, where the new factor's levels
#' are the union of the levels of the input factors.
#'
#' @param a factor
#' @param b factor
#'
#' @return factor
#' @export
#' @examples
#' fbind(iris$Species[c(1, 51, 101)], PlantGrowth$group[c(1, 11, 21)])
#' 
#' 
# lets load in our Rdata file
# After you have added documentation like above
document()

setwd("..")
load("./RD.200309.30.m.m3.p1.Rdata")

