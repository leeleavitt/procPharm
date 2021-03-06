% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelTools.R
\name{potassiumModeler}
\alias{potassiumModeler}
\title{Function to apply a general purpose potassium model which was trained on a potassium 
dose response experiment. This will apply the model to any potassium window
except the first k40 window, which has its own model to use}
\usage{
potassiumModeler(dat, modelLoc = NA)
}
\description{
Function to apply a general purpose potassium model which was trained on a potassium 
dose response experiment. This will apply the model to any potassium window
except the first k40 window, which has its own model to use
}
