% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logger.R
\name{logger}
\alias{logger}
\title{Function to track the various users data during analysis
The goal is to append a line to the log with}
\usage{
logger(functionName, timeInFunction, additionalInfo = NA, rate = F)
}
\description{
Needs to happen before the function goes
time1 <- proc.time()
Fill up your additional info...
additionalInfo <- choices
timeInFunction <- (proc.time() - time1)[3]
functionName <- as.character(match.call())[1]
additionalInfo <- choices
logger(functionName, timeInFunction, additionalInfo)
}
