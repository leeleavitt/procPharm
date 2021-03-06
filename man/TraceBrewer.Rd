% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tracePrep.R
\name{TraceBrewer}
\alias{TraceBrewer}
\title{This is our trace cleaning protocol}
\usage{
TraceBrewer(dat, bscore = F, blcPrep = T, verbose = F)
}
\arguments{
\item{dat}{is the RD.experiment to input}

\item{bscore}{logical, if true all scores will be replaced by bscore2 function}

\item{blcPrep}{logical, if this is true, an expensive baseline correction schema occurs. If false only some SNR occurs}
}
\description{
This is our trace cleaning protocol
}
