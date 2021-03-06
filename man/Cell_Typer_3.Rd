% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellTyping.R
\name{Cell_Typer_3}
\alias{Cell_Typer_3}
\title{Updated cell_typer which incorporates neural networks to distinguish n14 vs R13
Now only the proprioceptors and jagged need to be classified. Everything else is automated.}
\usage{
Cell_Typer_3(dat, modelDir = NA)
}
\arguments{
\item{dat}{RD.experiment}

\item{UL_classify}{boolean decision to classify large diameter cells}
}
\description{
Updated cell_typer which incorporates neural networks to distinguish n14 vs R13
Now only the proprioceptors and jagged need to be classified. Everything else is automated.
}
