% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coralSymbODEs.R
\name{solveCoral}
\alias{solveCoral}
\title{Solve coral-\emph{Symbiodinium} ODEs}
\usage{
solveCoral(
  times = c(0, 500),
  pars = defPars(),
  lambda = 5,
  method = "vode",
  ...
)
}
\arguments{
\item{times}{The time values at which output is desired. Default c(0,500).}

\item{pars}{The paramaters of the model, given as a list. See function \code{\link{defPars}}. Default defPars().}

\item{lambda}{The sensitivity of the runs. High values are more sensitive and small values are less sensitive. Default 5.}

\item{method}{The character method argument of ode() from deSolve desired. Default "vode".}

\item{...}{Any other arguments to be passed to ode().}
}
\value{
Matrix of values for fluxes, biomass, and host growth rate at explicitly desired time values.
}
\description{
Function solveCoral uses the deSolve package to generate values similar to coRal (Cunning). Is essentially a transcription of Ferdinand Pfab's Mathematica code to R.
}
\examples{
solveCoral()
solveCoral(times = seq(0,365,0.1), pars = defPars(), lambda = 10, atol = 0.01, rtol = 0.01)
}
\seealso{
\code{\link{defPars}}, \code{\link{initState}}
}
