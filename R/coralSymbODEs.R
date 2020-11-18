#' Solve coral-\emph{Symbiodinium} ODEs
#'
#' Function solveCoral uses the deSolve package to generate values similar to coRal (Cunning). Is essentially a transcription of Ferdinand Pfab's Mathematica code to R.
#'
#' @param times The time values at which output is desired. Default c(0,500).
#' @param pars The paramaters of the model, given as a list. See function \code{\link{defPars}}. Default defPars().
#' @param lambda The sensitivity of the runs. High values are more sensitive and small values are less sensitive. Default 5.
#' @param method The character method argument of ode() from deSolve desired. Default "vode".
#' @param ... Any other arguments to be passed to ode().
#' @return Matrix of values for fluxes, biomass, and host growth rate at explicitly desired time values.
#' @examples
#' solveCoral()
#' solveCoral(times = seq(0,365,0.1), pars = defPars(), lambda = 10, atol = 0.01, rtol = 0.01)
#' @seealso \code{\link{defPars}}, \code{\link{initState}}
#' @export
solveCoral <- function(times = c(0,500), pars = defPars(), lambda = 5, method = "vode", ...) {
  # Constant fluxes
  consts <- with(as.list(pars), {
    c(jX = mmk(X, KX, jXm),
      jN = mmk(N, KN, jNm),
      jHT = jHT0,
      rNH = jHT0 * nNH * sigmaNH,
      rNS = jST0 * nNS * sigmaNS)
  })
  # Solve system and return
  return(cbind(ode(y = initState(pars), times = times, func = coralODEs, parms = append(pars, c(lambda = lambda, consts)), method = method, ...),consts))
}  # End function solveCoral

# Helper methods
# ==============

# System of flux ODEs of the form dy.dt = lambda * (f(y_1, y_2, ...) - y)
coralODEs <- function(t, y, parameters) {
  return(with(as.list(c(y, parameters)), {
    list(c(# DEs
      eq(jHG, synth(jHGm, yC * rhoC * S / H + jX, (jN + nNX * jX + rNH) / nNH), lambda = lambda),  # jHG
      eq(rhoN, max(0, jN + nNX * jX + rNH - nNH * jHG / yC), lambda),  # rhoN
      eq(jeC, max(0, jX + rhoC * S / H - jHG / yC), lambda),  # jeC
      eq(jCO2, kCO2 * jeC, lambda),  # jCO2
      eq(jL, (1.26 + 1.39 * exp(-6.48 * S / H)) * L * astar, lambda),  # jL
      eq(rCH, sigmaCH * (jHT + (1 - yC) * jHG / yC), lambda),  # rCH
      eq(rCS, sigmaCS * (jST0 + (1-yC) * jSG / yC), lambda),  # rCS
      eq(jCP, synth(jCPm, yCL * jL, (jCO2 + rCH) * H / S + rCS) / cROS, lambda),  # X.jCP
      eq(jeL, max(0, jL - jCP / yCL), lambda),  # X.jeL
      eq(jNPQ, 1 / (1 / kNPQ + 1 / jeL), lambda),  # X.jNPQ
      eq(cROS, (1 + max(0, jeL - jNPQ) / kROS), lambda),  # X.cROS
      eq(jSG, synth(jSGm, yC * jCP, (rhoN * H / S + rNS) / nNS), lambda),  # X.jSG
      eq(rhoC, max(0, jCP - jSG / yC), lambda),  # X.rhoC
      eq(jST, jST0 * (1 + b * (cROS - 1)), lambda),  # X.jST
      (jHG - jHT) * H,  # H
      (jSG - jST) * S  # S
    ),
    c(dH.dt = (jHG - jHT) * H,dS.dt = (jSG - jST) * S))}))
}  # End function coralODEs

# Synthesis rate of a product given maximum rate m and substrates A and B.
synth <- function(m, A, B) {A * B * (A + B) * m / (A^2 * B + A * B^2 + A^2 * m + A * B * m + B^2 * m)}

# Uptake of a substrate A given half-saturation constant half and maximum rate max.
mmk <- function(A, half, max) {max * A / (A + half)}

# Standard form of the ODEs
eq <- function(x, aim, lambda) {lambda * (aim - x)} # dy.dt = lambda * (aim - y) | aim = f(y_1, y_2, ...)

#' Default parameters for solveCoral
#'
#' Helper function, based off of def_pars() from package coRal by Ross Cunning.
#' @examples defPars()
#' @return Named list of model parameters. Includes environment.
#' @seealso \code{\link{solveCoral}}
#' @export
defPars <- function() {
  return(c(
    jHT0=0.03,  # Host specific biomass turnover rate (d^-1)
    nNH=0.18,  # N:C ratio in host biomass (-)
    nNX=0.2,  # N:C ratio in prey biomass (-)
    sigmaNH=0.9,  # Proportion of host nitrogen turnover recycled (-)
    sigmaCH=0.1,  # Proportion of host carbon turnover recycled (-)
    jXm=0.13,  # Maximum specific host feeding rate (molX/CmolH/d)
    jNm=0.035,  # Maximum specific host DIN uptake rate (molN/CmolH/d)
    jHGm=1,  # Maximum specific host growth rate (CmolH/CmolH/d)
    jHG0 = 0.25,  # initial host growth rate
    jeC0 = 10,  # initial excess carbon flux
    kCO2=10,  # Rate of host CCM's (molCO2/molC/d)
    KN=1.5e-6,  # Half-saturation constant for host DIN uptake (molN/L)
    KX=1e-6,  # Half-saturation constant for host feeding (CmolX/L)
    initH=1,  # Initial host biomass (CmolH)
    yC=0.8,
    jST0=0.03,  # Symbiont specific biomass turnover rate (d^-1)
    nNS=0.13,  # N:C ratio in symbiont biomass (-)
    yCL=0.1,  # L:C ratio in fixed carbon (=quantum yield) (molC/mol ph)
    kNPQ=112,  # capacity of non-photochemical quenching (mol ph/CmolS/d)
    # calculated as 4x max. photochemical quenching (Gorbunov et al. 2001)
    kROS=80,  # amount of excess light beyond NPQ capacity (e.g., jeL-jNPQ) that doubles ROS production relative to baseline (mol ph/CmolS/d)
    cROS0 = 1,  # capacity for NPQ
    k=1,  # exponent on ROS production (-)
    astar=1.34,  # Symbiont specific cross-sectional area (m^2/C-molS)
    sigmaNS=0.9,  # Proportion of symbiont nitrogen turnover recylced (-)
    sigmaCS=0.9,  # Proportion of symbiont carbon turnover recycled (-)
    jCPm=2.8,  # Maximum specific photosynthate production rate (Cmol/CmolS/d)
    jSGm=0.25,  # Maximum specific symbiont growth rate (CmolS/CmolS/d)
    initS=1,  # Initial symbiont biomass (CmolS)
    b=5,  # Scaling parameter for bleaching response
    L=30,  # Environmental light (mol photons)
    N=1e-6,  # Environmental DIN (mol N)
    X=1e-7  # Environmental prey (mol X)
  ))
}  # End function defpars

#' Create initial state for solveCoral
#'
#' Helper function, based off of run_coral() from package coRal.
#' @param pars The named list of model parameters.
#' @param env The named list of environmental values.
#' @return A named numeric vector containing the initial flux and biomass values.
#' @seealso \code{\link{solveCoral}}, \code{\link{defPars}}
#' @examples initState(defPars())
#' @export
initState <- function(pars) {
  with(as.list(pars), {
    # Initial Host fluxes
    rhoN <- mmk(N, KN, jNm)
    jeC <- jeC0
    jCO2 <- kCO2 * jeC
    jHG <- jHG0
    rCH <- jHT0 * sigmaCH
    dH.Hdt <- jHGm
    H <- initH
    # Initial Symbiont fluxes
    jL <- L * astar
    jCP <- max(0, synth(jL * yCL, jCO2 * H / initS, jCPm))
    jeL <- max(0, jL - jCP / yCL)
    jNPQ <- kNPQ
    jSG <- jSGm/10
    rhoC <- jCP
    jST <- jST0
    rCS <- jST0 * sigmaCS
    cROS <- cROS0
    dS.Sdt <- jSGm
    S <- initS
    # Return named numeric vector
    return(c(# Initial flux Values
      jHG = jHG, rhoN = rhoN, jeC = jeC, jCO2 = jCO2, jL = jL, rCH = rCH,
      rCS = rCS, jCP = jCP, jeL = jeL, jNPQ = jNPQ, cROS = cROS, jSG = jSG,
      rhoC = rhoC, jST = jST,
      # Host and Symbiont biomass
      H = H,
      S = S))
  })
}  # End function initState

#' Wrapper function to coRal
#'
#' Has the same input and output as run_coral() from the coRal package
#' @export
coRal.solveCoralWrapper <- function(time, env, pars) {
  # time is equivalent to times
  # can only be used for constant env values at this point
  newPars <- append(pars, c(L = env[['L']][[1]],N = env[['N']][[1]],X = env[['X']][[1]]))
  run <- solveCoral(time, newPars)
  with(as.data.frame(run),{
    out <- data.frame(time,env$L,env$N,env$X,jN,rhoN,jeC,jCO2,jHG,rCH,
                      dH.Hdt=dH.dt/H,H,rNS,jL,jCP,jeL,jNPQ,jSG,rhoC,jST,rCS,
                      cROS,dS.Sdt=dS.dt/S,S)
    out
  })
}
