#'  Transcription of Ferdi's ODE model in Mathematica to R.
#'  Requires package deSolve, which is installed by install.packages("deSolve"), 
#'  and uses a modified version of coRal's def_pars() method.
library(deSolve)
source("defPars.R")

nsym <- 1 # Number of symbionts, limit to 1 for now
lambda <- 5 # coefficient in Ferdi's model
env <- list(L = 30, N = 1e-6, X = 1e-7) # Constant Environment

# Helper Methods
synth <- function(m, A, B) {A * B * (A + B) * m / (A^2 * B + A * B^2 + A^2 * m + A * B * m + B^2 * m)}
mmk <- function(A, half, max) {max * A / (A + half)}
eq <- function(x, aim) {lambda * (aim - x)} # dX.dt = lambda * (fX - X)

# Parameters
depars <- append(unlist(defPars(nsym = 1)), c(env))

# Vector of state variables and their initial values
destate <- {
    pars <- defPars(1)  # Parameters of the model
    # Initial Host fluxes
    rhoN <- mmk(env$N, pars$KN, pars$jNm)
    jeC <- 10
    jCO2 <- pars$kCO2 * jeC
    jHG <- 0.25
    rCH <- pars$jHT0 * pars$sigmaCH
    dH.Hdt <- pars$jHGm
    H <- pars$initH
    # Initial Symbiont fluxes
    jL <- env$L * pars$astar
    jCP <- max(0, synth(jL * pars$yCL, jCO2*H/pars$initS, pars$jCPm))
    jeL <- max(0, jL - jCP/pars$yCL)
    jNPQ <- pars$kNPQ
    jSG <- pars$jSGm/10
    rhoC <- jCP
    jST <- pars$jST0
    rCS <- pars$jST0 * pars$sigmaCS
    cROS <- 1
    dS.Sdt <- pars$jSGm
    S <- pars$initS
    c(
      # Initial Model Values
      X = c(
        jHG = jHG,
        rhoN = rhoN,
        jeC = jeC,
        jCO2 = jCO2,
        jL = jL,
        rCH = rCH,
        rCS = rCS,
        jCP = jCP,
        jeL = jeL,
        jNPQ = jNPQ,
        cROS = cROS,
        jSG = jSG,
        rhoC = rhoC,
        jST = jST
      ),
      # Host and Symbiont biomass
      H = H,
      S = S,
      # Initial Aims
      fX = c(
        jHG = jHG,
        rhoN = rhoN,
        jeC = jeC,
        jCO2 = jCO2,
        jL = jL,
        rCH = rCH,
        rCS = rCS,
        jCP = jCP,
        jeL = jeL,
        jNPQ = jNPQ,
        cROS = cROS,
        jSG = jSG,
        rhoC = rhoC,
        jST = jST
      )
    )}
# initd, initial rates of change, is only needed if using res argument of daspk()
# initd <- rep(0.1, 30)
# names(initd) <- paste0("d", names(destate))

# function for daspk() func argument, returns list, where first element is a
# vector containing right hand sides of DEs and the residual forms of the
# algebraic equations
coral <- function(t, y, parms) {
  # Constant values, including as states seemed to make the solver unhappy
  consts <- list(jX = mmk(parms$X, parms$KX, parms$jXm), 
                 jN = mmk(parms$N, parms$KN, parms$jNm),
                 jHT = parms$jHT0,
                 rNH = parms$jHT0 * parms$nNH * parms$sigmaNH,
                 rNS = parms$jST0 * parms$nNS * parms$sigmaNS)
  
  with(as.list(c(y, parms, consts)), {
    list(c(# DEs
         eq(X.jHG, fX.jHG),  # X.jHG
         eq(X.rhoN, fX.rhoN),  # X.rhoN
         eq(X.jeC, fX.jeC),  # X.jeC
         eq(X.jCO2, fX.jCO2),  # X.jCO2
         eq(X.jL, fX.jL),  # X.jL
         eq(X.rCH, fX.rCH),  # X.rCH
         eq(X.rCS, fX.rCS),  # X.rCS
         eq(X.jCP, fX.jCP),  # X.jCP
         eq(X.jeL, fX.jeL),  # X.jeL
         eq(X.jNPQ, fX.jNPQ),  # X.jNPQ
         eq(X.cROS, fX.cROS),  # X.cROS
         eq(X.jSG, fX.jSG),  # X.jSG
         eq(X.rhoC, fX.rhoC),  # X.rhoC
         eq(X.jST, fX.jST),  # X.jST
         (X.jHG - jHT) * H,  #X.jHG
         (X.jSG - X.jST) * S,  # X.jSG
         # Residuals, using consts jX, jN, jHT, rNH, rNS, and are of the form 0 = y - (calculation of y)
         fX.jHG - synth(jHGm, yC * fX.rhoC * S / H + jX, (jN + nNX * jX + rNH) / nNH),
         fX.rhoN - max(0, jN + nNX * jX + rNH - nNH * fX.jHG / yC),
         fX.jeC - max(0, jX + fX.rhoC * S / H - fX.jHG / yC),
         fX.jCO2 - kCO2 * fX.jeC,
         fX.jL - (1.26 + 1.39 * exp(-6.48 * S / H)) * L * astar,
         fX.rCH - sigmaCH * (jHT + (1 - yC) * fX.jHG / yC),
         fX.rCS - sigmaCS * (jST0 + (1-yC) * fX.jSG / yC),
         fX.jCP - synth(jCPm, yCL * fX.jL, (fX.jCO2 + fX.rCH) * H / S + fX.rCS) / fX.cROS,
         fX.jeL - max(0, fX.jL - fX.jCP / yCL),
         fX.jNPQ - 1 / (1 / kNPQ + 1 / fX.jeL),
         fX.cROS - (1 + max(0, fX.jeL - fX.jNPQ) / kROS),
         fX.jSG - synth(jSGm, yC * fX.jCP, (rhoN * H / S + rNS) / nNS),
         fX.rhoC - max(0, fX.jCP - fX.jSG / yC),
         fX.jST - jST0 * (1 + b * (fX.cROS - 1))
  ))})
}
# Mass matrix for DAE
mass <- diag(nrow = 30, ncol = 30)  # Create diagonal matrix of correct size
mass[17:30, 17:30] <- 0  # Set masses of algebraic expressions to 0

times <- seq(0, 100, 0.1)  # Sequence of times for run

# Run the model
out <- daspk(y = destate, times = times, func = coral, parms = depars, verbose = T, mass = mass)
plot(out[,1], out[,41]/out[,40], "l", lwd = 2, main = "S:H Biomass")
