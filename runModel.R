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
eq <- function(x, aim) {lambda * (aim - x)} # dy.dt = lambda * (aim - y) | aim = f(y_1, y_2, ...)

# Parameters
depars <- append(unlist(defPars(nsym = nsym)), c(env))

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
    c(# Initial Model Values
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
      jST = jST,
      # Host and Symbiont biomass
      H = H,
      S = S
    )}

# function for daspk() func argument, returns list, where first element is a
# vector containing right hand sides of DEs and the residual forms of the
# algebraic equations
coral <- function(t, y, parameters) {
  pars <- parameters
  # Constant values, including as states seemed to make the solver unhappy
  consts <- list(jX = mmk(pars$X, pars$KX, pars$jXm), 
                 jN = mmk(pars$N, pars$KN, pars$jNm),
                 jHT = pars$jHT0,
                 rNH = pars$jHT0 * pars$nNH * pars$sigmaNH,
                 rNS = pars$jST0 * pars$nNS * pars$sigmaNS)
  
  with(as.list(c(y, pars, consts)), {
    list(c(# DEs
         eq(jHG, synth(jHGm, yC * rhoC * S / H + jX, (jN + nNX * jX + rNH) / nNH)),  # jHG
         eq(rhoN, max(0, jN + nNX * jX + rNH - nNH * jHG / yC)),  # rhoN
         eq(jeC, max(0, jX + rhoC * S / H - jHG / yC)),  # jeC
         eq(jCO2, kCO2 * jeC),  # jCO2
         eq(jL, (1.26 + 1.39 * exp(-6.48 * S / H)) * L * astar),  # jL
         eq(rCH, sigmaCH * (jHT + (1 - yC) * jHG / yC)),  # rCH
         eq(rCS, sigmaCS * (jST0 + (1-yC) * jSG / yC)),  # rCS
         eq(jCP, synth(jCPm, yCL * jL, (jCO2 + rCH) * H / S + rCS) / cROS),  # X.jCP
         eq(jeL, max(0, jL - jCP / yCL)),  # X.jeL
         eq(jNPQ, 1 / (1 / kNPQ + 1 / jeL)),  # X.jNPQ
         eq(cROS, (1 + max(0, jeL - jNPQ) / kROS)),  # X.cROS
         eq(jSG, synth(jSGm, yC * jCP, (rhoN * H / S + rNS) / nNS)),  # X.jSG
         eq(rhoC, max(0, jCP - jSG / yC)),  # X.rhoC
         eq(jST, jST0 * (1 + b * (cROS - 1))),  # X.jST
         (jHG - jHT) * H,  # H
         (jSG - jST) * S  # S
  ))})
}

times <- seq(0, 365, 0.1)  # Sequence of times for run

# Run the model
system.time(outv <- ode(y = destate, times = times, func = coral, parms = depars, verbose = T))
system.time(out <- ode(y = destate, times = times, func = coral, parms = depars, verbose = F))
plot(out[,1], out[,17]/out[,16], "l", lwd = 2, main = "S:H Biomass")

#For verification against original package coRal
library(coRal)
system.time(master <- run_coral(times, init_env(times, c(30, 30, 0), c(1e-6, 1e-6, 0), c(1e-7, 1e-7, 0)), def_pars(1)))
plot(times, out[,17]/out[,16] - master$S/master$H, "l", lwd = 2)
