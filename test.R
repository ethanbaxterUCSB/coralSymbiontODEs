
{
  otherFunc <- function(a, b) {a + b}
  parameters <- c(myPar <- 2, X0 = c(1,1))
  state <- with(as.list(parameters), {c(X = c(jX = 1, hX = 1))})
  myFunc <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dX <- c(myPar, otherFunc(X.hX, X.jX))
      return(list(c(dX)))
    })
  }
  times <- seq(0, 100, 0.1)
  ode(y = state, times = times, func = myFunc, parms = parameters)
}
