"truncsurv" <- function( AMAX,  phi){
## If adults survive AMAX years beyond maturity, with a surv of phi
## up til then, at which point it's 0: what would the equiv surv of
## a population that had no cut-off death rate?

  # truncsurv ( 1, 0.9) # check: yep, should be ~0.5
  # truncsurv( 30, 0.9622) # 0.9458
  N <- phi ^ (0:AMAX)
  surv <- phi + 0*N
  surv[ length( surv)] <- 0
  N %**% surv / sum(N)
}
