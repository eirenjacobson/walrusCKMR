"truncsurv" <-
function( AMAX,  phi){
  # truncsurv ( 1, 0.9) # check: yep, should be ~0.5
  # truncsurv( 30, 0.9622) # 0.9458
  N <- phi ^ (0:AMAX)
  surv <- phi + 0*N
  surv[ length( surv)] <- 0
  N %**% surv / sum(N)
}
<bytecode: 0x000002212e243ac8>
