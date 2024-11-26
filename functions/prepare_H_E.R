prepare_H_E <- function( denv, Hbits){
  E <- list() # by pairwise comp type, eg MOP

  # H contains expected Hessian (d^2 lglk / dparams^2)
  # ie inverse covariance between parameters
  # Add it up from all lglk components
  
  # MOPs
  ncomp <- denv$n_comp_MOP_AYNL 
  H <- finfo_onetype(Hbits$DSP$DSP_MOP_EYEYNL, ncomp)
  E$n_MOP_EYEYNL <- ncomp * Hbits$Prkin$Pr_MOP_EYEYNL


  # Half sibs
  ncomp <- denv$n_comp_XmHSP_AY
  H <- H + finfo_onetype( Hbits$DSP$DSP_XmHSP_EY, ncomp)
  E$n_XmHSP_EY <- ncomp * Hbits$Prkin$Pr_XmHSP_EY

  # self
  ncomp <- denv$n_comp_selfP_YADY
  H <- H + finfo_onetype( Hbits$DSP$DSP_selfP_YADY, ncomp)
  E$n_selfP_YADY <- ncomp * Hbits$Prkin$Pr_selfP_YADY

return( list( H=H, E=E))
}


