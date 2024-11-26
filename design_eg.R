#######
## Here's some code for Design calculations
## Assumes that the following exists, eg from 'fit_walrus_ckmR_24_11.R'
# ptru, Vpar, denv (which contains sample size info)
# and various functions

# Quantities of interest:
stuff0 <- interesting_stuff( ptru)

# Just popdyn is quick (probs are slow) so don't parallel this
Dstuff0 <- numvbderiv( interesting_stuff, ptru, eps=1e-5)
Vstuff0 <- Dstuff0 %**% Vpar %**% t( Dstuff0)
SEstuff0 <- sqrt( diag( Vstuff0))
mpt0 <- make_precision_table( stuff0, SEstuff0, E_comp)

# Ready for import into Lyx--- which works pretty well
if( !dir.exists( 'results')){
  mkdir( results) # MVB 8/11--- tho might wanna put elsewhere anyway
}
write.csv( mpt0, file='results/out.csv')

# MVB 8/11: I have turned the block below into "unexecuted" rather than
# hashed-out, so that I can test it (even if we don't want it yet)
if( FALSE){
  ## Investigate contrib of CKMR data
  # Keep a bit of CKMR for identifiability
  H_comp_IMR <- list() 
  # ... the "comp" stands for "component", not "comparison"
  
  # This loops over the different lglk components
  # (in this case MOP, HSP, SelfP ... but I like writing generic code
  # although others don't always like reading it! ;)
  for( comptype in names( Hbits$DSP)){
    ncomp <- denv[[ sub( 'DSP', 'n_comp', comptype)]]
    if( !grepl( '(?i)self', comptype)){
      ncomp <- ncomp / 100 # downweight CKMR 
    }
    H_comp_IMR[[ comptype]] <- finfo_onetype( Hbits$DSP[[ comptype]], 
      ncomp)
  }

  # ... and combine
  H_IMR <- 0*H_comp_IMR[[1]]
  for( comptype in names( Hbits$DSP)){
    H_IMR <- H_IMR + H_comp_IMR[[ comptype]]
  }

  # ... with hindsight, I could just have applied the 
  # downweighting inside the last loop, without recalling finfo_onetype
  # Never Mind.

  Vnock <- solve( H_IMR)
  Vstuff_nock <- Dstuff0 %**% Vnock %**% t( Dstuff0)
  mpt_nock <- make_precision_table( stuff0, sqrt( diag( Vstuff_nock)), E_comp) 

  ## Historical & future samp sizes match the sims (I hope...)

  ## Let's tweak juve vs ad samp sizes in future...
  # More juves...
  mpt_mj <- make_precision_table_completely( 
      stuff0, Dstuff0, Hbits, denv, 
      adjust_future_sample_sel( denv, 2))
  Vmj <- attr( mpt_mj, 'Vpar')
  # sqrt( diag( Vmj)) is mostly _slightly_ worse than for Vpar
  # except juve surv, unsurprisingly! but improvement is small there

  # Less juves...
  mpt_lj <- make_precision_table_completely( 
      stuff0, Dstuff0, Hbits, denv, 
      adjust_future_sample_sel( denv, 0.5))
  Vlj <- attr( mpt_lj, 'Vpar')

  # Splunge the mpt's together...
  mpt_all <- rbind( mpt0[2,], mpt_mj[2,], mpt_lj[2,])
  Ecols <- which( startsWith( colnames( mpt0), 'E')) 
  mpt_all[ , Ecols] <- rbind( mpt0[1,Ecols], mpt_mj[1,Ecols], mpt_lj[1,Ecols])
}  
