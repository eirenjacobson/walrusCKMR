"doades" <- function( 
  m_SYEL,     # samp sizes
  Hbits,      # from get_Hbits
  reporteez,  # function that returns quants-of-interest
  lglk,       # prolly needed by reporteez, plus its env has "metadata"
  eps= 1e-5,  # for numvbderiv
  comp_wts= NULL
){
  ptru <- Hbits$PARS_FOR_H
  npars <- length( ptru)

  # Get param covar mat for this sampling scheme
  denv <- environment( lglk) # various useful things live here
  m_breakdown <- get_usable_SS( m_SYEL, denv$Devstage_A)
  
  # Other things required by make_n_comps:
  # IMO this is ugly code (why is it these 3 random-seeming quantities?!)
  # but getting it all generic and interdependent seems to 
  # decrease readability (at least in this case)
  ncomp_list <- make_n_comps( m_breakdown, 
      FIRST_PDYEAR= denv$FIRST_PDYEAR, 
      AMAT= denv$AMAT, 
      MAX_BGAP_HSPS= denv$MAX_BGAP_HSPS)
  
  H <- matrix( 0, npars, npars)
  
  # Loop over the different lglk components
  # (in this case MOP, HSP, SelfP) ... but I like writing generic code

  default_comp_wts <- rep( 1, length( Hbits$DSP))
  names( default_comp_wts) <- names( Hbits$DSP)
  
  if( !is.null( comp_wts)){ 
    # Allow user to specify just some of them, eg 
    # comp_wts= c( MOP=0)
    # Look for the right comp-type(s) inside name 
    # in the middle of the DSP name
    # dropping "s" for plural
    # I should know better than to write stuff like this...
    whicho <- match( sub( 's$', '', sprintf( '_%s', names( comp_wts))),
        sub( '^[^_]*(_[^_s]+)_.*', '\\1', names( Hbits$DSP)), 0)
    if( !all( whicho>0)){
stop( sprintf( "Don't understand which comp to weight: %s", 
        paste( names( comp_wts)[ whicho==0], collapse=', ')))}
    default_comp_wts[ whicho] <- comp_wts
  }
  comp_wts <- default_comp_wts
  
  for( comptype in names( Hbits$DSP)){
    this_ncomp <- ncomp_list[[ sub( 'DSP', 'n_comp', comptype)]]

stopifnot( my.all.equal( # dimensions had better match...
        unname( dimseq( Hbits$DSP[[ comptype]])[-1]), # 1st is params
        unname( dimseq( this_ncomp))
    ))
        
    this_H <- finfo_onetype( Hbits$DSP[[ comptype]], 
      this_ncomp)
    H <- H + this_H * comp_wts[ comptype]
  }
  # Other H-contributors, eg RE priors or general priors, should go here
  # H <- H + other_stuff
  
  Vpar <- solve( H)

  # Quants-of-interest: var via Delta-method
  # Just popdyn is quick (probs are slow) so don't parallel this
  stuff <- reporteez( ptru, lglk)
  
  Dstuff <- numvbderiv( reporteez, ptru, lglk=lglk, eps=1e-5)
  Vstuff <- Dstuff %**% Vpar %**% t( Dstuff)
  SEstuff <- sqrt( diag( Vstuff))

  SEstuff@V <- Vstuff
  
return( SEstuff)
}
