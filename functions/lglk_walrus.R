"lglk_walrus" <- function( 
    params, 
    want=c( 'lglk', 'probs', 'popdyn_only'),
    USE_NURSY_ONLY= FALSE 
){
  ##
  
  r"--{
USE_NURSY_ONLY: 
 TRUE => ignore non-nursy MOPs, which get prepped from sims by default as duplicates of the nursy ones
 NA => use some non-nursy and some nursy MOPs, if user has split dataset appropriately
 FALSE => only use nursy MOPs
}--"
  
  if( MAXAGERR>0){
    stop( "Ageing error not fully implemented yet")
  }
  
  ## First needs to have its envir set, to know about DATA
  # Pretty lazy partial test!
  stopifnot( exists( 'PDYEARS', environment( sys.function()), 
                     inherits=FALSE))
  
  ## Unpack parameters
  # get em one (or more) at a time
  # mitigates counting errors!
  iparam <- 0
  next_param <- function( n=1) {
    iparam <<- iparam + n
    unname( params[ seq( to=iparam, length=n)])
  }
  
  Nfad_ystart <- exp( next_param())
  RoI <- next_param()
  
  lgt_fadsurv <- next_param()
  diff_lgt_fjusurv <- next_param() # 0 means same surv as adult
  fadsurv <- inv.logit( lgt_fadsurv)
  fjusurv <- inv.logit( lgt_fadsurv + diff_lgt_fjusurv)
  
  # breeding cycle
  psicle <- offarray( inv.logit( next_param( 2)),
                      first=2, last=3)
  
  assign( 'last_params', params, environment( sys.function())) # debugging etc
  
  # Do what? Might return early...  
  want <- match.arg( want)
  
  ## General setup (boring housekeeping) if any
  # Well, there isn't any; all done in popdyn, or 
  # previously in add_data
  # so don't do it:
  # setup_walrus()

  ## Pop dyn
  popdyn_walrus()
  if( want=='popdyn_only'){
    return( returnList(
      Nfad_y,
      Nfju_y,
      Pbreedf,
      ppn_breedy,
      recip_ppn_breedy,
      fec_a,
      
      real_params= unlist( returnList( 
        Nfad_ystart,
        RoI,
        fadsurv,
        fjusurv,
        psicle
      ))
    ))}
  
 # movement_walrus() # even if only 1 region!
  
  ## Ideal probs (if we knew everything)
  # although these DO incorporate space
  ideal_MOP_walrus()
  ideal_XmHSP_walrus() # and ignore GGPs--- restrict comps instead
  ideal_selfP_walrus()
  
  ## Probs given observables
  age_blur_probs_walrus()
  
  # Housekeeping: useful to keep some stuff after function exits. 
  # Trust me, this works...
  ## WRONG VARS FOR WALRUS, OF COURSE!
  list2env( mget( cq( 
    fadsurv, psicle, Nfad_ystart, RoI, Nfad_y, fec_a)),
    envir=environment( sys.function()))
  
  returnables <- list()
  if( want=='probs'){
    # Only OBSERVABLE-conditioned ones, to avoid confusing despack
    returnables <- returnList(
     # Pr_MOP_EY, only applicable when max age error > 0?
      Pr_MOP_EYEYL,
      Pr_MOP_EYEYNL,
      Pr_XmHSP_EY,
      Pr_selfP_YADY
    )
    
    if( isT( USE_NURSY_ONLY)){
      returnables$Pr_MOP_EY <- NULL
    } else if( isF( USE_NURSY_ONLY)){
      returnables$Pr_MOP_EYN <- NULL
    }
  }
  # lglk (if there's any data)
  if( length( n_XmHSP_AY)>0){
    # May have some MOPs with Nursy, some without
    # Up to user to NOT duplicate!!
    lglk <- 0
    
    if( !isT( USE_NURSY_ONLY)){
      lglk <- lglk + ldpois( n_MOP_AYL, 
                             n_comp_MOP_AYL * Pr_MOP_EYEYL) # 
    }
    if( !isF( USE_NURSY_ONLY)){
      lglk <- lglk + ldpois( n_MOP_AYNL, 
                             n_comp_MOP_AYNL * Pr_MOP_EYEYNL) # 
    }
    lglk <- lglk + ldpois( n_XmHSP_AY, 
                           n_comp_XmHSP_AY * Pr_XmHSP_EY) # works
    
    lglk <- lglk + ldpois( n_selfP_YADY, 
                           n_comp_selfP_YADY * Pr_selfP_YADY) # works. maybe this should be EY, and I should make n_comps etc to match?
  } else {
    lglk <- 0
  }
  
  returnables$lglk <- c( lglk) # no dims please
  return( returnables)
}