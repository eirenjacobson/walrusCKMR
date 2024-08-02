"add_data" <- function( lglk= NULL,
                        indiv= NULL,
                        qk= NULL,
                        simfile= "WalrusSamples.RData",
                        simdir= "./simulation",
                        AMAX= 37, # and then you die
                        FIRST_PDYEAR= 2000,
                        SYEARS= if( missing( simfile)) 2010:2028,
                        AMAT= 6,
                        AMIN= 1, # calves not modelled 
                        MAX_WEANAGE= AMAT-1,
                        MAX_BGAP_HSPS= 2*AMAT + 2, # GGPs unlikely within this
                        FEC_ASYMP_AGE= 15,
                        MAXAGERR= 0, # bold!
                        Nursingness_fudge_ppn= 0.26, # ppn_breedy (this from ptru)
                        ... # just to allow ending all previous arglines with commas
){
  # Auto-trimming is good (remove calves)
  op <- options( offarray_table_warn_trim=FALSE)
  on.exit( options( op))
  
  if( ...length()){
    stop( sprintf( "Unknown args: %s", 
                   paste( ...names(), collapse=', ')))
  }
  
  # Some "constants", independent of sampling
  SEXES <- cq( F, M)
  NURSINGNESS <- cq( Yes, No)
  LETHALITY <- cq(LETHAL, NONLETHAL)
  
  # Feb: abolished distintion between "offspring" and "adult" ages.
  # Anything goes now, as long as it fits into PDYEARS
  AGES <- AMIN:AMAX
  
  # For self-recap: at 1st cap time, all adult ages are equiv
  SELF1AGES <- AMIN:AMAT # AMAT serves as "plus group" for selfs
  
  # Sim data: might be either be in simfile, or previously loaded
  # or might be non-existent--- in which case, will make all-zero
  # dataset waaaay below here
  if( is.null( indiv) && !is.null( simfile)){
    if( !file.exists( simfile)){
      simfile <- file.path( simdir, simfile)
    }
    
    what <- load( simfile) # indiv, or samples, or ....
    if( 'indiv' %not.in% what){
      # work out what its name was...
      for( i in what){
        indiv <- get( i)
        if( (indiv %is.a% 'data.frame') &&
            all( indiv %has.name% cq( Me, Mum, Dad, SampY)))
          break
      }
    } # now we should have indiv
    # last_raw_indiv <<- indiv # naughty code, but useful...
  } # loading of "raw" data
  
  if( !is.null( indiv)){ # any data?
    # MVB fix: very important!!!
    indiv <- within( indiv, {
      Dad[ Dad=='founder'] <- NA
      Mum[ Mum=='founder'] <- NA
    })
    
    samples <- indiv[!is.na(indiv$SampY),]
    set_SYEARS_etc( range( samples$SampY))
    
    if( is.null( qk)){
      recaps <- which( duplicated( samples$Me)) # EKJ NOTE this won't happen with sim data in main sim file
                                                # bc simulation overwrites capture years. AllSamples has this info.
      nrecaps <- length( recaps)
      nuniq <- nrow( samples) - nrecaps
      if( nrecaps>0){
        # Re-order to put them at the end; simpler later
        samples <- samples[ c( seq_len( nuniq) %except% recaps, recaps),]
      }
      # Normal quickin() wants entire sim input, including unsampled. But 
      # can do POPs & sibs (not GGPs) just from samples, as below
      # qk <- quickin( samples[ seq_len( nuniq),], max_gen = 2, 
      #     cousins. = FALSE, thiatics. = FALSE)
      qk <- quickin_only_popsib( samples[ seq_len( nuniq),])
      
      if( nrecaps>0){
        # Add in POPs & HSPs for the recap too
        me_again <- samples$Me[ -seq_len( nuniq)]
        me_first <- match( me_again, samples$Me)
        
        qk$POP <- add_recap_kin( qk$POP, me_first, me_again)
        qk$MatHSP <- add_recap_kin( qk$MatHSP, me_first, me_again)
      }
      
      # last_qk <<- qk
    }
    
    # convert to Bravo format
    # MVB: I've replaced %>% with |> since pipes are now base-R
    # ... NB careful, AgeLast is at DeathY
    
    extract.named( 
      rename.els( 
        samples[ cq( Me, Sex, SampY, BirthY, Lethality)],
        BirthY='B', SampY='Y', Lethality = 'L')
    ) # Me, Sex, Y, B, L now available
    
    if( samples %has.name% "Nursing"){
      Nursy <- samples$Nursing
      Nursingness_fudge_ppn <- NULL # didn't use it
    } else {
      # ... that field doesn't exist; fudge it
      # Using age-constant here; not ideal, but...
      Nursy <- ifelse( 
        (Sex=='F') & (Y-B >= AMAT) * 
          (runif( length( Y)) < Nursingness_fudge_ppn),
        'Yes', 'No')
      # Actually... this won't do! A modest number of actual MOPs might
      # be incompatible with that (25 / 390, actually).
    } 
    
    # Age at sampling:
    A <- Y - B
    
    # dply alternative
    #Samps <- samples |>
    #  select(Me, Sex, SampY, AgeLast) |>
    #  rename(Y = SampY, A = AgeLast)
    # extract.named( Samps)
    
    # find *indices* of matching pairs
    
    MOPs <- matrix( match( qk$POP, Me), 
                    nrow( qk$POP), 2)
    
    XmHSPs <- matrix( match( qk$MatHSP, Me), 
                      nrow( qk$MatHSP), 2)
    
    
    # Make sure in order. Even if DNAge imperfect, can prolly tell
    # which-is-which. I think later code assumes the birth-order is correct
    swappo <- B[ MOPs[,1]] > B[ MOPs[,2]]
    MOPs[ swappo,] <- MOPs[ swappo, 2:1]
    swappo <- B[ XmHSPs[,1]] > B[ XmHSPs[,2]]
    XmHSPs[ swappo,] <- XmHSPs[ swappo, 2:1]
    
    # Distro of age in samples (meant for ageing-error stuff)
    # weirdly, not linked to YEAR of sampling in code
    # but I guess that's QEQM assumption consistent!
    # Note this is a shortcut -- doesn't deconvolve est age to true age
    Pr_a_sampled <- offarray( table( A))
    Pr_a_sampled[] <- Pr_a_sampled / sum( Pr_a_sampled)
    
    # Offspring considered: only non-calfs born post 2000
    offposs <- (A>0) & (B >= FIRST_PDYEAR)
    # Parents: gotta mature by end of sampling
    # (this is just for efficiency)
    parposs <- (A>0) & (B < (max(PDYEARS) - AMAT))
    
    # Elim founders... shouldn't be any!!!
    print( sum( offposs & is.na( samples$Dad))) # 0, whew
    offposs <- offposs & !is.na( samples$Dad) # in case
    
    # Filter MOPs & XmHSPs based only on pop-dyn year-range *for now*
    MOPs <- MOPs[parposs[MOPs[,1]] & offposs[MOPs[,2]],]
    XmHSPs <- XmHSPs[offposs[ XmHSPs[,1]] & offposs[ XmHSPs[,2]],]
    
    if( !( samples %has.name% "Nursing")){
      r"--{
      Adjust Nursy so we don't get impossible gaps :/
      This is ghastly!!! Sims can easily add Nursy in future, but for now...
       - Can't be Nursy 2 years in a row
       - Must have been Nursy when that offspring was born
      }--"
      
      gap <- Y[ MOPs[,1]] - B[MOPs[,2]]
      Nursy[ MOPs[ abs( gap)==1,1] ] <- "No"
      Nursy[ MOPs[ gap==0,1] ] <- "Yes"
    }
    
    ## Self-recaps: done from separate file
    # NB that SS is _slightly_ different in the separate file
    # cos recaps are excluded from the original. But it's <2% so use orig
    # More consistent would be to use recapees *twice* in POP/HSP checks
    # but not worth fixing for now
    # Look for file in same folder as simfile
    allsimfile <- file.path( dirname( simfile),
                             'All' %&% basename( simfile))
    what <- load( allsimfile)        
    if( 'allsamples' %not.in% what){
      # work out what its name was...
      for( i in what){
        allsamples <- get( i)
        if(( allsamples %is.a% 'data.frame') &&
           all( allsamples %has.name% cq( Me, Mum, Dad, SampY)))
          break
      }
    } # now we should have allsamples
    
    allsamples <- rename.els( 
      allsamples[ cq( Me, Sex, SampY, BirthY, Lethality)],
      BirthY='B', SampY='Y', Lethality = 'L')
    allsamples$A <- with( allsamples, Y-B)
    # Females only, no calves
    allsamples <- allsamples %where% ((Sex=='F') & (A>0))
    
    # For Xtuple recaps (doesn't happen in "base data")
    # just use 1st and last...
    metab <- table( allsamples$Me)
    threecaps <- names( metab[ metab>2])
    if( length( threecaps)){
      scrungio <- which( allsamples$Me %in% threecaps)
      blurk <- allsamples$Me[ scrungio]
      firsto <- match( threecaps, blurk)
      lasto <- length( blurk) + 1 - match( threecaps, rev( blurk))
      middlo <- scrungio[ -c( firsto, lasto)]
      if( length( middlo)){
        allsamples <- allsamples[ -middlo,]
      }
    }
    
    # At most one recap each
    recaps <- names( metab[ metab==2])
    scrungio <- which( allsamples$Me %in% recaps)
    blurk <- allsamples$Me[ scrungio]
    firsto <- match( recaps, blurk)
    lasto <- length( blurk) + 1 - match( recaps, rev( blurk))
    
    selfPs <- cbind( scrungio[ firsto], scrungio[ lasto])
    
    ## SAMPLE SIZES: 
    # m_subset_xyz is samp size of subset-type samples with covars x,y,z
    # if subset omitted, it's all samples.
    # Now no attempt to pre-decide "poss par" and/or "poss off"
    # those will be incorped at the n_comp... stage
    # Data types will be live/dead, Male/Female, accompinfo/not
    # All as Y*A array
    # 'dimseq' arg in case ends-of-ranges don't have any samples
    

    m_SYANL <- offarray( table( S=Sex, Y, A, N=Nursy, L), 
                             dimseq=list( SEXES, SYEARS, AGES, NURSINGNESS, LETHALITY))
    
    m_F_YANL <- m_SYANL[ SLICE='F',,,,]
    
    m_F_YAL <- sumover( m_F_YANL, 'N')
    
    m_YA <- sumover(m_SYANL, cq( S, N, L)) 
    
    # For selfPs (noting comment above about slightly wrong SSize)
    m_F_YAselfL <- m_F_YAL[,SELF1AGES,]

    # ... aggregate all adults
    m_F_YAselfL[,AMAT,] <- sumover( m_F_YAL[,AMAT:AMAX,], 'A')
    m_F_YAself <- sumover(m_F_YAselfL, 'L')

    # note don't need to add a lethality dimension here
    # because this only applies to the second sample
    # and it doesn't matter bc make_n_comps 
    m_F_YD <- autoloop(
      Y=SYEARS, D=DEVSTAGES, SUMOVER=list( A=SELF1AGES),
      m_F_YAself[ Y, A] * (Devstage_A[ A]==D)
    )
    
  } else { # no sim data # EKJ didn't update this 
    set_SYEARS_etc( SYEARS) # perhaps from default
    MOPS <- mHSPs <- matrix( 0L, 0, 2)
    m_YA <- mM_live_YA <- 
      mF_live_YA <- offarray( 0, dimseq= list( SYEARS, AGES))
    Pr_a_sampled <- .....
  }
  
  ## Numbers of comparisons and kin-pairs
  
  r"--{
  Sometimes we will not use results of certain kin-comps. 
  EG POP when caught in same year (to avoid unweaned offspring being sampled near mother); 
  HSP when birth-gap is long enuf that GGP becomes likely. 
  These were _not_ necessarily filtered out when MOPs and HSPs were made earlier.
  }--"
  
  # create offarrays n_comps_<blah>
  extract.named( make_n_comps(
    m_YA,
    m_F_YANL,
    m_F_YAself, 
    m_F_YD,
    FIRST_PDYEAR,
    AMAT,
    MAX_BGAP_HSPS
  )) 
  
  
  ## MOPs: comps first. No same-capture-year. 
  # Condition on Nursy first, then just sum over that for no-nursy
  # number of pairs
  n_MOP_AYNL <- offarray(table(
    aj = A[MOPs[,2]],
    yj = Y[MOPs[,2]],
    ac = A[MOPs[,1]],
    yc = Y[MOPs[,1]],
    nc = Nursy[MOPs[,1]], 
    lc = L[MOPs[,1]]
  ), template = n_comp_MOP_AYNL) # use these dims 
  
  # Zero-out unusable cases, even if they exist
  n_MOP_AYNL[ MATSUB = whichoff( n_comp_MOP_AYNL, .==0)] <- 0
  n_MOP_AYNL[ MATSUB = whichoff( n_comp_MOP_AYNL, .==0)] <- 0
  n_MOP_AYL <- sumover( n_MOP_AYNL, 'nc')

  ## HSPs... make it compatible with offposs etc
  n_XmHSP_AY <- offarray( table(
    a1= A[ XmHSPs[,1]],
    y1= Y[ XmHSPs[,1]],
    a2= A[ XmHSPs[,2]],
    y2= Y[ XmHSPs[,2]]
  ), template = n_comp_XmHSP_AY) # use these dims
  
  # Zero-out unusable cases, even if they exist
  n_XmHSP_AY[ MATSUB = whichoff( n_comp_XmHSP_AY, .==0)] <- 0
  
  ## selfPs: comps first.
  # Indices are different here
  
  # number of self-recaps. 
  # MUST use A etc from allsamples NOT general A! 
  n_selfP_YADY <- with( allsamples,  # CRUCIAL!
                          offarray(table(
                            y1 = Y[selfPs[,1]],
                            aself1 = pmin( A[selfPs[,1]], AMAT),
                            d2 = Devstage_A[ A[ selfPs[,2]], NOOFF=TRUE], # NOOFF stops warnings
                            y2 = Y[selfPs[,2]]
                          ), template = n_comp_selfP_YADY)) # use these dims
  
  # Zero-out unusable cases, even if they exist
  n_selfP_YADY[ MATSUB = whichoff( n_comp_selfP_YADY, .==0)] <- 0
  
  ## Ageing error "data"--- assumed exact for now,
  # but (unfinished) code allows generality
  expand_by_maxagerr <- function( ARANGE)
    max( min( ARANGE)-MAXAGERR, 1) %upto% (max( ARANGE) + MAXAGERR)
  
  Pr_agerr_a <- offarray( 0, first=c( -(1+MAXAGERR), 1),
                          last=c( 1+MAXAGERR, AMAX))
  Pr_agerr_a[ 0,] <- 1 # no agerr for now, regardless of MAXAGERR
  
  # User only gets to observe age-blurred results. For now, no agerr
  # ... so we can just do this:
  # n_MOP_REY <- n_MOP_RAY etc etc etc
  for( nrey in ls( pattern='^n_.*_(R)?AY(N)?$')){
    assign( sub( 'AY', 'EY', nrey), get( nrey))
  }
  #  Selfs: gonna assume exact age up to AMAT
  
  env <- if( is.null( lglk)) .GlobalEnv else environment( lglk)
  env <- new.env( parent=env)
  e <- list2env( envir=env, returnList( 
    n_comp_MOP_AYL,
    n_MOP_AYL,
    n_comp_MOP_AYNL, 
    n_MOP_AYNL,
    n_comp_XmHSP_AY,
    n_XmHSP_AY,
    n_comp_selfP_YADY,
    n_selfP_YADY,
    
    m_SYANL,
    m_F_YAL,
    m_F_YANL,
    m_YA,
    m_F_YAself,
    m_F_YD,
    
    Pr_a_sampled, # samp age compo
    
    SEXES,
    LETHALITY,
    FIRST_PDYEAR,
    PDYEARS,
    MAX_Db= diff( range( PDYEARS)), # birth gap HSPs
    YSTART= min( PDYEARS),
    
    SYEARS,
    MAX_Dsyears= diff( range( SYEARS)), # recap gap selfPs
    
    MAXAGERR,
    Pr_agerr_a,
    AMIN,
    AMAX,
    AMAT,
    FEC_ASYMP_AGE,
    FECAGES= (AMAT-1):FEC_ASYMP_AGE,
    MAX_WEANAGE,
    MAX_BGAP_HSPS,
    MAX_popgap = 50, # who knows..?
    AGES,
    SELF1AGES,
    DEVSTAGES,
    Devstage_A,
    NURSINGNESS,
    
    # Slightly wider estimated-age-range
    EAGES= expand_by_maxagerr( AGES),
    SELFEAGES= SELF1AGES, # not allowing error here... yet...
    
    PCOHORTS,
    FIRST_PCOHORT= PCOHORTS[1],
    LAST_PCOHORT= tail( PCOHORTS, 1),
    
    MAX_Tsep= 100 # excessive, but too hard to precalc for now..!
  ))
  
  if( !is.null( lglk)){
    copy_lglk <- lglk # actually unnecessary, but clearer...
    environment( copy_lglk) <- e
    return( copy_lglk)
  } else {
    return( e)
  }
  
  
  
} # end function add_data
