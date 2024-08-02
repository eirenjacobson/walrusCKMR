"popdyn_walrus" <- function( nlocal=sys.parent()) mlocal({
  ## Population dynamics
  Nfad_y <- autoloop( Y=PDYEARS, {
    Nfad_ystart * exp( RoI * (Y-YSTART))
  })
  recip_Nfad_y <- 1/Nfad_y # I prefer calcing reciprocal once...
  
  # Juvenile version, assuming QEQM:
  # "Juves" are things included in Self comps (and probably as potential XHSPs and potential O's in MOPs, though those _could_ be something different I guess, at least if we had exact age)
  Nfju_y <- Nfad_y * 
    ((1+RoI-fadsurv) / (1-fjusurv)) *
    ((fjusurv ^ -(AMAT-AMIN)) - 1)
  Nf_yd <- offarray( 0, dimseq=list( PDYEARS, DEVSTAGES))
  Nf_yd[,'JuveF'] <- Nfju_y
  Nf_yd[,'AdF'] <- Nfad_y
  
  ## Breeding psicle
  # col= FROM, row= TO (Mark's version) so that T %*% oldprob = newprob
  # Birth happens when Mum is in state 2
  Tbreedf <- rbind( 
    c( 0, psicle[2], psicle[3]),
    c( 1, 0, 0),
    c( 0, 1-psicle[2], 1-psicle[3])
  )
  
  ## Ppn ad fem calving in any given year 
  
  # ... neglects ramping-up as maturity is reached
  # No longer going the eigenroute; get stationary eivec by...
  # ... solving eigeneq with eigenval 1
  M <- Tbreedf - diag( 3)
  M[3,] <- 1
  Pbreedf <- solve( M, c( 0,0,1)) # eqm distro of states
  # Whole (adult) popln, at any one time
  ppn_breedy <- Pbreedf[2] # tho *really* should be age-wted somehow I guess
  recip_ppn_breedy <- 1 / ppn_breedy  
  
  # ... and Do Something for TMB version! 
  # IIRC there's an explicit soln given form of Tbreedf. And if not,
  # then there's other ways, but they're harder.
  
  # Powers of the breed trans mat
  Pr_Bt_B0T <- offarray( 0, 
                         dimseq=list( Bt=1:3, B0=1:3, T=0:MAX_popgap))
  Pr_Bt_B0T[,,SLICE=0] <- diag( 3)
  for( T in 1:MAX_popgap){
    Pr_Bt_B0T[,,T] <- Tbreedf %**% Pr_Bt_B0T[,,SLICE=T-1]
  }
  
  ## Fec at age: start 2(?) yrs prior to AMAT in might-get-preg
  # max( FECAGES) serves as "mini plus group" here--- assumed const after
  # Construction is a bit ugly, coz gotta get offset right
  
  fec_a <- offarray( 
    c( Pr_Bt_B0T[ SLICE=2,SLICE=3,FECAGES-(AMAT-2)]),
    dimseq=list( FECAGES))
  fec_a <- fec_a * recip_ppn_breedy
  fec_a[ max( FECAGES)] <- 1 # force asymptote!
  
  ## Distro of rebreed prob, given mum still alive
  # This is not QUITE the same as Pr_w, since we forcibly
  # wean animals at MAX_WEANAGE, but mum might still not give birth
  # [sic] "This only applies at MAX_WEANAGE"
  Pr_breedyagain_Db <- Pr_Bt_B0T[ SLICE=2, SLICE=2, ]
  
  # FUTURE Bprob given state-at-sampling
  Pr_Breedyfuture_B0Tpos <- Pr_Bt_B0T[ SLICE=2,,]
  
  # Conditioned on NURSINGNESS (Yes => B0=2; No => B0 != 2)
  # I have just done this "the hard way" rather than via some
  # autoloop quasi-tensor-product thing
  Pr_Breedyfuture_NTpos <- offarray( 0, 
                                     dimseq= list( n=NURSINGNESS, T=0:MAX_popgap))
  Pr_Breedyfuture_NTpos[ 'Yes',] <- Pr_Breedyfuture_B0Tpos[ SLICE=2,]
  Pr_Breedyfuture_NTpos[ 'No',] <- (
    Pr_Breedyfuture_B0Tpos[ SLICE=1,] * Pbreedf[1] +
      Pr_Breedyfuture_B0Tpos[ SLICE=3,] * Pbreedf[3]
  ) / (Pbreedf[1] + Pbreedf[3])
  
  # PAST Bprob given state-at-sampling...
  # Also age-dept, cos I was @3 when aged 5 :/
  Pr_Breedypast_AB0Tneg <- autoloop(
    a=AGES, b0=1:3, t=0:MAX_popgap,{
      (a >= AMAT) * ((-t) >= AMAT - a) *
        Pr_Bt_B0T[ b0, 2, t] * 
        Pr_Bt_B0T[ 2, 3, pmax( (-t)-(AMAT-2-a), 0)] /
        Pr_Bt_B0T[ b0, 3, pmax( a-(AMAT-2), 0)] *
        1
    })
  
  # There are NANs, but AFAIK only where there's 0/0 and top should win!
  Pr_Breedypast_AB0Tneg[ VECSUB=
                           which( is.na( as.vector( Pr_Breedypast_AB0Tneg)))] <- 0
  
  Pr_Breedypast_ANTneg <- offarray( 0, 
                                    dimseq= list( a=AGES, n=NURSINGNESS, T=0:MAX_popgap))
  Pr_Breedypast_ANTneg[ ,'Yes',] <- Pr_Breedypast_AB0Tneg[ ,SLICE=2,]
  Pr_Breedypast_ANTneg[ ,'No',] <- (
    Pr_Breedypast_AB0Tneg[ ,SLICE=1,] * Pbreedf[1] +
      Pr_Breedypast_AB0Tneg[ ,SLICE=3,] * Pbreedf[3]
  ) / (Pbreedf[1] + Pbreedf[3])
  
  ## Weaning (ie offspring being INDEPENDENT of mum)
  r"--{
  This is ONLY relevant for the spatial case--- otherwise, it adds a lot of complications.
  
  All is implicitly conditional on mum surviving thru the weaning period, coz it's needed in XHSP cases where 2nd sib is born (by definition) after the first one weaned, and (by definition) with the mother alive; we take care of the possibility that mum dies before 2nd birthday (leading to Pr[sib]=0) separately.

  Consider 4 states in the weaning cycle, this time from an offspring-centric PoV. You start off in "Calf". Then when you are 1yo you transit to either "Mum preg" or "Mum resting". If "Mum preg", you will be "Weaned" the next year. If "Mum resting", you'll either stay there or go to "Mum preg". Once "Weaned", you stay there.
  
  This lets us calc the prob of being Weaned at-or-before (ie by) some age. Now, we actually need Pr[ *first* year of being weaned is w]. But we can get that just by differencing: Pr[ First weaned at w] = Pr[ weaned by w] - Pr[ weaned by (w-1)].
  
  Let's get into it!
}--"
  
  Tweano <- rbind(
    c(           0, 0,           0, 0), # TO 1yo
    c(   psicle[2], 0,   psicle[3], 0), # TO Mumpreg
    c( 1-psicle[2], 0, 1-psicle[3], 0), # TO Mumrest
    c(           0, 1,           0, 1)  # TO Weaned
  )
  
  Pr_Wstate_now <- c( 1, 0, 0, 0) # starting state of 1yo 
  Pr_weanedby <- rep( 0, MAX_WEANAGE) # start at age 1
  for( w in 1:MAX_WEANAGE){
    Pr_Wstate_now <- Tweano %**% Pr_Wstate_now
    Pr_weanedby[ w] <- Pr_Wstate_now[ 4]    
  }
  
  # Even the most devoted mum eventually kicks out the laziest child!
  # But NB4later:mum might not actually give birth then
  Pr_weanedby[ MAX_WEANAGE] <- 1
  
  # Pr weaning AT age
  Pr_w <- diff( c( 0, Pr_weanedby))
  
  r"--{
  If not weaned before that, then it's wean time for sure. # pmin() looks superfluous but avoids OOB. Given a sib is born at D later, clearly Mum did survive. Basic prob stuff:
  
  Pr[w|D,s] = Pr[w,D|s] / Pr[ D|s]
  Pr[w,D|s] = Pr[w|s] * Pr[D|w,s]
  Pr[D|w,s] is rebreedy over shorter interval
  
  denom = sum(numer)
  
  EXCEPT that if 1st sib is only weaned at MAX_WEANAGE, then mum may not actually have given birth then. So, find Pr[w=MAX_WEANAGE] by subtraction afterwards. But calc full array first "nominally", to get dims right
  }--"
  
  Pr_w_Db <- autoloop( w=1:MAX_WEANAGE, Db=0:MAX_Db,
                       (Db >= w) *    # MUST be weaned by Db, by assumption
                         Pr_w[ w] *     # this specific age
                         Pr_breedyagain_Db[ pmax( 0, Db-w)] / # then get breedy again
                         (Pr_breedyagain_Db[ Db] + (Pr_breedyagain_Db[Db]==0))
                       # ... NB denom looks weird but is just to guard against 0/0
  )
  
  # Get those MAX_WEANAGE probs by subtraction, coz 1st sib 
  # must be weaned by birth of 2nd
  Pr_w_Db[ MAX_WEANAGE,] <- 1 - colSums( Pr_w_Db[ 1:(MAX_WEANAGE-1),])
  
  
  ## Survival of adults; used in XHSP code
  Pr_fadsurv_t <- autoloop(
    t= 0:MAX_Tsep, # I *think* MAX_Tsep is limit...
    fadsurv ^ t
  )
  
  r"---{    
  Age-specific version (ie from given starting age). Only needed for starting ages up to AMAT; for animals known to be adult (eg parents), it's the same as Pr_fadsurv_t. Allows different survs for JUVEF and ADF. Currently used in ideal_selfP; should be used in ideal_MOP too.
  }---"
  
  # Create then loop is easiest here...  
  Pr_fsurv_ta <- offarray( 1, 
                           dimseq=list( t= 0:MAX_Dsyears, a= SELF1AGES)
  )
  
  for( tt in 1 %upto% MAX_Dsyears){
    Pr_fsurv_ta[ tt,] <- Pr_fsurv_ta[ tt-1,] * 
      ifelse( SELF1AGES + tt-1 >= AMAT, fadsurv, fjusurv)
  }
  
  ## Ageing error
  r"--{
  We are not actually using this; maybe leave it for now? But then I'd have to change the code :(
  
  Assumes Pr_a_sampled is an INPUT without error (from measured sample age compos). 
  But we are not using sampled ages for lglk on age compo per se, despite unselective-within-devstage assumption reqd for SelfP. 
  Conservative because of latter. Liberal because of no-error. Simple because it is!
  But still not used...
  }--"
  
  agerr_indices <- list(
    a= AGES,
    e= EAGES,
    y= SYEARS) # y unused 4NOW ; for when age compo changes btwn years
  
  # this is the numerator of Bayes' theorem
  Pr_a_ey <- autoloop( 
    indices= agerr_indices,
    (abs( a-e) <= MAXAGERR) *
      Pr_agerr_a[ clippo( e-a, -MAXAGERR, MAXAGERR), a] *
      Pr_a_sampled[ a]
  )
  
  # ... must sum to 1
  cumpr <- sumover( Pr_a_ey, 'a')
  # this is the denominator of Bayes' theorem
  recip_cumpr <- 1/cumpr # minimize divisions by doing once only (faster)
  Pr_a_ey <- autoloop(
    indices= agerr_indices,
    Pr_a_ey[a,e,y] * recip_cumpr[ e,y]
  )
  
  
  ## Movement stuff done separately
})