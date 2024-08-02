"ideal_MOP_walrus" <- function( nlocal=sys.parent()) mlocal({
  r"---{
  Standard CK prob. If yc<bj, mother would have to survive post her sampling, and could have been immature then, so Pr_survival needs to take account of her age then.
  
  TRO is *not* adjusted to allow for ramp-up in fec around maturity. Small beer.
  
  PCOHORTS is potential-parent's birth-year. It's a wide range. We do not need the pop dyn to operate then, however; it is enough to know she was alive then.
  }---"
  
  Pr_MOP_byby_nonlethal <- autoloop( 
    bj= PDYEARS, yj = SYEARS,
    bc= PCOHORTS, yc= SYEARS,{ 
      aj = yj - bj;
      (yc > bc) * # you can't be sampled before you exist!
        (yc - bc <= AMAX) * # coz then you die      
        (bc+AMAT <= bj) * # gotta be mature
        fec_a[ clippo( bj-bc, FECAGES)] * # cow's ERO (cf avg adult)
        recip_Nfad_y[ bj] *  # competition
        # How many more years must Mum live, starting from sampling?
        # Noting that if juv > 0 mum must have also survived its first year  
        Pr_fsurv_ta[ pmax( 0, (bj+ifelse( aj>0, 1, 0))-yc), 
                     clippo( yc-bc, SELF1AGES)]
    })
  
  Pr_MOP_byby_lethal <- autoloop( 
    bj= PDYEARS, yj = SYEARS,
    bc= PCOHORTS, yc= SYEARS, {
      aj = yj - bj;
      (yc > bc) * # you can't be sampled before you exist!
        (yc - bc <= AMAX) * # coz then you die      
        (bc+AMAT <= bj) * # gotta be mature
        fec_a[ clippo( bj-bc, FECAGES)] * # cow's ERO (cf avg adult)
        recip_Nfad_y[ bj] *   # competition
        # lethal sample must not be before juv birth year & 
        # if juv is > 0, pot mom can't have died in its birth year
        ifelse((bj > yc) | (aj > 0 & yc == bj), 0, 1) 
    })
  
  Pr_MOP_bybyl <- offarray( 0, dimseq=list( bj=PDYEARS, yj=SYEARS, bc=PCOHORTS, yc=SYEARS, lc=LETHALITY))
  Pr_MOP_bybyl[,,,,SLICE='LETHAL'] <- Pr_MOP_byby_lethal
  Pr_MOP_bybyl[,,,,SLICE='NONLETHAL'] <- Pr_MOP_byby_nonlethal
  
  # Eiren's translation (not run)
  if(FALSE){
    for (bj in PDYEARS){
      for (yj in SYEARS){
        for(bc in PCOHORTS){
          for(yc in SYEARS)
            aj = yj - bj # age of juvenile when sampled
          ac_bj = ifelse(bj - bc <= 0, 1, bj - bc) # age of cow when j is born
          Pr_MOP_bby_lethal[bj, yj, bc, yc] <- 
            ifelse(yc > bc, 1, 0) *
            ifelse(yc - bc <= AMAX, 1, 0) *
            ifelse(bc + AMAT <= bj, 1, 0) *
            fec_a[ac_bj] * # note this isn't great bc will end up with neg indices
            recip_Nfad_y[ bj] *
            # lethal sample must not be before juv birth year & 
            # if juv is > 0, pot mom can't have died in its birth year
            ifelse((bj > yc) | (aj > 0 & yc == bj), 0, 1) 
        }
      }
    }
  } # end if FALSE
  
  r"---{
    ... below: in that Pr_fsurv_ta[], the pmin( AMAT,...) is there coz surv is constant once adult, so the array index can stop at AMAT. The pmax() is there cos if O's birth is before potential M was sampled, she was certainly not-yet-dead back at O's birth (first line picks up that she must still have been born & mature).
  }---"
  
  # Version conditioned on adult's Nursingness at sampling
  # (ie with/without calf). I have basically run out of letters...
  # NON-SPATIAL ONLY!
  
  # Separate calcs for ad-samp before/after ju-birth
  # Easiest to compute both inside one autoloop
  # zeroing if it doesn't apply
  
  # NB these "look" superficially different from Pr_MOP_bby
  # only becoz that has fec_a which is normed to 1 for "adults"
  # but here we work with "absolute" fecundity, and must 
  # incorp ppn_breedy explicitly
  
  TPOSNEG <- autoloop( 
    bj= PDYEARS, yj = SYEARS,
    bc= PCOHORTS, yc= SYEARS, 
    nc= NURSINGNESS, lc = LETHALITY, {
      acclip <- clippo( yc - bc, AMIN, AMAX); # for past only
      aj = yj - bj;
      
      # NOTE this will only apply to nonlethal samples
      TPOS <- # potential juv is born in the future
        (lc != "LETHAL") * # TPOS impossible if cow sampling is lethal
        (yc < bj) *  # birth must be IN FUTURE (after sampling of adult) 
        (yc > bc) * # born before sampled!
        (yc - bc <= AMAX) * # coz then you die
        (bc+AMAT <= bj) * # must be mature...
        # ... and in right Bstate, given w/wo calf at samp
        Pr_Breedyfuture_NTpos[ nc, max(bj-yc, 1)] * 
        recip_Nfad_y[ bj] * recip_ppn_breedy * # competition
        # ... and survival from yc til j's birth 
        Pr_fsurv_ta[ pmax( 0, (bj + ifelse(aj>0, 1, 0))-yc), 
                     # adjust so cow must live through juv's first year
                     clippo( yc-bc, SELF1AGES)]
      
      TNEG <- # in the past
        (yc > bc) * # born before sampled!
        (yc - bc <= AMAX) * # coz then you die        
        (bc+AMAT <= bj) * # must have been mature...
        # ... and in right Bstate, given w/wo calf at samp
        Pr_Breedypast_ANTneg[ acclip,  nc, abs( bj-yc)] *
        recip_Nfad_y[ bj] * recip_ppn_breedy * # competition
        1* # no survival required
        ifelse((bj > yc) | ((bj == yc) & (lc == "LETHAL") & (aj > 0)), 0, 1)
      # only intersted in past cases in TNEG
      # if the birth year of the pot juv is the same as the year the pot cow
      # was sampled lethally and the pot juv is > 0, that is impossible
      
      returnList( TPOS, TNEG)
    })
  # Exactly one of past or future case will apply...
  # ... other will be 0, by construction
  Pr_MOP_bybynl <- with( TPOSNEG, TPOS + TNEG)
  
  # Checking code (to be run manually from debugger)
  if( FALSE){ 
    r"--{
     Marginal MOP prob should be *between* cond on Yes & on No coz wted average. I think equality (to either one) is only possible if all are 0, or at least one is 0 anyway.

     Use of FEC_ASYMP_AGE will slightly break that coz it's an approx, so for testing make FEC_ASYMP_AGE=AMAX-1 (-1 "just becoz").
     }--"
    
    # Find offending cases...
    glump <- Pr_MOP_bby
    glump[] <- pmax( 
      c( Pr_MOP_bbyn[,,,SLICE='Yes']), 
      c( Pr_MOP_bbyn[,,,SLICE='No'])
    )
    erk <- whichoff( glump, .==0 & Pr_MOP_bby>0)
    nrow( erk) # 0; OK, have fixed it now!
    
    # summary( erk[,3] - erk[,2]) # aha..!
    # All v. old; should have auto-died before sampling
    
    reverk <- whichoff( glump, .>0 & Pr_MOP_bby==0) # empty. Phew.
  }
})