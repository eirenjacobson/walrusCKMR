"age_blur_probs_walrus" <- function( nlocal=sys.parent()) mlocal({
  r"---{
  Assumes PRE-KNOWN matrix of ageing "error" (uncertainty). MAXAGERR is greatest poss +/- deviation between truth & actual. "E" stands for "estimated age", in EAGES; might be wider than true age range.
  
  Gonna special-case MAXAGERR==0 coz code should be much faster (and that's the base case anyway).
  
  If MAXAGERR>0, this should be re-written to loop over AGERR=-MAXAGERR:+MAXAGERR only, since that will much quicker than running over the full age range; for most measured ages, most _real_ ages are not possible (you can't really be 20 if DNAage says 30).
}---"
  
  #     abs( B[ HSPs[,1]] - B[ HSPs[,2]]) <= MAX_BGAP_HSPS,] 
  
  eyHSP <- list( 
    e1= EAGES, y1= SYEARS,
    e2= EAGES, y2= SYEARS)
  
  if( MAXAGERR==0){ # e1==a1, e2==a2
    Pr_XmHSP_EY <- autoloop( 
      indices= eyHSP, {
        b1 <- y1-e1
        b1clip <- pmax( b1, FIRST_PDYEAR)
        b2 <- y2-e2
        b2clip <- pmax( b2, FIRST_PDYEAR)
        
        (b1 >= b1clip) * (b2 >= b2clip) *
          Pr_XmHSP_by[b1clip, y1, b2clip, y2] 
      })
  } else {
    Pr_XmHSP_EY <- autoloop( 
      indices= eyHSP, SUMOVER= list( a1= AGES, a2= AGES), {
        # Anti OOB
        b1 <- y1-a1
        b1clip <- pmax( b1, FIRST_PDYEAR)
        b2 <- y2-a2
        b2clip <- pmax( b2, FIRST_PDYEAR)
        
        (b1 >= b1clip) * (b2 >= b2clip) *
          Pr_XmHSP_by[b1clip, y1, b2clip, y2] * 
          Pr_a_ey[ a1, e1, y1] * Pr_a_ey[ a2, e2, y2]
      })
  }
  
  ## MOPs: 
  r"---{
  NB for efficiency, birth-year of Parent is pushed up to min( PCOHORTS), since any parent born prior to that will have reached asymp fec by the time the offspring is born
  }---"

  eyeylMOP <- list(
    ej= EAGES, yj= SYEARS,
    ec= EAGES, yc= SYEARS, lc=LETHALITY) 
  
  if( MAXAGERR==0){ # ej==aj, ec==ac
    Pr_MOP_EYEYL <- autoloop(
      indices= eyeylMOP, {
        bj <- yj-ej
        bjclip <- pmax( bj, FIRST_PDYEAR)
        
        bc <- yc-ec
        bcclip <- clippo( bc, FIRST_PCOHORT, LAST_PCOHORT)
        
        (bj==bjclip) * (bc <= LAST_PCOHORT) *
          Pr_MOP_bybyl[bjclip, yj, bcclip, yc, lc]
      })
    
    # Nursy case:
    # By now, I'm just ignoring spatial; cut to chase!
    Pr_MOP_EYEYNL <- autoloop(
      indices= list(
        ej= EAGES, yj= SYEARS,
        ec= EAGES, yc= SYEARS, nc=NURSINGNESS, lc = LETHALITY),
      {
        bj <- yj-ej
        bjclip <- pmax( bj, FIRST_PDYEAR)
        
        bc <- yc-ec
        bcclip <- clippo( bc, FIRST_PCOHORT, LAST_PCOHORT)
        
        (bj==bjclip) * (bc <= LAST_PCOHORT) *
          Pr_MOP_bybynl[ bjclip, yj, bcclip, yc, nc, lc] *
          1
      })
  } else { # AGE ERROR
    # Incomplete (no Nursy) and inefficient anyway (bc sums over all possible ages, not just nearby ages)
    Pr_MOP_EY <- autoloop(
      indices= eyeylMOP, SUMOVER= list( aj= AGES, ac= AGES), {
        bj <- yj-aj
        bjclip <- pmax( bj, FIRST_PDYEAR)
        
        bc <- yc-ac
        bcclip <- clippo( bc, FIRST_PCOHORT, LAST_PCOHORT)
        
        (bj==bjclip) * (bc <= LAST_PCOHORT) *
          Pr_MOP_bybyl[bjclip, yj, bcclip, yc, lc] *
          Pr_a_ey[ aj, ej, yj] * Pr_a_ey[ ac, ec, yc]
      })
    
    Pr_MOP_eyeyl <- autoloop(
      indices= list(ej = AGES, yj = SYEARS, ec = AGES, yc = SYEARS, lc = LETHALITY), 
      SUMOVER= list( aj= AGES, ac= AGES), {
        bj <- yj-aj
        bjclip <- pmax( bj, FIRST_PDYEAR)
        
        bc <- yc-ac
        bcclip <- clippo( bc, FIRST_PCOHORT, LAST_PCOHORT)
        
        (bj==bjclip) * (bc <= LAST_PCOHORT) *
          Pr_MOP_bybyl[bjclip, yj, bcclip, yc, lc] *
          Pr_a_ey[ aj, ej, yj] * Pr_a_ey[ ac, ec, yc]
      })
  }
  
  ## NB NB: For Nursingness variant, I've ignored space and agerr
  # and done the calc directly in ideal_MOP_prob()
  # straight to MOP_EYN
  
  eySelf <- list( 
    e1= SELFEAGES, y1= SYEARS,
    e2= SELFEAGES, y2= SYEARS)
  
  # abs( y2-y1) below since prob may as well be computed for impossibles---
  # but gotta avoid OOB---
  # easier than deciding which are possible!
  # NB two statements here, for clarity. autoloop() returns only the last one, 
  # unless you specifically ask it to return several things.
  
  Pr_selfP_EY <- autoloop( 
    indices=eySelf, SUMOVER= list( a1= SELF1AGES, a2= SELF1AGES), {
      d2 <- ifelse( a2 >= AMAT, "AdF", "JuveF")
      Pr_selfP_YADY[y1, pmin( a1, AMAT), d2, y2] * 
        Pr_a_ey[ a1, e1, y1] * Pr_a_ey[ a2, e2, y2]
    })
  
}) # end age_blur_probs_walrus