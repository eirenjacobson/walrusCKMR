"ideal_XmHSP_walrus" <- function( nlocal=sys.parent()) mlocal({
  r"---{
  Non-spatial first. We know Mum was in state 2 when first was born. What's the prob she was again in state 2 when second was born? And who else would have been? Those are pre-calced in movement_walrus.
  Let's go with the "asymmetric" version here, where b1<b2 by design. Since it might get calced for all (b1,b2), I use abs() to head off OOBery
  }---"
  
  
  Pr_XmHSP_bb <- autoloop(
    b1= PDYEARS, b2= PDYEARS,
    (b1 < b2) *   # non-symmn
      Pr_fadsurv_t[ abs( b2-b1)] * # avoid OOB
      Pr_breedyagain_Db[ abs( b2-b1)] * 
      recip_Nfad_y[ pmax( b1, b2)] *
      recip_ppn_breedy
  )
  
  
  byHSP <- list( 
    b1= PDYEARS, y1= SYEARS,
    b2= PDYEARS, y2= SYEARS)
  
  r"---{
  Spatial: the number of movement-years separating the two *would* just be the sum of their ages, plus inter-birth interval... except that #1 may stay with mum for several years til weaning. But weaning must be no later than birth of #2. So we need a weany-sum, like with MOPs.
  }---"
  
  Pr_XmHSP_by <- autoloop( 
    indices= byHSP, SUMOVER= list( w1=1:MAX_WEANAGE),
    Pr_XmHSP_bb[ b1, b2] *
      Pr_w_Db[ w1, abs( b2-b1)] #* # anti OOB
      #Pr_rr_t[ r1, r2, clippo( 
     #   (y2-b2) + # #2 since birth
     #     tsep( y1, b2, y1+w1), # mum @ birth of #2
     #   0)] # ...I think!
  )
})