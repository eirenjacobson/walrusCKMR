"ideal_selfP_walrus" <- function( nlocal=sys.parent()) mlocal({
  r"---{  
Condition on *1st* sample's age, only; and on Development-stage (d) of 2nd sample (Juve or Ad). Assume unselective sampling *within* Devstage, but not across Devstages.

Also, assume first sample is taken before second, and non-lethally..! (Data org should exclude lethal 1st samples from comparisons)
}---"
  
  # Slightly different set of indices cf MOP
  # SELF1AGES is all juve ages, plus AMAT which serves as "plus group"
  # DEVSTAGES is (JuveF, AdF)
  # Devstage_A maps a true age into DEVSTAGE
  
  yady_12 <- list( 
    y1= SYEARS, a1= SELF1AGES,
    d2= DEVSTAGES, y2= SYEARS)
  
  # abs( y2-y1) below since prob may as well be computed for impossibles---
  # but gotta avoid OOB---
  # easier than deciding which are possible!
  Pr_selfP_YADY <- autoloop( 
    indices= yady_12,
    (y2 > y1) * # no same-year (also must censor in data)
     # Pr_rr_t[ r1, r2, abs( y2-y1)] * # movement
     # recip_Pr_r[ r1] * recip_Pr_r[ r2] *
      ( Devstage_A[ pmin( a1 + abs(y2-y1), AMAT)]==d2) * # right Devstage?
      Pr_fsurv_ta[ abs( y2-y1), a1] * # survival
      1/Nf_yd[ y2, d2] *
      1 # eases debugging if all lines end in "*"
  )
  
})