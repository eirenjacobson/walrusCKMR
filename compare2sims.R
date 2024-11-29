#### Compare *totals* of kin-pairs by type
# There are more checkables than just totals, but it's a start

compcheck <- array( 0, 
    c( 3, 3, length( suffixes)),
    dimnames=list( 
      cq( MOP, XmHSP, SelfP),
      cq( Exp, Obs, P),
      suffixes
    )
  )

# Loop over sim datasets, loading...
sim_lglks <- list() # we will keep this for later

for( i in seq_along( suffixes)){
  # Lglk fun for *this* dataset:
  # ... and a place to keep the data summary
  sim_lglks[[i]] <- add_data( lglk_walrus,
      simfile = sprintf( 'WalrusSamples_%s.RData', suffixes[i]),
      YSTART = denv$YSTART) # presumably 2015
  denvi <- environment( sim_lglks[[ i]]) # where stuff lives

  compcheck[ 'MOP', 'Obs', i] <- sum(denvi$n_MOP_EYEYL)
  compcheck[ 'XmHSP', 'Obs', i] <- sum(denvi$n_XmHSP_EYEY)
  compcheck[ 'SelfP', 'Obs', i] <- sum(denvi$n_selfP_EYDY) 
} # for suffix (ie sim dataset)

save( sim_lglks, file='sim_lglks.RData')

# Expecteds:
compcheck[ 'MOP', 'Exp', ] <- sum(E$n_MOP_EYEYL) # same per sim
compcheck[ 'XmHSP', 'Exp', ] <- sum(E$n_XmHSP_EYEY)
compcheck[ 'SelfP', 'Exp', ] <- sum(E$n_selfP_EYDY) 

# All probs at once:
compcheck[ , 'P', ] <- ppois( compcheck[,'Obs',], compcheck[,'Exp',])

#### Now check lglk derivs at true pars:
Dlglk0 <- matrix( 0, length( suffixes), length( ptru))

for( i in seq_along( suffixes)){
  Dlglk0[ i,] <- numvbderiv( sim_lglks[[i]], ptru, want='just_lglk')
}

### Birth gap example--- dont' wanna lose this code
# NB it's also OK to compare vs _one_ sim dataset, like so
# but best is prolly to average the obs across sims
if( FALSE){
  # Birth-gap check: 
  onn <- as.data.frame( denv$n_XmHSP_EYEY)
  onn$bgap <- with( onn, (y2-a2) - (y1-a1))
  obs_mean_bgap <- with( onn, sum( bgap * response)) / sum( onn$response)
  enn <- as.data.frame( E$n_XmHSP_EYEY)
  enn$bgap <- with( enn, (y2-a2) - (y1-a1))
  exp_mean_bgap <- with( enn, sum( bgap * response)) / sum( enn$response)
}


# MVB version of loading the Nfad (and any other) results into a matrix
# Though Nfad & RoI (and anything else) shouldn't vary much, should they..?
sim_summ <- NULL
etemp <- new.env()
for (i in seq_along(suffixes)){
  load( sprintf( './simulation/Nfad_RoI_%s.RData', suffixes[ i]), 
      envir=etemp)
  # I happen to know it's called 'out'
  
  if( i==1) {
    sim_summ <- matrix( 0, length( suffixes), length( etemp$out))
    dimnames( sim_summ) <- list( suffixes, names( etemp$out))
  }
  
  sim_summ[i,] <- e$out
  
  #  suffix <- suffixes[i]
  #  load(paste0("./simulation/Nfad_RoI_", suffix, ".RData"))
}
save( sim_summ, file='sim_summ.RData')  


