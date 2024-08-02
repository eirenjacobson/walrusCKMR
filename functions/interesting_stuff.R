"interesting_stuff" <- function( pars, ..., data_env=denv){
  lglk <- lglk_walrus 
  environment( lglk) <- data_env # no copy penalty on envirs
  raw_stuff <- lglk( pars, ..., want='popdyn_only')
  
  # Only *scalar* interesting things are allowed
  # Vectors are not interesting. End of.
  
  # Make it easy to change the list of calculees, and 
  # ... avoid having to explicitly refer to raw-stuff everywhere
  e <- new.env( parent=list2env( raw_stuff)) 
  evalq( envir=e, {
    # Well, these are not *that* interesting, really...  
    Nfad_2005 <- Nfad_y[ 2005]
    Nfad_2015 <- Nfad_y[ 2015]
    Nfad_2025 <- Nfad_y[ 2025]
    Nfju_2025 <- Nfju_y[ 2025]
    rel_ad_15_25 <- Nfad_y[ 2025]/Nfad_y[ 2015] - 1
    extract.named( real_params[ cq( fadsurv, fjusurv)])
    ppn_breedy <- ppn_breedy
  })
  
  # Must return a vector of interesting scalars!
  stuff <- unlist( as.list( e))
  
  # Sorting by names helps a bit... abunds first
  iscap <- substring( names( stuff), 1, 1) %in% LETTERS
  o <- order( !iscap, names( stuff))
  stuff <- stuff[ o]
  return( stuff)  
}