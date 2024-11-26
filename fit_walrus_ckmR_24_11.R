## Walrus design & testing on simulated data
library( mvbutils)
library( fastmatch) # for quickin
library( offarray) # now should auto-turn-off BC
library( kinsimmer)
library( despack)
library( doParallel)

r"--{
  MVB: only calc the Hessian stuff _once_ (unless sampling setup changes). We don't actually need _data_ for that, just "true" params (or working guess), number of comps, and lglk function. The _expected_ Hessian, which is what finfo_onetype() and so on calculate, does not depend on observed numbers of kin-pairs.
  
  Comparisons with "actual" simulated MOPs etc come later, in a separate script. However, we do need _one_ simulated dataset here, coz it has sample size info (n_comp_<BLAH>)
}--"

# load all of the functions in the functions folder
# functionfiles <- paste0("./functions/", list.files("./functions"))
functionfiles <- list.files( './functions', full=TRUE) # with path
# library(purrr)
# map(functionfiles, source)
# Without purrr:
lapply( functionfiles, source) 

# Eiren's suffixes here (can these just be deduced from simulation folder?)
suffixes <- c("Sc00_Sd20241023", "Sc00_Sd20241121", "Sc00_Sd20839572", 
              "Sc00_Sd29571492", "Sc00_Sd76937593", "Sc00_Sd92759173",
              "Sc00_Sd41850183", "Sc00_Sd38519472", "Sc00_Sd35719375",
              "Sc00_Sd57394720") 

# The "generic" lglk function needs to be told what the data are
# You will need to set 'simfile=...' arg; default is for MVB
lglk_with_data <- add_data( lglk_walrus,
    simfile = sprintf( 'WalrusSamples_%s.RData', suffixes[1]))
    # paste0("WalrusSamples_", suffix, ".RData")) 
denv <- environment( lglk_with_data) # where stuff lives


## True parameter values: mostly from Eiren's sim notes
# but Nfad2000 and RoI are stored on file
print( load( sprintf( './simulation/Nfad_RoI_%s.RData', suffixes[ 1])))
# variable called 'out'

ptru <- c(
  log_Nfad_y0 = log( out[ 'Nfad_2000']), # in AD2000
  RoI= out[ 'RoI'], # bugger-all
  lgt_fadsurv= logit( 0.9622), # FWIW that's 3.237
  diff_lgt_fjusurv=  logit( 0.9) - logit( 0.9622), # 0 => same as adult
  lpsi= logit( c( 0.1, 0.5)), # Pr[preg in 2nd year], Pr[preg from y>2]
  # No spatial params now
  NULL # so all previous lines can end with comma
)

lglk0 <- lglk_with_data( ptru)
# For list of prob arrays:
# probars <- lglk_with_data( ptru, want='probs')
# For popdyn goodies:
popdyn0 <- lglk_with_data( ptru, want='popdyn_only')
probs0 <- lglk_with_data( ptru, want='probs')

# Set up for parallel numderiv
ncores <- parallel::detectCores( logical=FALSE)-2 # leave spare...
ncores <- pmin( ncores, 8) # for now...
CLUSTO <- makeCluster( ncores)
# For some reason, that gives problems 1st time (windows security)
# and then is OK 2nd time..?
registerDoParallel( CLUSTO, ncores, nocompile=TRUE)

# Cluster may need explicit exports, unless 
# ... all funs are in current env
#findo <- find( 'lglk_walrus')
#if(  findo != '.GlobalEnv'){
  export_funs <- rownames( 
      mvbutils::foodweb( where='.GlobalEnv', prune='lglk_walrus', 
          ancestors=FALSE, descendents=TRUE)$funmat)
#} else export_funs <- character()

#export_funs <- "lglk_walrus"

# For parallel, need a new copy of lglk_with_data
# that can find all the data, using list not envir
# Various orta's don't work AFAICS; this does.

parallel_data_list <- c( as.list( denv), 
    mget( export_funs, inherits=TRUE))

parallel_get_probs <- function( params, DATALIST, ...){
  extract.named( DATALIST) # including lglk_walrus
  environment( lglk_walrus) <- environment() # knows stuff in DATALIST
  lglk_walrus( params, want='probs', ...) # also returns lglk
}

# just checking that parallel setup works
checko <- parallel_get_probs( ptru, parallel_data_list)

# Get bits needed to compute Hessian
Hbits <- get_Hbits( 
  ptru, parallel_get_probs,
  numderiv_fun= numvbderiv_parallel,
  Pargs=list( DATALIST=parallel_data_list, USE_NURSY_ONLY=FALSE),
  Dargs=list( 
      eps=1e-5, 
      FOREACH_ARGS=list( 
          .packages= c( 'mvbutils', 'offarray')
        ),
      PARALLEL=TRUE
    )
)    

stopImplicitCluster()
stopCluster( CLUSTO)
rm( CLUSTO)

## Now get expected Hessian and numbers of kin-pairs
# Those shouldn't vary between simulations (unless samp sizes change)
# but we need the "ncomps" from *a* simulated dataset
extract.named( prepare_H_E( denv, Hbits)) 
# ... now H & E exist
Vpar <- solve( H) # Param covariance matrix

## MVB: OK that's the "do it once" stuff 
# Now there's two separate paths: 
# 1. Compare obs & exp etc in simulations:  in "compare2sims.R"
# 2. Design calculations: in "design_eg.R"




