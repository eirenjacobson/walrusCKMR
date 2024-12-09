## Walrus design & testing on simulated data
library( mvbutils)
library( atease) # use @ to access attributes--- lazy is good sometimes
library( fastmatch) # for quickin
library( offarray) # now should auto-turn-off BC
library( kinsimmer)
library( despack)
library( doParallel)

r"--{
  MVB: only calc the Hessian ingredients (Hbits) once, unless popdyn changes. We don't actually need _data_ for that, just "true" params (or working guess), number of comps, and lglk function. The _expected_ Hessian, which is what finfo_onetype() and so on calculate, does not depend on observed numbers of kin-pairs.
  
  Comparisons with "actual" simulated MOPs etc come later, in a separate script. However, we do need _one_ simulated dataset here, coz it has sample size info (n_comp_<BLAH>)
}--"

# load all of the functions in the functions folder
functionfiles <- list.files( './functions', full=TRUE) # with path
# Probably MVB only: don't accidentally overwrite!
if( any( sub( '[.][Rr]$', '', basename( functionfiles)) %in% 
    lsall( .GlobalEnv) )){ 
  KABOOM <- yes.no( 'Reload & overwrite functions..?')
  if( KABOOM){
    lapply( functionfiles, source) 
  }
}

# Suffixes from simulation folder:
suffixes <- dir( 'simulation', patt='^WalrusSamples.*RData') |>
    xsub( 'WalrusSamples_', '') |>
    xsub( '[.]RData', '')

# The "generic" lglk function needs to be told what the data are
# You will need to set 'simfile=...' arg; default is for MVB

# For some reason, suffixes[1] only has sampling thru 2017 (I guess it is
# all the _known_ extant data) so I am running with suffixes[4] here
# coz it goes thru 2027
# by which time there will be a LOT of pairs :) !!!

TEST_SUFFIX <- 4 # choose a "scenario" that works nicely...

lglk_with_data <- add_data( lglk_walrus,
    simfile = sprintf( 'WalrusSamples_%s.RData', suffixes[ TEST_SUFFIX]),
    YSTART= 2015,      #  more stable parametrization (no math difference)
    SYEARS= 2013:2027)
    # SYEARS explicit, rather than inferred from simfile; make sure 
    # ... it's the same in all cases
denv <- environment( lglk_with_data) # where stuff lives
r"--{
  There will prolly be a warning from add_data() about zapping (ie ignoring) actual observed kin.  Now check eg denv$zap_MOP if that was mentioned in the warning. When I checked, all zappings were sensible; ie, such pairs could reasonably exist, but we have deliberately decided not to use such comparisons, in make_n_comps()
}--"


## True parameter values: mostly from Eiren's sim notes
# but Nfad2000 and RoI are stored on file
print( load( sprintf( './simulation/Nfad_RoI_%s.RData', 
    suffixes[ TEST_SUFFIX])))
# variable called 'out', which has "real" abund & RoI param values

# Adults in simulation all die @~37yo, ie ~30yrs after maturity.
# Survival before that is 0.9622, so average adult survival is lower:
effective_adult_surv <- truncsurv( 30, 0.9622) # 0.9458

# NB NB: CHANGED this 
Nfad_2015 <- out[ 'Nfad_2000'] * exp( 15*out['RoI'])
ptru <- c(
  log_Nfad_y0 = log( Nfad_2015), # NB add_data( ... YSTART=2015)
  RoI= out[ 'RoI'], # bugger-all
  lgt_fadsurv= logit( effective_adult_surv),
  diff_lgt_fjusurv=  logit( 0.9) - logit( 0.9622), # 0 => same as adult
  lpsi= logit( c( 0.1, 0.5)), # Pr[preg in 2nd year], Pr[preg from y>2]
  # No spatial params now
  NULL # so all previous lines can end with comma
)

lglk0 <- lglk_with_data( ptru, want='just_lglk')
# Dlglk0 <- numvbderiv( lglk_with_data, ptru, want='just_lglk')

# For list of prob arrays:
# probars <- lglk_with_data( ptru, want='probs')
# For popdyn goodies:
popdyn0 <- lglk_with_data( ptru, want='popdyn')
probs0 <- lglk_with_data( ptru, want='probs')

# Set up for parallel numderiv
ncores <- parallel::detectCores( logical=FALSE)-2 # leave spare...
ncores <- pmin( ncores, 8) # for now...
CLUSTO <- makeCluster( ncores)
# For some reason, that gives problems 1st time (windows security)
# and then is OK 2nd time..?
registerDoParallel( CLUSTO, ncores, nocompile=TRUE)

# Cluster may need explicit exports--- might not if all functions
# are in .GlobalEnv, but I'm not sure. Maybe it needs them anyway.
# I'm not sure. I don't care! Just export them explicitly anyway!
# mvbutils::fo
export_funs <- rownames( mvbutils::foodweb( 
    where='.GlobalEnv', prune='lglk_walrus', 
    plotting=FALSE, ancestors=FALSE, descendents=TRUE)$funmat
  )

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
    Pargs=list( DATALIST=parallel_data_list),
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
# but we need the "ncomps" from *one* simulated dataset
extract.named( prepare_H_E( denv, Hbits)) 
# ... now H & E exist
Vpar <- solve( H) # Param covariance matrix

## MVB: OK that's the "do it once" stuff 
# Now there's two separate paths: 
# 1. Compare obs & exp etc in simulations:  in "compare2sims.R"
# 2. Design calculations: in "design_eg.R"




