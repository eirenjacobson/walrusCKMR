## Walrus design & testing on simulated data
library( mvbutils)
library( fastmatch) # for quickin
library( offarray) # now should auto-turn-off BC
library( kinsimmer)
library( despack)
library( doParallel)

library(purrr)

# load all of the functions in the functions folder
functionfiles <- paste0("./functions/", list.files("./functions"))
map(functionfiles, source)

# Use the simulations as a guide to future samp size etc
# also set up for actual testing on data1
# The "generic" lglk function needs to be told what the data are
# You will need to set 'simfile=...' arg; default is for MVB
lglk_with_data <- add_data( lglk_walrus, 
                            simfile = "WalrusSamples_RealisticNoLethality.RData") 
denv <- environment( lglk_with_data) # where stuff lives

## These come from Eiren's sim notes
ptru <- c(
  log_Nfad_y0 = log(68418), # in AD2000
  RoI= -0.006515221, # bugger all (3% over 80years =...)
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

# Hbits_Nursy <- get_Hbits( 
#   ptru, parallel_get_probs,
#   numderiv_fun= numvbderiv_parallel,
#   Pargs=list( DATALIST=parallel_data_list, USE_NURSY_ONLY=TRUE),
#   Dargs=list( 
#       eps=1e-5, 
#       FOREACH_ARGS=list( 
#           .packages= c( 'mvbutils', 'offarray')
#         ),
#       PARALLEL=TRUE
#     )
# )    

stopImplicitCluster()
stopCluster( CLUSTO)
rm( CLUSTO)

# Want this for all comp-types:
# H$MOP <- finfo_onetype( Hbits$DSP$DSP_MOP_REY, denv$n_comp_MOP_REY)
# fisher info gives you the informational value of a single comparison of type X
# also Ekin
# Via loop below

E_comp <- H_comp <- list()

# mops
ncomp <- denv[["n_comp_MOP_AYNL"]] 
# H_comp contains inverse covariance between parameters
H_comp[["DSP_MOP_EYEYNL"]] <- finfo_onetype(Hbits$DSP[["DSP_MOP_EYEYNL"]], ncomp)
# E_comp is expected number of kin pairs 
E_comp[["MOP_EYEYNL"]] <- ncomp * Hbits$Prkin[["Pr_MOP_EYEYNL"]]

# check expected versus observed number of kin pairs
sum(E_comp$MOP_EYEYNL)
sum(denv$n_MOP_AYL)

# half sibs
ncomp <- denv[["n_comp_XmHSP_AY"]]
H_comp[["DSP_XmHSP_EY"]] <- finfo_onetype( Hbits$DSP[[ "DSP_XmHSP_EY"]], ncomp)
E_comp[["XmHSP_EY"]] <- ncomp * Hbits$Prkin[["Pr_XmHSP_EY"]]

# check expected versus observed number of kin pairs
sum(E_comp$XmHSP_EY) 
sum(denv$n_XmHSP_AY) 

# self
ncomp <- denv[["n_comp_selfP_YADY"]]
H_comp[["DSP_selfP_YADY"]] <- finfo_onetype( Hbits$DSP[[ "DSP_selfP_YADY"]], ncomp)
E_comp[["selfP_YADY"]] <- ncomp * Hbits$Prkin[["Pr_selfP_YADY"]]

sum(E_comp$selfP_YADY) 
sum(denv$n_selfP_YADY) 

# EKJ did this bit manually (above) bc names not consistent re. AY or EY 
# for( comptype in names( Hbits$DSP)){
#   ncomp <- denv[[ sub( 'DSP', 'n_comp', comptype)]]
#   H_comp[[ comptype]] <- finfo_onetype( Hbits$DSP[[ comptype]], ncomp)
#   E_comp[[ sub( 'DSP_', '', comptype)]] <- ncomp * 
#       Hbits$Prkin[[ sub( 'DSP', 'Pr', comptype)]]
# }

 # Birth-gap check: 
 onn <- as.data.frame( denv$n_XmHSP_AY)
 onn$bgap <- with( onn, (y2-a2) - (y1-a1))
 obs_mean_bgap <- with( onn, sum( bgap * response)) / sum( onn$response)
 enn <- as.data.frame( E_comp$XmHSP_EY)
 enn$bgap <- with( enn, (y2-a2) - (y1-a1))
 exp_mean_bgap <- with( enn, sum( bgap * response)) / sum( enn$response)

# Add up Hessian bits:
H <- 0*H_comp[[1]]
for( comptype in names( Hbits$DSP)[2:4]){
  H <- H + H_comp[[ comptype]]
}

# What would the approx MLE be?
# Shift is Dlglk (score) "times" inv(Hess) ; plus/minus!
Vpar <- Hinv <- solve( H)
parshift <- Hinv %**% Hbits$Dnonprob$lglk
parshift / sqrt( diag( Vpar)) # in SDs
# actually not that bad... all well within 2SD
# Should this be +shift or -shift? "Yes" ;)

# Would like Hbits$Dnonprob$lglk on average should be zero (across sims)
# relative to the actual Hessian (expect to have high variances)
# SAVE THESE AS OUTPUT so they can be compared across simulations
# can tell us which parameters might be biased

# Well, we can try it empirically:
lglk_plus_shift <- lglk_with_data( ptru + parshift)
lglk_minus_shift <- lglk_with_data( ptru - parshift)
lglk0

# ... plus shift is definitly better! "Only" by 10 units so not tooo drastic...
 popdyn_shift <- lglk_with_data( ptru + parshift, want='popdyn')
 sumEshift <- numeric()
 forEshift <- lglk_with_data( ptru + parshift, want='probs')
 for( comptype in names( Hbits$DSP)){
   ncomp <- denv[[ sub( 'DSP', 'n_comp', comptype)]]
   sumEshift[ sub( 'DSP_', '', comptype)] <- sum( 
       ncomp * forEshift[[ sub( 'DSP', 'Pr', comptype)]]
     )
 }
 rbind( sumEshift, sumbo)
 
 
 # mops
 ncomp <- denv[["n_comp_MOP_AYL"]] 
 # H_comp contains inverse covariance between parameters
 sumEshift[["DSP_MOP_EYEYNL"]] <- sum(ncomp * forEshift$Pr_MOP_EYEYL)

 ncomp <- denv[["n_comp_XmHSP_AY"]] 
 # H_comp contains inverse covariance between parameters
 sumEshift[["DSP_XmHSP_EY"]] <- sum(ncomp * forEshift$Pr_XmHSP_EY)
 
 ncomp <- denv[["n_comp_selfP_YADY"]] 
 # H_comp contains inverse covariance between parameters
 sumEshift[["n_comp_selfP_YADY"]] <- sum(ncomp * forEshift$Pr_selfP_YADY)
 
 
# sumEshift is now quite reasonable cf obs kinpair tots
# and params have not changed that much (juve surv mainly)
# 
# # Using Nursy status:
# Hbits_Nursy <- get_Hbits( 
#   ptru, parallel_get_probs,
#   numderiv_fun= numvbderiv_parallel,
#   Pargs=list( DATALIST=parallel_data_list, USE_NURSY_ONLY=TRUE),
#   Dargs=list( 
#       eps=1e-5, 
#       FOREACH_ARGS=list( 
#           .packages= c( 'mvbutils', 'offarray')
#         ),
#       PARALLEL=TRUE
#     )
# )    
# 
# H_comp_Nursy <- list()
# for( comptype in names( Hbits_Nursy$DSP)){
#   ncomp <- denv[[ sub( 'DSP', 'n_comp', comptype)]]
#   H_comp_Nursy[[ comptype]] <- finfo_onetype( 
#       Hbits_Nursy$DSP[[ comptype]], ncomp)
# }
# 
# # Add up Hessian bits:
# H_Nursy <- 0*H_comp_Nursy[[1]]
# for( comptype in names( Hbits_Nursy$DSP)){
#   H_Nursy <- H_Nursy + H_comp_Nursy[[ comptype]]
# }


# Quantities of interest:
stuff0 <- interesting_stuff( ptru)

# Just popdyn is quick (probs are slow) so don't parallel this
Dstuff0 <- numvbderiv( interesting_stuff, ptru, eps=1e-5)
Vstuff0 <- Dstuff0 %**% Vpar %**% t( Dstuff0)
SEstuff0 <- sqrt( diag( Vstuff0))
mpt0 <- make_precision_table( stuff0, SEstuff0, E_comp)

# Ready for import into Lyx--- which works pretty well
if( !dir.exists( 'results')){
  mkdir( results) # MVB 8/11--- tho might wanna put elsewhere anyway
}
write.csv( mpt0, file='results/out.csv')

# MVB 8/11: I have turned the block below into "unexecuted" rather than
# hashed-out, so that I can test it (even if we don't want it yet)
if( FALSE){
  ## Investigate contrib of CKMR data
  # Keep a bit of CKMR for identifiability
  H_comp_IMR <- list()
  for( comptype in names( Hbits$DSP)){
    ncomp <- denv[[ sub( 'DSP', 'n_comp', comptype)]]
    if( !grepl( '(?i)self', comptype)){
      ncomp <- ncomp / 100
    }
    H_comp_IMR[[ comptype]] <- finfo_onetype( Hbits$DSP[[ comptype]], 
      ncomp)
  }

  # ... and combine
  H_IMR <- 0*H_comp_IMR[[1]]
  for( comptype in names( Hbits$DSP)){
    H_IMR <- H_IMR + H_comp_IMR[[ comptype]]
  }

  # ... with hindsight, I could just have applied the 
  # downweighting inside the last loop, without recalling finfo_onetype
  # Never Mind.

  Vnock <- solve( H_IMR)
  Vstuff_nock <- Dstuff0 %**% Vnock %**% t( Dstuff0)
  mpt_nock <- make_precision_table( stuff0, sqrt( diag( Vstuff_nock)), E_comp) 

  ## Historical & future samp sizes match the sims (I hope...)

  ## Let's tweak juve vs ad samp sizes in future...
  # More juves...
  mpt_mj <- make_precision_table_completely( 
      stuff0, Dstuff0, Hbits, denv, 
      adjust_future_sample_sel( denv, 2))
  Vmj <- attr( mpt_mj, 'Vpar')
  # sqrt( diag( Vmj)) is mostly _slightly_ worse than for Vpar
  # except juve surv, unsurprisingly! but improvement is small there

  # Less juves...
  mpt_lj <- make_precision_table_completely( 
      stuff0, Dstuff0, Hbits, denv, 
      adjust_future_sample_sel( denv, 0.5))
  Vlj <- attr( mpt_lj, 'Vpar')

  # Splunge the mpt's together...
  mpt_all <- rbind( mpt0[2,], mpt_mj[2,], mpt_lj[2,])
  Ecols <- which( startsWith( colnames( mpt0), 'E')) 
  mpt_all[ , Ecols] <- rbind( mpt0[1,Ecols], mpt_mj[1,Ecols], mpt_lj[1,Ecols])
}  


