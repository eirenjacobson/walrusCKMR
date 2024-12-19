
library( mvbutils)
library( atease) # use @ to access attributes--- lazy is good sometimes
library( fastmatch) # for quickin
library( offarray) # now should auto-turn-off BC
library( kinsimmer)
library( despack)
library( doParallel)

functionfiles <- list.files( './functions', full=TRUE) # with path
lapply( functionfiles, source) 

# Suffixes from simulation folder:
suffixes <- dir( 'simulation', patt='^WalrusSamples.*RData') |>
  xsub( 'WalrusSamples_', '') |>
  xsub( '[.]RData', '')

# only want d1-3 and s1-7
# d0 is "null" for testing only
# s0 is all seven years of sampling, for testing only
design_sc <- which(substr(suffixes, 12, 19) == "20241023" & 
                     substr(suffixes, 2,2) != 0 & substr(suffixes, 8,8) != 0) # length should be 42

# set up comparisons

compdesign <- array(NA, 
                    c( 3, 2, 7, 2, 2, 7, 3),
                    dimnames=list( 
                      D = cq(1, 2, 3),
                      L = cq( 1, 2),
                      S = cq(1, 2, 3, 4, 5, 6, 7),
                      CKMR = cq(Yes, No),
                      SelfP = cq(Yes, No),
                      Type = cq(Nfad_2015, Nfad_2020, Nfad_2025, fadsurv, fjusurv, ppn_breedy, rel_ad_15_25),
                      Value = cq(True, Est, SE)
                    )
)

# set up giant dataframe to hold results

design_df <- expand.grid(ID = suffixes[design_sc], CKMR = c("Yes", "No"), 
                         Value = c("TRUE", "Est", "SE", "CV"))
design_df <- cbind.data.frame(design_df, "D" = NA, "L" = NA, "S" = NA,
                              "Nfad_2015" = NA, "Nfad_2020" = NA, "Nfad_2025" = NA, 
                              "fadsurv" = NA, "fjusurv" = NA,
                              "ppn_breedy" = NA, "rel_ad_15_25" = NA)
samplesizes <- data.frame(ID = suffixes[design_sc], "Samples" = NA, "AllSamples" = NA, 
                          "MOPs" = NA, "XmHSPs" = NA, "SelfP" = NA)
 
for (i in 1:length(design_sc)){

  TEST_SUFFIX <- design_sc[i] # choose a "scenario" 
  
  # load data to check sample sizes
  
  # current values
  
  d <- as.numeric(substr(suffixes[TEST_SUFFIX], 2, 2))
  l <- as.numeric(substr(suffixes[TEST_SUFFIX], 5, 5))
  s <- as.numeric(substr(suffixes[TEST_SUFFIX], 8, 8))
  
  design_df[which(design_df$ID == suffixes[TEST_SUFFIX]), "D"] <- d
  design_df[which(design_df$ID == suffixes[TEST_SUFFIX]), "L"] <- l
  design_df[which(design_df$ID == suffixes[TEST_SUFFIX]), "S"] <- s
  
  # always load the corresponding file for S0, because that is the only one saved from each sim
  load( paste0("./simulation/Nfad_RoI_D", d, "_L", l, "_S0_Sd20241023.RData"))
  
  # what are the demo pars
  di <- substr(suffixes[TEST_SUFFIX], 1, 2)
  if(di == "D1"){phi1 <- 0.99; phi2 <- 0.55; effective_adult_surv <- truncsurvsenesce(phi1, phi2) }else
    if(di == "D2"){phi1 <- 0.985; phi2 <- 0.5; effective_adult_surv <- truncsurvsenesce(phi1, phi2)}else
      if(di == "D3"){phi1 <- 0.99; phi2 <- 0.6; effective_adult_surv <- truncsurvsenesce(phi1, phi2)}else
      {phi <- 0.9622; effective_adult_surv <- truncsurv(phi=phi)}
  
  lglk_with_data <- add_data( lglk_walrus,
                              simfile = sprintf( 'WalrusSamples_%s.RData', suffixes[ TEST_SUFFIX]),
                              YSTART= 2015,      #  more stable parametrization (no math difference)
                              SYEARS= 2013:2028)
  # SYEARS explicit, rather than inferred from simfile; make sure 
  # ... it's the same in all cases
  denv <- environment( lglk_with_data) # where stuff lives
  
  # keep track of sample sizes
  samplesizes[i,]$Samples <- denv$nsamples
  samplesizes[i,]$AllSamples <- denv$nallsamples
  samplesizes[i,]$MOPs <- sum(denv$n_MOP_EYEYL)
  samplesizes[i,]$XmHSPs <- sum(denv$n_XmHSP_EYEY)
  samplesizes[i,]$SelfP <- sum(denv$n_selfP_EYDY)
  
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
  
  source("./scripts/design_eg_ekj.R")
}

save(compdesign, file = "./results/compdesign.RData")
save(design_df, file = "./results/design_df.RData")
save(samplesizes, file = "./results/samplesizes.RData")

