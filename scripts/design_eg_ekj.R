#######
## Here's some code for Design calculations
## Assumes that the following exists, eg from 'fit_walrus_ckmR_24_11.R'
# ptru, Vpar, denv (which contains sample size info)
# and various functions

# current values

d <- 1
l <- 2
s <- 3

# Ready for import into Lyx--- which works pretty well
if( !dir.exists( 'results')){
  mkdir( 'results') # MVB 8/11--- tho might wanna put elsewhere anyway
}

compdesign <- array(NA, 
                    c( 3, 3, 7, 2, 2, 6, 3),
                    dimnames=list( 
                      D = cq( 1, 2, 3),
                      L = cq( 1, 2),
                      S = cq( 1, 2, 3, 4, 5, 6, 7),
                      CKMR = cq(Yes, No),
                      SelfP = cq(Yes, No),
                      Type = cq(Nfad_2015, Nfad_2025, fadsurv, fjusurv, ppn_breedy, rel_ad_15_25),
                      Value = cq(True, Est, SE)
                      
                    )
)

compdesign[d, l, s, "Yes", "Yes", , "Est"] <- interesting_stuff(ptru, lglk_with_data)


# Need a SS scenario
denv <- environment( lglk_with_data)
m_SYEL_0 <- denv$m_SYEL

basic0 <- doades( 
  m_SYEL_0,          # could be anything
  Hbits,             # *assume* this p'tic Hbits matches this lglk_with_data
  interesting_stuff, # look at function code
  lglk_with_data
)

compdesign[d,l,s,"Yes", "Yes", ,"SE"] <- as.numeric(basic0)

basic0_noCK <- doades( 
  m_SYEL_0,          # could be anything
  Hbits,             # *assume* this p'tic Hbits matches this lglk_with_data
  interesting_stuff, # look at function code
  lglk_with_data,
  comp_wts= c( MOP=0.01, XmHSP=0.01)
)

compdesign[d,l,s,"No", "Yes", ,"SE"]  <- as.numeric(basic0_noCK)

save(compdesign, file = "./results/compdesign.RData")

tout <- round(compdesign[d,l,s,,"Yes",,"SE"], digits = 2)

write.csv( tout, file=paste0("results/compdesign_D", d, "_L", l, "_S", s, ".csv"))

cvs <- compdesign[1, 2, 1, "Yes", "Yes",, "SE"]/compdesign[1, 2, 3, "Yes", "Yes",, "Est"]


