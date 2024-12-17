#######
## Here's some code for Design calculations
## Assumes that the following exists, eg from 'fit_walrus_ckmR_24_11.R'
# ptru, Vpar, denv (which contains sample size info)
# and various functions

compdesign[d, l, s, "Yes", "Yes", , "Est"] <- as.numeric(interesting_stuff(ptru, lglk_with_data))
design_df[which(design_df$ID == suffixes[TEST_SUFFIX] & design_df$CKMR == "Yes" & design_df$Value == "Est"),7:13] <- as.numeric(interesting_stuff(ptru, lglk_with_data))

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
design_df[which(design_df$ID == suffixes[TEST_SUFFIX] & design_df$CKMR == "Yes" & design_df$Value == "SE"),7:13] <- as.numeric(basic0)

basic0_noCK <- doades( 
  m_SYEL_0,          # could be anything
  Hbits,             # *assume* this p'tic Hbits matches this lglk_with_data
  interesting_stuff, # look at function code
  lglk_with_data,
  comp_wts= c( MOP=0.01, XmHSP=0.01)
)

compdesign[d,l,s,"No", "Yes", ,"SE"]  <- as.numeric(basic0_noCK)
design_df[which(design_df$ID == suffixes[TEST_SUFFIX] & design_df$CKMR == "No" & design_df$Value == "SE"),7:13] <- as.numeric(basic0)

tout <- round(compdesign[d,l,s,,"Yes",,"SE"], digits = 2)

write.csv( tout, file=paste0("results/compdesign_D", d, "_L", l, "_S", s, ".csv"))

cvs <- compdesign[1, 2, 1, "Yes", "Yes",, "SE"]/compdesign[1, 2, 3, "Yes", "Yes",, "Est"]
design_df[which(design_df$ID == suffixes[TEST_SUFFIX] & design_df$CKMR == "Yes" & design_df$Value == "CV"),7:13] <- compdesign[1, 2, 1, "Yes", "Yes",, "SE"]/compdesign[1, 2, 3, "Yes", "Yes",, "Est"]


write.csv( cvs, file=paste0("results/cvs_D", d, "_L", l, "_S", s, ".csv"))

