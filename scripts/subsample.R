
subsample <- function(suffix, samplingscenario){
  
  load(paste0("./simulation/WalrusSamples_", suffix, ".RData"))
  load(paste0("./simulation/AllWalrusSamples_", suffix, ".RData"))
  
  if(samplingscenario == 1){ # leave out 2026 and 2027
    
    samples <- samples[-which(samples$SampY %in% c(2026, 2027)),]
    allsamples <- allsamples[-which(allsamples$SampY %in% c(2026, 2027)),]
    
    save(samples, file = paste0("./simulation/WalrusSamples_", 
                                substr(suffix, 1,6), "S", samplingscenario,
                                substr(suffix, 9, 19), ".RData"))
    save(allsamples, file = paste0("./simulation/AllWalrusSamples_", 
                                substr(suffix, 1,6), "S", samplingscenario, 
                                substr(suffix, 9, 19), ".RData"))
    
  }
  
  if(samplingscenario == 2){ # leave out 2027
    samples <- samples[-which(samples$SampY == 2027),]
    allsamples <- allsamples[-which(allsamples$SampY == 2027),]
    
    save(samples, file = paste0("./simulation/WalrusSamples_", 
                                substr(suffix, 1,6), "S", samplingscenario,
                                substr(suffix, 9, 19), ".RData"))
    save(allsamples, file = paste0("./simulation/AllWalrusSamples_", 
                                   substr(suffix, 1,6), "S", samplingscenario,
                                   substr(suffix, 9, 19), ".RData"))
    
  }
} # end function