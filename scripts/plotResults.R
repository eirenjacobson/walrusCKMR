
library(ggplot2)
library(dplyr)
library(tidyr)

# check expected versus observed comparison types

suffixes <- c("Sc00_Sd20241023", "Sc00_Sd20241121", "Sc00_Sd20839572", 
              "Sc00_Sd29571492", "Sc00_Sd76937593", "Sc00_Sd92759173",
              "Sc00_Sd41850183", "Sc00_Sd38519472", "Sc00_Sd35719375",
              "Sc00_Sd57394720") 


for (i in 1:length(suffixes)){
  
  load(paste0("./results/compcheck_", suffixes[i], ".RData"))
  
  if(i == 1){allcomparisons <- compcheck} else {
  
  allcomparisons <- rbind.data.frame(allcomparisons, compcheck)
  
  }}

ggplot(allcomparisons) +
  geom_boxplot(aes(x=Type, y = P)) +
  theme_bw()

all_Dnonprob_lglk <- data.frame(matrix(rep(NA, 6*10), nrow = 10, ncol = 6))
names(all_Dnonprob_lglk) <- names(Dnonprob_lglk)

for (i in 1:length(suffixes)){
  
  load(paste0("./results/Dnonprob_lglk_", suffixes[i], ".RData"))
  
    all_Dnonprob_lglk[i,] <- Dnonprob_lglk
    
}

all_Dnonprob_lglk %>%
  pivot_longer(names_to = "Type", values_to = "Value", cols = 1:6) %>%
  filter(Type != "RoI") %>%
  ggplot() +
  geom_boxplot(aes(x=Type, y = Value)) +
  theme_bw()


