library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(PNWColors)

load("./results/design_df_filled.RData")

lhdf <- filter(design_df, Value == "Est") %>%
  select(ID, CKMR, Value, D, L, S, fadsurv, fjusurv, ppn_breedy) %>%
  pivot_longer(cols = c(fadsurv, fjusurv, ppn_breedy), names_to = "Par", 
               values_to = "Est") %>%
  filter(Par == "fadsurv") %>%
  group_by(CKMR, L) %>%
  summarize(Mean = round(mean(Est, na.rm=TRUE), digits = 2)) 

cvdf <- na.omit(filter(design_df, Value == "CV")) %>%
  select(ID, CKMR, Value, D, L, S, fadsurv, fjusurv, ppn_breedy) %>%
  pivot_longer(cols = c(fadsurv, fjusurv, ppn_breedy), names_to = "Par", 
               values_to = "CV") %>%
  mutate(CKMR = as.factor(CKMR)) %>%
  mutate(L = factor(L, levels=2:1, labels=c("Yes", "No")))


comparecvs <-  na.omit(filter(design_df, Value == "CV")) %>%
  select(CKMR, Value, D, L, S, fadsurv, fjusurv, ppn_breedy) %>%
  pivot_longer(cols = c(fadsurv, fjusurv, ppn_breedy), names_to = "Par", 
               values_to = "CV") %>%
  filter(Par == "ppn_breedy") %>%
  group_by(CKMR) %>%
#  pivot_wider(names_from = CKMR, values_from = CV) %>%
#  mutate(Gain = Yes - No) %>%
  group_by(D) %>%
  summarize(mean(CV))

ggplot(cvdf) +
  geom_point(aes(x=S, y = Est))+
  facet_grid(D ~ Par)


cvdf %>% filter(Par == "fadsurv") %>% 
  group_by(D) %>%
  ggplot() + 
  geom_point(aes(x=S,y=Est,color = L, shape = CKMR)) + 
  facet_wrap(~D)

lhappendix <- na.omit(filter(design_df, Value == "CV")) %>%
  select(D, L, S, CKMR, fadsurv, fjusurv, ppn_breedy) %>%
  mutate(L = factor(L, levels=2:1, labels=c("Yes", "No"))) %>%
  rename("Lethal Samples" = L, "Sampling Scenario" = S,
         "Demographic Scenario" = D) %>%
  mutate(fadsurv = round(fadsurv, digits = 2),
         fjusurv = round(fjusurv, digits = 2),
         ppn_breedy = round(ppn_breedy, digits = 2)) %>%
  rename("Adult Female Survival" = fadsurv, 
         "Juvenile Female Survival" = fjusurv,
         "P. Adult Female in State 2" = ppn_breedy)

write.csv(lhappendix, file = "./results/lhappendix.csv")

Nappendix <- na.omit(filter(design_df, Value == "CV")) %>%
  select(D, L, S, CKMR, Nfad_2015, Nfad_2020, Nfad_2025) %>%
  mutate(L = factor(L, levels=2:1, labels=c("Yes", "No"))) %>%
  rename("Lethal Samples" = L, "Sampling Scenario" = S,
         "Demographic Scenario" = D) %>%
  mutate(Nfad_2015 = round(Nfad_2015, digits = 2),
         Nfad_2020 = round(Nfad_2020, digits = 2),
         Nfad_2025 = round(Nfad_2025, digits = 2)) %>%
  rename("2015 Adult Females" = Nfad_2015, 
         "2020 Adult Females" = Nfad_2020,
         "2025 Adult Females" = Nfad_2025)

write.csv(Nappendix, file = "./results/Nappendix.csv")
