library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(PNWColors)

load("./results/design_df_filled.RData")

estdf <- filter(design_df, Value == "Est", CKMR == "Yes") %>%
  select(ID, CKMR, Value, D, L, S, Nfad_2015, Nfad_2020, Nfad_2025) %>%
  pivot_longer(cols = c(Nfad_2015, Nfad_2020, Nfad_2025), names_to = "Year", values_to = "Est") %>%
  mutate(Year = substr(Year, 6, 9))

cvdf <- na.omit(filter(design_df, Value == "CV")) %>%
  select(ID, CKMR, Value, D, L, S, Nfad_2015, Nfad_2020, Nfad_2025) %>%
  pivot_longer(cols = c(Nfad_2015, Nfad_2020, Nfad_2025), names_to = "Year", values_to = "CV") %>%
  mutate(Year = substr(Year, 6, 9)) %>%
  mutate(CKMR = as.factor(CKMR)) %>%
  mutate(L = factor(L, levels=2:1, labels=c("Yes", "No")))

sedf <- na.omit(filter(design_df, Value == "SE")) %>%
  select(ID, CKMR, Value, D, L, S, Nfad_2015, Nfad_2020, Nfad_2025) %>%
  pivot_longer(cols = c(Nfad_2015, Nfad_2020, Nfad_2025), names_to = "Year", values_to = "SE") %>%
  mutate(Year = substr(Year, 6, 9))

alldf <- left_join(estdf, sedf, by = c("ID", "CKMR", "D", "L", "S", "Year")) %>%
  mutate(LCI = Est - 1.96*SE) %>%
  mutate(UCI = Est + 1.96*SE)


ggplot() +
  geom_point(data = filter(estdf, D == 1), aes(x=Year, y = Est))+
  geom_line(data = filter(estdf, D == 1), aes(x=Year, y = Est))+
  geom_errorbar(data = filter(alldf, D == 1), aes(x=Year, ymin = LCI, ymax = UCI)) +
  facet_grid(~S)

# check mean gain in CV with ICKMR

cvdf %>% 
  pivot_wider(names_from = CKMR, values_from = CV) %>%
  mutate(Difference = Yes - No) %>%
  mutate(Dval = substr(ID, 1, 2)) %>%
  group_by(Dval) %>%
  summarize(MeanD = mean(Difference))

# check mean gain in CV with lethal

cvdf %>% 
  select(-ID) %>%
  pivot_wider(names_from = L, values_from = CV) %>%
  mutate(Difference = Yes - No) %>%
  group_by(interaction(D, CKMR)) %>%
  summarize(MeanD = mean(Difference))

############### Plot of CVs on N across years, sampling, demo scen, lethal, ckmr

sampling_names <- list(
  '1'="+ 100% '25-'26",
  '2'="+ 100% '25-'27",
  '3'="+ 100% '25-'28",
  '4'="+ 75% '25-'26",
  '5'="+ 75% through '27",
  '6'="+ 75% through '28",
  '7'="100% '23-'25"
)

demo_names <- list(
  '1' = "Stable Population",
  '2' = "Decreasing Population",
  '3' = "Increasing Population"
)

plot_labeller <- function(variable,value){
  if (variable=='S') {
    return(sampling_names[value])
  } else {
    return(demo_names[value])
  }
}

pal <- pnw_palette("Starfish", 2)

ggplot() +
  geom_point(data = cvdf, aes(x=Year, y = CV, color = L, shape = CKMR))+
  facet_grid(D~S, labeller = plot_labeller) +
  labs(colour = "Lethal Samples", shape="CKMR", y="Expected CV # Adult Females") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey") +
  theme_bw() +
  scale_color_manual(values = c("#105b58", "#ba7999")) +
  theme(axis.text.x = element_text(angle = -70,vjust = 0.3, hjust=0)) +
  theme(legend.position = "bottom")

ggsave(plot = last_plot(), file = "./figures/MegaResults.png", 
       width = 8.5, height = 5.5, units = "in")

##################################### Plot of sampling effort v cv

effort <- data.frame("S" = c(1:7),
                     "Effort" = c(3, 4, 5, 2.5, 3.25, 4, 3 ))

cveffort <- left_join(cvdf, effort, by = "S")

plot_labeller <- function(variable,value){
  if (variable=='S') {
    return(sampling_names[value])
  } else {
    return(demo_names[value])
  }
}

ggplot(filter(cveffort, Year == 2025)) + 
  geom_point(aes(x=Effort, y = CV, color = as.factor(L), shape = as.factor(CKMR))) +
  facet_wrap(~D, nrow = 3, labeller = plot_labeller) +
  labs(colour = "Lethal Samples", shape="CKMR", y="Expected CV") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("#105b58", "#ba7999")) +
  
  ylab("Expected CV for # Adult Females in 2025")+
  xlab("Sampling Effort (Years)")+
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(plot = last_plot(), file = "./figures/EffortVCV.png", 
       width = 4, height = 5, units = "in")


# compare equal effort scenarios

p1 <- cveffort[which(cveffort$Effort==3),] %>%
  group_by(CKMR, D, L, Year) %>%
  select(CKMR, Value, D, L, S, Year, CV) %>%
  pivot_wider(values_from = CV, names_from = S) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1, color = "grey") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey") +
  geom_point(aes(x=`1`, y=`7`,color = as.factor(L), shape = as.factor(CKMR))) +
  labs(colour = "Lethal Samples", shape="CKMR", y="Expected CV") +
  scale_color_manual(values = c("#105b58", "#ba7999")) +
  
  facet_wrap(~Year)+
  theme_bw()+
  labs(x = "Expected CV on # Adult Females under S1",
       y = "Expected CV on # Adult Females under S7",
       title = "Three Years of Sampling Effort") +
  theme(legend.position = "none")



p2 <- cveffort[which(cveffort$Effort==4),] %>%
  group_by(CKMR, D, L, Year) %>%
  select(CKMR, Value, D, L, S, Year, CV) %>%
  pivot_wider(values_from = CV, names_from = S) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1, color = "grey") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey") +
  geom_point(aes(x=`2`, y=`6`, color = as.factor(L), shape = as.factor(CKMR))) +
  scale_color_manual(values = c("#105b58", "#ba7999")) +
  
  labs(colour = "Lethal Samples", shape="CKMR", y="Expected CV") +
  facet_wrap(~Year)+
  theme_bw()+
  labs(x = "Expected CV on # Adult Females under S2",
       y = "Expected CV on # Adult Females under S6",
       title = "Four Years of Sampling Effort") +
  theme(legend.position = "bottom")

p <- ggarrange(p1, p2, nrow = 2)

ggsave(plot = p, file = "./figures/CVsSameEffort.png", 
       width = 7, height = 7, units = "in")
