library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(PNWColors)

load("./results/design_df_filled.RData")

# take out scenario 9
design_df <- design_df[-which(design_df$S == "S9"),]

# filter to keep L1 and L3
design_df <- design_df %>% filter(L == "L1" | L == "L3")

estdf <- filter(design_df, Value == "Est", CKMR == "Yes") %>%
  select(ID, CKMR, Value, D, L, S, Nfad_2015, Nfad_2020, Nfad_2025) %>%
  pivot_longer(cols = c(Nfad_2015, Nfad_2020, Nfad_2025), names_to = "Year", values_to = "Est") %>%
  mutate(Year = substr(Year, 6, 9))

cvdf <- na.omit(filter(design_df, Value == "CV")) %>%
  select(ID, CKMR, Value, D, L, S, Nfad_2015, Nfad_2020, Nfad_2025) %>%
  pivot_longer(cols = c(Nfad_2015, Nfad_2020, Nfad_2025), names_to = "Year", values_to = "CV") %>%
  mutate(Year = substr(Year, 6, 9)) %>%
  mutate(CKMR = as.factor(CKMR)) %>%
  mutate(L = factor(L, levels=paste0("L", 4:1), labels=c("1000", "500", "100", "0")))

sedf <- na.omit(filter(design_df, Value == "SE")) %>%
  select(ID, CKMR, Value, D, L, S, Nfad_2015, Nfad_2020, Nfad_2025) %>%
  pivot_longer(cols = c(Nfad_2015, Nfad_2020, Nfad_2025), names_to = "Year", values_to = "SE") %>%
  mutate(Year = substr(Year, 6, 9))

alldf <- left_join(estdf, sedf, by = c("ID", "CKMR", "D", "L", "S", "Year")) %>%
  mutate(LCI = Est - 1.96*SE) %>%
  mutate(UCI = Est + 1.96*SE)

# RESULTS TEXT 

# check mean gain in CV with CKMR
# mean decrease in CV on Nfad w/ w/o CKMR, by Demo
cvdf %>% 
  pivot_wider(names_from = CKMR, values_from = CV) %>%
  mutate(Difference = Yes - No) %>%
  mutate(Dval = substr(ID, 1, 2)) %>%
  group_by(Dval) %>%
  summarize(MeanD = mean(Difference)) %>%
  mutate(MeanD = round(MeanD, digits = 2))

# check mean gain in CV with lethal

cvdf %>% 
  select(-ID) %>%
  pivot_wider(names_from = L, values_from = CV) %>%
  mutate(Difference = `500` - `0`) %>%
  #  group_by(interaction(D, CKMR)) %>%
  group_by(CKMR) %>%
  summarize(MeanD = round(mean(Difference), digits = 4))

###




ggplot() +
  geom_point(data = filter(estdf, D == "D1"), aes(x=Year, y = Est))+
  geom_line(data = filter(estdf, D == "D1"), aes(x=Year, y = Est))+
  geom_errorbar(data = filter(alldf, D == "D1"), aes(x=Year, ymin = LCI, ymax = UCI)) +
  facet_grid(~S)

############### Plot of CVs on N across years, sampling, demo scen, lethal, ckmr

sampling_names <- list(
  '1'="+ 100% '25",
  '2'="+ 100% '25-'26",
  '3'="+ 100% '25-'27",
  '4'="+ 100% '25-'28",
  '5'="+ 75% '25",
  '6'="+ 75% '25-'26",
  '7'="+ 75% '25-'27",
  '8'="+ 75% '25-'28")

demo_names <- list(
  '1' = "Stationary Population",
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
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey") +
  geom_point(data = cvdf, aes(x=Year, y = CV, color = L, shape = CKMR))+
  facet_grid(D~S, labeller = plot_labeller) +
  labs(colour = "Lethal Samples", shape="CKMR", y="Expected CV # Adult Females") +
  theme_bw() +
  scale_color_manual(values = c("#105b58", "#ba7999")) +
  theme(axis.text.x = element_text(angle = -70,vjust = 0.3, hjust=0)) +
  theme(legend.position = "bottom")

ggsave(plot = last_plot(), file = "./figures/MegaResults.png", 
       width = 8.5, height = 5.5, units = "in")

#### FIgure 2 with subset of scenarios


sampling_ids <- list(
  '1'="S1",
  '2'="S2",
  '3'="S3",
  '4'="S4",
  '5'="S5",
  '6'="S6",
  '7'="S7",
  '8'="S8")

plot_labeller <- function(variable,value){
  if (variable=='S') {
    return(sampling_names[value])
  } else {
    return(sampling_ids[value])
  }
}

sub.df <- cvdf %>%
  filter(D=="D1" & S %in% c("S1", "S2", "S3", "S4") ) %>%
  mutate(Group = factor(paste0(L, CKMR), levels = c("0No", "0Yes", "500No", "500Yes"))) %>%
  mutate(S2 = S)

ggplot() +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey") +
  geom_point(data = sub.df, 
                           aes(x=Year, y = CV, color = Group, shape = Group, group = Group))+
  facet_wrap(S2~S, labeller = plot_labeller, nrow = 1) +
  labs(colour = "Lethal Samples", shape="CKMR", y="Expected CV # Adult Females") +
  theme(axis.text.x = element_text(angle = -70,vjust = 0.3, hjust=0)) +
  scale_colour_manual(name = NULL,
                      labels = c("IMR w/o Lethal", "ICKMR w/o Lethal", "IMR w/ Lethal", "ICKMR w/ Lethal"),
                      values = c(c("#ba7999", "#ba7999"), c("#105b58", "#105b58"))) +   
  scale_shape_manual(name = NULL,
                     labels = c("IMR w/o Lethal", "ICKMR w/o Lethal", "IMR w/ Lethal", "ICKMR w/ Lethal"),
                     values = c(19, 17, 19, 17))+
#  theme(legend.position = "bottom") +
  ylim(c(0, NA))+
  theme_bw()

ggsave(plot = last_plot(), file = "./figures/ResultsS1toS4Stationary.png", 
       width = 8.5, height = 4, units = "in")

##################################### Plot of sampling effort v cv

sampling_names <- list(
  '2015'="2015",
  '2020'="2020",
  '2025'="2025")

demo_names <- list(
  '1' = "Stationary Population",
  '2' = "Decreasing Population",
  '3' = "Increasing Population"
)

plot_labeller <- function(variable,value){
  if (variable=='Year') {
    return(sampling_names[value])
  } else {
    return(demo_names[value])
  }
}

effort <- data.frame("S" = c("S1", "S2", "S3", "S4", 
                             "S5", "S6", "S7", "S8", "S9"),
                     "Effort" = c(2, 3, 4, 5, 1.75, 2.5, 3.25, 4, 3))

cveffort <- cvdf %>%
  left_join(effort, by = "S") %>%
  mutate(Group = factor(paste0(L, CKMR), levels = c("0No", "0Yes", "500No", "500Yes")))
  

ggplot(filter(cveffort, D == "D1" & Year == 2025)) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey") +
  geom_point(aes(x=Effort, y = CV, color = Group, shape = Group, group = Group)) +
  scale_colour_manual(name = NULL,
                      labels = c("IMR w/o Lethal", "ICKMR w/o Lethal", "IMR w/ Lethal", "ICKMR w/ Lethal"),
                      values = c(c("#ba7999", "#ba7999"), c("#105b58", "#105b58"))) +   
  scale_shape_manual(name = NULL,
                     labels = c("IMR w/o Lethal", "ICKMR w/o Lethal", "IMR w/ Lethal", "ICKMR w/ Lethal"),
                     values = c(19, 17, 19, 17))+
#  facet_grid(D~Year, labeller = plot_labeller) +
  labs(colour = "Lethal Samples", shape="CKMR", y="Expected CV") +
  ylab("Expected CV for # Adult Females in 2025")+
  xlab("Total Future Survey Effort (Years)")+
  ylim(c(0, NA))+
  theme_bw() 

ggsave(plot = last_plot(), file = "./figures/EffortVCV.png", 
       width = 8.5, height = 4, units = "in")


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
  geom_point(aes(x=`S3`, y=`S8`, color = as.factor(L), shape = as.factor(CKMR))) +
  scale_color_manual(values = c("#105b58", "#ba7999")) +
  labs(colour = "Lethal Samples", shape="CKMR", y="Expected CV") +
  facet_wrap(~Year)+
  theme_bw()+
  labs(x = "Expected CV on # Adult Females under S2",
       y = "Expected CV on # Adult Females under S6") +
  theme(legend.position = "bottom")

ggsave(plot = p2, file = "./figures/CVw4YrsEff.png", 
        width = 6, height = 4, units = "in")

#p <- ggarrange(p1, p2, nrow = 2)

#ggsave(plot = p, file = "./figures/CVsSameEffort.png", 
#       width = 7, height = 7, units = "in")
