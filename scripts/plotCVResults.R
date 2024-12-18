library(tidyr)
library(dplyr)
library(ggplot2)

estdf <- filter(design_df, Value == "Est", CKMR == "Yes") %>%
  select(ID, CKMR, Value, D, L, S, Nfad_2015, Nfad_2020, Nfad_2025) %>%
  pivot_longer(cols = c(Nfad_2015, Nfad_2020, Nfad_2025), names_to = "Year", values_to = "Est") %>%
  mutate(Year = substr(Year, 6, 9))

cvdf <- na.omit(filter(design_df, Value == "CV")) %>%
  select(ID, CKMR, Value, D, L, S, Nfad_2015, Nfad_2020, Nfad_2025) %>%
  pivot_longer(cols = c(Nfad_2015, Nfad_2020, Nfad_2025), names_to = "Year", values_to = "CV") %>%
  mutate(Year = substr(Year, 6, 9))

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

ggplot() +
  geom_point(data = cvdf, aes(x=Year, y = CV, color = as.factor(L), shape = as.factor(CKMR)))+
  facet_grid(D~S, labeller = plot_labeller) +
  ylab("Expected CV")+
  scale_color_discrete(name = "Lethal Samples", labels = c("No", "Yes"))+
  scale_shape_discrete(name = "CKMR", labels = c("Yes", "No")) +
  geom_hline(yintercept = 0.15, linetype = "dashed", color = "grey") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -70,vjust = 0.3, hjust=0))

ggsave(plot = last_plot(), file = "./figures/MegaResults.png", 
       width = 10, height = 6, units = "in")

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
  scale_color_discrete(name = "Lethal Samples", labels = c("No", "Yes"))+
  scale_shape_discrete(name = "CKMR", labels = c("Yes", "No")) +
  geom_hline(yintercept = 0.15, linetype = "dashed", color = "grey") +
  ylab("Expected CV for #Adult Females in 2025")+
  xlab("Sampling Effort (Years)")+
  theme_bw()


ggsave(plot = last_plot(), file = "./figures/EffortVCV.png", 
       width = 8, height = 6, units = "in")
