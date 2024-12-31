
kinship_types <- list(
  'MOPs' = "Mother-Offspring Pairs",
  'XmHSPs' = "Half-Sibling Pairs",
  'SelfP' = "Self Recaptures"
)

plot_labeller <- function(variable,value){

    return(kinship_types[value])

}

samplesizes %>%
  pivot_longer(cols = c(MOPs, XmHSPs, SelfP), names_to = "Type", values_to = "Number") %>%
  ggplot() +
  geom_point(aes(x=AllSamples, y = Number, group = Type)) +
  labs(x= "Number of Samples Collected", y = "Number of Kin Pairs")+
  ylim(c(0, 1000))+
  facet_wrap(~Type, nrow = 3, labeller = plot_labeller) +
  theme_bw()

ggsave(plot = last_plot(), filename = "figures/samplesizes.png", 
       width = 6, height = 8, units = "in")
