library(knitr)
library(kableExtra)
library(dplyr)

# Table A6
LH_expected_CVs <- read.csv("../results/lhappendix.csv") %>%
  select(-X) %>%
  arrange(Sampling.Scenario, Demographic.Scenario, Lethal.Samples, CKMR) %>%
  select(Sampling.Scenario, Demographic.Scenario, Lethal.Samples, CKMR,
         Adult.Female.Survival, Juvenile.Female.Survival,
         P..Adult.Female.in.State.2)


LHcoln <- c("Sampling Scenario", "Demographic Scenario", "Lethal Samples",
            "CKMR", "$\\phi_A$", "$\\phi_J$", "$\\psi_2$")



# Table A7
N_expected_CVs <- read.csv("../results/Nappendix.csv") %>%
  select(-X) %>%
  arrange(Sampling.Scenario, Demographic.Scenario, Lethal.Samples, CKMR) %>%
  select(Sampling.Scenario, Demographic.Scenario, Lethal.Samples, CKMR,
         X2015.Adult.Females, X2020.Adult.Females, X2025.Adult.Females)

Ncoln <- c("Sampling Scenario", "Demographic Scenario", "Lethal Samples",
          "CKMR", "$N_{2015,A}$", "$N_{2020,A}$", "$N_{2025,A}$")



tableset <- function(data, set, coln){

  data_sub <- subset(data, Sampling.Scenario %in% paste0("S", set))

  # just label first instance of sampling/demo scenarios
  data_sub$Demographic.Scenario[duplicated(
    paste0(data_sub$Sampling.Scenario, data_sub$Demographic.Scenario))] <- ""
  data_sub$Sampling.Scenario[duplicated(data_sub$Sampling.Scenario)] <- ""

  rownames(data_sub) <- NULL
  xx <- kable(data_sub, format="latex", col.names=coln, escape=FALSE,
              longtable=TRUE, booktabs=TRUE, linesep="") %>%
    column_spec(c(1:3), width="2cm")

  print(xx)
}
