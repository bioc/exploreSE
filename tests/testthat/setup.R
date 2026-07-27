library(DeeDeeExperiment)
library(airway)

data(airway)
airway <- DeeDeeExperiment(airway)

data("baseline_de", package = "exploreSE")
data("controlled_de", package = "exploreSE")

se_with_de <- addDEA(airway, baseline)
se_with_de <- addDEA(se_with_de, cell_controlled)

data("baseline_up_go", package = "exploreSE")
data("baseline_dn_go", package = "exploreSE")

data("controlled_up_go", package = "exploreSE")
data("controlled_dn_go", package = "exploreSE")


se_with_go <- DeeDeeExperiment::addFEA(se_with_de, gos_baseline_up, "baseline")
se_with_go <- DeeDeeExperiment::addFEA(se_with_go, gos_baseline_dn, "baseline")
se_with_go <- DeeDeeExperiment::addFEA(
  se_with_go,
  controlled_up_go,
  "controlled"
)
se_with_go <- DeeDeeExperiment::addFEA(
  se_with_go,
  controlled_dn_go,
  "bascontrolledeline"
)
