library(DeeDeeExperiment)
library(airway)

data(airway)
airway <- DESeq2::estimateSizeFactors(DESeq2::DESeqDataSet(
  airway,
  design = ~dex
))
airway <- DeeDeeExperiment(airway)

data("baseline_de", package = "exploreSE")
data("controlled_de", package = "exploreSE")

se_with_de <- addDEA(airway, baseline)
se_with_de <- addDEA(se_with_de, cell_controlled)

data("baseline_up_go", package = "exploreSE")
data("baseline_down_go", package = "exploreSE")

data("controlled_up_go", package = "exploreSE")
data("controlled_down_go", package = "exploreSE")


se_with_go <- DeeDeeExperiment::addFEA(se_with_de, baseline_up_go, "baseline")
se_with_go <- DeeDeeExperiment::addFEA(se_with_go, baseline_down_go, "baseline")
se_with_go <- DeeDeeExperiment::addFEA(
  se_with_go,
  controlled_up_go,
  "controlled"
)
se_with_go <- DeeDeeExperiment::addFEA(
  se_with_go,
  controlled_down_go,
  "controlled"
)
