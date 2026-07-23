library(DESeq2)
library(DeeDeeExperiment)
library(airway)

data(airway)

airway <- airway[rowData(airway)$gene_biotype == "protein_coding", ]

airway <- DESeqDataSet(airway, design = ~dex)
airway <- DESeq(airway)
baseline <- results(airway)

design(airway) <- ~ cell + dex
airway <- DESeq(airway)
cell_controlled <- results(airway)

dd <- DeeDeeExperiment(airway)

se_with_de <- addDEA(dd, baseline)
se_with_de <- addDEA(se_with_de, cell_controlled)

se_with_go <- get.gos(
  obj = se_with_de,
  NAME = "baseline",
  gene_type = "ENSEMBL"
)
se_with_go <- get.gos(
  obj = se_with_go,
  NAME = "cell_controlled",
  gene_type = "ENSEMBL"
)

# save(se_with_de, file = "data/se_with_des.RData", compress = "xz")
# save(se_with_go, file = "data/se_with_gos.RData", compress = "xz")
