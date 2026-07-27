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

universe_baseline <- baseline %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::pull(gene)

up_baseline <- baseline %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::filter(padj < 0.05, log2FoldChange > 1) %>%
  dplyr::pull(gene)
dn_baseline <- baseline %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::filter(padj < 0.05, log2FoldChange < -1) %>%
  dplyr::pull(gene)


baseline_up_go <- clusterProfiler::enrichGO(
  up_baseline,
  org.Hs.eg.db::org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  universe = universe_baseline
)


baseline_dn_go <- clusterProfiler::enrichGO(
  dn_baseline,
  org.Hs.eg.db::org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  universe = universe_baseline
)


universe_controlled <- cell_controlled %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::pull(gene)

up_controlled <- cell_controlled %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::filter(padj < 0.05, log2FoldChange > 1) %>%
  dplyr::pull(gene)
dn_controlled <- cell_controlled %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::filter(padj < 0.05, log2FoldChange < -1) %>%
  dplyr::pull(gene)


controlled_up_go <- clusterProfiler::enrichGO(
  up_controlled,
  org.Hs.eg.db::org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  universe = universe_controlled
)
controlled_dn_go <- clusterProfiler::enrichGO(
  dn_controlled,
  org.Hs.eg.db::org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  universe = universe_controlled
)
# save(baseline, file = "data/baseline_de.RData", compress = "xz")
# save(cell_controlled, file = "data/controlled_de.RData", compress = "xz")

# save(baseline_up_go, file = "data/baseline_up_go.RData", compress = "xz")
# save(baseline_dn_go, file = "data/baseline_dn_go.RData", compress = "xz")

# save(controlled_up_go, file = "data/controlled_up_go.RData", compress = "xz")
# save(controlled_dn_go, file = "data/controlled_dn_go.RData", compress = "xz")
