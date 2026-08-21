
<!-- README.md is generated from README.Rmd. Please edit that file -->

# exploreSE

<!-- badges: start -->

<!-- badges: end -->

exploreSE is package that provides an interactive Shiny-based user
interface for exploring transcriptional data and analysis results stored
in
[summarizedExperiment](https://www.bioconductor.org/packages/release/bioc/html/tidySummarizedExperiment.html)
or
[DeeDeeExperiment](https://bioconductor.org/packages//release/bioc/html/DeeDeeExperiment.html)
format. The aim is to facilitate easy comparison between different model
approaches on a single data set.

The package can be installed via the `BiocManager` by starting R and
running:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("exploreSE")
```

Once installed, the packages can be accessed through the following bit
of code:

``` r
library(exploreSE)
```

To explore the results of differential expression analysis, we need to
organise them in the structure and you will learn how to do that in this
vignette. If you want to explore on your own, you can run the following
code and start exploring the example data.

``` r
example(exploreSE, ask = FALSE)
```

# Setting up the data

For the purposes of this example, we’ll be using the `airway` data. It
is part of the
[airway](https://www.bioconductor.org/packages/release/data/experiment/html/airway.html)
package, containing RNA-Seq data from different airway smooth muscle
cell lines either untreated or treated with dexamethasone.

``` r
library(airway)
data(airway)
airway
#> class: RangedSummarizedExperiment 
#> dim: 63677 8 
#> metadata(1): ''
#> assays(1): counts
#> rownames(63677): ENSG00000000003 ENSG00000000005 ... ENSG00000273492
#>   ENSG00000273493
#> rowData names(10): gene_id gene_name ... seq_coord_system symbol
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(9): SampleName cell ... Sample BioSample
```

Out of the box, the `airway` object is `summarizedExperiment` object,
containing eight samples from four cell lines. This information is
available in `coldata(airway)`.

``` r
colData(airway)
#> DataFrame with 8 rows and 9 columns
#>            SampleName     cell      dex    albut        Run avgLength
#>              <factor> <factor> <factor> <factor>   <factor> <integer>
#> SRR1039508 GSM1275862  N61311     untrt    untrt SRR1039508       126
#> SRR1039509 GSM1275863  N61311     trt      untrt SRR1039509       126
#> SRR1039512 GSM1275866  N052611    untrt    untrt SRR1039512       126
#> SRR1039513 GSM1275867  N052611    trt      untrt SRR1039513        87
#> SRR1039516 GSM1275870  N080611    untrt    untrt SRR1039516       120
#> SRR1039517 GSM1275871  N080611    trt      untrt SRR1039517       126
#> SRR1039520 GSM1275874  N061011    untrt    untrt SRR1039520       101
#> SRR1039521 GSM1275875  N061011    trt      untrt SRR1039521        98
#>            Experiment    Sample    BioSample
#>              <factor>  <factor>     <factor>
#> SRR1039508  SRX384345 SRS508568 SAMN02422669
#> SRR1039509  SRX384346 SRS508567 SAMN02422675
#> SRR1039512  SRX384349 SRS508571 SAMN02422678
#> SRR1039513  SRX384350 SRS508572 SAMN02422670
#> SRR1039516  SRX384353 SRS508575 SAMN02422682
#> SRR1039517  SRX384354 SRS508576 SAMN02422673
#> SRR1039520  SRX384357 SRS508579 SAMN02422683
#> SRR1039521  SRX384358 SRS508580 SAMN02422677
```

In this vignette, we will subset `airway` to only contain protein-coding
genes. The required information is stored in its `rowData`, accessible
through `rowData(airway)`.

``` r
rowData(airway)
#> DataFrame with 63677 rows and 10 columns
#>                         gene_id     gene_name  entrezid   gene_biotype
#>                     <character>   <character> <integer>    <character>
#> ENSG00000000003 ENSG00000000003        TSPAN6        NA protein_coding
#> ENSG00000000005 ENSG00000000005          TNMD        NA protein_coding
#> ENSG00000000419 ENSG00000000419          DPM1        NA protein_coding
#> ENSG00000000457 ENSG00000000457         SCYL3        NA protein_coding
#> ENSG00000000460 ENSG00000000460      C1orf112        NA protein_coding
#> ...                         ...           ...       ...            ...
#> ENSG00000273489 ENSG00000273489 RP11-180C16.1        NA      antisense
#> ENSG00000273490 ENSG00000273490        TSEN34        NA protein_coding
#> ENSG00000273491 ENSG00000273491  RP11-138A9.2        NA        lincRNA
#> ENSG00000273492 ENSG00000273492    AP000230.1        NA        lincRNA
#> ENSG00000273493 ENSG00000273493  RP11-80H18.4        NA        lincRNA
#>                 gene_seq_start gene_seq_end              seq_name seq_strand
#>                      <integer>    <integer>           <character>  <integer>
#> ENSG00000000003       99883667     99894988                     X         -1
#> ENSG00000000005       99839799     99854882                     X          1
#> ENSG00000000419       49551404     49575092                    20         -1
#> ENSG00000000457      169818772    169863408                     1         -1
#> ENSG00000000460      169631245    169823221                     1          1
#> ...                        ...          ...                   ...        ...
#> ENSG00000273489      131178723    131182453                     7         -1
#> ENSG00000273490       54693789     54697585 HSCHR19LRC_LRC_J_CTG1          1
#> ENSG00000273491      130600118    130603315          HG1308_PATCH          1
#> ENSG00000273492       27543189     27589700                    21          1
#> ENSG00000273493       58315692     58315845                     3          1
#>                 seq_coord_system        symbol
#>                        <integer>   <character>
#> ENSG00000000003               NA        TSPAN6
#> ENSG00000000005               NA          TNMD
#> ENSG00000000419               NA          DPM1
#> ENSG00000000457               NA         SCYL3
#> ENSG00000000460               NA      C1orf112
#> ...                          ...           ...
#> ENSG00000273489               NA RP11-180C16.1
#> ENSG00000273490               NA        TSEN34
#> ENSG00000273491               NA  RP11-138A9.2
#> ENSG00000273492               NA    AP000230.1
#> ENSG00000273493               NA  RP11-80H18.4
```

``` r
dim(airway)
#> [1] 63677     8
airway <- airway[rowData(airway)$gene_biotype == "protein_coding", ]
dim(airway)
#> [1] 22810     8
```

Let’s start our anaylsis. Using the
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
package, we will convert airway into an `DESeqDataSet` and build a
simple model based the dexamethasone stimulation stored in the `dex`
variable.

``` r
library(DESeq2)
airway <- DESeqDataSet(airway, design = ~dex)
airway <- DESeq(airway)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
baseline <- results(airway)
```

In addition, we’ll build a second model, where we control for the effect
of the cell line, stored in the `cell` variable:

``` r
design(airway) <- ~ cell + dex
airway <- DESeq(airway)
#> using pre-existing size factors
#> estimating dispersions
#> found already estimated dispersions, replacing these
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
cell_controlled <- results(airway)
```

Using the
[DeeDeeExperiment](https://bioconductor.org/packages//release/bioc/html/DeeDeeExperiment.html)
package, we will store these results in the their respective slots.

``` r
library(DeeDeeExperiment)
#> Loading required package: SingleCellExperiment
airway <- DeeDeeExperiment(airway)
airway <- addDEA(airway, baseline)
airway <- addDEA(airway, cell_controlled)
```

Let’s take a look:

``` r
airway
#> class: DeeDeeExperiment 
#> dim: 22810 8 
#> metadata(3): '' version singlecontrast
#> assays(4): counts mu H cooks
#> rownames(22810): ENSG00000000003 ENSG00000000005 ... ENSG00000273482
#>   ENSG00000273490
#> rowData names(50): gene_id gene_name ... cell_controlled_pvalue
#>   cell_controlled_padj
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(10): SampleName cell ... BioSample sizeFactor
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> dea(2): baseline, cell_controlled 
#> fea(0):
```

In addition to differential expression anaylsis, we can also explore
some biological data mining, in this case GO ORA. Here we use the
`get.gos()` function from the `exploreSE` package, but it just populates
the FEA slot of the DeeDeeExperiment.

``` r
airway <- get.gos(obj = airway, NAME = "baseline", gene_type = "ENSEMBL")
#> 
#> 
#> Found 5112 gene sets in `enrichResult` object, of which 57 are significant.
#> Converting for usage within the DeeDeeExperiment framework...
#> ✔ Renamed FEA entries: "up_go" to "baseline_up_go"
#> Found 5112 gene sets in `enrichResult` object, of which 77 are significant.
#> Converting for usage within the DeeDeeExperiment framework...
#> ✔ Renamed FEA entries: "dn_go" to "baseline_down_go"
airway <- get.gos(obj = airway, NAME = "cell_controlled", gene_type = "ENSEMBL")
#> Found 5079 gene sets in `enrichResult` object, of which 118 are significant.
#> Converting for usage within the DeeDeeExperiment framework...
#> ✔ Renamed FEA entries: "up_go" to "cell_controlled_up_go"
#> Found 5079 gene sets in `enrichResult` object, of which 32 are significant.
#> Converting for usage within the DeeDeeExperiment framework...
#> ✔ Renamed FEA entries: "dn_go" to "cell_controlled_down_go"
```

``` r
airway
#> class: DeeDeeExperiment 
#> dim: 22810 8 
#> metadata(3): '' version singlecontrast
#> assays(4): counts mu H cooks
#> rownames(22810): ENSG00000000003 ENSG00000000005 ... ENSG00000273482
#>   ENSG00000273490
#> rowData names(50): gene_id gene_name ... cell_controlled_pvalue
#>   cell_controlled_padj
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(10): SampleName cell ... BioSample sizeFactor
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> dea(2): baseline, cell_controlled 
#> fea(4): baseline_up_go, baseline_down_go, cell_controlled_up_go, cell_controlled_down_go
```

Now that our data is ready, we can explore the data.

# Launching the App

The simplest way to launch the app is through a call to the
`exploreSE()` function. Without any arguments, the app opens and you can
load in any .RDS file. Alternatively, you can supply the `file`
argument, defining the path to the file, or the `object` argument when
you have the object already loaded.

``` r
app <- exploreSE(object = airway)
shiny::runApp(app, port = 1234)
```
