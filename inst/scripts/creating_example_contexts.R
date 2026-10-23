# Creates the example context lists:
#   inst/extdata/example_contexts.csv - one context per column, for the upload
#   data/example_contexts.RData       - the same contexts as a named list, for
#                                       exploreSE(context = ...)
# The csv is the source for the R object, so both always contain the same genes.
# Run from the package root.

# dex_response -----
# well-known glucocorticoid response genes, compiled from:
# - human dermal fibroblasts: doi:10.3390/biom11050619
# - human airway smooth muscle cells (the airway data set itself):
#   doi:10.1371/journal.pone.0099625
# - mouse keratinocytes: doi:10.1016/j.mce.2015.05.015
dex_response <- c(
  "FKBP5", "TSC22D3", "DUSP1", "PER1", "KLF15", "DDIT4",
  "ERRFI1", "NFKBIA", "MT2A", "ZFP36", "CRISPLD2"
)

# cardiovascular genetic risk factors -----
# curated from the publications below; duplicated genes within a list removed

# aortic valve stenosis, high priority genes
# Small et al. 2025, doi:10.1038/s41588-025-02417-6
aortic_valve_stenosis_high <- c(
  "ALDH1A2", "CDK8", "HAS1", "LDLR", "LPL", "PALLD", "RNF213", "TPPP", "TRIOBP",
  "ABO", "ACAN", "ACTR2", "ADRB1", "AFAP1", "ALPL", "ANGPTL4", "ATP13A3",
  "CDH2", "CHD4", "CMKLR1", "COL4A2", "COMP", "ERG", "FLNB", "GPC6", "HEATR1",
  "LTBP4", "MTAP", "MUC4", "PALMD", "PCSK6", "PCSK9", "PODXL", "PRMT5", "PSRC1",
  "RGS14", "SDK2", "SLC39A8", "STRN", "SVIL", "TNS1", "TSEN15", "ADAMTS5",
  "ADAMTS7", "ANKH", "APOB", "ARHGAP24", "ATP2B1", "BCAR1", "BCAT1", "CCDC34",
  "CEP83", "DLST", "DST", "DYSF", "ELMO1", "ENPP1", "FAM107B", "FAM118A",
  "FBN2", "FNDC3B", "HEY2", "HNF1A", "IL6", "IL6R", "JMJD1C", "KLF2", "LYPLAL1",
  "MEGF9", "MRC2", "NAV1", "NFATC4", "NOTCH4", "ODF2L", "PDGFA", "PDGFD",
  "PDGFRB", "PLCE1", "PLEC", "RAP1GAP2", "RBM38", "RNF144A", "SEMA3C", "SMAD9",
  "UMPS", "VCAN", "ZFPM1", "ZNF219", "ZPR1"
)

# aortic valve stenosis, high and medium priority genes
# Small et al. 2025, doi:10.1038/s41588-025-02417-6
aortic_valve_stenosis_high_medium <- c(
  "ALDH1A2", "CDK8", "HAS1", "LDLR", "LPL", "PALLD", "RNF213", "TPPP", "TRIOBP",
  "ABO", "ACAN", "ACTR2", "ADRB1", "AFAP1", "ALPL", "ANGPTL4", "ATP13A3",
  "CDH2", "CHD4", "CMKLR1", "COL4A2", "COMP", "ERG", "FLNB", "GPC6", "HEATR1",
  "LTBP4", "MTAP", "MUC4", "PALMD", "PCSK6", "PCSK9", "PODXL", "PRMT5", "PSRC1",
  "RGS14", "SDK2", "SLC39A8", "STRN", "SVIL", "TNS1", "TSEN15", "ADAMTS5",
  "ADAMTS7", "ANKH", "APOB", "ARHGAP24", "ATP2B1", "BCAR1", "BCAT1", "CCDC34",
  "CEP83", "DLST", "DST", "DYSF", "ELMO1", "ENPP1", "FAM107B", "FAM118A",
  "FBN2", "FNDC3B", "HEY2", "HNF1A", "IL6", "IL6R", "JMJD1C", "KLF2", "LYPLAL1",
  "MEGF9", "MRC2", "NAV1", "NFATC4", "NOTCH4", "ODF2L", "PDGFA", "PDGFD",
  "PDGFRB", "PLCE1", "PLEC", "RAP1GAP2", "RBM38", "RNF144A", "SEMA3C", "SMAD9",
  "UMPS", "VCAN", "ZFPM1", "ZNF219", "ZPR1", "ABCG8", "ABL1", "ADAMTS7P1",
  "ADTRP", "AKAP13", "ANKS3", "ANLN", "APCS", "APOH", "ARHGEF26", "ARID1B",
  "ASCC2", "ASXL2", "ATG7", "BAHCC1", "BRCA1", "C11orf58", "CACNA1C", "CCDC188",
  "CDC27", "CDK6", "CDKN1A", "CEP120", "CEP170", "CERS2", "CFAP95", "CFHR2",
  "CHST3", "CIMAP1A", "CNTN2", "COBL", "COLEC11", "CREB5", "CTAGE1", "CUBN",
  "CYP24A1", "CYP2D6", "DAAM1", "DAP", "DDAH1", "DGKB", "DLC1", "DLG5", "EBF1",
  "EBF2", "EIF5AL1", "FADS1", "FAM110B", "FAM168A", "FAM171B", "FARP1", "FAT4",
  "FER", "FES", "FGF23", "FN1", "FTO", "GATA2", "GCKR", "GLI3", "GLIS1", "HAO1",
  "HAUS6", "HHEX", "HMGA2", "HMGCR", "ISL1", "JAG1", "KCNJ2", "KYAT3", "LEF1",
  "LMO4", "LOC100130987", "LOC102723966", "LOC105372043", "LPA", "LRMDA",
  "LRP12", "LRRFIP2", "LYSMD4", "LYZL1", "MAD2L1", "MAF", "MAML2", "MAML3",
  "MAP3K7CL", "MARK3", "MECOM", "MFHAS1", "MFSD6", "MTA1", "MTA3", "MYC",
  "MYEOV", "MYLK4", "MYOCD", "NDP", "NF1", "NFATC1", "NKX2-5", "NKX2-6", "NMBR",
  "NOL4L", "NYAP2", "OR8A1", "PDE3A", "PDE7B", "PDGFRA", "PHACTR1", "PHEX",
  "PHF2", "PIP5K1B", "PLEKHM2", "POGLUT3", "POLR1F", "PPP1R3C", "PPTC7",
  "PRDM16", "PRRX1", "PTK2", "RASSF6", "RNF144B", "RNFT1", "RPS23", "RPS6KA3",
  "S1PR2", "SARM1", "SCAF1", "SCARB1", "SGCB", "SGCG", "SHROOM3", "SIPA1L1",
  "SLC25A26", "SLC35F1", "SLC66A3", "SLC6A15", "SMAD3", "SMC2", "SMG6", "SOCS3",
  "SORBS1", "SORL1", "SPDYE17", "SPRY2", "SREBF1", "ST3GAL6", "STAG1", "STARD9",
  "STN1", "SUFU", "SYN3", "TBX20", "TBX4", "TMED2", "TMEM100", "TMEM18",
  "TNFSF12", "TNS2", "TP53BP2", "TP53INP2", "TRIB1", "TRIM59", "TTC32",
  "UGT2B15", "WARS2", "WASF2", "WWTR1", "ZEB2", "ZFPM2", "ZNF512B", "ZNF827"
)

# bicuspid aortic valve
# Theriault et al. 2026, doi:10.1161/CIRCULATIONAHA.125.074752
bicuspid_aortic_valve <- c(
  "WNT4", "PALMD", "RCSD1", "ARPC5", "STRN", "SLC8A1", "ACTR2", "ZEB2", "ERBB4",
  "PDZRN3", "TMEM44", "MUC4", "LEF1", "FAT4", "MAML3", "PALLD", "SUB1", "OSMR",
  "RPS23", "PRDM6", "MYLK4", "HDGFL1", "VTA1", "LPA", "DEFB136", "ADAMTSL1",
  "CPM", "ATP2B1", "TBX3", "FGF9", "DAOA", "DPF3", "SS18", "KANK2", "RBM38",
  "ARFGAP3", "NOTCH1", "SMAD6", "GATA4", "GATA5", "GATA6", "ADAMTS19", "ROBO4",
  "FBN2", "MAT2A", "MIB1"
)

# thoracic aortic aneurysm
# Klarin et al. 2023, doi:10.1038/s41588-023-01420-z
# Murdock et al. 2026 (medRxiv preprint), doi:10.64898/2026.09.16.26363254
thoracic_aortic_aneurysm <- c(
  "SPSB1", "AC092684.1", "OSR1", "COL6A3", "ULK4", "ZNF827", "PRDM6",
  "CTD-2337A12.1", "TBX20", "ELN", "SCAI", "TCF7L2", "PLCE1", "ADAMTS8", "LRP1",
  "AL354828.1", "DLEU1", "FBN1", "THSD4", "CDH13", "SMG6", "ACTA2", "SMAD3",
  "TDGF1", "VDAC2", "COL3A1", "TDRD6", "FNDC3B", "ZSCAN22", "PCDH1", "SART3",
  "USP54", "ROCK1", "PKD1", "MUS81", "TGFBR2", "ARHGEF12", "RNMT", "INTS8",
  "PARVG", "MOGS", "URM1", "UBA6", "RPL7", "BAG1", "KRT24", "CHSY1", "SLFN11",
  "MYLK", "KRT79", "DCBLD2", "PSAP", "DDX19A", "RGS20", "RNF39", "ARFGEF1",
  "OXCT1", "ARHGAP23", "CCDC80", "RUSC2", "FLNA", "STK33", "DENND2A", "FGF1",
  "E2F6", "GAS2", "COL16A1", "MYCBP2", "KIT", "NAP1L4", "F10", "IMPG2", "TRMU",
  "KRT86", "HDAC5", "USP34", "DRAM1", "NELL1", "ENPP1", "NRXN1", "KCNS2",
  "TAF3", "TMEM63A", "RPL18", "GPHN", "GABBR1", "CLEC16A", "MTHFD2", "SNX33",
  "RXRA", "RPUSD2", "ITCH", "PSMD13", "MEIOC", "CREM", "HHIP", "SPAG9", "CHD6",
  "VCAN", "KMT5B"
)

# abdominal aortic aneurysm
# Roychowdhury et al. 2023, doi:10.1038/s41588-023-01510-y
abdominal_aortic_aneurysm <- c(
  "SPSB1", "MFAP2", "PCSK9", "TGFBR3", "CELSR2", "ECM1", "IL6R", "DNM3",
  "SMYD2", "TGFB2", "HLX", "GDF7", "MCFD2", "ANTXR1", "TGFBRAP1", "SH3RF3",
  "IL1F10", "IFIH1", "MAP3K20", "IRS1", "NICN1", "ZBTB38", "LRRC34", "LDB2",
  "AC058822.1", "SPP1", "ADH1C", "NFKB1", "ZNF827", "TERT", "FST", "C5orf67",
  "HMGCR", "LHFPL2", "VCAN", "SLCO4C1", "FER", "PRDM6", "CDKN1A", "VEGFA",
  "COL21A1", "TUBE1", "CITED2", "IGF2R", "THBS2", "SUGCT", "CDK6", "PINX1",
  "LPL", "KCTD9", "EPHX2", "FBXO32", "TRIB1", "AL359922.1", "MUSK", "DAB2IP",
  "ZNF365", "PLAU", "TSPAN14", "LIPA", "PLCE1", "STN1", "TCF7L2", "ODF3",
  "NUCB2", "BDNF", "DGKZ", "KCTD14", "NOX4", "YAP1", "RDX", "ZPR1", "ADAMTS8",
  "AEBP2", "HDAC7", "LRP1", "MON2", "SCARB1", "FGF9", "COL4A1", "F7", "DAAM1",
  "EIF2B2", "HHIPL1", "FBN1", "ADAM10", "RORA", "FBXL22", "SMAD3", "ARIH1",
  "LOXL1", "CHRNA5", "MCTP2", "CETP", "WWP2", "ZFHX3", "CFDP1", "CMIP",
  "CRISPLD2", "RAB34", "TUBD1", "MRC2", "ABCA6", "RBBP8", "SMAD4", "PLEKHJ1",
  "ADAMTS10", "LDLR", "LTBP4", "PLAUR", "APOE", "PRKD2", "JAG1", "PLTP",
  "PREX1", "ZBTB46", "ERG", "TXNRD2", "TTC28", "CCDC157", "H6PD", "TDRD15",
  "PDGFRA", "SLC22A2", "INHBA", "JMJD1C", "MMP12", "TGFB3", "CTAGE1", "PSMG1"
)

# increased ascending aorta diameter
# Pirruccello et al. 2022, doi:10.1038/s41588-021-00962-4
ascending_aorta_diameter <- c(
  "VTI1A", "CACNB2", "ARID5B", "JMJD1C", "NOC3L", "ENTPD1", "NCAM1", "ADAMTS8",
  "ABCC8", "DCDC1", "ANO1", "C12orf75", "MED13L", "RERG", "ABCC9", "USP15",
  "LLPH", "GALNT4", "CRADD", "FGF9", "DLEU1", "ASB2", "DISP2", "FBN1", "THSD4",
  "GNAO1", "WWP2", "HYDIN", "CDH13", "CBFA2T3", "MAP2K4", "SMG6", "ATAD5",
  "COPRS", "CD68", "PIEZO2", "CCNE1", "HMCN1", "RYR2", "SF3A3", "EDN2", "FGGY",
  "SPSB1", "SLX4IP", "SLC24A3", "KCNE2", "TNRC6B", "ZEB2", "MBD5", "FIGN",
  "OSR1", "COL6A3", "GATA2", "FGD5", "MASP1", "ULK4", "FLNB", "LRIG1", "HAND2",
  "FGF5", "PRDM6", "KCNMB1", "ATP6AP1L", "PCSK1", "ADTRP", "EDN1", "AIG1",
  "ESR1", "CDKN1A", "COL21A1", "TBX20", "ELN", "SEMA3D", "ANGPT1", "HAS2",
  "FBXO32", "TRAPPC9", "SLC25A37", "PI15", "MSRA", "LPAR1", "ENPP2"
)

# increased descending aorta diameter
# Pirruccello et al. 2022, doi:10.1038/s41588-021-00962-4
descending_aorta_diameter <- c(
  "SVIL", "PLCE1", "PCSK7", "ADAMTS8", "ALDH2", "STAT6", "ARHGEF40", "MYH6",
  "NFATC4", "OTUB2", "DISP2", "FBN1", "LOXL1", "ADAMTS7", "WWP2", "TBC1D16",
  "CTIF", "DOT1L", "ACTN4", "SPTBN4", "FBXO46", "NAV1", "PKN2", "JAG1",
  "CCDC141", "GDF7", "COL6A3", "ZBTB20", "MASP1", "ULK4", "FLNB", "LMCD1",
  "FIP1L1", "AFAP1", "FGF5", "PDLIM5", "FER", "NKX2-5", "MICA", "CDKN1A",
  "NT5E", "ZNF292", "TWIST1", "TBX20", "MAL2", "PI15", "FHL1"
)

# coronary artery disease, high priority genes
# Aragam et al. 2022, doi:10.1038/s41588-022-01233-6
coronary_artery_disease_high <- c(
  "LDLR", "NOS3", "PCSK9", "ANGPTL4", "LIPA", "LPL", "ABCG8", "APOB", "APOE",
  "CCM2", "CDK8", "FN1", "HNF1A", "PNPLA3", "REST", "SCARB1", "ACVRL1", "APOC3",
  "ARHGAP42", "ARHGEF26", "CETP", "COL4A1", "CXCR4", "EDNRA", "FES", "GUCY1A3",
  "IL6R", "LAMB2", "LIPC", "LOX", "MC4R", "MECOM", "MYH11", "NME7", "PHACTR1",
  "PLG", "RRP1B", "SEMA3C", "SGCD", "SREBF1", "SWAP70", "TBX20", "TCF21",
  "TGFB1", "AAK1", "ADAMTS7", "AFAP1L2", "ARNTL", "ARVCF", "ATP2B1", "ATXN2",
  "BCAR1", "BCAS3", "BMP1", "C1S", "CABIN1", "CDH13", "CENPW", "CFTR", "CXCL12",
  "DAB2IP", "DOCK8", "F10", "FBN2", "FER", "FGD6", "FHL3", "FLT1", "GEM",
  "GIGYF2", "GPSM1", "HTRA1", "ITGA1", "JCAD", "LIPG", "LMAN1", "LOXL4", "MAFB",
  "MAP1S", "MMP13", "MRAS", "MYO9B", "NBEAL1", "NCOA6", "NGF", "NRP1", "OSM",
  "PLCE1", "PLCG2", "PLTP", "PPAP2B", "PRDM16", "PTRF", "RRBP1", "SCUBE1",
  "SEMA5A", "SERPINA1", "SHROOM3", "SLC2A12", "SMAD3", "SMAD7", "SMG6", "SYN3",
  "TBX3", "TFAP2B", "TNS1", "TRIB1", "TRIM5", "UHRF1BP1", "UNC13D", "VEGFA",
  "VPS11", "WT1", "XPO6", "ZKSCAN1", "ZNF652", "ZNF827", "B9D2", "CCDC30",
  "EIF2B2", "LEPRE1", "LPA", "MLH3"
)

# coronary artery disease, high and medium priority genes
# Aragam et al. 2022, doi:10.1038/s41588-022-01233-6
coronary_artery_disease_high_medium <- c(
  "LDLR", "NOS3", "PCSK9", "ANGPTL4", "LIPA", "LPL", "ABCG8", "APOB", "APOE",
  "CCM2", "CDK8", "FN1", "HNF1A", "PNPLA3", "REST", "SCARB1", "ACVRL1", "APOC3",
  "ARHGAP42", "ARHGEF26", "CETP", "COL4A1", "CXCR4", "EDNRA", "FES", "GUCY1A3",
  "IL6R", "LAMB2", "LIPC", "LOX", "MC4R", "MECOM", "MYH11", "NME7", "PHACTR1",
  "PLG", "RRP1B", "SEMA3C", "SGCD", "SREBF1", "SWAP70", "TBX20", "TCF21",
  "TGFB1", "AAK1", "ADAMTS7", "AFAP1L2", "ARNTL", "ARVCF", "ATP2B1", "ATXN2",
  "BCAR1", "BCAS3", "BMP1", "C1S", "CABIN1", "CDH13", "CENPW", "CFTR", "CXCL12",
  "DAB2IP", "DOCK8", "F10", "FBN2", "FER", "FGD6", "FHL3", "FLT1", "GEM",
  "GIGYF2", "GPSM1", "HTRA1", "ITGA1", "JCAD", "LIPG", "LMAN1", "LOXL4", "MAFB",
  "MAP1S", "MMP13", "MRAS", "MYO9B", "NBEAL1", "NCOA6", "NGF", "NRP1", "OSM",
  "PLCE1", "PLCG2", "PLTP", "PPAP2B", "PRDM16", "PTRF", "RRBP1", "SCUBE1",
  "SEMA5A", "SERPINA1", "SHROOM3", "SLC2A12", "SMAD3", "SMAD7", "SMG6", "SYN3",
  "TBX3", "TFAP2B", "TNS1", "TRIB1", "TRIM5", "UHRF1BP1", "UNC13D", "VEGFA",
  "VPS11", "WT1", "XPO6", "ZKSCAN1", "ZNF652", "ZNF827", "B9D2", "CCDC30",
  "EIF2B2", "LEPRE1", "LPA", "MLH3", "ABHD2", "ABO", "ACER2", "ACVR2A",
  "ARHGAP15", "ARHGAP21", "ARHGEF12", "ARID1A", "ATP1B3", "C5", "CAMK1D", "CEL",
  "CELSR2", "CMIP", "CNNM2", "CNPY2", "CRYAB", "CTAGE1", "CYP17A1", "DAGLB",
  "DMXL1", "DNAH10", "DOCK5", "DOCK9", "DYNLRB2", "FAM46A", "FDX1", "FGD5",
  "FGF1", "FHL5", "FIGN", "FNDC3B", "FOXC1", "GATA6", "GGCX", "GNPDA2", "HBS1L",
  "HDAC9", "HDGFL1", "HHAT", "HHIPL1", "HMGCR", "IRS1", "ITGB3", "KCNK5",
  "KIAA0040", "KLF4", "MAD1L1", "MAD2L1", "MAGI3", "MAP3K1", "MFGE8", "MIA3",
  "MORF4L", "MTAP", "N4BP2L2", "OSR1", "PALLD", "PARP12", "PAWR", "PDE1A",
  "PDE1C", "PDE3A", "PDE5A", "PDGFD", "PECAM1", "PID1", "PIP4K2B", "PLEKHG1",
  "POLG", "PPP2R3A", "PRKCE", "PTK7", "RAC1", "RASGEF1B", "RDX", "RGL3",
  "SAP130", "SERPINH1", "SIPA1", "SKI", "SKP2", "SLC22A1", "SLC5A3", "SORT1",
  "ST3GAL4", "STAT2", "TARID", "TCF7L2", "TIPARP", "TSPAN14", "TWIST1", "UMPS",
  "WIPI1", "ZC3HC1", "ZEB2", "ZNF781"
)

# genes associated with familial/monogenic aortic diseases
# TODO: source to be added
monogenic_aortic_disease <- c(
  "FBN1", "TGFBR1", "TGFBR2", "SMAD3", "TGFB2", "TGFB3", "SMAD2", "COL3A1",
  "SKI", "SLC2A10", "BGN", "EFEMP2", "FBLN5", "LTBP4", "LTBP2", "LTBP3", "PKD1",
  "PKD2", "FLNA", "COL1A2", "ACTA2", "MYH11", "MYLK", "PRKG1", "LOX", "MAT2A",
  "FOXE3", "THSD4", "MFAP5", "NOTCH1", "SMAD6", "GATA4", "GATA5", "GATA6",
  "ADAMTS19", "ROBO4", "FBN2", "MIB1", "ELN", "LDLR", "APOB", "PCSK9", "LPA",
  "NKX2-5"
)

# write the csv, padding the shorter columns with empty cells -----
contexts <- list(
  dex_response = dex_response,
  aortic_valve_stenosis_high = aortic_valve_stenosis_high,
  aortic_valve_stenosis_high_medium = aortic_valve_stenosis_high_medium,
  bicuspid_aortic_valve = bicuspid_aortic_valve,
  thoracic_aortic_aneurysm = thoracic_aortic_aneurysm,
  abdominal_aortic_aneurysm = abdominal_aortic_aneurysm,
  ascending_aorta_diameter = ascending_aorta_diameter,
  descending_aorta_diameter = descending_aorta_diameter,
  coronary_artery_disease_high = coronary_artery_disease_high,
  coronary_artery_disease_high_medium = coronary_artery_disease_high_medium,
  monogenic_aortic_disease = monogenic_aortic_disease
)
n_rows <- max(lengths(contexts))
context_table <- as.data.frame(
  lapply(contexts, \(x) c(x, rep("", n_rows - length(x))))
)
utils::write.csv(
  context_table,
  "inst/extdata/example_contexts.csv",
  row.names = FALSE
)

# build the R object from the csv -----
example_contexts <- utils::read.csv(
  "inst/extdata/example_contexts.csv",
  check.names = FALSE,
  colClasses = "character"
) |>
  as.list() |>
  lapply(\(x) x[!is.na(x) & x != ""])

save(
  example_contexts,
  file = "data/example_contexts.RData",
  compress = "xz"
)
