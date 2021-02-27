
library(celda)

sce <- readRDS("../Data/sce.rds")

de1 <- differentialExpression(sce,
    useAssay = "decontXcounts",
    c1 = c(15, 16, 17, 18, 19, 20))

de1[Gene %in% c("CD3D", "CD3E", "CD3G"), ]
#    Gene Pvalue  Log2_FC FDR
# 1: CD3D      0 6.200828   0
# 2: CD3E      0 5.565626   0
# 3: CD3G      0 4.046326   0

de2 <- differentialExpression(sce,
    useAssay = "decontXcounts",
    c1 = c(17, 18, 19),
    c2 = c(15, 16, 20))

de2[Gene %in% c("CD8A", "CD8B"), ]
#    Gene        Pvalue  Log2_FC           FDR
# 1: CD8B 5.910581e-176 3.287508 4.978778e-172
# 2: CD8A 7.534401e-169 3.024783 5.077282e-165

de3 <- differentialExpression(sce,
    useAssay = "decontXcounts",
    c1 = c(17),
    c2 = c(18, 19))

de3[Gene %in% c("CCR7"), ]
#    Gene       Pvalue  Log2_FC          FDR
# 1: CCR7 3.415604e-88 3.279663 3.110415e-85

de4 <- differentialExpression(sce,
    useAssay = "decontXcounts",
    c1 = c(18),
    c2 = c(17, 19))

de4[Gene %in% c("GNLY", "KLRG1", "GZMA", "GZMH"), ]
#     Gene        Pvalue  Log2_FC           FDR
# 1:  GZMA 5.602638e-104 4.657639 3.775506e-100
# 2:  GZMH  5.513208e-88 3.599870  2.322025e-84
# 3: KLRG1  1.439579e-48 3.012203  1.054460e-45
# 4:  GNLY  2.251357e-34 2.259016  8.719222e-32

de5 <- differentialExpression(sce,
    useAssay = "decontXcounts",
    c1 = c(1),
    c2 = c(2, 3, 4))

plasmac <- which(celdaClusters(sce) == 1)
#[1] 2904

de5[Gene %in% c("IGHG1", "IGHG3", "IGLC2"), ]
#     Gene       Pvalue  Log2_FC          FDR
# 1: IGHG1 8.906256e-11 1.583893 2.326259e-08
# 2: IGHG3 1.030868e-10 1.483946 2.651455e-08
# 3: IGLC2 5.638171e-04 3.802954 1.860652e-02

counts(sce)[c("CD79A", "CD79B", "CD19", "IGHG1", "IGHG3", "IGLC2", "IGLC3"),
    plasmac]
# CD79A CD79B  CD19 IGHG1 IGHG3 IGLC2
#    10     8     2  1982  1409  8029

sort(counts(sce)[, plasmac], decreasing = TRUE)[1:20]
# IGLC2  IGHG1  IGLC3  IGHG3 MALAT1    B2M  RPL41  RPS18  RPLP1   RPS2  RPL13
#  8029   1982   1480   1409    580    399    359    339    339    269    259
# HLA-B TMSB10  RPL10  RPL21   RPS6  RPS19 RPL18A RPL13A   RPS3
#   239    234    232    220    207    205    199    193    191

sum(counts(sce)[c("IGHG1", "IGHG3", "IGLC2", "IGLC3"), plasmac]) /
    sum(counts(sce)[, plasmac])
#[1] 0.2662923




