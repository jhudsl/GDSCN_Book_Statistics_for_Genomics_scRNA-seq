# (PART\*) ANNOTATION {-}

# Marker gene detection

We have managed to identify clusters of related cells in the data. However, these clusters aren't very useful until we can identify the biology meaning of each group. This is where functional annotation comes in.

The real art, and the greatest challenge, in an scRNA-seq analysis comes when interpreting the results. Up to this point (cleaning and clustering the data), the analysis and computation has been straightforward. Figuring out the biological state that each cluster represents, on the other hand, is more difficult, as it requires applying prior biological knowledge to the dataset. 

Thanks to previous research, we know many _marker genes_, or genes can be used to identify particular cell types. These genes are differentially expressed across cell types, and by examining the expression profiles of multiple marker genes across all the clusters, we can assign particular cell type identities to each cluster. 

## Calculating and ranking effect size summary statistics

We begin by comparing each pair of clusters and calculating scores for expression differences between the two for each gene. We have multiple options for the statistics used to compare expression values.



  * **AUC** (area under the curve) is the probability that a randomly chosen observation from cluster A is greater than a randomly chosen observation from cluster B. This statistic is a way to quantify how well we can distinguish between two distributions (clusters) in a pairwise comparison. An AUC of 1 means all values in cluster A are greater than any value from cluster B and suggests upregulation. An AUC of of 0.5 means the two clusters are indistinguishable from each other, while an AUC of 0 suggests the marker gene observations in cluster A are downregulated compared to those in cluster B. 
  
  * **Cohen's _d_** is a standardized log-fold change, and can be thought of as the number of standard deviations that separate the two groups. Positive values of Cohen's _d_ suggest that our cluster of interest (cluster A) are upregulated compared to cluster B, while negative values suggest the marker gene observations in cluster A are downregulated compared to cluster B.
  
  * **log-fold change (logFC)** is a measure of whether there is a difference in expression between clusters. Keep in mind that these values ignore the magnitude of the change. As with the others, positive values indicate upregulation in the cluster of interest (cluster A), while negative values indicate downregulation.

For each of these statistics, `scoreMarkers` calculates mean, median, minimum value (min), maximum value (max), and minimum rank (rank; the smallest rank of each gene across all pairwise comparisons). For most of these measures, a larger number indicates upregulation. For minimum rank, however, a small value means the gene is one of the top upregulated genes.

AUC or Cohen’s _d_ are effective regardless of the magnitude of the expression values and thus are good choices for general marker detection. The log-fold change in the detected proportion is specifically useful for identifying binary changes in expression.

For this exercise, we're going to focus on upregulated markers, since those are particularly useful for identifying cell types in a heterogeneous population like the Zeisel dataset. We use the `findMarkers` command, which quickly identifies potential marker genes.


```r
# identify those genes which are upregulated in some clusters compared to others
#markers <- findMarkers(sce.zeisel.qc, direction = "up")
#marker.set <- markers[["1"]]
#head(marker.set, 10)
```

This dataframe shows us the in log-fold expression change for each potential marker gene between cluster 1 and every other cluster.

## Comparing gene expression levels across clusters

Once we've identified potential marker genes, we can use a heatmap to compare gene expression in each cell between clusters. 


```r
# pull the top most upregulated markers from cluster 1 (compared to the rest of the clusters) and look at their expression in all clusters
#top.markers <- rownames(marker.set)[marker.set$Top <= 10]
#plotHeatmap(sce.zeisel.PCA.1000, features = top.markers, order_columns_by = "label")
```

In this heatmap, clusters are on the horizontal, while the top upregulated genes in cluster 1 are on the vertical. The magnitude of the log-fold expression change is indicated by color of each cell. 

We can also create a heatmap that shows the mean log-fold change of cluster 1 cells compared to the mean of each other cluster. This can simplify the heatmap and is useful when dealing with many clusters.




```r
# AnVIL::install("pheatmap")
#library(pheatmap)

#this heatmap lets us compare the average expression of the gene within a cluster compared to the other clusters
#logFCs <- getMarkerEffects(marker.set[1:50,])
#pheatmap(logFCs, breaks = seq(-5, 5, length.out = 101))
```

Here we see that three genes are generally upregulated in Cluster 1 compared to the other clusters: _Gad1_, _Gad2_, and _Slc6a1_. This is where prior biological knowledge comes in handy, as both _Gad1_ and _Slc6a1_ are known interneuron markers (Zeng et al. 2012).

::: {.fyi}
QUESTION
1. Are there any groups or patterns you see in the second heatmap that look interesting?
:::

# Cell type annotation

Unless you are already an expert in neuronal cell expression, you probably didn't know _Gad1_ and _Slc6a1_ are known interneuron markers until you were told. Luckily, we have several high-quality references databases that can be used for annotating scRNA-seq datasets. Some of these references are useful for identifying cell types (the Zeisel dataset we have been using is a well-known reference for identifying neuronal cell types). Others, such as the Gene Ontology (GO) or Kyoto Encyclopedia of Genes and Genomes (KEGG) collections, are useful for identifying biological processes associated with each cluster.

## Annotating cell types

There are three basic strategies for annotating datasets: match the expression profile of each individual cell to the expression profile of cells from a reference dataset (the reference dataset approach); identify sets of marker genes highly expressed in each cell and match to gene sets from known cell types (the gene set approach); or perform a gene-set enrichment analysis on the marker genes that define each cluster. 

The Zeisel dataset we have been working with has actually been annotated all this time.


```r
# calculate the top marker genes assigned to each cell type (level1class) in the Zeisel dataset
wilcox.z <- pairwiseWilcox(sce.zeisel.qc, sce.zeisel.qc$level1class, 
    lfc = 1, direction = "up")
markers.z <- getTopMarkers(wilcox.z$statistics, wilcox.z$pairs,
    pairwise = FALSE, n = 50)

# look at how many cell-type categories there are, as well as how many cells assigned to each category
lengths(markers.z)
```

```
## astrocytes_ependymal    endothelial-mural         interneurons 
##                   78                   83                  119 
##            microglia     oligodendrocytes        pyramidal CA1 
##                   69                   80                  124 
##         pyramidal SS 
##                  147
```

This dataset has grouped the individual cells into 7 different categories of neuronal subtypes. We can use the gene sets that define these categories to annotate a second brain scRNA-seq dataset from Tasic et al. (2016). (We will not go through all the data cleaning and clustering steps for this new dataset - having made it this far, we trust you can do that!)




```r
#load the new dataset
sce.tasic <- TasicBrainData()
```

We first create the gene set lists using the `GSEABase` package. The `AUCell` package identifies and ranks marker sets highly expressed in each cell using an area under the curve (AUC) approach. We can assign cell type identity in the Tasic dataset by taking the marker set with the highest AUC as the label for the cell.


```r
# AnVIL::install("GSEABase")
library(GSEABase)

#create a dataset that contains just the information about the Zeisel cell-type categories and the marker genes that define each cell-type
all.sets <- lapply(names(markers.z), function(x) {
    GeneSet(markers.z[[x]], setName = x)        
})
all.sets <- GeneSetCollection(all.sets)

# AnVIL::install("AUCell")
library(AUCell)

# rank genes by expression values within each cell
rankings <- AUCell_buildRankings(counts(sce.tasic),
    plotStats = FALSE, verbose = FALSE)

# calculate AUC for each previously-defined marker set (from Zeisel) in the Tasic data
cell.aucs <- AUCell_calcAUC(all.sets, rankings)
results <- t(assay(cell.aucs))
```

After assigning cell type identities to each cluster, a researcher should always verify that the identities make sense. Since the Tasic dataset has also been annotated, we can compare our annotation to the researcher-provided annotation as a sort of sanity check.


```r
# assign cell type identity in the Tasic dataset by assumig the marker set with the top AUC is the proper label
new.labels <- colnames(results)[max.col(results)]

# compare our annotations to the annotations provided by the Tasic dataset
tab <- table(new.labels, sce.tasic$broad_type)
tab
```

```
##                       
## new.labels             Astrocyte Endothelial Cell GABA-ergic Neuron
##   astrocytes_ependymal        43                2                 0
##   endothelial-mural            0               27                 0
##   interneurons                 0                0               759
##   microglia                    0                0                 0
##   oligodendrocytes             0                0                 1
##   pyramidal SS                 0                0                 1
##                       
## new.labels             Glutamatergic Neuron Microglia Oligodendrocyte
##   astrocytes_ependymal                    0         0               0
##   endothelial-mural                       0         0               0
##   interneurons                            2         0               0
##   microglia                               0        22               0
##   oligodendrocytes                        0         0              38
##   pyramidal SS                          810         0               0
##                       
## new.labels             Oligodendrocyte Precursor Cell Unclassified
##   astrocytes_ependymal                             19            4
##   endothelial-mural                                 0            2
##   interneurons                                      0           15
##   microglia                                         0            1
##   oligodendrocytes                                  3            0
##   pyramidal SS                                      0           60
```

As you can see, the labels applied by the researchers and by our annotation mostly match. (Pyramidal SS nerves are primarily glutamatergic, so although the categories are labeled differently, those two categories do indeed match!) The major exception is the oligodendrocyte precursor cells, which our annotation called astrocytes. However, this mismatch isn't as concerning as you might think, once you know that both astrocytes and oligodendrocytes come from the same precursor cell lineage.

::: {.fyi}
QUESTION
1. Should we be concerned when 1 or 2 cells are assigned a different cell type identity by our annotation than by the researchers? Why or why not?
:::


```r
sessionInfo()
```

```
## R version 4.1.3 (2022-03-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.4 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] AUCell_1.16.0               GSEABase_1.56.0            
##  [3] graph_1.72.0                annotate_1.72.0            
##  [5] XML_3.99-0.9                AnnotationDbi_1.56.2       
##  [7] BiocSingular_1.10.0         scran_1.22.1               
##  [9] scater_1.22.0               ggplot2_3.3.5              
## [11] scuttle_1.4.0               scRNAseq_2.8.0             
## [13] SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0
## [15] Biobase_2.54.0              GenomicRanges_1.46.1       
## [17] GenomeInfoDb_1.30.1         IRanges_2.28.0             
## [19] S4Vectors_0.32.4            BiocGenerics_0.40.0        
## [21] MatrixGenerics_1.6.0        matrixStats_0.61.0         
## 
## loaded via a namespace (and not attached):
##   [1] AnnotationHub_3.2.2           BiocFileCache_2.2.1          
##   [3] igraph_1.3.1                  lazyeval_0.2.2               
##   [5] BiocParallel_1.28.3           digest_0.6.29                
##   [7] ensembldb_2.18.4              htmltools_0.5.2              
##   [9] viridis_0.6.2                 fansi_1.0.3                  
##  [11] magrittr_2.0.3                memoise_2.0.1                
##  [13] ScaledMatrix_1.2.0            cluster_2.1.2                
##  [15] limma_3.50.3                  Biostrings_2.62.0            
##  [17] R.utils_2.11.0                prettyunits_1.1.1            
##  [19] colorspace_2.0-3              blob_1.2.3                   
##  [21] rappdirs_0.3.3                ggrepel_0.9.1                
##  [23] xfun_0.26                     dplyr_1.0.8                  
##  [25] crayon_1.5.1                  RCurl_1.98-1.6               
##  [27] jsonlite_1.8.0                glue_1.6.2                   
##  [29] gtable_0.3.0                  zlibbioc_1.40.0              
##  [31] XVector_0.34.0                DelayedArray_0.20.0          
##  [33] scales_1.2.0                  DBI_1.1.2                    
##  [35] edgeR_3.36.0                  Rcpp_1.0.8.3                 
##  [37] viridisLite_0.4.0             xtable_1.8-4                 
##  [39] progress_1.2.2                dqrng_0.3.0                  
##  [41] bit_4.0.4                     rsvd_1.0.5                   
##  [43] metapod_1.2.0                 httr_1.4.2                   
##  [45] ellipsis_0.3.2                R.methodsS3_1.8.1            
##  [47] pkgconfig_2.0.3               sass_0.4.1                   
##  [49] dbplyr_2.1.1                  locfit_1.5-9.5               
##  [51] utf8_1.2.2                    tidyselect_1.1.2             
##  [53] rlang_1.0.2                   later_1.3.0                  
##  [55] munsell_0.5.0                 BiocVersion_3.14.0           
##  [57] tools_4.1.3                   cachem_1.0.6                 
##  [59] cli_3.2.0                     generics_0.1.2               
##  [61] RSQLite_2.2.12                ExperimentHub_2.2.1          
##  [63] evaluate_0.15                 stringr_1.4.0                
##  [65] fastmap_1.1.0                 yaml_2.3.5                   
##  [67] knitr_1.33                    bit64_4.0.5                  
##  [69] purrr_0.3.4                   KEGGREST_1.34.0              
##  [71] AnnotationFilter_1.18.0       sparseMatrixStats_1.6.0      
##  [73] mime_0.12                     R.oo_1.24.0                  
##  [75] xml2_1.3.3                    biomaRt_2.50.3               
##  [77] compiler_4.1.3                beeswarm_0.4.0               
##  [79] filelock_1.0.2                curl_4.3.2                   
##  [81] png_0.1-7                     interactiveDisplayBase_1.32.0
##  [83] statmod_1.4.36                tibble_3.1.6                 
##  [85] bslib_0.3.1                   stringi_1.7.6                
##  [87] GenomicFeatures_1.46.5        lattice_0.20-45              
##  [89] bluster_1.4.0                 ProtGenerics_1.26.0          
##  [91] Matrix_1.4-0                  vctrs_0.4.1                  
##  [93] pillar_1.7.0                  lifecycle_1.0.1              
##  [95] BiocManager_1.30.16           jquerylib_0.1.4              
##  [97] BiocNeighbors_1.12.0          data.table_1.14.2            
##  [99] bitops_1.0-7                  irlba_2.3.5                  
## [101] httpuv_1.6.5                  rtracklayer_1.54.0           
## [103] R6_2.5.1                      BiocIO_1.4.0                 
## [105] bookdown_0.24                 promises_1.2.0.1             
## [107] gridExtra_2.3                 vipor_0.4.5                  
## [109] assertthat_0.2.1              rjson_0.2.21                 
## [111] withr_2.5.0                   GenomicAlignments_1.30.0     
## [113] Rsamtools_2.10.0              GenomeInfoDbData_1.2.7       
## [115] parallel_4.1.3                hms_1.1.1                    
## [117] grid_4.1.3                    beachmat_2.10.0              
## [119] rmarkdown_2.10                DelayedMatrixStats_1.16.0    
## [121] Rtsne_0.16                    shiny_1.7.1                  
## [123] ggbeeswarm_0.6.0              restfulr_0.0.13
```
