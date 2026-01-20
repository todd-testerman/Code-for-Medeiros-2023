Medeiros et al., 2023
================
Todd Testerman
2023-08-22

Load packages

``` r
library(pacman)
pacman::p_load(decontam, phyloseq, data.table, ggplot2, BiocManager, qiime2R, DESeq2, tidyverse, RColorBrewer, viridis, vegan, pheatmap, patchwork, ggpubr, lme4, nlme, microViz, rstatix, ANCOMBC, microbiome)
```

Import QIIME2 objects and metadata

``` r
phy = qza_to_phyloseq("table.qza", "rooted-tree.qza", "taxonomy.qza","Alzheimers_metadata.txt",tmp = "C:/tmp")
```

Agglomerate at phylum level, normalize, melt object, subset to just
Bacteroidetes and Firmicutes, compile listof group comparisons, build
plot and save

``` r
phy_phylum_glom = tax_glom(phy, "Phylum")
phy_normalize = transform_sample_counts(phy_phylum_glom, function(x) (x / sum(x)))
melted_phy = psmelt(phy_normalize)
melted_phy_subset = subset(melted_phy, Phylum %in% c("Bacteroidetes", "Firmicutes"))
my_comparisons = list(c("Baseline","Control"),c("Control", "Probiotic"), c("Baseline", "Probiotic"))
p2 <- ggplot(melted_phy_subset, aes(x=Cohort, y=Abundance)) + 
  geom_boxplot(outlier.shape = NA) + geom_dotplot(binaxis='y', stackdir='center', dotsize = 0) + theme_bw() + facet_wrap(~Phylum) + geom_point(aes(colour = Cohort), size = 5) + stat_compare_means(comparisons = my_comparisons, method = "t.test", label.x = 1.25, label = "p.signif") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = 20)) + ylab("Relative Abundance")  + guides(colour = guide_legend(override.aes = list(size=10))) + scale_color_discrete(labels = c('Baseline', 'ADC', 'ADP'))
p2
```

    ## Bin width defaults to 1/30 of the range of the data. Pick better value with `binwidth`.

![](rmakrdown_for_github_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggsave("phylum_comp.png", height = 10, width = 12)
```

    ## Bin width defaults to 1/30 of the range of the data. Pick better value with `binwidth`.

``` r
ggsave("phylum_comp.tiff", height = 10, width = 12, compression = "lzw")
```

    ## Bin width defaults to 1/30 of the range of the data. Pick better value with `binwidth`.

Calculate exact p-values for group comparisons

``` r
melted_phy_subset_bacteroidetes = subset(melted_phy, Phylum %in% "Bacteroidetes")
bacteroidetes_pvalue = compare_means(Abundance ~ Cohort, p.adjust.method = "bonferroni", data = melted_phy_subset_bacteroidetes, method = "t.test")
bacteroidetes_pvalue
```

    ## # A tibble: 3 x 8
    ##   .y.       group1    group2        p p.adj p.format p.signif method
    ##   <chr>     <chr>     <chr>     <dbl> <dbl> <chr>    <chr>    <chr> 
    ## 1 Abundance Probiotic Baseline 0.0153 0.046 0.015    *        T-test
    ## 2 Abundance Probiotic Control  0.0730 0.22  0.073    ns       T-test
    ## 3 Abundance Baseline  Control  0.695  1     0.695    ns       T-test

``` r
melted_phy_subset_firmicutes = subset(melted_phy, Phylum %in% "Firmicutes")
firmicutes_pvalue = compare_means(Abundance ~ Cohort, p.adjust.method = "bonferroni", data = melted_phy_subset_firmicutes, method = "t.test")
firmicutes_pvalue
```

    ## # A tibble: 3 x 8
    ##   .y.       group1   group2        p p.adj p.format p.signif method
    ##   <chr>     <chr>    <chr>     <dbl> <dbl> <chr>    <chr>    <chr> 
    ## 1 Abundance Baseline Control   0.308  0.92 0.31     ns       T-test
    ## 2 Abundance Baseline Probiotic 0.103  0.31 0.10     ns       T-test
    ## 3 Abundance Control  Probiotic 0.886  1    0.89     ns       T-test

Shapiro Normality Test to determine if data is normal - *it is normal
for Bacteroidetes and Firmicutes* - *use t-test*

``` r
melted_phy_subset_bacteroidetes_distinct = melted_phy_subset_bacteroidetes %>% distinct(Sample, .keep_all = T)
melted_phy_subset_bacteroidetes_distinct_df = as_data_frame(melted_phy_subset_bacteroidetes_distinct)
```

    ## Warning: `as_data_frame()` was deprecated in tibble 2.0.0.
    ## Please use `as_tibble()` instead.
    ## The signature and semantics have changed, see `?as_tibble`.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

``` r
melted_phy_subset_bacteroidetes_distinct_df %>% group_by(Cohort) %>% shapiro_test(Abundance)
```

    ## # A tibble: 3 x 4
    ##   Cohort    variable  statistic     p
    ##   <chr>     <chr>         <dbl> <dbl>
    ## 1 Baseline  Abundance     0.937 0.515
    ## 2 Control   Abundance     0.910 0.469
    ## 3 Probiotic Abundance     0.976 0.910

``` r
#Firmicutes
melted_phy_subset_firmicutes_distinct = melted_phy_subset_firmicutes %>% distinct(Sample, .keep_all = T)
melted_phy_subset_firmicutes_distinct_df = as_data_frame(melted_phy_subset_firmicutes_distinct)

melted_phy_subset_firmicutes_distinct_df %>% group_by(Cohort) %>% shapiro_test(Abundance)
```

    ## # A tibble: 3 x 4
    ##   Cohort    variable  statistic     p
    ##   <chr>     <chr>         <dbl> <dbl>
    ## 1 Baseline  Abundance     0.933 0.477
    ## 2 Control   Abundance     0.913 0.484
    ## 3 Probiotic Abundance     0.982 0.946

``` r
sessionInfo()
```

    ## R version 4.1.0 (2021-05-18)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19045)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] microbiome_1.14.0           ANCOMBC_1.2.0              
    ##  [3] rstatix_0.7.0               microViz_0.7.6             
    ##  [5] nlme_3.1-152                lme4_1.1-29                
    ##  [7] Matrix_1.3-4                ggpubr_0.4.0               
    ##  [9] patchwork_1.1.1             pheatmap_1.0.12            
    ## [11] vegan_2.6-2                 lattice_0.20-44            
    ## [13] permute_0.9-7               viridis_0.6.2              
    ## [15] viridisLite_0.4.0           RColorBrewer_1.1-3         
    ## [17] forcats_0.5.1               stringr_1.4.0              
    ## [19] dplyr_1.0.9                 purrr_0.3.4                
    ## [21] readr_2.1.2                 tidyr_1.2.0                
    ## [23] tibble_3.1.7                tidyverse_1.3.1            
    ## [25] DESeq2_1.32.0               SummarizedExperiment_1.22.0
    ## [27] Biobase_2.52.0              MatrixGenerics_1.4.0       
    ## [29] matrixStats_0.62.0          GenomicRanges_1.44.0       
    ## [31] GenomeInfoDb_1.28.0         IRanges_2.26.0             
    ## [33] S4Vectors_0.30.0            BiocGenerics_0.38.0        
    ## [35] qiime2R_0.99.6              BiocManager_1.30.17        
    ## [37] ggplot2_3.3.6               data.table_1.14.2          
    ## [39] phyloseq_1.36.0             decontam_1.12.0            
    ## [41] pacman_0.5.1               
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.4.0           backports_1.4.1        Hmisc_4.7-0           
    ##   [4] plyr_1.8.7             igraph_1.3.1           splines_4.1.0         
    ##   [7] BiocParallel_1.26.0    digest_0.6.29          foreach_1.5.2         
    ##  [10] htmltools_0.5.2        fansi_1.0.3            magrittr_2.0.3        
    ##  [13] checkmate_2.1.0        memoise_2.0.1          cluster_2.1.2         
    ##  [16] tzdb_0.3.0             Biostrings_2.60.0      annotate_1.70.0       
    ##  [19] modelr_0.1.8           jpeg_0.1-9             colorspace_2.0-3      
    ##  [22] rvest_1.0.2            blob_1.2.3             rbibutils_2.2.8       
    ##  [25] haven_2.5.0            xfun_0.30              crayon_1.5.1          
    ##  [28] RCurl_1.98-1.6         jsonlite_1.8.0         genefilter_1.74.0     
    ##  [31] survival_3.3-1         iterators_1.0.14       ape_5.6-2             
    ##  [34] glue_1.6.2             gtable_0.3.0           zlibbioc_1.38.0       
    ##  [37] XVector_0.32.0         DelayedArray_0.18.0    car_3.0-13            
    ##  [40] Rhdf5lib_1.14.0        abind_1.4-5            scales_1.2.0          
    ##  [43] DBI_1.1.2              Rcpp_1.0.8.3           xtable_1.8-4          
    ##  [46] htmlTable_2.4.0        foreign_0.8-82         bit_4.0.4             
    ##  [49] Formula_1.2-4          DT_0.22                truncnorm_1.0-8       
    ##  [52] htmlwidgets_1.5.4      httr_1.4.3             ellipsis_0.3.2        
    ##  [55] farver_2.1.0           pkgconfig_2.0.3        NADA_1.6-1.1          
    ##  [58] XML_3.99-0.9           nnet_7.3-17            dbplyr_2.1.1          
    ##  [61] locfit_1.5-9.5         utf8_1.2.2             labeling_0.4.2        
    ##  [64] tidyselect_1.1.2       rlang_1.0.2            reshape2_1.4.4        
    ##  [67] AnnotationDbi_1.54.0   munsell_0.5.0          cellranger_1.1.0      
    ##  [70] tools_4.1.0            cachem_1.0.6           cli_3.3.0             
    ##  [73] generics_0.1.2         RSQLite_2.2.14         ade4_1.7-19           
    ##  [76] broom_0.8.0            evaluate_0.15          biomformat_1.20.0     
    ##  [79] fastmap_1.1.0          yaml_2.3.5             fs_1.5.2              
    ##  [82] knitr_1.39             bit64_4.0.5            KEGGREST_1.32.0       
    ##  [85] xml2_1.3.3             compiler_4.1.0         rstudioapi_0.13       
    ##  [88] png_0.1-7              ggsignif_0.6.3         zCompositions_1.4.0-1 
    ##  [91] reprex_2.0.1           geneplotter_1.70.0     stringi_1.7.6         
    ##  [94] highr_0.9              nloptr_2.0.1           multtest_2.48.0       
    ##  [97] vctrs_0.4.1            pillar_1.7.0           lifecycle_1.0.1       
    ## [100] rhdf5filters_1.4.0     Rdpack_2.3             bitops_1.0-7          
    ## [103] R6_2.5.1               latticeExtra_0.6-29    gridExtra_2.3         
    ## [106] codetools_0.2-18       boot_1.3-28            MASS_7.3-54           
    ## [109] assertthat_0.2.1       rhdf5_2.36.0           withr_2.5.0           
    ## [112] GenomeInfoDbData_1.2.6 mgcv_1.8-36            hms_1.1.1             
    ## [115] grid_4.1.0             rpart_4.1.16           minqa_1.2.4           
    ## [118] rmarkdown_2.14         carData_3.0-5          Rtsne_0.16            
    ## [121] lubridate_1.8.0        base64enc_0.1-3
