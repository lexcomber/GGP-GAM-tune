# GGP-GAM-tune
Deeper investigation into Spatially Varying Coefficient models using GAMs with GP smooths

```{r}
> sessionInfo()
R version 4.3.2 (2023-10-31)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.5

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

time zone: Europe/London
tzcode source: internal

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] stgam_0.0.1.3     spdep_1.3-1       spData_2.3.0      GWmodel_2.3-1     Rcpp_1.0.13       sp_2.1-4          robustbase_0.99-1 doParallel_1.0.17
 [9] iterators_1.0.14  foreach_1.5.2     mgcv_1.9-1        nlme_3.1-164      cols4all_0.8      cowplot_1.1.3     sf_1.0-18         spmoran_0.2.2.9  
[17] lubridate_1.9.3   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4       purrr_1.0.2       readr_2.1.5       tidyr_1.3.1       tibble_3.2.1     
[25] ggplot2_3.5.1     tidyverse_2.0.0  

loaded via a namespace (and not attached):
 [1] DBI_1.2.3             deldir_2.0-2          s2_1.1.6              permute_0.9-7         sandwich_3.1-0        rlang_1.1.4          
 [7] magrittr_2.0.3        multcomp_1.4-25       e1071_1.7-14          compiler_4.3.2        png_0.1-8             vctrs_0.6.5          
[13] maps_3.4.2            pkgconfig_2.0.3       wk_0.9.2              crayon_1.5.3          labeling_0.4.3        utf8_1.2.4           
[19] tzdb_0.4.0            spacesXYZ_1.3-0       xfun_0.48             LearnBayes_2.15.1     cluster_2.1.6         R6_2.5.1             
[25] stringi_1.8.4         RColorBrewer_1.1-3    boot_1.3-28.1         knitr_1.48            fields_15.2           zoo_1.8-12           
[31] FNN_1.1.3.2           Matrix_1.6-4          splines_4.3.2         timechange_0.3.0      tidyselect_1.2.1      abind_1.4-8          
[37] vegan_2.6-4           codetools_0.2-19      lattice_0.22-5        intervals_0.15.4      withr_3.0.1           rARPACK_0.11-0       
[43] coda_0.19-4           survival_3.5-7        units_0.8-5           proxy_0.4-27          xts_0.13.1            pillar_1.9.0         
[49] KernSmooth_2.23-22    generics_0.1.3        spacetime_1.3-1       hms_1.1.3             munsell_0.5.1         scales_1.3.0         
[55] class_7.3-22          glue_1.8.0            spatialreg_1.3-1      tools_4.3.2           RSpectra_0.16-1       mvtnorm_1.2-4        
[61] dotCall64_1.1-1       grid_4.3.2            colorspace_2.1-1      cli_3.6.3             spam_2.10-0           fansi_1.0.6          
[67] ggthemes_5.0.0        viridisLite_0.4.2     gtable_0.3.5          DEoptimR_1.1-3        classInt_0.4-10       TH.data_1.1-2        
[73] farver_2.1.2          lifecycle_1.0.4       microbenchmark_1.4.10 MASS_7.3-60.0.1      
```
