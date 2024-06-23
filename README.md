# scGeneXpress : R package to detect expressed genes

scGeneXpress is a reference-based method to detect expressed genes in single cell RNAseq data. The method quantify gene expression levels in three categories including low, medium and high, to classifying whether genes are expressed in a population or not.

The gene quantification algorithm is a single-cell adaptation of RefBool [1] which relies on bulk RNA-seq data.

*[1] RefBool: a reference-based algorithm for discretizing gene expression data, S. Jung et al. - Bioinformatics, 2017*


# Guidelines


## Install scGeneXpress

```{r}
#Install devtools
install.packages("devtools")
library(devtools)

#Install scGeneXpress from GitHub
devtools::install_github("saschajung/scGeneXpress")
```

## Load the necesarry libraries

```{r}
#Load scGeneXpress
library(scGeneXpress)

#Load dependencies
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  require(foreach)
  require(parallel)
  require(doParallel)
  require(reshape2)
  require(tictoc)
  require(progressr)
  require(scuttle)
  require(Matrix)
  require(sctransform)
  require(LaplacesDemon)
  library(Hmisc)
})
```

## Identity genes identification

```{r}
#Load your raw scRNA-seq matrix data: cells in columns, genes in rows
# mtx =
#Load your metadata: two columns (1) Cell ID, (2) Cluster name
# meta =
#Specify a directory to save the results
# d = "~/Desktop/"
#Prefix to give for the results file
# resName = "Example"
#Provide the link to the background to use (cell type, subtype or phenotype) 
# pb = "~/Desktop/backgrounds/human/"

#Run scGeneXpress
run_scGeneXpress(data = mtx,
            metadata = meta,
            org="human", #change if mouse
            dir = d,
            file.name = resName,
            precBack = pb,
            discretize = T,
            sig.frames = T,
            ncores=2,
            fixseed=1234)
```

## R Session

```{r}
R version 4.0.3 (2020-10-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /mnt/irisgpfs/users/cbarlier/libraries/miniconda3/envs/r4-base-clone9/lib/libopenblasp-r0.3.12.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] sctransform_0.3.2           qlcMatrix_0.9.7            
 [3] sparsesvd_0.2               slam_0.1-48                
 [5] quminorm_0.1.0              scater_1.18.6              
 [7] ggplot2_3.3.5               SingleCellExperiment_1.12.0
 [9] SummarizedExperiment_1.20.0 Biobase_2.50.0             
[11] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
[13] IRanges_2.24.1              S4Vectors_0.28.1           
[15] BiocGenerics_0.36.1         MatrixGenerics_1.2.1       
[17] matrixStats_0.59.0          Matrix_1.3-4               
[19] stringr_1.4.0               data.table_1.14.0          
[21] progressr_0.9.0             tictoc_1.0.1               
[23] reshape2_1.4.4              doParallel_1.0.16          
[25] iterators_1.0.13            foreach_1.5.1              

loaded via a namespace (and not attached):
 [1] viridis_0.6.1             BiocSingular_1.6.0       
 [3] VGAM_1.1-5                viridisLite_0.4.0        
 [5] splines_4.0.3             DelayedMatrixStats_1.12.3
 [7] scuttle_1.0.4             sads_0.4.2               
 [9] assertthat_0.2.1          GenomeInfoDbData_1.2.4   
[11] vipor_0.4.5               globals_0.14.0           
[13] numDeriv_2016.8-1.1       pillar_1.6.1             
[15] lattice_0.20-41           glue_1.4.2               
[17] beachmat_2.6.4            bbmle_1.0.24             
[19] digest_0.6.27             XVector_0.30.0           
[21] colorspace_2.0-2          plyr_1.8.6               
[23] pkgconfig_2.0.3           listenv_0.8.0            
[25] zlibbioc_1.36.0           purrr_0.3.4              
[27] mvtnorm_1.1-2             scales_1.1.1             
[29] BiocParallel_1.24.1       tibble_3.1.2             
[31] docopt_0.7.1              generics_0.1.0           
[33] ellipsis_0.3.2            withr_2.4.2              
[35] magrittr_2.0.1            crayon_1.4.1             
[37] parallelly_1.26.1         future_1.21.0            
[39] fansi_0.5.0               MASS_7.3-54              
[41] beeswarm_0.4.0            poilog_0.4               
[43] tools_4.0.3               GUILDS_1.3               
[45] lifecycle_1.0.0           munsell_0.5.0            
[47] DelayedArray_0.16.3       irlba_2.3.3              
[49] compiler_4.0.3            rsvd_1.0.5               
[51] rlang_0.4.11              grid_4.0.3               
[53] RCurl_1.98-1.3            BiocNeighbors_1.8.2      
[55] bitops_1.0-7              gtable_0.3.0             
[57] codetools_0.2-18          DBI_1.1.1                
[59] R6_2.5.0                  gridExtra_2.3            
[61] dplyr_1.0.7               bdsmatrix_1.3-4          
[63] future.apply_1.7.0        utf8_1.2.1               
[65] stringi_1.6.2             ggbeeswarm_0.6.0         
[67] Rcpp_1.0.7                vctrs_0.3.8              
[69] tidyselect_1.1.1          sparseMatrixStats_1.2.1 
```
