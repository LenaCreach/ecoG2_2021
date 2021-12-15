R Notebook
================

``` r
library("knitr")
library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
```

    ## Loading required package: ggplot2

    ## Loading required package: gridExtra

    ## Loading required package: dada2

    ## Loading required package: Rcpp

    ## Loading required package: phyloseq

    ## Loading required package: DECIPHER

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

    ## Loading required package: parallel

    ## Loading required package: phangorn

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     complement

    ##   ggplot2 gridExtra     dada2  phyloseq  DECIPHER  phangorn 
    ##      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE

``` r
set.seed(100)
miseq_path <- "./MiSeq_SOP"
list.files(miseq_path)
```

    ## character(0)

``` r
set.seed(100)
miseq_path <- "/home/rstudio/MiSeq_SOP"
list.files(miseq_path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

#Sort ensures forward/reverse reads are in same order

``` r
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
fnFs[1:3]
```

    ## [1] "/home/rstudio/MiSeq_SOP/F3D0_S188_L001_R1_001.fastq"  
    ## [2] "/home/rstudio/MiSeq_SOP/F3D1_S189_L001_R1_001.fastq"  
    ## [3] "/home/rstudio/MiSeq_SOP/F3D141_S207_L001_R1_001.fastq"

``` r
plotQualityProfile(fnFs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](cc1_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
plotQualityProfile(fnRs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](cc1_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))
```

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN = 0,maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress=TRUE, multithread = TRUE) #On Windows set multithread=FALSE
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

``` r
derepFs <- derepFastq(filtFs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D0_F_filt.fastq.gz

    ## Encountered 1979 unique sequences from 7113 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D1_F_filt.fastq.gz

    ## Encountered 1639 unique sequences from 5299 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D141_F_filt.fastq.gz

    ## Encountered 1477 unique sequences from 5463 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D142_F_filt.fastq.gz

    ## Encountered 904 unique sequences from 2914 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D143_F_filt.fastq.gz

    ## Encountered 939 unique sequences from 2941 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D144_F_filt.fastq.gz

    ## Encountered 1267 unique sequences from 4312 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D145_F_filt.fastq.gz

    ## Encountered 1756 unique sequences from 6741 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D146_F_filt.fastq.gz

    ## Encountered 1438 unique sequences from 4560 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D147_F_filt.fastq.gz

    ## Encountered 3590 unique sequences from 15637 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D148_F_filt.fastq.gz

    ## Encountered 2762 unique sequences from 11413 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D149_F_filt.fastq.gz

    ## Encountered 3021 unique sequences from 12017 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D150_F_filt.fastq.gz

    ## Encountered 1566 unique sequences from 5032 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D2_F_filt.fastq.gz

    ## Encountered 3707 unique sequences from 18075 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D3_F_filt.fastq.gz

    ## Encountered 1479 unique sequences from 6250 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D5_F_filt.fastq.gz

    ## Encountered 1195 unique sequences from 4052 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D6_F_filt.fastq.gz

    ## Encountered 1832 unique sequences from 7369 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D7_F_filt.fastq.gz

    ## Encountered 1183 unique sequences from 4765 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D8_F_filt.fastq.gz

    ## Encountered 1382 unique sequences from 4871 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D9_F_filt.fastq.gz

    ## Encountered 1709 unique sequences from 6504 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/Mock_F_filt.fastq.gz

    ## Encountered 897 unique sequences from 4314 total sequences read.

``` r
derepRs <- derepFastq(filtRs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D0_R_filt.fastq.gz

    ## Encountered 1660 unique sequences from 7113 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D1_R_filt.fastq.gz

    ## Encountered 1349 unique sequences from 5299 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D141_R_filt.fastq.gz

    ## Encountered 1335 unique sequences from 5463 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D142_R_filt.fastq.gz

    ## Encountered 853 unique sequences from 2914 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D143_R_filt.fastq.gz

    ## Encountered 880 unique sequences from 2941 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D144_R_filt.fastq.gz

    ## Encountered 1286 unique sequences from 4312 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D145_R_filt.fastq.gz

    ## Encountered 1803 unique sequences from 6741 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D146_R_filt.fastq.gz

    ## Encountered 1265 unique sequences from 4560 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D147_R_filt.fastq.gz

    ## Encountered 3414 unique sequences from 15637 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D148_R_filt.fastq.gz

    ## Encountered 2522 unique sequences from 11413 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D149_R_filt.fastq.gz

    ## Encountered 2771 unique sequences from 12017 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D150_R_filt.fastq.gz

    ## Encountered 1415 unique sequences from 5032 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D2_R_filt.fastq.gz

    ## Encountered 3290 unique sequences from 18075 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D3_R_filt.fastq.gz

    ## Encountered 1390 unique sequences from 6250 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D5_R_filt.fastq.gz

    ## Encountered 1134 unique sequences from 4052 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D6_R_filt.fastq.gz

    ## Encountered 1635 unique sequences from 7369 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D7_R_filt.fastq.gz

    ## Encountered 1084 unique sequences from 4765 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D8_R_filt.fastq.gz

    ## Encountered 1161 unique sequences from 4871 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D9_R_filt.fastq.gz

    ## Encountered 1502 unique sequences from 6504 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/Mock_R_filt.fastq.gz

    ## Encountered 732 unique sequences from 4314 total sequences read.

``` r
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames
```

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
plotErrors(errF)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](cc1_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
plotErrors(errR)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](cc1_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

``` r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

``` r
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

``` r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
```

``` r
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  85 186   5   2

``` r
seqtabNoC <- removeBimeraDenovo(seqtabAll)
```

``` bash
cd ~
wget  https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
```

    ## --2021-12-15 08:50:00--  https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 137283333 (131M) [application/octet-stream]
    ## Saving to: ‘silva_nr99_v138.1_train_set.fa.gz.1’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0%  917K 2m26s
    ##     50K .......... .......... .......... .......... ..........  0%  880K 2m29s
    ##    100K .......... .......... .......... .......... ..........  0%  866K 2m31s
    ##    150K .......... .......... .......... .......... ..........  0% 1.58M 2m14s
    ##    200K .......... .......... .......... .......... ..........  0% 1.08M 2m11s
    ##    250K .......... .......... .......... .......... ..........  0% 1.71M 2m2s
    ##    300K .......... .......... .......... .......... ..........  0% 1.74M 1m55s
    ##    350K .......... .......... .......... .......... ..........  0% 1.64M 1m51s
    ##    400K .......... .......... .......... .......... ..........  0% 1.31M 1m50s
    ##    450K .......... .......... .......... .......... ..........  0% 1.78M 1m46s
    ##    500K .......... .......... .......... .......... ..........  0% 1.27M 1m46s
    ##    550K .......... .......... .......... .......... ..........  0%  952K 1m48s
    ##    600K .......... .......... .......... .......... ..........  0%  347K 2m10s
    ##    650K .......... .......... .......... .......... ..........  0% 1.24M 2m8s
    ##    700K .......... .......... .......... .......... ..........  0% 3.03M 2m2s
    ##    750K .......... .......... .......... .......... ..........  0% 2.57M 1m58s
    ##    800K .......... .......... .......... .......... ..........  0% 3.03M 1m53s
    ##    850K .......... .......... .......... .......... ..........  0% 3.43M 1m49s
    ##    900K .......... .......... .......... .......... ..........  0% 1.64M 1m47s
    ##    950K .......... .......... .......... .......... ..........  0% 3.73M 1m44s
    ##   1000K .......... .......... .......... .......... ..........  0% 13.2M 99s
    ##   1050K .......... .......... .......... .......... ..........  0% 2.31M 97s
    ##   1100K .......... .......... .......... .......... ..........  0% 5.18M 94s
    ##   1150K .......... .......... .......... .......... ..........  0% 3.69M 92s
    ##   1200K .......... .......... .......... .......... ..........  0% 3.97M 89s
    ##   1250K .......... .......... .......... .......... ..........  0% 3.25M 87s
    ##   1300K .......... .......... .......... .......... ..........  1% 4.83M 85s
    ##   1350K .......... .......... .......... .......... ..........  1% 1.50M 85s
    ##   1400K .......... .......... .......... .......... ..........  1% 3.15M 83s
    ##   1450K .......... .......... .......... .......... ..........  1% 2.13M 83s
    ##   1500K .......... .......... .......... .......... ..........  1% 2.54M 82s
    ##   1550K .......... .......... .......... .......... ..........  1% 10.3M 79s
    ##   1600K .......... .......... .......... .......... ..........  1% 2.09M 79s
    ##   1650K .......... .......... .......... .......... ..........  1% 4.10M 77s
    ##   1700K .......... .......... .......... .......... ..........  1% 3.14M 76s
    ##   1750K .......... .......... .......... .......... ..........  1% 3.58M 75s
    ##   1800K .......... .......... .......... .......... ..........  1% 2.12M 75s
    ##   1850K .......... .......... .......... .......... ..........  1% 1.30M 75s
    ##   1900K .......... .......... .......... .......... ..........  1% 1.79M 75s
    ##   1950K .......... .......... .......... .......... ..........  1% 2.56M 75s
    ##   2000K .......... .......... .......... .......... ..........  1% 7.99M 73s
    ##   2050K .......... .......... .......... .......... ..........  1% 1.91M 73s
    ##   2100K .......... .......... .......... .......... ..........  1% 3.16M 72s
    ##   2150K .......... .......... .......... .......... ..........  1% 2.99M 72s
    ##   2200K .......... .......... .......... .......... ..........  1% 2.27M 71s
    ##   2250K .......... .......... .......... .......... ..........  1% 2.30M 71s
    ##   2300K .......... .......... .......... .......... ..........  1%  930K 72s
    ##   2350K .......... .......... .......... .......... ..........  1%  488K 76s
    ##   2400K .......... .......... .......... .......... ..........  1% 1.38M 77s
    ##   2450K .......... .......... .......... .......... ..........  1% 2.12M 76s
    ##   2500K .......... .......... .......... .......... ..........  1% 2.58M 76s
    ##   2550K .......... .......... .......... .......... ..........  1% 1.71M 76s
    ##   2600K .......... .......... .......... .......... ..........  1% 1.79M 76s
    ##   2650K .......... .......... .......... .......... ..........  2% 65.1M 74s
    ##   2700K .......... .......... .......... .......... ..........  2% 1.82M 74s
    ##   2750K .......... .......... .......... .......... ..........  2%  823K 76s
    ##   2800K .......... .......... .......... .......... ..........  2% 5.87M 75s
    ##   2850K .......... .......... .......... .......... ..........  2% 3.83M 74s
    ##   2900K .......... .......... .......... .......... ..........  2% 7.07M 73s
    ##   2950K .......... .......... .......... .......... ..........  2% 2.15M 73s
    ##   3000K .......... .......... .......... .......... ..........  2% 1.44M 73s
    ##   3050K .......... .......... .......... .......... ..........  2% 2.48M 73s
    ##   3100K .......... .......... .......... .......... ..........  2% 2.48M 72s
    ##   3150K .......... .......... .......... .......... ..........  2% 2.34M 72s
    ##   3200K .......... .......... .......... .......... ..........  2% 3.85M 71s
    ##   3250K .......... .......... .......... .......... ..........  2% 75.8M 70s
    ##   3300K .......... .......... .......... .......... ..........  2% 3.39M 70s
    ##   3350K .......... .......... .......... .......... ..........  2% 7.62M 69s
    ##   3400K .......... .......... .......... .......... ..........  2% 7.56M 68s
    ##   3450K .......... .......... .......... .......... ..........  2% 4.22M 68s
    ##   3500K .......... .......... .......... .......... ..........  2% 96.6M 67s
    ##   3550K .......... .......... .......... .......... ..........  2% 3.65M 66s
    ##   3600K .......... .......... .......... .......... ..........  2% 9.24M 65s
    ##   3650K .......... .......... .......... .......... ..........  2% 6.20M 65s
    ##   3700K .......... .......... .......... .......... ..........  2% 11.3M 64s
    ##   3750K .......... .......... .......... .......... ..........  2% 2.68M 64s
    ##   3800K .......... .......... .......... .......... ..........  2% 86.9M 63s
    ##   3850K .......... .......... .......... .......... ..........  2% 4.51M 62s
    ##   3900K .......... .......... .......... .......... ..........  2% 2.19M 62s
    ##   3950K .......... .......... .......... .......... ..........  2% 5.27M 62s
    ##   4000K .......... .......... .......... .......... ..........  3% 11.6M 61s
    ##   4050K .......... .......... .......... .......... ..........  3% 89.9M 61s
    ##   4100K .......... .......... .......... .......... ..........  3% 2.39M 60s
    ##   4150K .......... .......... .......... .......... ..........  3% 3.99M 60s
    ##   4200K .......... .......... .......... .......... ..........  3% 2.51M 60s
    ##   4250K .......... .......... .......... .......... ..........  3% 81.0M 59s
    ##   4300K .......... .......... .......... .......... ..........  3% 1.29M 60s
    ##   4350K .......... .......... .......... .......... ..........  3% 2.28M 60s
    ##   4400K .......... .......... .......... .......... ..........  3% 1.65M 60s
    ##   4450K .......... .......... .......... .......... ..........  3% 70.2M 59s
    ##   4500K .......... .......... .......... .......... ..........  3% 1.68M 59s
    ##   4550K .......... .......... .......... .......... ..........  3% 2.71M 59s
    ##   4600K .......... .......... .......... .......... ..........  3% 3.90M 59s
    ##   4650K .......... .......... .......... .......... ..........  3% 2.24M 59s
    ##   4700K .......... .......... .......... .......... ..........  3% 32.4M 58s
    ##   4750K .......... .......... .......... .......... ..........  3% 1.59M 58s
    ##   4800K .......... .......... .......... .......... ..........  3% 1.66M 58s
    ##   4850K .......... .......... .......... .......... ..........  3% 2.59M 58s
    ##   4900K .......... .......... .......... .......... ..........  3% 82.0M 58s
    ##   4950K .......... .......... .......... .......... ..........  3% 4.59M 57s
    ##   5000K .......... .......... .......... .......... ..........  3% 87.8M 57s
    ##   5050K .......... .......... .......... .......... ..........  3% 5.38M 56s
    ##   5100K .......... .......... .......... .......... ..........  3% 5.39M 56s
    ##   5150K .......... .......... .......... .......... ..........  3% 3.74M 56s
    ##   5200K .......... .......... .......... .......... ..........  3% 85.8M 55s
    ##   5250K .......... .......... .......... .......... ..........  3% 3.25M 55s
    ##   5300K .......... .......... .......... .......... ..........  3% 90.6M 55s
    ##   5350K .......... .......... .......... .......... ..........  4% 1.94M 55s
    ##   5400K .......... .......... .......... .......... ..........  4% 2.55M 55s
    ##   5450K .......... .......... .......... .......... ..........  4% 76.6M 54s
    ##   5500K .......... .......... .......... .......... ..........  4% 4.77M 54s
    ##   5550K .......... .......... .......... .......... ..........  4% 11.8M 53s
    ##   5600K .......... .......... .......... .......... ..........  4% 63.8M 53s
    ##   5650K .......... .......... .......... .......... ..........  4% 3.92M 53s
    ##   5700K .......... .......... .......... .......... ..........  4% 4.14M 53s
    ##   5750K .......... .......... .......... .......... ..........  4% 10.7M 52s
    ##   5800K .......... .......... .......... .......... ..........  4% 12.5M 52s
    ##   5850K .......... .......... .......... .......... ..........  4% 3.28M 52s
    ##   5900K .......... .......... .......... .......... ..........  4% 4.74M 51s
    ##   5950K .......... .......... .......... .......... ..........  4% 59.4M 51s
    ##   6000K .......... .......... .......... .......... ..........  4% 8.74M 51s
    ##   6050K .......... .......... .......... .......... ..........  4% 6.04M 50s
    ##   6100K .......... .......... .......... .......... ..........  4% 75.5M 50s
    ##   6150K .......... .......... .......... .......... ..........  4% 8.15M 50s
    ##   6200K .......... .......... .......... .......... ..........  4% 81.3M 49s
    ##   6250K .......... .......... .......... .......... ..........  4% 2.76M 49s
    ##   6300K .......... .......... .......... .......... ..........  4% 4.12M 49s
    ##   6350K .......... .......... .......... .......... ..........  4% 3.50M 49s
    ##   6400K .......... .......... .......... .......... ..........  4% 2.77M 49s
    ##   6450K .......... .......... .......... .......... ..........  4% 93.3M 49s
    ##   6500K .......... .......... .......... .......... ..........  4% 5.01M 48s
    ##   6550K .......... .......... .......... .......... ..........  4% 81.5M 48s
    ##   6600K .......... .......... .......... .......... ..........  4% 7.02M 48s
    ##   6650K .......... .......... .......... .......... ..........  4% 76.3M 47s
    ##   6700K .......... .......... .......... .......... ..........  5% 4.71M 47s
    ##   6750K .......... .......... .......... .......... ..........  5% 1.60M 47s
    ##   6800K .......... .......... .......... .......... ..........  5% 78.6M 47s
    ##   6850K .......... .......... .......... .......... ..........  5%  801K 48s
    ##   6900K .......... .......... .......... .......... ..........  5% 79.8M 47s
    ##   6950K .......... .......... .......... .......... ..........  5% 1.16M 48s
    ##   7000K .......... .......... .......... .......... ..........  5% 6.92M 48s
    ##   7050K .......... .......... .......... .......... ..........  5% 1.22M 48s
    ##   7100K .......... .......... .......... .......... ..........  5% 7.45M 48s
    ##   7150K .......... .......... .......... .......... ..........  5% 1.63M 48s
    ##   7200K .......... .......... .......... .......... ..........  5% 74.4M 48s
    ##   7250K .......... .......... .......... .......... ..........  5% 2.64M 48s
    ##   7300K .......... .......... .......... .......... ..........  5% 2.18M 48s
    ##   7350K .......... .......... .......... .......... ..........  5% 9.43M 47s
    ##   7400K .......... .......... .......... .......... ..........  5% 81.0M 47s
    ##   7450K .......... .......... .......... .......... ..........  5% 1.94M 47s
    ##   7500K .......... .......... .......... .......... ..........  5% 75.2M 47s
    ##   7550K .......... .......... .......... .......... ..........  5% 2.54M 47s
    ##   7600K .......... .......... .......... .......... ..........  5% 2.25M 47s
    ##   7650K .......... .......... .......... .......... ..........  5% 69.4M 47s
    ##   7700K .......... .......... .......... .......... ..........  5% 2.55M 47s
    ##   7750K .......... .......... .......... .......... ..........  5% 66.5M 46s
    ##   7800K .......... .......... .......... .......... ..........  5% 3.64M 46s
    ##   7850K .......... .......... .......... .......... ..........  5% 12.1M 46s
    ##   7900K .......... .......... .......... .......... ..........  5% 3.50M 46s
    ##   7950K .......... .......... .......... .......... ..........  5% 5.30M 46s
    ##   8000K .......... .......... .......... .......... ..........  6% 1.47M 46s
    ##   8050K .......... .......... .......... .......... ..........  6% 65.3M 46s
    ##   8100K .......... .......... .......... .......... ..........  6% 1.50M 46s
    ##   8150K .......... .......... .......... .......... ..........  6% 49.1M 45s
    ##   8200K .......... .......... .......... .......... ..........  6% 4.50M 45s
    ##   8250K .......... .......... .......... .......... ..........  6% 2.39M 45s
    ##   8300K .......... .......... .......... .......... ..........  6% 75.1M 45s
    ##   8350K .......... .......... .......... .......... ..........  6% 1.40M 45s
    ##   8400K .......... .......... .......... .......... ..........  6% 1.14M 46s
    ##   8450K .......... .......... .......... .......... ..........  6% 52.8M 45s
    ##   8500K .......... .......... .......... .......... ..........  6% 69.7M 45s
    ##   8550K .......... .......... .......... .......... ..........  6%  998K 46s
    ##   8600K .......... .......... .......... .......... ..........  6% 63.4M 45s
    ##   8650K .......... .......... .......... .......... ..........  6% 1.76M 45s
    ##   8700K .......... .......... .......... .......... ..........  6% 73.3M 45s
    ##   8750K .......... .......... .......... .......... ..........  6% 2.14M 45s
    ##   8800K .......... .......... .......... .......... ..........  6% 73.2M 45s
    ##   8850K .......... .......... .......... .......... ..........  6% 2.83M 45s
    ##   8900K .......... .......... .......... .......... ..........  6% 68.8M 45s
    ##   8950K .......... .......... .......... .......... ..........  6% 3.03M 45s
    ##   9000K .......... .......... .......... .......... ..........  6% 58.0M 44s
    ##   9050K .......... .......... .......... .......... ..........  6% 76.8M 44s
    ##   9100K .......... .......... .......... .......... ..........  6% 1.12M 44s
    ##   9150K .......... .......... .......... .......... ..........  6% 60.1M 44s
    ##   9200K .......... .......... .......... .......... ..........  6% 1.50M 44s
    ##   9250K .......... .......... .......... .......... ..........  6% 61.7M 44s
    ##   9300K .......... .......... .......... .......... ..........  6% 3.96M 44s
    ##   9350K .......... .......... .......... .......... ..........  7% 69.1M 44s
    ##   9400K .......... .......... .......... .......... ..........  7% 4.05M 44s
    ##   9450K .......... .......... .......... .......... ..........  7% 75.6M 43s
    ##   9500K .......... .......... .......... .......... ..........  7% 2.64M 43s
    ##   9550K .......... .......... .......... .......... ..........  7% 67.5M 43s
    ##   9600K .......... .......... .......... .......... ..........  7% 1.97M 43s
    ##   9650K .......... .......... .......... .......... ..........  7% 71.5M 43s
    ##   9700K .......... .......... .......... .......... ..........  7% 79.0M 43s
    ##   9750K .......... .......... .......... .......... ..........  7% 1.50M 43s
    ##   9800K .......... .......... .......... .......... ..........  7% 91.7M 43s
    ##   9850K .......... .......... .......... .......... ..........  7% 2.03M 43s
    ##   9900K .......... .......... .......... .......... ..........  7% 80.9M 43s
    ##   9950K .......... .......... .......... .......... ..........  7% 1.63M 43s
    ##  10000K .......... .......... .......... .......... ..........  7% 63.7M 43s
    ##  10050K .......... .......... .......... .......... ..........  7% 1.52M 43s
    ##  10100K .......... .......... .......... .......... ..........  7% 71.0M 43s
    ##  10150K .......... .......... .......... .......... ..........  7% 3.97M 42s
    ##  10200K .......... .......... .......... .......... ..........  7% 2.43M 42s
    ##  10250K .......... .......... .......... .......... ..........  7% 77.1M 42s
    ##  10300K .......... .......... .......... .......... ..........  7%  936K 43s
    ##  10350K .......... .......... .......... .......... ..........  7% 50.1M 42s
    ##  10400K .......... .......... .......... .......... ..........  7% 3.55M 42s
    ##  10450K .......... .......... .......... .......... ..........  7% 9.23M 42s
    ##  10500K .......... .......... .......... .......... ..........  7% 10.4M 42s
    ##  10550K .......... .......... .......... .......... ..........  7% 14.4M 42s
    ##  10600K .......... .......... .......... .......... ..........  7% 38.1M 42s
    ##  10650K .......... .......... .......... .......... ..........  7% 17.1M 42s
    ##  10700K .......... .......... .......... .......... ..........  8% 23.0M 41s
    ##  10750K .......... .......... .......... .......... ..........  8% 21.1M 41s
    ##  10800K .......... .......... .......... .......... ..........  8% 25.7M 41s
    ##  10850K .......... .......... .......... .......... ..........  8% 3.46M 41s
    ##  10900K .......... .......... .......... .......... ..........  8% 28.3M 41s
    ##  10950K .......... .......... .......... .......... ..........  8% 68.3M 41s
    ##  11000K .......... .......... .......... .......... ..........  8% 21.4M 40s
    ##  11050K .......... .......... .......... .......... ..........  8% 31.5M 40s
    ##  11100K .......... .......... .......... .......... ..........  8% 5.39M 40s
    ##  11150K .......... .......... .......... .......... ..........  8% 36.4M 40s
    ##  11200K .......... .......... .......... .......... ..........  8% 2.89M 40s
    ##  11250K .......... .......... .......... .......... ..........  8% 19.2M 40s
    ##  11300K .......... .......... .......... .......... ..........  8% 4.05M 40s
    ##  11350K .......... .......... .......... .......... ..........  8% 60.0M 40s
    ##  11400K .......... .......... .......... .......... ..........  8% 84.8M 39s
    ##  11450K .......... .......... .......... .......... ..........  8% 1.84M 39s
    ##  11500K .......... .......... .......... .......... ..........  8% 84.9M 39s
    ##  11550K .......... .......... .......... .......... ..........  8% 2.71M 39s
    ##  11600K .......... .......... .......... .......... ..........  8% 87.2M 39s
    ##  11650K .......... .......... .......... .......... ..........  8% 2.15M 39s
    ##  11700K .......... .......... .......... .......... ..........  8% 66.3M 39s
    ##  11750K .......... .......... .......... .......... ..........  8% 78.3M 39s
    ##  11800K .......... .......... .......... .......... ..........  8% 1.62M 39s
    ##  11850K .......... .......... .......... .......... ..........  8% 73.0M 39s
    ##  11900K .......... .......... .......... .......... ..........  8% 80.0M 39s
    ##  11950K .......... .......... .......... .......... ..........  8% 1.21M 39s
    ##  12000K .......... .......... .......... .......... ..........  8% 77.3M 39s
    ##  12050K .......... .......... .......... .......... ..........  9% 1.53M 39s
    ##  12100K .......... .......... .......... .......... ..........  9% 80.0M 39s
    ##  12150K .......... .......... .......... .......... ..........  9% 79.1M 38s
    ##  12200K .......... .......... .......... .......... ..........  9% 1.53M 39s
    ##  12250K .......... .......... .......... .......... ..........  9% 83.1M 38s
    ##  12300K .......... .......... .......... .......... ..........  9% 2.75M 38s
    ##  12350K .......... .......... .......... .......... ..........  9% 32.2M 38s
    ##  12400K .......... .......... .......... .......... ..........  9% 78.6M 38s
    ##  12450K .......... .......... .......... .......... ..........  9% 14.3M 38s
    ##  12500K .......... .......... .......... .......... ..........  9% 74.4M 38s
    ##  12550K .......... .......... .......... .......... ..........  9% 45.5M 38s
    ##  12600K .......... .......... .......... .......... ..........  9% 17.2M 38s
    ##  12650K .......... .......... .......... .......... ..........  9% 76.1M 37s
    ##  12700K .......... .......... .......... .......... ..........  9% 28.8M 37s
    ##  12750K .......... .......... .......... .......... ..........  9% 18.6M 37s
    ##  12800K .......... .......... .......... .......... ..........  9% 44.5M 37s
    ##  12850K .......... .......... .......... .......... ..........  9% 20.4M 37s
    ##  12900K .......... .......... .......... .......... ..........  9% 5.69M 37s
    ##  12950K .......... .......... .......... .......... ..........  9% 4.13M 37s
    ##  13000K .......... .......... .......... .......... ..........  9%  103M 37s
    ##  13050K .......... .......... .......... .......... ..........  9% 11.0M 36s
    ##  13100K .......... .......... .......... .......... ..........  9% 3.66M 36s
    ##  13150K .......... .......... .......... .......... ..........  9% 36.7M 36s
    ##  13200K .......... .......... .......... .......... ..........  9% 10.5M 36s
    ##  13250K .......... .......... .......... .......... ..........  9% 65.7M 36s
    ##  13300K .......... .......... .......... .......... ..........  9% 96.4M 36s
    ##  13350K .......... .......... .......... .......... ..........  9% 10.1M 36s
    ##  13400K .......... .......... .......... .......... .......... 10% 87.0M 36s
    ##  13450K .......... .......... .......... .......... .......... 10%  113M 35s
    ##  13500K .......... .......... .......... .......... .......... 10% 4.30M 35s
    ##  13550K .......... .......... .......... .......... .......... 10% 79.5M 35s
    ##  13600K .......... .......... .......... .......... .......... 10% 13.5M 35s
    ##  13650K .......... .......... .......... .......... .......... 10% 91.2M 35s
    ##  13700K .......... .......... .......... .......... .......... 10% 26.4M 35s
    ##  13750K .......... .......... .......... .......... .......... 10% 28.2M 35s
    ##  13800K .......... .......... .......... .......... .......... 10%  100M 35s
    ##  13850K .......... .......... .......... .......... .......... 10% 16.5M 35s
    ##  13900K .......... .......... .......... .......... .......... 10% 72.5M 34s
    ##  13950K .......... .......... .......... .......... .......... 10% 87.4M 34s
    ##  14000K .......... .......... .......... .......... .......... 10% 19.7M 34s
    ##  14050K .......... .......... .......... .......... .......... 10% 29.0M 34s
    ##  14100K .......... .......... .......... .......... .......... 10% 99.4M 34s
    ##  14150K .......... .......... .......... .......... .......... 10% 30.6M 34s
    ##  14200K .......... .......... .......... .......... .......... 10% 16.0M 34s
    ##  14250K .......... .......... .......... .......... .......... 10% 44.9M 34s
    ##  14300K .......... .......... .......... .......... .......... 10% 13.8M 33s
    ##  14350K .......... .......... .......... .......... .......... 10% 50.9M 33s
    ##  14400K .......... .......... .......... .......... .......... 10%  103M 33s
    ##  14450K .......... .......... .......... .......... .......... 10% 18.3M 33s
    ##  14500K .......... .......... .......... .......... .......... 10% 26.8M 33s
    ##  14550K .......... .......... .......... .......... .......... 10% 66.4M 33s
    ##  14600K .......... .......... .......... .......... .......... 10% 7.33M 33s
    ##  14650K .......... .......... .......... .......... .......... 10% 21.0M 33s
    ##  14700K .......... .......... .......... .......... .......... 11% 26.5M 33s
    ##  14750K .......... .......... .......... .......... .......... 11% 21.4M 32s
    ##  14800K .......... .......... .......... .......... .......... 11% 26.8M 32s
    ##  14850K .......... .......... .......... .......... .......... 11% 81.5M 32s
    ##  14900K .......... .......... .......... .......... .......... 11% 5.29M 32s
    ##  14950K .......... .......... .......... .......... .......... 11% 84.7M 32s
    ##  15000K .......... .......... .......... .......... .......... 11% 8.91M 32s
    ##  15050K .......... .......... .......... .......... .......... 11% 80.3M 32s
    ##  15100K .......... .......... .......... .......... .......... 11%  102M 32s
    ##  15150K .......... .......... .......... .......... .......... 11% 13.5M 32s
    ##  15200K .......... .......... .......... .......... .......... 11% 97.6M 32s
    ##  15250K .......... .......... .......... .......... .......... 11% 12.7M 31s
    ##  15300K .......... .......... .......... .......... .......... 11% 83.8M 31s
    ##  15350K .......... .......... .......... .......... .......... 11% 83.8M 31s
    ##  15400K .......... .......... .......... .......... .......... 11% 7.98M 31s
    ##  15450K .......... .......... .......... .......... .......... 11% 87.5M 31s
    ##  15500K .......... .......... .......... .......... .......... 11% 5.06M 31s
    ##  15550K .......... .......... .......... .......... .......... 11% 27.0M 31s
    ##  15600K .......... .......... .......... .......... .......... 11% 91.8M 31s
    ##  15650K .......... .......... .......... .......... .......... 11% 8.36M 31s
    ##  15700K .......... .......... .......... .......... .......... 11% 80.8M 31s
    ##  15750K .......... .......... .......... .......... .......... 11% 82.6M 31s
    ##  15800K .......... .......... .......... .......... .......... 11% 11.0M 30s
    ##  15850K .......... .......... .......... .......... .......... 11% 93.6M 30s
    ##  15900K .......... .......... .......... .......... .......... 11% 11.9M 30s
    ##  15950K .......... .......... .......... .......... .......... 11% 79.9M 30s
    ##  16000K .......... .......... .......... .......... .......... 11% 7.50M 30s
    ##  16050K .......... .......... .......... .......... .......... 12% 67.7M 30s
    ##  16100K .......... .......... .......... .......... .......... 12% 91.5M 30s
    ##  16150K .......... .......... .......... .......... .......... 12% 18.4M 30s
    ##  16200K .......... .......... .......... .......... .......... 12% 74.1M 30s
    ##  16250K .......... .......... .......... .......... .......... 12%  108M 30s
    ##  16300K .......... .......... .......... .......... .......... 12% 27.9M 30s
    ##  16350K .......... .......... .......... .......... .......... 12% 35.7M 29s
    ##  16400K .......... .......... .......... .......... .......... 12% 25.1M 29s
    ##  16450K .......... .......... .......... .......... .......... 12% 41.9M 29s
    ##  16500K .......... .......... .......... .......... .......... 12% 83.8M 29s
    ##  16550K .......... .......... .......... .......... .......... 12% 21.7M 29s
    ##  16600K .......... .......... .......... .......... .......... 12% 38.9M 29s
    ##  16650K .......... .......... .......... .......... .......... 12% 20.5M 29s
    ##  16700K .......... .......... .......... .......... .......... 12% 79.1M 29s
    ##  16750K .......... .......... .......... .......... .......... 12% 35.4M 29s
    ##  16800K .......... .......... .......... .......... .......... 12% 11.3M 29s
    ##  16850K .......... .......... .......... .......... .......... 12% 84.2M 29s
    ##  16900K .......... .......... .......... .......... .......... 12% 17.9M 28s
    ##  16950K .......... .......... .......... .......... .......... 12% 70.1M 28s
    ##  17000K .......... .......... .......... .......... .......... 12% 75.8M 28s
    ##  17050K .......... .......... .......... .......... .......... 12% 13.6M 28s
    ##  17100K .......... .......... .......... .......... .......... 12% 86.9M 28s
    ##  17150K .......... .......... .......... .......... .......... 12% 15.8M 28s
    ##  17200K .......... .......... .......... .......... .......... 12% 77.3M 28s
    ##  17250K .......... .......... .......... .......... .......... 12% 91.9M 28s
    ##  17300K .......... .......... .......... .......... .......... 12% 16.6M 28s
    ##  17350K .......... .......... .......... .......... .......... 12% 3.73M 28s
    ##  17400K .......... .......... .......... .......... .......... 13% 81.4M 28s
    ##  17450K .......... .......... .......... .......... .......... 13% 84.3M 28s
    ##  17500K .......... .......... .......... .......... .......... 13% 27.7M 28s
    ##  17550K .......... .......... .......... .......... .......... 13% 17.2M 27s
    ##  17600K .......... .......... .......... .......... .......... 13% 75.0M 27s
    ##  17650K .......... .......... .......... .......... .......... 13% 15.0M 27s
    ##  17700K .......... .......... .......... .......... .......... 13% 76.1M 27s
    ##  17750K .......... .......... .......... .......... .......... 13% 76.8M 27s
    ##  17800K .......... .......... .......... .......... .......... 13% 15.2M 27s
    ##  17850K .......... .......... .......... .......... .......... 13% 81.9M 27s
    ##  17900K .......... .......... .......... .......... .......... 13% 88.6M 27s
    ##  17950K .......... .......... .......... .......... .......... 13% 18.0M 27s
    ##  18000K .......... .......... .......... .......... .......... 13% 67.6M 27s
    ##  18050K .......... .......... .......... .......... .......... 13% 88.4M 27s
    ##  18100K .......... .......... .......... .......... .......... 13% 23.4M 27s
    ##  18150K .......... .......... .......... .......... .......... 13% 55.8M 27s
    ##  18200K .......... .......... .......... .......... .......... 13% 87.3M 26s
    ##  18250K .......... .......... .......... .......... .......... 13% 15.1M 26s
    ##  18300K .......... .......... .......... .......... .......... 13% 84.0M 26s
    ##  18350K .......... .......... .......... .......... .......... 13% 15.0M 26s
    ##  18400K .......... .......... .......... .......... .......... 13% 51.1M 26s
    ##  18450K .......... .......... .......... .......... .......... 13% 75.7M 26s
    ##  18500K .......... .......... .......... .......... .......... 13% 17.0M 26s
    ##  18550K .......... .......... .......... .......... .......... 13% 64.7M 26s
    ##  18600K .......... .......... .......... .......... .......... 13% 76.6M 26s
    ##  18650K .......... .......... .......... .......... .......... 13% 16.7M 26s
    ##  18700K .......... .......... .......... .......... .......... 13% 74.3M 26s
    ##  18750K .......... .......... .......... .......... .......... 14% 6.18M 26s
    ##  18800K .......... .......... .......... .......... .......... 14% 74.4M 26s
    ##  18850K .......... .......... .......... .......... .......... 14% 83.5M 26s
    ##  18900K .......... .......... .......... .......... .......... 14% 14.8M 26s
    ##  18950K .......... .......... .......... .......... .......... 14% 70.1M 25s
    ##  19000K .......... .......... .......... .......... .......... 14% 17.0M 25s
    ##  19050K .......... .......... .......... .......... .......... 14% 58.9M 25s
    ##  19100K .......... .......... .......... .......... .......... 14% 77.0M 25s
    ##  19150K .......... .......... .......... .......... .......... 14% 20.1M 25s
    ##  19200K .......... .......... .......... .......... .......... 14% 5.28M 25s
    ##  19250K .......... .......... .......... .......... .......... 14% 74.2M 25s
    ##  19300K .......... .......... .......... .......... .......... 14% 83.7M 25s
    ##  19350K .......... .......... .......... .......... .......... 14% 4.26M 25s
    ##  19400K .......... .......... .......... .......... .......... 14% 64.7M 25s
    ##  19450K .......... .......... .......... .......... .......... 14% 82.3M 25s
    ##  19500K .......... .......... .......... .......... .......... 14% 23.9M 25s
    ##  19550K .......... .......... .......... .......... .......... 14% 47.9M 25s
    ##  19600K .......... .......... .......... .......... .......... 14% 84.0M 25s
    ##  19650K .......... .......... .......... .......... .......... 14% 23.0M 25s
    ##  19700K .......... .......... .......... .......... .......... 14% 67.7M 25s
    ##  19750K .......... .......... .......... .......... .......... 14% 19.7M 24s
    ##  19800K .......... .......... .......... .......... .......... 14% 73.5M 24s
    ##  19850K .......... .......... .......... .......... .......... 14% 75.5M 24s
    ##  19900K .......... .......... .......... .......... .......... 14% 11.3M 24s
    ##  19950K .......... .......... .......... .......... .......... 14% 66.4M 24s
    ##  20000K .......... .......... .......... .......... .......... 14% 15.6M 24s
    ##  20050K .......... .......... .......... .......... .......... 14% 65.1M 24s
    ##  20100K .......... .......... .......... .......... .......... 15% 81.4M 24s
    ##  20150K .......... .......... .......... .......... .......... 15% 15.8M 24s
    ##  20200K .......... .......... .......... .......... .......... 15% 82.3M 24s
    ##  20250K .......... .......... .......... .......... .......... 15% 86.3M 24s
    ##  20300K .......... .......... .......... .......... .......... 15% 16.1M 24s
    ##  20350K .......... .......... .......... .......... .......... 15% 69.5M 24s
    ##  20400K .......... .......... .......... .......... .......... 15% 16.4M 24s
    ##  20450K .......... .......... .......... .......... .......... 15% 79.0M 24s
    ##  20500K .......... .......... .......... .......... .......... 15% 83.5M 24s
    ##  20550K .......... .......... .......... .......... .......... 15% 16.3M 24s
    ##  20600K .......... .......... .......... .......... .......... 15% 78.9M 23s
    ##  20650K .......... .......... .......... .......... .......... 15% 17.3M 23s
    ##  20700K .......... .......... .......... .......... .......... 15% 63.9M 23s
    ##  20750K .......... .......... .......... .......... .......... 15% 64.1M 23s
    ##  20800K .......... .......... .......... .......... .......... 15% 36.2M 23s
    ##  20850K .......... .......... .......... .......... .......... 15% 24.7M 23s
    ##  20900K .......... .......... .......... .......... .......... 15% 79.0M 23s
    ##  20950K .......... .......... .......... .......... .......... 15% 68.6M 23s
    ##  21000K .......... .......... .......... .......... .......... 15% 23.3M 23s
    ##  21050K .......... .......... .......... .......... .......... 15% 71.6M 23s
    ##  21100K .......... .......... .......... .......... .......... 15% 76.4M 23s
    ##  21150K .......... .......... .......... .......... .......... 15% 44.9M 23s
    ##  21200K .......... .......... .......... .......... .......... 15% 26.4M 23s
    ##  21250K .......... .......... .......... .......... .......... 15% 53.0M 23s
    ##  21300K .......... .......... .......... .......... .......... 15% 37.6M 23s
    ##  21350K .......... .......... .......... .......... .......... 15% 24.3M 23s
    ##  21400K .......... .......... .......... .......... .......... 15% 29.9M 23s
    ##  21450K .......... .......... .......... .......... .......... 16% 89.3M 22s
    ##  21500K .......... .......... .......... .......... .......... 16% 29.7M 22s
    ##  21550K .......... .......... .......... .......... .......... 16% 31.0M 22s
    ##  21600K .......... .......... .......... .......... .......... 16% 31.6M 22s
    ##  21650K .......... .......... .......... .......... .......... 16% 28.6M 22s
    ##  21700K .......... .......... .......... .......... .......... 16% 81.6M 22s
    ##  21750K .......... .......... .......... .......... .......... 16% 27.4M 22s
    ##  21800K .......... .......... .......... .......... .......... 16% 33.5M 22s
    ##  21850K .......... .......... .......... .......... .......... 16% 30.5M 22s
    ##  21900K .......... .......... .......... .......... .......... 16% 76.9M 22s
    ##  21950K .......... .......... .......... .......... .......... 16% 25.3M 22s
    ##  22000K .......... .......... .......... .......... .......... 16% 28.9M 22s
    ##  22050K .......... .......... .......... .......... .......... 16% 34.0M 22s
    ##  22100K .......... .......... .......... .......... .......... 16% 26.3M 22s
    ##  22150K .......... .......... .......... .......... .......... 16% 67.5M 22s
    ##  22200K .......... .......... .......... .......... .......... 16% 35.7M 22s
    ##  22250K .......... .......... .......... .......... .......... 16% 32.8M 22s
    ##  22300K .......... .......... .......... .......... .......... 16% 32.7M 22s
    ##  22350K .......... .......... .......... .......... .......... 16% 41.0M 22s
    ##  22400K .......... .......... .......... .......... .......... 16% 71.8M 21s
    ##  22450K .......... .......... .......... .......... .......... 16% 24.7M 21s
    ##  22500K .......... .......... .......... .......... .......... 16% 42.6M 21s
    ##  22550K .......... .......... .......... .......... .......... 16% 24.4M 21s
    ##  22600K .......... .......... .......... .......... .......... 16% 86.4M 21s
    ##  22650K .......... .......... .......... .......... .......... 16% 35.0M 21s
    ##  22700K .......... .......... .......... .......... .......... 16% 27.7M 21s
    ##  22750K .......... .......... .......... .......... .......... 17% 32.4M 21s
    ##  22800K .......... .......... .......... .......... .......... 17% 28.7M 21s
    ##  22850K .......... .......... .......... .......... .......... 17% 73.1M 21s
    ##  22900K .......... .......... .......... .......... .......... 17% 39.5M 21s
    ##  22950K .......... .......... .......... .......... .......... 17% 23.7M 21s
    ##  23000K .......... .......... .......... .......... .......... 17% 38.5M 21s
    ##  23050K .......... .......... .......... .......... .......... 17% 89.0M 21s
    ##  23100K .......... .......... .......... .......... .......... 17% 24.2M 21s
    ##  23150K .......... .......... .......... .......... .......... 17% 60.1M 21s
    ##  23200K .......... .......... .......... .......... .......... 17% 16.5M 21s
    ##  23250K .......... .......... .......... .......... .......... 17% 71.6M 21s
    ##  23300K .......... .......... .......... .......... .......... 17% 92.9M 21s
    ##  23350K .......... .......... .......... .......... .......... 17% 12.2M 21s
    ##  23400K .......... .......... .......... .......... .......... 17% 84.6M 20s
    ##  23450K .......... .......... .......... .......... .......... 17% 14.4M 20s
    ##  23500K .......... .......... .......... .......... .......... 17% 76.7M 20s
    ##  23550K .......... .......... .......... .......... .......... 17% 12.0M 20s
    ##  23600K .......... .......... .......... .......... .......... 17% 66.9M 20s
    ##  23650K .......... .......... .......... .......... .......... 17% 84.3M 20s
    ##  23700K .......... .......... .......... .......... .......... 17% 84.9M 20s
    ##  23750K .......... .......... .......... .......... .......... 17% 25.4M 20s
    ##  23800K .......... .......... .......... .......... .......... 17% 74.2M 20s
    ##  23850K .......... .......... .......... .......... .......... 17% 64.7M 20s
    ##  23900K .......... .......... .......... .......... .......... 17% 24.5M 20s
    ##  23950K .......... .......... .......... .......... .......... 17% 37.0M 20s
    ##  24000K .......... .......... .......... .......... .......... 17% 23.8M 20s
    ##  24050K .......... .......... .......... .......... .......... 17% 84.7M 20s
    ##  24100K .......... .......... .......... .......... .......... 18% 42.2M 20s
    ##  24150K .......... .......... .......... .......... .......... 18% 20.8M 20s
    ##  24200K .......... .......... .......... .......... .......... 18% 34.9M 20s
    ##  24250K .......... .......... .......... .......... .......... 18% 81.1M 20s
    ##  24300K .......... .......... .......... .......... .......... 18% 27.7M 20s
    ##  24350K .......... .......... .......... .......... .......... 18% 28.9M 20s
    ##  24400K .......... .......... .......... .......... .......... 18% 26.8M 20s
    ##  24450K .......... .......... .......... .......... .......... 18% 37.0M 20s
    ##  24500K .......... .......... .......... .......... .......... 18% 77.4M 20s
    ##  24550K .......... .......... .......... .......... .......... 18% 22.3M 19s
    ##  24600K .......... .......... .......... .......... .......... 18% 39.4M 19s
    ##  24650K .......... .......... .......... .......... .......... 18% 23.2M 19s
    ##  24700K .......... .......... .......... .......... .......... 18% 89.4M 19s
    ##  24750K .......... .......... .......... .......... .......... 18% 34.3M 19s
    ##  24800K .......... .......... .......... .......... .......... 18% 19.0M 19s
    ##  24850K .......... .......... .......... .......... .......... 18% 49.0M 19s
    ##  24900K .......... .......... .......... .......... .......... 18% 22.6M 19s
    ##  24950K .......... .......... .......... .......... .......... 18% 61.1M 19s
    ##  25000K .......... .......... .......... .......... .......... 18% 36.2M 19s
    ##  25050K .......... .......... .......... .......... .......... 18% 24.4M 19s
    ##  25100K .......... .......... .......... .......... .......... 18% 33.3M 19s
    ##  25150K .......... .......... .......... .......... .......... 18% 23.6M 19s
    ##  25200K .......... .......... .......... .......... .......... 18% 81.4M 19s
    ##  25250K .......... .......... .......... .......... .......... 18% 35.0M 19s
    ##  25300K .......... .......... .......... .......... .......... 18% 23.3M 19s
    ##  25350K .......... .......... .......... .......... .......... 18% 26.5M 19s
    ##  25400K .......... .......... .......... .......... .......... 18% 94.9M 19s
    ##  25450K .......... .......... .......... .......... .......... 19% 30.3M 19s
    ##  25500K .......... .......... .......... .......... .......... 19% 29.7M 19s
    ##  25550K .......... .......... .......... .......... .......... 19% 29.1M 19s
    ##  25600K .......... .......... .......... .......... .......... 19% 23.7M 19s
    ##  25650K .......... .......... .......... .......... .......... 19% 62.1M 19s
    ##  25700K .......... .......... .......... .......... .......... 19% 89.2M 19s
    ##  25750K .......... .......... .......... .......... .......... 19% 15.6M 19s
    ##  25800K .......... .......... .......... .......... .......... 19% 95.9M 18s
    ##  25850K .......... .......... .......... .......... .......... 19% 12.4M 18s
    ##  25900K .......... .......... .......... .......... .......... 19% 81.4M 18s
    ##  25950K .......... .......... .......... .......... .......... 19% 80.3M 18s
    ##  26000K .......... .......... .......... .......... .......... 19% 95.5M 18s
    ##  26050K .......... .......... .......... .......... .......... 19% 7.64M 18s
    ##  26100K .......... .......... .......... .......... .......... 19% 71.5M 18s
    ##  26150K .......... .......... .......... .......... .......... 19% 77.8M 18s
    ##  26200K .......... .......... .......... .......... .......... 19% 87.8M 18s
    ##  26250K .......... .......... .......... .......... .......... 19% 22.8M 18s
    ##  26300K .......... .......... .......... .......... .......... 19% 78.4M 18s
    ##  26350K .......... .......... .......... .......... .......... 19% 61.2M 18s
    ##  26400K .......... .......... .......... .......... .......... 19% 90.6M 18s
    ##  26450K .......... .......... .......... .......... .......... 19% 29.3M 18s
    ##  26500K .......... .......... .......... .......... .......... 19% 22.3M 18s
    ##  26550K .......... .......... .......... .......... .......... 19% 78.4M 18s
    ##  26600K .......... .......... .......... .......... .......... 19% 99.2M 18s
    ##  26650K .......... .......... .......... .......... .......... 19% 81.7M 18s
    ##  26700K .......... .......... .......... .......... .......... 19% 20.9M 18s
    ##  26750K .......... .......... .......... .......... .......... 19% 65.5M 18s
    ##  26800K .......... .......... .......... .......... .......... 20% 94.7M 18s
    ##  26850K .......... .......... .......... .......... .......... 20% 91.4M 18s
    ##  26900K .......... .......... .......... .......... .......... 20% 45.9M 18s
    ##  26950K .......... .......... .......... .......... .......... 20% 74.0M 18s
    ##  27000K .......... .......... .......... .......... .......... 20% 40.3M 18s
    ##  27050K .......... .......... .......... .......... .......... 20% 71.5M 18s
    ##  27100K .......... .......... .......... .......... .......... 20% 24.5M 18s
    ##  27150K .......... .......... .......... .......... .......... 20% 73.4M 17s
    ##  27200K .......... .......... .......... .......... .......... 20% 31.4M 17s
    ##  27250K .......... .......... .......... .......... .......... 20% 55.0M 17s
    ##  27300K .......... .......... .......... .......... .......... 20% 83.4M 17s
    ##  27350K .......... .......... .......... .......... .......... 20% 69.8M 17s
    ##  27400K .......... .......... .......... .......... .......... 20% 28.3M 17s
    ##  27450K .......... .......... .......... .......... .......... 20% 37.8M 17s
    ##  27500K .......... .......... .......... .......... .......... 20% 66.9M 17s
    ##  27550K .......... .......... .......... .......... .......... 20% 74.7M 17s
    ##  27600K .......... .......... .......... .......... .......... 20% 43.8M 17s
    ##  27650K .......... .......... .......... .......... .......... 20% 29.4M 17s
    ##  27700K .......... .......... .......... .......... .......... 20% 91.0M 17s
    ##  27750K .......... .......... .......... .......... .......... 20% 43.4M 17s
    ##  27800K .......... .......... .......... .......... .......... 20% 91.4M 17s
    ##  27850K .......... .......... .......... .......... .......... 20% 24.0M 17s
    ##  27900K .......... .......... .......... .......... .......... 20% 90.7M 17s
    ##  27950K .......... .......... .......... .......... .......... 20% 62.5M 17s
    ##  28000K .......... .......... .......... .......... .......... 20% 59.9M 17s
    ##  28050K .......... .......... .......... .......... .......... 20% 30.7M 17s
    ##  28100K .......... .......... .......... .......... .......... 20% 89.1M 17s
    ##  28150K .......... .......... .......... .......... .......... 21% 69.1M 17s
    ##  28200K .......... .......... .......... .......... .......... 21% 94.3M 17s
    ##  28250K .......... .......... .......... .......... .......... 21% 13.2M 17s
    ##  28300K .......... .......... .......... .......... .......... 21% 85.5M 17s
    ##  28350K .......... .......... .......... .......... .......... 21% 79.6M 17s
    ##  28400K .......... .......... .......... .......... .......... 21% 18.7M 17s
    ##  28450K .......... .......... .......... .......... .......... 21% 68.9M 17s
    ##  28500K .......... .......... .......... .......... .......... 21% 77.9M 17s
    ##  28550K .......... .......... .......... .......... .......... 21% 86.8M 17s
    ##  28600K .......... .......... .......... .......... .......... 21% 17.9M 16s
    ##  28650K .......... .......... .......... .......... .......... 21% 69.3M 16s
    ##  28700K .......... .......... .......... .......... .......... 21% 27.5M 16s
    ##  28750K .......... .......... .......... .......... .......... 21% 63.7M 16s
    ##  28800K .......... .......... .......... .......... .......... 21% 57.4M 16s
    ##  28850K .......... .......... .......... .......... .......... 21% 68.8M 16s
    ##  28900K .......... .......... .......... .......... .......... 21% 71.8M 16s
    ##  28950K .......... .......... .......... .......... .......... 21% 86.0M 16s
    ##  29000K .......... .......... .......... .......... .......... 21% 26.6M 16s
    ##  29050K .......... .......... .......... .......... .......... 21% 75.9M 16s
    ##  29100K .......... .......... .......... .......... .......... 21% 71.2M 16s
    ##  29150K .......... .......... .......... .......... .......... 21% 53.7M 16s
    ##  29200K .......... .......... .......... .......... .......... 21% 19.7M 16s
    ##  29250K .......... .......... .......... .......... .......... 21% 78.8M 16s
    ##  29300K .......... .......... .......... .......... .......... 21% 94.6M 16s
    ##  29350K .......... .......... .......... .......... .......... 21% 15.2M 16s
    ##  29400K .......... .......... .......... .......... .......... 21% 84.6M 16s
    ##  29450K .......... .......... .......... .......... .......... 22% 85.7M 16s
    ##  29500K .......... .......... .......... .......... .......... 22% 99.2M 16s
    ##  29550K .......... .......... .......... .......... .......... 22% 21.1M 16s
    ##  29600K .......... .......... .......... .......... .......... 22% 84.4M 16s
    ##  29650K .......... .......... .......... .......... .......... 22% 78.7M 16s
    ##  29700K .......... .......... .......... .......... .......... 22% 83.3M 16s
    ##  29750K .......... .......... .......... .......... .......... 22% 36.0M 16s
    ##  29800K .......... .......... .......... .......... .......... 22% 35.8M 16s
    ##  29850K .......... .......... .......... .......... .......... 22% 43.4M 16s
    ##  29900K .......... .......... .......... .......... .......... 22% 88.0M 16s
    ##  29950K .......... .......... .......... .......... .......... 22% 68.5M 16s
    ##  30000K .......... .......... .......... .......... .......... 22% 38.2M 16s
    ##  30050K .......... .......... .......... .......... .......... 22% 67.1M 16s
    ##  30100K .......... .......... .......... .......... .......... 22% 59.4M 16s
    ##  30150K .......... .......... .......... .......... .......... 22% 82.1M 16s
    ##  30200K .......... .......... .......... .......... .......... 22% 52.8M 15s
    ##  30250K .......... .......... .......... .......... .......... 22% 28.9M 15s
    ##  30300K .......... .......... .......... .......... .......... 22% 89.3M 15s
    ##  30350K .......... .......... .......... .......... .......... 22% 64.3M 15s
    ##  30400K .......... .......... .......... .......... .......... 22% 76.5M 15s
    ##  30450K .......... .......... .......... .......... .......... 22% 32.0M 15s
    ##  30500K .......... .......... .......... .......... .......... 22% 87.0M 15s
    ##  30550K .......... .......... .......... .......... .......... 22% 40.2M 15s
    ##  30600K .......... .......... .......... .......... .......... 22% 97.2M 15s
    ##  30650K .......... .......... .......... .......... .......... 22% 17.6M 15s
    ##  30700K .......... .......... .......... .......... .......... 22% 79.0M 15s
    ##  30750K .......... .......... .......... .......... .......... 22% 82.5M 15s
    ##  30800K .......... .......... .......... .......... .......... 23% 27.3M 15s
    ##  30850K .......... .......... .......... .......... .......... 23% 70.4M 15s
    ##  30900K .......... .......... .......... .......... .......... 23% 88.0M 15s
    ##  30950K .......... .......... .......... .......... .......... 23% 71.5M 15s
    ##  31000K .......... .......... .......... .......... .......... 23% 31.2M 15s
    ##  31050K .......... .......... .......... .......... .......... 23% 79.7M 15s
    ##  31100K .......... .......... .......... .......... .......... 23% 44.7M 15s
    ##  31150K .......... .......... .......... .......... .......... 23% 64.5M 15s
    ##  31200K .......... .......... .......... .......... .......... 23% 61.6M 15s
    ##  31250K .......... .......... .......... .......... .......... 23% 22.0M 15s
    ##  31300K .......... .......... .......... .......... .......... 23% 57.7M 15s
    ##  31350K .......... .......... .......... .......... .......... 23% 76.7M 15s
    ##  31400K .......... .......... .......... .......... .......... 23% 92.6M 15s
    ##  31450K .......... .......... .......... .......... .......... 23% 17.7M 15s
    ##  31500K .......... .......... .......... .......... .......... 23% 71.2M 15s
    ##  31550K .......... .......... .......... .......... .......... 23% 72.7M 15s
    ##  31600K .......... .......... .......... .......... .......... 23%  105M 15s
    ##  31650K .......... .......... .......... .......... .......... 23% 27.5M 15s
    ##  31700K .......... .......... .......... .......... .......... 23% 71.7M 15s
    ##  31750K .......... .......... .......... .......... .......... 23% 57.8M 15s
    ##  31800K .......... .......... .......... .......... .......... 23% 91.2M 15s
    ##  31850K .......... .......... .......... .......... .......... 23% 42.6M 15s
    ##  31900K .......... .......... .......... .......... .......... 23% 27.0M 15s
    ##  31950K .......... .......... .......... .......... .......... 23% 79.1M 15s
    ##  32000K .......... .......... .......... .......... .......... 23% 87.1M 14s
    ##  32050K .......... .......... .......... .......... .......... 23%  106M 14s
    ##  32100K .......... .......... .......... .......... .......... 23% 17.6M 14s
    ##  32150K .......... .......... .......... .......... .......... 24% 81.7M 14s
    ##  32200K .......... .......... .......... .......... .......... 24% 34.2M 14s
    ##  32250K .......... .......... .......... .......... .......... 24% 69.6M 14s
    ##  32300K .......... .......... .......... .......... .......... 24% 5.55M 14s
    ##  32350K .......... .......... .......... .......... .......... 24% 81.5M 14s
    ##  32400K .......... .......... .......... .......... .......... 24% 99.8M 14s
    ##  32450K .......... .......... .......... .......... .......... 24% 17.9M 14s
    ##  32500K .......... .......... .......... .......... .......... 24% 65.6M 14s
    ##  32550K .......... .......... .......... .......... .......... 24% 80.9M 14s
    ##  32600K .......... .......... .......... .......... .......... 24%  102M 14s
    ##  32650K .......... .......... .......... .......... .......... 24% 97.6M 14s
    ##  32700K .......... .......... .......... .......... .......... 24% 40.8M 14s
    ##  32750K .......... .......... .......... .......... .......... 24% 45.4M 14s
    ##  32800K .......... .......... .......... .......... .......... 24% 87.4M 14s
    ##  32850K .......... .......... .......... .......... .......... 24% 4.95M 14s
    ##  32900K .......... .......... .......... .......... .......... 24% 75.1M 14s
    ##  32950K .......... .......... .......... .......... .......... 24% 84.3M 14s
    ##  33000K .......... .......... .......... .......... .......... 24%  108M 14s
    ##  33050K .......... .......... .......... .......... .......... 24% 29.4M 14s
    ##  33100K .......... .......... .......... .......... .......... 24% 81.7M 14s
    ##  33150K .......... .......... .......... .......... .......... 24% 73.5M 14s
    ##  33200K .......... .......... .......... .......... .......... 24%  104M 14s
    ##  33250K .......... .......... .......... .......... .......... 24% 18.7M 14s
    ##  33300K .......... .......... .......... .......... .......... 24% 79.5M 14s
    ##  33350K .......... .......... .......... .......... .......... 24% 82.9M 14s
    ##  33400K .......... .......... .......... .......... .......... 24%  103M 14s
    ##  33450K .......... .......... .......... .......... .......... 24% 7.52M 14s
    ##  33500K .......... .......... .......... .......... .......... 25% 86.0M 14s
    ##  33550K .......... .......... .......... .......... .......... 25% 73.6M 14s
    ##  33600K .......... .......... .......... .......... .......... 25%  111M 14s
    ##  33650K .......... .......... .......... .......... .......... 25% 15.3M 14s
    ##  33700K .......... .......... .......... .......... .......... 25% 89.7M 14s
    ##  33750K .......... .......... .......... .......... .......... 25% 94.9M 14s
    ##  33800K .......... .......... .......... .......... .......... 25% 3.47M 14s
    ##  33850K .......... .......... .......... .......... .......... 25% 74.2M 14s
    ##  33900K .......... .......... .......... .......... .......... 25% 94.7M 14s
    ##  33950K .......... .......... .......... .......... .......... 25% 96.3M 14s
    ##  34000K .......... .......... .......... .......... .......... 25% 10.7M 14s
    ##  34050K .......... .......... .......... .......... .......... 25% 77.7M 14s
    ##  34100K .......... .......... .......... .......... .......... 25% 89.2M 14s
    ##  34150K .......... .......... .......... .......... .......... 25% 97.1M 14s
    ##  34200K .......... .......... .......... .......... .......... 25%  111M 13s
    ##  34250K .......... .......... .......... .......... .......... 25% 20.5M 13s
    ##  34300K .......... .......... .......... .......... .......... 25%  104M 13s
    ##  34350K .......... .......... .......... .......... .......... 25% 83.6M 13s
    ##  34400K .......... .......... .......... .......... .......... 25%  108M 13s
    ##  34450K .......... .......... .......... .......... .......... 25% 15.8M 13s
    ##  34500K .......... .......... .......... .......... .......... 25% 97.3M 13s
    ##  34550K .......... .......... .......... .......... .......... 25% 97.9M 13s
    ##  34600K .......... .......... .......... .......... .......... 25%  106M 13s
    ##  34650K .......... .......... .......... .......... .......... 25% 16.3M 13s
    ##  34700K .......... .......... .......... .......... .......... 25%  107M 13s
    ##  34750K .......... .......... .......... .......... .......... 25% 90.2M 13s
    ##  34800K .......... .......... .......... .......... .......... 25% 18.4M 13s
    ##  34850K .......... .......... .......... .......... .......... 26% 86.4M 13s
    ##  34900K .......... .......... .......... .......... .......... 26% 89.8M 13s
    ##  34950K .......... .......... .......... .......... .......... 26% 94.6M 13s
    ##  35000K .......... .......... .......... .......... .......... 26% 23.0M 13s
    ##  35050K .......... .......... .......... .......... .......... 26% 86.1M 13s
    ##  35100K .......... .......... .......... .......... .......... 26% 87.3M 13s
    ##  35150K .......... .......... .......... .......... .......... 26% 86.3M 13s
    ##  35200K .......... .......... .......... .......... .......... 26% 13.0M 13s
    ##  35250K .......... .......... .......... .......... .......... 26% 80.9M 13s
    ##  35300K .......... .......... .......... .......... .......... 26% 91.7M 13s
    ##  35350K .......... .......... .......... .......... .......... 26% 98.1M 13s
    ##  35400K .......... .......... .......... .......... .......... 26% 7.09M 13s
    ##  35450K .......... .......... .......... .......... .......... 26% 96.7M 13s
    ##  35500K .......... .......... .......... .......... .......... 26% 80.4M 13s
    ##  35550K .......... .......... .......... .......... .......... 26% 87.4M 13s
    ##  35600K .......... .......... .......... .......... .......... 26% 6.40M 13s
    ##  35650K .......... .......... .......... .......... .......... 26%  103M 13s
    ##  35700K .......... .......... .......... .......... .......... 26%  118M 13s
    ##  35750K .......... .......... .......... .......... .......... 26%  102M 13s
    ##  35800K .......... .......... .......... .......... .......... 26% 2.67M 13s
    ##  35850K .......... .......... .......... .......... .......... 26%  111M 13s
    ##  35900K .......... .......... .......... .......... .......... 26%  118M 13s
    ##  35950K .......... .......... .......... .......... .......... 26%  100M 13s
    ##  36000K .......... .......... .......... .......... .......... 26%  119M 13s
    ##  36050K .......... .......... .......... .......... .......... 26%  115M 13s
    ##  36100K .......... .......... .......... .......... .......... 26% 9.19M 13s
    ##  36150K .......... .......... .......... .......... .......... 27% 77.7M 13s
    ##  36200K .......... .......... .......... .......... .......... 27%  110M 13s
    ##  36250K .......... .......... .......... .......... .......... 27%  124M 13s
    ##  36300K .......... .......... .......... .......... .......... 27%  123M 13s
    ##  36350K .......... .......... .......... .......... .......... 27% 5.78M 13s
    ##  36400K .......... .......... .......... .......... .......... 27% 86.9M 13s
    ##  36450K .......... .......... .......... .......... .......... 27%  101M 13s
    ##  36500K .......... .......... .......... .......... .......... 27%  112M 13s
    ##  36550K .......... .......... .......... .......... .......... 27% 3.73M 13s
    ##  36600K .......... .......... .......... .......... .......... 27%  103M 13s
    ##  36650K .......... .......... .......... .......... .......... 27% 80.1M 13s
    ##  36700K .......... .......... .......... .......... .......... 27% 95.5M 13s
    ##  36750K .......... .......... .......... .......... .......... 27% 17.0M 12s
    ##  36800K .......... .......... .......... .......... .......... 27% 86.3M 12s
    ##  36850K .......... .......... .......... .......... .......... 27%  121M 12s
    ##  36900K .......... .......... .......... .......... .......... 27%  103M 12s
    ##  36950K .......... .......... .......... .......... .......... 27%  102M 12s
    ##  37000K .......... .......... .......... .......... .......... 27% 12.0M 12s
    ##  37050K .......... .......... .......... .......... .......... 27% 66.4M 12s
    ##  37100K .......... .......... .......... .......... .......... 27% 96.4M 12s
    ##  37150K .......... .......... .......... .......... .......... 27% 4.00M 12s
    ##  37200K .......... .......... .......... .......... .......... 27% 39.9M 12s
    ##  37250K .......... .......... .......... .......... .......... 27% 92.8M 12s
    ##  37300K .......... .......... .......... .......... .......... 27% 82.5M 12s
    ##  37350K .......... .......... .......... .......... .......... 27% 3.11M 12s
    ##  37400K .......... .......... .......... .......... .......... 27% 71.8M 12s
    ##  37450K .......... .......... .......... .......... .......... 27% 87.0M 12s
    ##  37500K .......... .......... .......... .......... .......... 28% 97.0M 12s
    ##  37550K .......... .......... .......... .......... .......... 28%  104M 12s
    ##  37600K .......... .......... .......... .......... .......... 28% 15.0M 12s
    ##  37650K .......... .......... .......... .......... .......... 28%  115M 12s
    ##  37700K .......... .......... .......... .......... .......... 28% 83.9M 12s
    ##  37750K .......... .......... .......... .......... .......... 28%  105M 12s
    ##  37800K .......... .......... .......... .......... .......... 28% 6.34M 12s
    ##  37850K .......... .......... .......... .......... .......... 28% 73.0M 12s
    ##  37900K .......... .......... .......... .......... .......... 28%  116M 12s
    ##  37950K .......... .......... .......... .......... .......... 28% 92.0M 12s
    ##  38000K .......... .......... .......... .......... .......... 28% 5.91M 12s
    ##  38050K .......... .......... .......... .......... .......... 28%  111M 12s
    ##  38100K .......... .......... .......... .......... .......... 28%  103M 12s
    ##  38150K .......... .......... .......... .......... .......... 28%  104M 12s
    ##  38200K .......... .......... .......... .......... .......... 28% 6.52M 12s
    ##  38250K .......... .......... .......... .......... .......... 28%  100M 12s
    ##  38300K .......... .......... .......... .......... .......... 28%  104M 12s
    ##  38350K .......... .......... .......... .......... .......... 28% 19.0M 12s
    ##  38400K .......... .......... .......... .......... .......... 28% 89.9M 12s
    ##  38450K .......... .......... .......... .......... .......... 28%  110M 12s
    ##  38500K .......... .......... .......... .......... .......... 28%  117M 12s
    ##  38550K .......... .......... .......... .......... .......... 28% 15.0M 12s
    ##  38600K .......... .......... .......... .......... .......... 28% 79.1M 12s
    ##  38650K .......... .......... .......... .......... .......... 28% 83.2M 12s
    ##  38700K .......... .......... .......... .......... .......... 28% 91.9M 12s
    ##  38750K .......... .......... .......... .......... .......... 28% 28.5M 12s
    ##  38800K .......... .......... .......... .......... .......... 28% 46.4M 12s
    ##  38850K .......... .......... .......... .......... .......... 29% 95.3M 12s
    ##  38900K .......... .......... .......... .......... .......... 29% 16.4M 12s
    ##  38950K .......... .......... .......... .......... .......... 29% 11.0M 12s
    ##  39000K .......... .......... .......... .......... .......... 29%  116M 12s
    ##  39050K .......... .......... .......... .......... .......... 29%  136M 12s
    ##  39100K .......... .......... .......... .......... .......... 29%  143M 12s
    ##  39150K .......... .......... .......... .......... .......... 29% 8.14M 12s
    ##  39200K .......... .......... .......... .......... .......... 29%  108M 12s
    ##  39250K .......... .......... .......... .......... .......... 29%  132M 12s
    ##  39300K .......... .......... .......... .......... .......... 29%  134M 12s
    ##  39350K .......... .......... .......... .......... .......... 29% 9.48M 12s
    ##  39400K .......... .......... .......... .......... .......... 29% 39.6M 12s
    ##  39450K .......... .......... .......... .......... .......... 29%  125M 12s
    ##  39500K .......... .......... .......... .......... .......... 29%  137M 12s
    ##  39550K .......... .......... .......... .......... .......... 29% 17.5M 12s
    ##  39600K .......... .......... .......... .......... .......... 29%  110M 12s
    ##  39650K .......... .......... .......... .......... .......... 29%  121M 12s
    ##  39700K .......... .......... .......... .......... .......... 29% 8.09M 11s
    ##  39750K .......... .......... .......... .......... .......... 29% 76.3M 11s
    ##  39800K .......... .......... .......... .......... .......... 29%  114M 11s
    ##  39850K .......... .......... .......... .......... .......... 29%  122M 11s
    ##  39900K .......... .......... .......... .......... .......... 29%  124M 11s
    ##  39950K .......... .......... .......... .......... .......... 29% 33.1M 11s
    ##  40000K .......... .......... .......... .......... .......... 29% 99.7M 11s
    ##  40050K .......... .......... .......... .......... .......... 29% 50.7M 11s
    ##  40100K .......... .......... .......... .......... .......... 29% 89.1M 11s
    ##  40150K .......... .......... .......... .......... .......... 29% 27.1M 11s
    ##  40200K .......... .......... .......... .......... .......... 30%  121M 11s
    ##  40250K .......... .......... .......... .......... .......... 30% 23.7M 11s
    ##  40300K .......... .......... .......... .......... .......... 30%  112M 11s
    ##  40350K .......... .......... .......... .......... .......... 30% 56.3M 11s
    ##  40400K .......... .......... .......... .......... .......... 30% 19.3M 11s
    ##  40450K .......... .......... .......... .......... .......... 30% 91.8M 11s
    ##  40500K .......... .......... .......... .......... .......... 30%  123M 11s
    ##  40550K .......... .......... .......... .......... .......... 30% 90.0M 11s
    ##  40600K .......... .......... .......... .......... .......... 30% 18.4M 11s
    ##  40650K .......... .......... .......... .......... .......... 30% 93.9M 11s
    ##  40700K .......... .......... .......... .......... .......... 30%  108M 11s
    ##  40750K .......... .......... .......... .......... .......... 30% 64.7M 11s
    ##  40800K .......... .......... .......... .......... .......... 30% 31.5M 11s
    ##  40850K .......... .......... .......... .......... .......... 30% 49.0M 11s
    ##  40900K .......... .......... .......... .......... .......... 30% 74.9M 11s
    ##  40950K .......... .......... .......... .......... .......... 30% 89.9M 11s
    ##  41000K .......... .......... .......... .......... .......... 30% 45.8M 11s
    ##  41050K .......... .......... .......... .......... .......... 30% 46.9M 11s
    ##  41100K .......... .......... .......... .......... .......... 30% 45.7M 11s
    ##  41150K .......... .......... .......... .......... .......... 30% 39.8M 11s
    ##  41200K .......... .......... .......... .......... .......... 30%  107M 11s
    ##  41250K .......... .......... .......... .......... .......... 30% 84.8M 11s
    ##  41300K .......... .......... .......... .......... .......... 30% 34.9M 11s
    ##  41350K .......... .......... .......... .......... .......... 30% 45.7M 11s
    ##  41400K .......... .......... .......... .......... .......... 30%  131M 11s
    ##  41450K .......... .......... .......... .......... .......... 30% 91.7M 11s
    ##  41500K .......... .......... .......... .......... .......... 30% 31.2M 11s
    ##  41550K .......... .......... .......... .......... .......... 31% 28.3M 11s
    ##  41600K .......... .......... .......... .......... .......... 31%  131M 11s
    ##  41650K .......... .......... .......... .......... .......... 31%  129M 11s
    ##  41700K .......... .......... .......... .......... .......... 31% 40.6M 11s
    ##  41750K .......... .......... .......... .......... .......... 31% 31.2M 11s
    ##  41800K .......... .......... .......... .......... .......... 31%  112M 11s
    ##  41850K .......... .......... .......... .......... .......... 31%  138M 11s
    ##  41900K .......... .......... .......... .......... .......... 31% 45.2M 11s
    ##  41950K .......... .......... .......... .......... .......... 31% 33.7M 11s
    ##  42000K .......... .......... .......... .......... .......... 31%  106M 11s
    ##  42050K .......... .......... .......... .......... .......... 31% 31.5M 11s
    ##  42100K .......... .......... .......... .......... .......... 31%  116M 11s
    ##  42150K .......... .......... .......... .......... .......... 31% 21.3M 11s
    ##  42200K .......... .......... .......... .......... .......... 31%  127M 11s
    ##  42250K .......... .......... .......... .......... .......... 31% 74.8M 11s
    ##  42300K .......... .......... .......... .......... .......... 31%  110M 11s
    ##  42350K .......... .......... .......... .......... .......... 31% 24.7M 11s
    ##  42400K .......... .......... .......... .......... .......... 31%  142M 11s
    ##  42450K .......... .......... .......... .......... .......... 31% 75.2M 11s
    ##  42500K .......... .......... .......... .......... .......... 31% 23.6M 11s
    ##  42550K .......... .......... .......... .......... .......... 31% 93.2M 11s
    ##  42600K .......... .......... .......... .......... .......... 31%  118M 11s
    ##  42650K .......... .......... .......... .......... .......... 31%  128M 10s
    ##  42700K .......... .......... .......... .......... .......... 31% 24.1M 10s
    ##  42750K .......... .......... .......... .......... .......... 31% 31.3M 10s
    ##  42800K .......... .......... .......... .......... .......... 31%  126M 10s
    ##  42850K .......... .......... .......... .......... .......... 31%  135M 10s
    ##  42900K .......... .......... .......... .......... .......... 32% 36.6M 10s
    ##  42950K .......... .......... .......... .......... .......... 32% 30.3M 10s
    ##  43000K .......... .......... .......... .......... .......... 32%  134M 10s
    ##  43050K .......... .......... .......... .......... .......... 32%  134M 10s
    ##  43100K .......... .......... .......... .......... .......... 32% 63.5M 10s
    ##  43150K .......... .......... .......... .......... .......... 32% 19.6M 10s
    ##  43200K .......... .......... .......... .......... .......... 32%  113M 10s
    ##  43250K .......... .......... .......... .......... .......... 32%  131M 10s
    ##  43300K .......... .......... .......... .......... .......... 32%  126M 10s
    ##  43350K .......... .......... .......... .......... .......... 32% 23.1M 10s
    ##  43400K .......... .......... .......... .......... .......... 32% 97.0M 10s
    ##  43450K .......... .......... .......... .......... .......... 32% 79.5M 10s
    ##  43500K .......... .......... .......... .......... .......... 32%  102M 10s
    ##  43550K .......... .......... .......... .......... .......... 32% 25.0M 10s
    ##  43600K .......... .......... .......... .......... .......... 32%  132M 10s
    ##  43650K .......... .......... .......... .......... .......... 32% 19.5M 10s
    ##  43700K .......... .......... .......... .......... .......... 32% 18.1M 10s
    ##  43750K .......... .......... .......... .......... .......... 32% 61.4M 10s
    ##  43800K .......... .......... .......... .......... .......... 32% 74.9M 10s
    ##  43850K .......... .......... .......... .......... .......... 32% 6.80M 10s
    ##  43900K .......... .......... .......... .......... .......... 32% 68.1M 10s
    ##  43950K .......... .......... .......... .......... .......... 32% 59.3M 10s
    ##  44000K .......... .......... .......... .......... .......... 32% 79.2M 10s
    ##  44050K .......... .......... .......... .......... .......... 32% 16.2M 10s
    ##  44100K .......... .......... .......... .......... .......... 32% 78.1M 10s
    ##  44150K .......... .......... .......... .......... .......... 32% 67.1M 10s
    ##  44200K .......... .......... .......... .......... .......... 33% 63.7M 10s
    ##  44250K .......... .......... .......... .......... .......... 33% 13.5M 10s
    ##  44300K .......... .......... .......... .......... .......... 33% 57.9M 10s
    ##  44350K .......... .......... .......... .......... .......... 33% 52.4M 10s
    ##  44400K .......... .......... .......... .......... .......... 33% 58.6M 10s
    ##  44450K .......... .......... .......... .......... .......... 33% 66.0M 10s
    ##  44500K .......... .......... .......... .......... .......... 33% 75.2M 10s
    ##  44550K .......... .......... .......... .......... .......... 33% 75.0M 10s
    ##  44600K .......... .......... .......... .......... .......... 33% 65.9M 10s
    ##  44650K .......... .......... .......... .......... .......... 33% 65.3M 10s
    ##  44700K .......... .......... .......... .......... .......... 33% 85.1M 10s
    ##  44750K .......... .......... .......... .......... .......... 33% 51.4M 10s
    ##  44800K .......... .......... .......... .......... .......... 33% 37.3M 10s
    ##  44850K .......... .......... .......... .......... .......... 33% 24.0M 10s
    ##  44900K .......... .......... .......... .......... .......... 33% 58.0M 10s
    ##  44950K .......... .......... .......... .......... .......... 33% 59.6M 10s
    ##  45000K .......... .......... .......... .......... .......... 33% 52.1M 10s
    ##  45050K .......... .......... .......... .......... .......... 33% 53.6M 10s
    ##  45100K .......... .......... .......... .......... .......... 33% 60.7M 10s
    ##  45150K .......... .......... .......... .......... .......... 33% 5.65M 10s
    ##  45200K .......... .......... .......... .......... .......... 33% 57.8M 10s
    ##  45250K .......... .......... .......... .......... .......... 33% 53.4M 10s
    ##  45300K .......... .......... .......... .......... .......... 33% 62.3M 10s
    ##  45350K .......... .......... .......... .......... .......... 33% 54.0M 10s
    ##  45400K .......... .......... .......... .......... .......... 33% 60.4M 10s
    ##  45450K .......... .......... .......... .......... .......... 33% 23.9M 10s
    ##  45500K .......... .......... .......... .......... .......... 33% 83.5M 10s
    ##  45550K .......... .......... .......... .......... .......... 34% 69.2M 10s
    ##  45600K .......... .......... .......... .......... .......... 34% 87.4M 10s
    ##  45650K .......... .......... .......... .......... .......... 34% 90.0M 10s
    ##  45700K .......... .......... .......... .......... .......... 34% 44.0M 10s
    ##  45750K .......... .......... .......... .......... .......... 34% 67.5M 10s
    ##  45800K .......... .......... .......... .......... .......... 34% 73.5M 10s
    ##  45850K .......... .......... .......... .......... .......... 34% 57.4M 10s
    ##  45900K .......... .......... .......... .......... .......... 34% 88.7M 10s
    ##  45950K .......... .......... .......... .......... .......... 34% 45.3M 10s
    ##  46000K .......... .......... .......... .......... .......... 34% 31.8M 10s
    ##  46050K .......... .......... .......... .......... .......... 34% 70.5M 10s
    ##  46100K .......... .......... .......... .......... .......... 34% 95.1M 9s
    ##  46150K .......... .......... .......... .......... .......... 34% 78.5M 9s
    ##  46200K .......... .......... .......... .......... .......... 34% 93.8M 9s
    ##  46250K .......... .......... .......... .......... .......... 34% 23.0M 9s
    ##  46300K .......... .......... .......... .......... .......... 34% 51.1M 9s
    ##  46350K .......... .......... .......... .......... .......... 34% 17.1M 9s
    ##  46400K .......... .......... .......... .......... .......... 34% 50.3M 9s
    ##  46450K .......... .......... .......... .......... .......... 34% 56.1M 9s
    ##  46500K .......... .......... .......... .......... .......... 34% 54.8M 9s
    ##  46550K .......... .......... .......... .......... .......... 34% 50.9M 9s
    ##  46600K .......... .......... .......... .......... .......... 34% 15.1M 9s
    ##  46650K .......... .......... .......... .......... .......... 34% 55.9M 9s
    ##  46700K .......... .......... .......... .......... .......... 34% 55.4M 9s
    ##  46750K .......... .......... .......... .......... .......... 34% 56.7M 9s
    ##  46800K .......... .......... .......... .......... .......... 34% 80.9M 9s
    ##  46850K .......... .......... .......... .......... .......... 34% 83.8M 9s
    ##  46900K .......... .......... .......... .......... .......... 35% 29.9M 9s
    ##  46950K .......... .......... .......... .......... .......... 35% 66.9M 9s
    ##  47000K .......... .......... .......... .......... .......... 35% 74.6M 9s
    ##  47050K .......... .......... .......... .......... .......... 35% 63.4M 9s
    ##  47100K .......... .......... .......... .......... .......... 35% 63.7M 9s
    ##  47150K .......... .......... .......... .......... .......... 35% 39.6M 9s
    ##  47200K .......... .......... .......... .......... .......... 35% 63.3M 9s
    ##  47250K .......... .......... .......... .......... .......... 35% 63.6M 9s
    ##  47300K .......... .......... .......... .......... .......... 35% 87.4M 9s
    ##  47350K .......... .......... .......... .......... .......... 35% 44.4M 9s
    ##  47400K .......... .......... .......... .......... .......... 35% 66.2M 9s
    ##  47450K .......... .......... .......... .......... .......... 35% 81.7M 9s
    ##  47500K .......... .......... .......... .......... .......... 35% 86.8M 9s
    ##  47550K .......... .......... .......... .......... .......... 35% 77.8M 9s
    ##  47600K .......... .......... .......... .......... .......... 35% 16.5M 9s
    ##  47650K .......... .......... .......... .......... .......... 35% 71.7M 9s
    ##  47700K .......... .......... .......... .......... .......... 35% 67.7M 9s
    ##  47750K .......... .......... .......... .......... .......... 35% 77.1M 9s
    ##  47800K .......... .......... .......... .......... .......... 35% 88.5M 9s
    ##  47850K .......... .......... .......... .......... .......... 35% 14.7M 9s
    ##  47900K .......... .......... .......... .......... .......... 35% 61.6M 9s
    ##  47950K .......... .......... .......... .......... .......... 35% 64.9M 9s
    ##  48000K .......... .......... .......... .......... .......... 35% 85.3M 9s
    ##  48050K .......... .......... .......... .......... .......... 35% 83.1M 9s
    ##  48100K .......... .......... .......... .......... .......... 35% 41.5M 9s
    ##  48150K .......... .......... .......... .......... .......... 35% 56.2M 9s
    ##  48200K .......... .......... .......... .......... .......... 35% 62.0M 9s
    ##  48250K .......... .......... .......... .......... .......... 36% 54.5M 9s
    ##  48300K .......... .......... .......... .......... .......... 36% 77.8M 9s
    ##  48350K .......... .......... .......... .......... .......... 36% 58.7M 9s
    ##  48400K .......... .......... .......... .......... .......... 36% 30.9M 9s
    ##  48450K .......... .......... .......... .......... .......... 36% 73.0M 9s
    ##  48500K .......... .......... .......... .......... .......... 36% 16.0M 9s
    ##  48550K .......... .......... .......... .......... .......... 36% 57.5M 9s
    ##  48600K .......... .......... .......... .......... .......... 36% 77.2M 9s
    ##  48650K .......... .......... .......... .......... .......... 36% 68.8M 9s
    ##  48700K .......... .......... .......... .......... .......... 36% 86.6M 9s
    ##  48750K .......... .......... .......... .......... .......... 36% 44.8M 9s
    ##  48800K .......... .......... .......... .......... .......... 36% 76.4M 9s
    ##  48850K .......... .......... .......... .......... .......... 36% 47.2M 9s
    ##  48900K .......... .......... .......... .......... .......... 36% 61.3M 9s
    ##  48950K .......... .......... .......... .......... .......... 36% 75.5M 9s
    ##  49000K .......... .......... .......... .......... .......... 36% 51.5M 9s
    ##  49050K .......... .......... .......... .......... .......... 36% 20.9M 9s
    ##  49100K .......... .......... .......... .......... .......... 36% 64.4M 9s
    ##  49150K .......... .......... .......... .......... .......... 36% 57.0M 9s
    ##  49200K .......... .......... .......... .......... .......... 36% 87.3M 9s
    ##  49250K .......... .......... .......... .......... .......... 36% 91.9M 9s
    ##  49300K .......... .......... .......... .......... .......... 36% 46.6M 9s
    ##  49350K .......... .......... .......... .......... .......... 36% 68.9M 9s
    ##  49400K .......... .......... .......... .......... .......... 36% 79.1M 9s
    ##  49450K .......... .......... .......... .......... .......... 36% 28.6M 9s
    ##  49500K .......... .......... .......... .......... .......... 36% 88.6M 9s
    ##  49550K .......... .......... .......... .......... .......... 36% 46.9M 9s
    ##  49600K .......... .......... .......... .......... .......... 37% 15.6M 9s
    ##  49650K .......... .......... .......... .......... .......... 37% 65.2M 9s
    ##  49700K .......... .......... .......... .......... .......... 37% 74.6M 9s
    ##  49750K .......... .......... .......... .......... .......... 37% 65.6M 9s
    ##  49800K .......... .......... .......... .......... .......... 37% 94.3M 9s
    ##  49850K .......... .......... .......... .......... .......... 37% 99.7M 9s
    ##  49900K .......... .......... .......... .......... .......... 37% 37.1M 9s
    ##  49950K .......... .......... .......... .......... .......... 37% 25.1M 9s
    ##  50000K .......... .......... .......... .......... .......... 37% 82.6M 8s
    ##  50050K .......... .......... .......... .......... .......... 37% 94.9M 8s
    ##  50100K .......... .......... .......... .......... .......... 37% 94.0M 8s
    ##  50150K .......... .......... .......... .......... .......... 37% 83.3M 8s
    ##  50200K .......... .......... .......... .......... .......... 37% 22.1M 8s
    ##  50250K .......... .......... .......... .......... .......... 37% 66.5M 8s
    ##  50300K .......... .......... .......... .......... .......... 37% 96.6M 8s
    ##  50350K .......... .......... .......... .......... .......... 37% 80.9M 8s
    ##  50400K .......... .......... .......... .......... .......... 37% 94.6M 8s
    ##  50450K .......... .......... .......... .......... .......... 37% 28.1M 8s
    ##  50500K .......... .......... .......... .......... .......... 37% 80.5M 8s
    ##  50550K .......... .......... .......... .......... .......... 37% 69.5M 8s
    ##  50600K .......... .......... .......... .......... .......... 37% 79.4M 8s
    ##  50650K .......... .......... .......... .......... .......... 37% 95.1M 8s
    ##  50700K .......... .......... .......... .......... .......... 37% 73.7M 8s
    ##  50750K .......... .......... .......... .......... .......... 37% 24.1M 8s
    ##  50800K .......... .......... .......... .......... .......... 37% 81.5M 8s
    ##  50850K .......... .......... .......... .......... .......... 37% 97.1M 8s
    ##  50900K .......... .......... .......... .......... .......... 38%  102M 8s
    ##  50950K .......... .......... .......... .......... .......... 38% 83.7M 8s
    ##  51000K .......... .......... .......... .......... .......... 38% 98.0M 8s
    ##  51050K .......... .......... .......... .......... .......... 38% 32.9M 8s
    ##  51100K .......... .......... .......... .......... .......... 38% 18.4M 8s
    ##  51150K .......... .......... .......... .......... .......... 38% 84.0M 8s
    ##  51200K .......... .......... .......... .......... .......... 38% 96.4M 8s
    ##  51250K .......... .......... .......... .......... .......... 38% 64.0M 8s
    ##  51300K .......... .......... .......... .......... .......... 38% 56.4M 8s
    ##  51350K .......... .......... .......... .......... .......... 38% 36.6M 8s
    ##  51400K .......... .......... .......... .......... .......... 38% 44.9M 8s
    ##  51450K .......... .......... .......... .......... .......... 38% 52.0M 8s
    ##  51500K .......... .......... .......... .......... .......... 38% 58.5M 8s
    ##  51550K .......... .......... .......... .......... .......... 38% 46.3M 8s
    ##  51600K .......... .......... .......... .......... .......... 38% 86.3M 8s
    ##  51650K .......... .......... .......... .......... .......... 38% 31.9M 8s
    ##  51700K .......... .......... .......... .......... .......... 38%  104M 8s
    ##  51750K .......... .......... .......... .......... .......... 38% 90.9M 8s
    ##  51800K .......... .......... .......... .......... .......... 38% 48.3M 8s
    ##  51850K .......... .......... .......... .......... .......... 38% 90.6M 8s
    ##  51900K .......... .......... .......... .......... .......... 38% 43.4M 8s
    ##  51950K .......... .......... .......... .......... .......... 38% 65.0M 8s
    ##  52000K .......... .......... .......... .......... .......... 38% 97.5M 8s
    ##  52050K .......... .......... .......... .......... .......... 38% 91.2M 8s
    ##  52100K .......... .......... .......... .......... .......... 38% 67.3M 8s
    ##  52150K .......... .......... .......... .......... .......... 38% 95.4M 8s
    ##  52200K .......... .......... .......... .......... .......... 38% 37.6M 8s
    ##  52250K .......... .......... .......... .......... .......... 39% 94.3M 8s
    ##  52300K .......... .......... .......... .......... .......... 39% 92.1M 8s
    ##  52350K .......... .......... .......... .......... .......... 39% 34.2M 8s
    ##  52400K .......... .......... .......... .......... .......... 39% 40.9M 8s
    ##  52450K .......... .......... .......... .......... .......... 39% 83.2M 8s
    ##  52500K .......... .......... .......... .......... .......... 39%  103M 8s
    ##  52550K .......... .......... .......... .......... .......... 39% 80.2M 8s
    ##  52600K .......... .......... .......... .......... .......... 39% 38.6M 8s
    ##  52650K .......... .......... .......... .......... .......... 39%  117M 8s
    ##  52700K .......... .......... .......... .......... .......... 39% 43.7M 8s
    ##  52750K .......... .......... .......... .......... .......... 39% 88.1M 8s
    ##  52800K .......... .......... .......... .......... .......... 39% 93.2M 8s
    ##  52850K .......... .......... .......... .......... .......... 39% 23.2M 8s
    ##  52900K .......... .......... .......... .......... .......... 39% 92.2M 8s
    ##  52950K .......... .......... .......... .......... .......... 39% 76.2M 8s
    ##  53000K .......... .......... .......... .......... .......... 39%  100M 8s
    ##  53050K .......... .......... .......... .......... .......... 39% 92.8M 8s
    ##  53100K .......... .......... .......... .......... .......... 39%  118M 8s
    ##  53150K .......... .......... .......... .......... .......... 39% 16.9M 8s
    ##  53200K .......... .......... .......... .......... .......... 39% 92.8M 8s
    ##  53250K .......... .......... .......... .......... .......... 39%  119M 8s
    ##  53300K .......... .......... .......... .......... .......... 39% 15.0M 8s
    ##  53350K .......... .......... .......... .......... .......... 39% 77.9M 8s
    ##  53400K .......... .......... .......... .......... .......... 39%  125M 8s
    ##  53450K .......... .......... .......... .......... .......... 39% 71.9M 8s
    ##  53500K .......... .......... .......... .......... .......... 39% 99.7M 8s
    ##  53550K .......... .......... .......... .......... .......... 39% 27.3M 8s
    ##  53600K .......... .......... .......... .......... .......... 40% 88.9M 8s
    ##  53650K .......... .......... .......... .......... .......... 40%  116M 8s
    ##  53700K .......... .......... .......... .......... .......... 40% 65.7M 8s
    ##  53750K .......... .......... .......... .......... .......... 40%  100M 8s
    ##  53800K .......... .......... .......... .......... .......... 40% 73.4M 8s
    ##  53850K .......... .......... .......... .......... .......... 40% 86.2M 8s
    ##  53900K .......... .......... .......... .......... .......... 40% 48.5M 8s
    ##  53950K .......... .......... .......... .......... .......... 40% 24.0M 8s
    ##  54000K .......... .......... .......... .......... .......... 40% 93.3M 8s
    ##  54050K .......... .......... .......... .......... .......... 40%  127M 8s
    ##  54100K .......... .......... .......... .......... .......... 40%  127M 8s
    ##  54150K .......... .......... .......... .......... .......... 40% 88.3M 8s
    ##  54200K .......... .......... .......... .......... .......... 40% 20.5M 8s
    ##  54250K .......... .......... .......... .......... .......... 40% 93.0M 8s
    ##  54300K .......... .......... .......... .......... .......... 40%  102M 8s
    ##  54350K .......... .......... .......... .......... .......... 40% 85.6M 8s
    ##  54400K .......... .......... .......... .......... .......... 40%  127M 8s
    ##  54450K .......... .......... .......... .......... .......... 40% 21.3M 8s
    ##  54500K .......... .......... .......... .......... .......... 40% 89.9M 7s
    ##  54550K .......... .......... .......... .......... .......... 40% 96.0M 7s
    ##  54600K .......... .......... .......... .......... .......... 40% 83.2M 7s
    ##  54650K .......... .......... .......... .......... .......... 40% 48.2M 7s
    ##  54700K .......... .......... .......... .......... .......... 40% 35.6M 7s
    ##  54750K .......... .......... .......... .......... .......... 40% 82.8M 7s
    ##  54800K .......... .......... .......... .......... .......... 40%  110M 7s
    ##  54850K .......... .......... .......... .......... .......... 40% 97.2M 7s
    ##  54900K .......... .......... .......... .......... .......... 40% 93.6M 7s
    ##  54950K .......... .......... .......... .......... .......... 41% 11.5M 7s
    ##  55000K .......... .......... .......... .......... .......... 41% 59.7M 7s
    ##  55050K .......... .......... .......... .......... .......... 41%  100M 7s
    ##  55100K .......... .......... .......... .......... .......... 41% 67.9M 7s
    ##  55150K .......... .......... .......... .......... .......... 41% 80.3M 7s
    ##  55200K .......... .......... .......... .......... .......... 41%  125M 7s
    ##  55250K .......... .......... .......... .......... .......... 41% 27.0M 7s
    ##  55300K .......... .......... .......... .......... .......... 41% 58.0M 7s
    ##  55350K .......... .......... .......... .......... .......... 41% 36.7M 7s
    ##  55400K .......... .......... .......... .......... .......... 41% 94.7M 7s
    ##  55450K .......... .......... .......... .......... .......... 41%  116M 7s
    ##  55500K .......... .......... .......... .......... .......... 41% 89.8M 7s
    ##  55550K .......... .......... .......... .......... .......... 41% 46.2M 7s
    ##  55600K .......... .......... .......... .......... .......... 41% 29.4M 7s
    ##  55650K .......... .......... .......... .......... .......... 41%  129M 7s
    ##  55700K .......... .......... .......... .......... .......... 41%  114M 7s
    ##  55750K .......... .......... .......... .......... .......... 41%  113M 7s
    ##  55800K .......... .......... .......... .......... .......... 41% 90.6M 7s
    ##  55850K .......... .......... .......... .......... .......... 41% 12.0M 7s
    ##  55900K .......... .......... .......... .......... .......... 41% 96.9M 7s
    ##  55950K .......... .......... .......... .......... .......... 41% 49.0M 7s
    ##  56000K .......... .......... .......... .......... .......... 41% 93.7M 7s
    ##  56050K .......... .......... .......... .......... .......... 41%  120M 7s
    ##  56100K .......... .......... .......... .......... .......... 41% 38.3M 7s
    ##  56150K .......... .......... .......... .......... .......... 41% 57.6M 7s
    ##  56200K .......... .......... .......... .......... .......... 41%  112M 7s
    ##  56250K .......... .......... .......... .......... .......... 41% 12.8M 7s
    ##  56300K .......... .......... .......... .......... .......... 42%  117M 7s
    ##  56350K .......... .......... .......... .......... .......... 42% 49.5M 7s
    ##  56400K .......... .......... .......... .......... .......... 42% 85.5M 7s
    ##  56450K .......... .......... .......... .......... .......... 42%  130M 7s
    ##  56500K .......... .......... .......... .......... .......... 42% 15.7M 7s
    ##  56550K .......... .......... .......... .......... .......... 42% 87.4M 7s
    ##  56600K .......... .......... .......... .......... .......... 42% 67.6M 7s
    ##  56650K .......... .......... .......... .......... .......... 42% 90.0M 7s
    ##  56700K .......... .......... .......... .......... .......... 42%  123M 7s
    ##  56750K .......... .......... .......... .......... .......... 42% 29.0M 7s
    ##  56800K .......... .......... .......... .......... .......... 42% 21.8M 7s
    ##  56850K .......... .......... .......... .......... .......... 42% 67.7M 7s
    ##  56900K .......... .......... .......... .......... .......... 42% 74.2M 7s
    ##  56950K .......... .......... .......... .......... .......... 42% 94.5M 7s
    ##  57000K .......... .......... .......... .......... .......... 42%  108M 7s
    ##  57050K .......... .......... .......... .......... .......... 42% 39.0M 7s
    ##  57100K .......... .......... .......... .......... .......... 42%  119M 7s
    ##  57150K .......... .......... .......... .......... .......... 42% 36.9M 7s
    ##  57200K .......... .......... .......... .......... .......... 42% 27.4M 7s
    ##  57250K .......... .......... .......... .......... .......... 42%  135M 7s
    ##  57300K .......... .......... .......... .......... .......... 42% 48.3M 7s
    ##  57350K .......... .......... .......... .......... .......... 42%  115M 7s
    ##  57400K .......... .......... .......... .......... .......... 42% 34.7M 7s
    ##  57450K .......... .......... .......... .......... .......... 42%  107M 7s
    ##  57500K .......... .......... .......... .......... .......... 42% 69.4M 7s
    ##  57550K .......... .......... .......... .......... .......... 42% 51.8M 7s
    ##  57600K .......... .......... .......... .......... .......... 43% 30.8M 7s
    ##  57650K .......... .......... .......... .......... .......... 43%  109M 7s
    ##  57700K .......... .......... .......... .......... .......... 43%  110M 7s
    ##  57750K .......... .......... .......... .......... .......... 43% 78.2M 7s
    ##  57800K .......... .......... .......... .......... .......... 43% 37.1M 7s
    ##  57850K .......... .......... .......... .......... .......... 43%  127M 7s
    ##  57900K .......... .......... .......... .......... .......... 43% 47.5M 7s
    ##  57950K .......... .......... .......... .......... .......... 43% 53.0M 7s
    ##  58000K .......... .......... .......... .......... .......... 43% 95.6M 7s
    ##  58050K .......... .......... .......... .......... .......... 43% 59.0M 7s
    ##  58100K .......... .......... .......... .......... .......... 43% 99.1M 7s
    ##  58150K .......... .......... .......... .......... .......... 43% 32.8M 7s
    ##  58200K .......... .......... .......... .......... .......... 43% 64.4M 7s
    ##  58250K .......... .......... .......... .......... .......... 43%  107M 7s
    ##  58300K .......... .......... .......... .......... .......... 43% 98.5M 7s
    ##  58350K .......... .......... .......... .......... .......... 43% 66.2M 7s
    ##  58400K .......... .......... .......... .......... .......... 43% 67.6M 7s
    ##  58450K .......... .......... .......... .......... .......... 43% 37.2M 7s
    ##  58500K .......... .......... .......... .......... .......... 43% 41.5M 7s
    ##  58550K .......... .......... .......... .......... .......... 43% 88.7M 7s
    ##  58600K .......... .......... .......... .......... .......... 43%  125M 7s
    ##  58650K .......... .......... .......... .......... .......... 43%  116M 7s
    ##  58700K .......... .......... .......... .......... .......... 43% 13.7M 7s
    ##  58750K .......... .......... .......... .......... .......... 43% 99.2M 7s
    ##  58800K .......... .......... .......... .......... .......... 43% 95.5M 7s
    ##  58850K .......... .......... .......... .......... .......... 43% 94.4M 7s
    ##  58900K .......... .......... .......... .......... .......... 43% 99.2M 7s
    ##  58950K .......... .......... .......... .......... .......... 44%  104M 7s
    ##  59000K .......... .......... .......... .......... .......... 44% 16.8M 7s
    ##  59050K .......... .......... .......... .......... .......... 44% 34.3M 7s
    ##  59100K .......... .......... .......... .......... .......... 44% 53.6M 7s
    ##  59150K .......... .......... .......... .......... .......... 44% 73.5M 7s
    ##  59200K .......... .......... .......... .......... .......... 44%  117M 7s
    ##  59250K .......... .......... .......... .......... .......... 44%  128M 7s
    ##  59300K .......... .......... .......... .......... .......... 44% 13.1M 7s
    ##  59350K .......... .......... .......... .......... .......... 44% 33.7M 7s
    ##  59400K .......... .......... .......... .......... .......... 44% 98.9M 7s
    ##  59450K .......... .......... .......... .......... .......... 44% 36.4M 7s
    ##  59500K .......... .......... .......... .......... .......... 44%  100M 7s
    ##  59550K .......... .......... .......... .......... .......... 44% 84.1M 7s
    ##  59600K .......... .......... .......... .......... .......... 44% 45.1M 7s
    ##  59650K .......... .......... .......... .......... .......... 44% 16.1M 7s
    ##  59700K .......... .......... .......... .......... .......... 44% 70.0M 7s
    ##  59750K .......... .......... .......... .......... .......... 44%  110M 7s
    ##  59800K .......... .......... .......... .......... .......... 44%  124M 7s
    ##  59850K .......... .......... .......... .......... .......... 44%  132M 6s
    ##  59900K .......... .......... .......... .......... .......... 44%  124M 6s
    ##  59950K .......... .......... .......... .......... .......... 44% 23.0M 6s
    ##  60000K .......... .......... .......... .......... .......... 44% 48.7M 6s
    ##  60050K .......... .......... .......... .......... .......... 44% 95.3M 6s
    ##  60100K .......... .......... .......... .......... .......... 44%  101M 6s
    ##  60150K .......... .......... .......... .......... .......... 44%  103M 6s
    ##  60200K .......... .......... .......... .......... .......... 44%  105M 6s
    ##  60250K .......... .......... .......... .......... .......... 44% 52.7M 6s
    ##  60300K .......... .......... .......... .......... .......... 45% 61.8M 6s
    ##  60350K .......... .......... .......... .......... .......... 45% 40.9M 6s
    ##  60400K .......... .......... .......... .......... .......... 45% 72.2M 6s
    ##  60450K .......... .......... .......... .......... .......... 45%  101M 6s
    ##  60500K .......... .......... .......... .......... .......... 45% 68.5M 6s
    ##  60550K .......... .......... .......... .......... .......... 45% 41.2M 6s
    ##  60600K .......... .......... .......... .......... .......... 45% 37.7M 6s
    ##  60650K .......... .......... .......... .......... .......... 45%  101M 6s
    ##  60700K .......... .......... .......... .......... .......... 45%  133M 6s
    ##  60750K .......... .......... .......... .......... .......... 45% 50.8M 6s
    ##  60800K .......... .......... .......... .......... .......... 45% 18.2M 6s
    ##  60850K .......... .......... .......... .......... .......... 45%  116M 6s
    ##  60900K .......... .......... .......... .......... .......... 45%  107M 6s
    ##  60950K .......... .......... .......... .......... .......... 45% 98.2M 6s
    ##  61000K .......... .......... .......... .......... .......... 45% 86.7M 6s
    ##  61050K .......... .......... .......... .......... .......... 45% 18.6M 6s
    ##  61100K .......... .......... .......... .......... .......... 45% 79.1M 6s
    ##  61150K .......... .......... .......... .......... .......... 45% 54.2M 6s
    ##  61200K .......... .......... .......... .......... .......... 45%  109M 6s
    ##  61250K .......... .......... .......... .......... .......... 45% 89.6M 6s
    ##  61300K .......... .......... .......... .......... .......... 45% 30.1M 6s
    ##  61350K .......... .......... .......... .......... .......... 45% 56.3M 6s
    ##  61400K .......... .......... .......... .......... .......... 45% 72.9M 6s
    ##  61450K .......... .......... .......... .......... .......... 45% 83.1M 6s
    ##  61500K .......... .......... .......... .......... .......... 45% 90.5M 6s
    ##  61550K .......... .......... .......... .......... .......... 45% 73.5M 6s
    ##  61600K .......... .......... .......... .......... .......... 45% 75.1M 6s
    ##  61650K .......... .......... .......... .......... .......... 46% 73.7M 6s
    ##  61700K .......... .......... .......... .......... .......... 46% 85.0M 6s
    ##  61750K .......... .......... .......... .......... .......... 46% 79.4M 6s
    ##  61800K .......... .......... .......... .......... .......... 46% 72.9M 6s
    ##  61850K .......... .......... .......... .......... .......... 46% 72.3M 6s
    ##  61900K .......... .......... .......... .......... .......... 46% 32.1M 6s
    ##  61950K .......... .......... .......... .......... .......... 46% 41.7M 6s
    ##  62000K .......... .......... .......... .......... .......... 46% 84.2M 6s
    ##  62050K .......... .......... .......... .......... .......... 46% 92.0M 6s
    ##  62100K .......... .......... .......... .......... .......... 46% 92.1M 6s
    ##  62150K .......... .......... .......... .......... .......... 46% 64.4M 6s
    ##  62200K .......... .......... .......... .......... .......... 46% 87.4M 6s
    ##  62250K .......... .......... .......... .......... .......... 46% 94.7M 6s
    ##  62300K .......... .......... .......... .......... .......... 46% 96.3M 6s
    ##  62350K .......... .......... .......... .......... .......... 46% 83.3M 6s
    ##  62400K .......... .......... .......... .......... .......... 46% 29.8M 6s
    ##  62450K .......... .......... .......... .......... .......... 46% 92.8M 6s
    ##  62500K .......... .......... .......... .......... .......... 46% 76.0M 6s
    ##  62550K .......... .......... .......... .......... .......... 46% 68.5M 6s
    ##  62600K .......... .......... .......... .......... .......... 46% 96.2M 6s
    ##  62650K .......... .......... .......... .......... .......... 46% 90.5M 6s
    ##  62700K .......... .......... .......... .......... .......... 46% 20.0M 6s
    ##  62750K .......... .......... .......... .......... .......... 46% 71.2M 6s
    ##  62800K .......... .......... .......... .......... .......... 46% 98.4M 6s
    ##  62850K .......... .......... .......... .......... .......... 46% 95.9M 6s
    ##  62900K .......... .......... .......... .......... .......... 46%  102M 6s
    ##  62950K .......... .......... .......... .......... .......... 46% 15.1M 6s
    ##  63000K .......... .......... .......... .......... .......... 47% 39.6M 6s
    ##  63050K .......... .......... .......... .......... .......... 47% 77.3M 6s
    ##  63100K .......... .......... .......... .......... .......... 47%  107M 6s
    ##  63150K .......... .......... .......... .......... .......... 47% 86.8M 6s
    ##  63200K .......... .......... .......... .......... .......... 47% 53.8M 6s
    ##  63250K .......... .......... .......... .......... .......... 47% 70.1M 6s
    ##  63300K .......... .......... .......... .......... .......... 47% 83.2M 6s
    ##  63350K .......... .......... .......... .......... .......... 47% 56.0M 6s
    ##  63400K .......... .......... .......... .......... .......... 47% 76.2M 6s
    ##  63450K .......... .......... .......... .......... .......... 47%  108M 6s
    ##  63500K .......... .......... .......... .......... .......... 47% 32.1M 6s
    ##  63550K .......... .......... .......... .......... .......... 47% 89.0M 6s
    ##  63600K .......... .......... .......... .......... .......... 47%  106M 6s
    ##  63650K .......... .......... .......... .......... .......... 47% 87.0M 6s
    ##  63700K .......... .......... .......... .......... .......... 47%  103M 6s
    ##  63750K .......... .......... .......... .......... .......... 47% 33.1M 6s
    ##  63800K .......... .......... .......... .......... .......... 47% 89.2M 6s
    ##  63850K .......... .......... .......... .......... .......... 47% 87.8M 6s
    ##  63900K .......... .......... .......... .......... .......... 47% 82.2M 6s
    ##  63950K .......... .......... .......... .......... .......... 47% 93.4M 6s
    ##  64000K .......... .......... .......... .......... .......... 47% 64.3M 6s
    ##  64050K .......... .......... .......... .......... .......... 47% 23.9M 6s
    ##  64100K .......... .......... .......... .......... .......... 47% 83.6M 6s
    ##  64150K .......... .......... .......... .......... .......... 47% 91.7M 6s
    ##  64200K .......... .......... .......... .......... .......... 47% 51.9M 6s
    ##  64250K .......... .......... .......... .......... .......... 47%  102M 6s
    ##  64300K .......... .......... .......... .......... .......... 47% 86.4M 6s
    ##  64350K .......... .......... .......... .......... .......... 48% 53.2M 6s
    ##  64400K .......... .......... .......... .......... .......... 48% 80.0M 6s
    ##  64450K .......... .......... .......... .......... .......... 48% 41.6M 6s
    ##  64500K .......... .......... .......... .......... .......... 48% 90.0M 6s
    ##  64550K .......... .......... .......... .......... .......... 48% 87.6M 6s
    ##  64600K .......... .......... .......... .......... .......... 48% 81.3M 6s
    ##  64650K .......... .......... .......... .......... .......... 48%  116M 6s
    ##  64700K .......... .......... .......... .......... .......... 48% 38.7M 6s
    ##  64750K .......... .......... .......... .......... .......... 48% 57.3M 6s
    ##  64800K .......... .......... .......... .......... .......... 48% 89.0M 6s
    ##  64850K .......... .......... .......... .......... .......... 48% 89.7M 6s
    ##  64900K .......... .......... .......... .......... .......... 48% 96.4M 6s
    ##  64950K .......... .......... .......... .......... .......... 48% 45.3M 6s
    ##  65000K .......... .......... .......... .......... .......... 48% 64.1M 6s
    ##  65050K .......... .......... .......... .......... .......... 48% 23.9M 6s
    ##  65100K .......... .......... .......... .......... .......... 48% 92.4M 6s
    ##  65150K .......... .......... .......... .......... .......... 48% 76.1M 6s
    ##  65200K .......... .......... .......... .......... .......... 48%  115M 6s
    ##  65250K .......... .......... .......... .......... .......... 48% 87.8M 6s
    ##  65300K .......... .......... .......... .......... .......... 48% 47.4M 6s
    ##  65350K .......... .......... .......... .......... .......... 48% 68.6M 6s
    ##  65400K .......... .......... .......... .......... .......... 48% 76.0M 6s
    ##  65450K .......... .......... .......... .......... .......... 48% 25.4M 6s
    ##  65500K .......... .......... .......... .......... .......... 48%  101M 6s
    ##  65550K .......... .......... .......... .......... .......... 48% 20.4M 6s
    ##  65600K .......... .......... .......... .......... .......... 48%  103M 6s
    ##  65650K .......... .......... .......... .......... .......... 49% 97.2M 6s
    ##  65700K .......... .......... .......... .......... .......... 49%  107M 6s
    ##  65750K .......... .......... .......... .......... .......... 49% 95.2M 6s
    ##  65800K .......... .......... .......... .......... .......... 49% 17.4M 6s
    ##  65850K .......... .......... .......... .......... .......... 49% 45.5M 6s
    ##  65900K .......... .......... .......... .......... .......... 49% 77.9M 6s
    ##  65950K .......... .......... .......... .......... .......... 49% 66.8M 6s
    ##  66000K .......... .......... .......... .......... .......... 49%  100M 6s
    ##  66050K .......... .......... .......... .......... .......... 49%  117M 5s
    ##  66100K .......... .......... .......... .......... .......... 49% 74.1M 5s
    ##  66150K .......... .......... .......... .......... .......... 49% 68.2M 5s
    ##  66200K .......... .......... .......... .......... .......... 49% 62.8M 5s
    ##  66250K .......... .......... .......... .......... .......... 49% 38.1M 5s
    ##  66300K .......... .......... .......... .......... .......... 49% 78.7M 5s
    ##  66350K .......... .......... .......... .......... .......... 49% 72.6M 5s
    ##  66400K .......... .......... .......... .......... .......... 49% 93.7M 5s
    ##  66450K .......... .......... .......... .......... .......... 49% 73.1M 5s
    ##  66500K .......... .......... .......... .......... .......... 49% 13.7M 5s
    ##  66550K .......... .......... .......... .......... .......... 49% 65.1M 5s
    ##  66600K .......... .......... .......... .......... .......... 49% 76.7M 5s
    ##  66650K .......... .......... .......... .......... .......... 49% 90.6M 5s
    ##  66700K .......... .......... .......... .......... .......... 49%  112M 5s
    ##  66750K .......... .......... .......... .......... .......... 49% 43.2M 5s
    ##  66800K .......... .......... .......... .......... .......... 49% 65.0M 5s
    ##  66850K .......... .......... .......... .......... .......... 49%  105M 5s
    ##  66900K .......... .......... .......... .......... .......... 49% 33.0M 5s
    ##  66950K .......... .......... .......... .......... .......... 49% 95.5M 5s
    ##  67000K .......... .......... .......... .......... .......... 50%  112M 5s
    ##  67050K .......... .......... .......... .......... .......... 50%  115M 5s
    ##  67100K .......... .......... .......... .......... .......... 50%  112M 5s
    ##  67150K .......... .......... .......... .......... .......... 50% 12.7M 5s
    ##  67200K .......... .......... .......... .......... .......... 50% 82.1M 5s
    ##  67250K .......... .......... .......... .......... .......... 50% 92.6M 5s
    ##  67300K .......... .......... .......... .......... .......... 50%  112M 5s
    ##  67350K .......... .......... .......... .......... .......... 50%  102M 5s
    ##  67400K .......... .......... .......... .......... .......... 50% 13.6M 5s
    ##  67450K .......... .......... .......... .......... .......... 50% 83.9M 5s
    ##  67500K .......... .......... .......... .......... .......... 50% 92.7M 5s
    ##  67550K .......... .......... .......... .......... .......... 50% 80.4M 5s
    ##  67600K .......... .......... .......... .......... .......... 50%  116M 5s
    ##  67650K .......... .......... .......... .......... .......... 50%  115M 5s
    ##  67700K .......... .......... .......... .......... .......... 50% 76.8M 5s
    ##  67750K .......... .......... .......... .......... .......... 50% 66.0M 5s
    ##  67800K .......... .......... .......... .......... .......... 50% 98.1M 5s
    ##  67850K .......... .......... .......... .......... .......... 50% 99.9M 5s
    ##  67900K .......... .......... .......... .......... .......... 50%  116M 5s
    ##  67950K .......... .......... .......... .......... .......... 50% 93.7M 5s
    ##  68000K .......... .......... .......... .......... .......... 50% 45.9M 5s
    ##  68050K .......... .......... .......... .......... .......... 50% 48.6M 5s
    ##  68100K .......... .......... .......... .......... .......... 50% 35.8M 5s
    ##  68150K .......... .......... .......... .......... .......... 50% 72.4M 5s
    ##  68200K .......... .......... .......... .......... .......... 50%  118M 5s
    ##  68250K .......... .......... .......... .......... .......... 50% 50.8M 5s
    ##  68300K .......... .......... .......... .......... .......... 50%  114M 5s
    ##  68350K .......... .......... .......... .......... .......... 51% 29.7M 5s
    ##  68400K .......... .......... .......... .......... .......... 51% 99.2M 5s
    ##  68450K .......... .......... .......... .......... .......... 51% 92.0M 5s
    ##  68500K .......... .......... .......... .......... .......... 51%  102M 5s
    ##  68550K .......... .......... .......... .......... .......... 51% 92.5M 5s
    ##  68600K .......... .......... .......... .......... .......... 51%  101M 5s
    ##  68650K .......... .......... .......... .......... .......... 51% 6.56M 5s
    ##  68700K .......... .......... .......... .......... .......... 51%  107M 5s
    ##  68750K .......... .......... .......... .......... .......... 51% 44.5M 5s
    ##  68800K .......... .......... .......... .......... .......... 51% 91.8M 5s
    ##  68850K .......... .......... .......... .......... .......... 51%  116M 5s
    ##  68900K .......... .......... .......... .......... .......... 51%  122M 5s
    ##  68950K .......... .......... .......... .......... .......... 51% 21.8M 5s
    ##  69000K .......... .......... .......... .......... .......... 51%  118M 5s
    ##  69050K .......... .......... .......... .......... .......... 51% 18.4M 5s
    ##  69100K .......... .......... .......... .......... .......... 51%  103M 5s
    ##  69150K .......... .......... .......... .......... .......... 51%  105M 5s
    ##  69200K .......... .......... .......... .......... .......... 51% 8.45M 5s
    ##  69250K .......... .......... .......... .......... .......... 51%  104M 5s
    ##  69300K .......... .......... .......... .......... .......... 51% 70.2M 5s
    ##  69350K .......... .......... .......... .......... .......... 51% 78.0M 5s
    ##  69400K .......... .......... .......... .......... .......... 51%  121M 5s
    ##  69450K .......... .......... .......... .......... .......... 51% 21.6M 5s
    ##  69500K .......... .......... .......... .......... .......... 51% 94.7M 5s
    ##  69550K .......... .......... .......... .......... .......... 51% 5.99M 5s
    ##  69600K .......... .......... .......... .......... .......... 51% 93.1M 5s
    ##  69650K .......... .......... .......... .......... .......... 51%  122M 5s
    ##  69700K .......... .......... .......... .......... .......... 52% 23.9M 5s
    ##  69750K .......... .......... .......... .......... .......... 52% 99.5M 5s
    ##  69800K .......... .......... .......... .......... .......... 52%  120M 5s
    ##  69850K .......... .......... .......... .......... .......... 52% 27.6M 5s
    ##  69900K .......... .......... .......... .......... .......... 52%  119M 5s
    ##  69950K .......... .......... .......... .......... .......... 52% 8.68M 5s
    ##  70000K .......... .......... .......... .......... .......... 52% 95.2M 5s
    ##  70050K .......... .......... .......... .......... .......... 52% 93.6M 5s
    ##  70100K .......... .......... .......... .......... .......... 52%  127M 5s
    ##  70150K .......... .......... .......... .......... .......... 52%  112M 5s
    ##  70200K .......... .......... .......... .......... .......... 52%  128M 5s
    ##  70250K .......... .......... .......... .......... .......... 52% 2.58M 5s
    ##  70300K .......... .......... .......... .......... .......... 52% 98.5M 5s
    ##  70350K .......... .......... .......... .......... .......... 52% 91.2M 5s
    ##  70400K .......... .......... .......... .......... .......... 52%  103M 5s
    ##  70450K .......... .......... .......... .......... .......... 52%  132M 5s
    ##  70500K .......... .......... .......... .......... .......... 52%  133M 5s
    ##  70550K .......... .......... .......... .......... .......... 52% 16.2M 5s
    ##  70600K .......... .......... .......... .......... .......... 52% 4.51M 5s
    ##  70650K .......... .......... .......... .......... .......... 52%  112M 5s
    ##  70700K .......... .......... .......... .......... .......... 52% 37.4M 5s
    ##  70750K .......... .......... .......... .......... .......... 52% 47.2M 5s
    ##  70800K .......... .......... .......... .......... .......... 52%  116M 5s
    ##  70850K .......... .......... .......... .......... .......... 52%  129M 5s
    ##  70900K .......... .......... .......... .......... .......... 52% 37.2M 5s
    ##  70950K .......... .......... .......... .......... .......... 52% 28.3M 5s
    ##  71000K .......... .......... .......... .......... .......... 52% 32.5M 5s
    ##  71050K .......... .......... .......... .......... .......... 53% 16.5M 5s
    ##  71100K .......... .......... .......... .......... .......... 53% 19.9M 5s
    ##  71150K .......... .......... .......... .......... .......... 53% 15.7M 5s
    ##  71200K .......... .......... .......... .......... .......... 53% 20.1M 5s
    ##  71250K .......... .......... .......... .......... .......... 53% 20.5M 5s
    ##  71300K .......... .......... .......... .......... .......... 53% 21.4M 5s
    ##  71350K .......... .......... .......... .......... .......... 53% 18.2M 5s
    ##  71400K .......... .......... .......... .......... .......... 53% 20.8M 5s
    ##  71450K .......... .......... .......... .......... .......... 53% 21.2M 5s
    ##  71500K .......... .......... .......... .......... .......... 53% 21.3M 5s
    ##  71550K .......... .......... .......... .......... .......... 53% 18.4M 5s
    ##  71600K .......... .......... .......... .......... .......... 53% 22.8M 5s
    ##  71650K .......... .......... .......... .......... .......... 53% 22.7M 5s
    ##  71700K .......... .......... .......... .......... .......... 53% 22.4M 5s
    ##  71750K .......... .......... .......... .......... .......... 53% 20.8M 5s
    ##  71800K .......... .......... .......... .......... .......... 53% 24.5M 5s
    ##  71850K .......... .......... .......... .......... .......... 53% 23.4M 5s
    ##  71900K .......... .......... .......... .......... .......... 53% 22.6M 5s
    ##  71950K .......... .......... .......... .......... .......... 53% 20.4M 5s
    ##  72000K .......... .......... .......... .......... .......... 53% 22.1M 5s
    ##  72050K .......... .......... .......... .......... .......... 53% 24.6M 5s
    ##  72100K .......... .......... .......... .......... .......... 53% 25.2M 5s
    ##  72150K .......... .......... .......... .......... .......... 53% 22.9M 5s
    ##  72200K .......... .......... .......... .......... .......... 53% 24.9M 5s
    ##  72250K .......... .......... .......... .......... .......... 53% 27.0M 5s
    ##  72300K .......... .......... .......... .......... .......... 53% 27.5M 5s
    ##  72350K .......... .......... .......... .......... .......... 54% 23.9M 5s
    ##  72400K .......... .......... .......... .......... .......... 54% 27.0M 5s
    ##  72450K .......... .......... .......... .......... .......... 54% 28.0M 5s
    ##  72500K .......... .......... .......... .......... .......... 54% 27.9M 5s
    ##  72550K .......... .......... .......... .......... .......... 54% 25.2M 5s
    ##  72600K .......... .......... .......... .......... .......... 54% 28.5M 5s
    ##  72650K .......... .......... .......... .......... .......... 54% 29.6M 5s
    ##  72700K .......... .......... .......... .......... .......... 54% 29.2M 5s
    ##  72750K .......... .......... .......... .......... .......... 54% 23.9M 5s
    ##  72800K .......... .......... .......... .......... .......... 54% 29.0M 5s
    ##  72850K .......... .......... .......... .......... .......... 54% 28.1M 5s
    ##  72900K .......... .......... .......... .......... .......... 54% 30.0M 5s
    ##  72950K .......... .......... .......... .......... .......... 54% 27.2M 5s
    ##  73000K .......... .......... .......... .......... .......... 54% 29.6M 5s
    ##  73050K .......... .......... .......... .......... .......... 54% 30.9M 5s
    ##  73100K .......... .......... .......... .......... .......... 54% 30.3M 5s
    ##  73150K .......... .......... .......... .......... .......... 54% 27.3M 5s
    ##  73200K .......... .......... .......... .......... .......... 54% 32.3M 5s
    ##  73250K .......... .......... .......... .......... .......... 54% 32.4M 5s
    ##  73300K .......... .......... .......... .......... .......... 54% 31.2M 5s
    ##  73350K .......... .......... .......... .......... .......... 54% 27.3M 5s
    ##  73400K .......... .......... .......... .......... .......... 54% 30.6M 5s
    ##  73450K .......... .......... .......... .......... .......... 54% 30.2M 5s
    ##  73500K .......... .......... .......... .......... .......... 54% 33.9M 5s
    ##  73550K .......... .......... .......... .......... .......... 54% 28.4M 5s
    ##  73600K .......... .......... .......... .......... .......... 54% 31.9M 5s
    ##  73650K .......... .......... .......... .......... .......... 54% 34.1M 5s
    ##  73700K .......... .......... .......... .......... .......... 55% 31.6M 5s
    ##  73750K .......... .......... .......... .......... .......... 55% 29.4M 5s
    ##  73800K .......... .......... .......... .......... .......... 55% 33.5M 5s
    ##  73850K .......... .......... .......... .......... .......... 55% 35.4M 5s
    ##  73900K .......... .......... .......... .......... .......... 55% 33.4M 5s
    ##  73950K .......... .......... .......... .......... .......... 55% 30.1M 5s
    ##  74000K .......... .......... .......... .......... .......... 55% 34.5M 5s
    ##  74050K .......... .......... .......... .......... .......... 55% 34.6M 5s
    ##  74100K .......... .......... .......... .......... .......... 55% 35.0M 5s
    ##  74150K .......... .......... .......... .......... .......... 55% 29.3M 5s
    ##  74200K .......... .......... .......... .......... .......... 55% 33.8M 5s
    ##  74250K .......... .......... .......... .......... .......... 55% 37.2M 5s
    ##  74300K .......... .......... .......... .......... .......... 55% 34.7M 5s
    ##  74350K .......... .......... .......... .......... .......... 55% 30.5M 4s
    ##  74400K .......... .......... .......... .......... .......... 55% 36.5M 4s
    ##  74450K .......... .......... .......... .......... .......... 55% 34.7M 4s
    ##  74500K .......... .......... .......... .......... .......... 55% 37.3M 4s
    ##  74550K .......... .......... .......... .......... .......... 55% 31.2M 4s
    ##  74600K .......... .......... .......... .......... .......... 55% 37.1M 4s
    ##  74650K .......... .......... .......... .......... .......... 55% 37.8M 4s
    ##  74700K .......... .......... .......... .......... .......... 55% 36.5M 4s
    ##  74750K .......... .......... .......... .......... .......... 55% 31.2M 4s
    ##  74800K .......... .......... .......... .......... .......... 55% 37.5M 4s
    ##  74850K .......... .......... .......... .......... .......... 55% 36.2M 4s
    ##  74900K .......... .......... .......... .......... .......... 55% 36.2M 4s
    ##  74950K .......... .......... .......... .......... .......... 55% 31.8M 4s
    ##  75000K .......... .......... .......... .......... .......... 55% 36.7M 4s
    ##  75050K .......... .......... .......... .......... .......... 56% 37.2M 4s
    ##  75100K .......... .......... .......... .......... .......... 56% 33.7M 4s
    ##  75150K .......... .......... .......... .......... .......... 56% 29.6M 4s
    ##  75200K .......... .......... .......... .......... .......... 56% 36.5M 4s
    ##  75250K .......... .......... .......... .......... .......... 56% 35.7M 4s
    ##  75300K .......... .......... .......... .......... .......... 56% 35.8M 4s
    ##  75350K .......... .......... .......... .......... .......... 56% 32.3M 4s
    ##  75400K .......... .......... .......... .......... .......... 56% 36.8M 4s
    ##  75450K .......... .......... .......... .......... .......... 56% 37.4M 4s
    ##  75500K .......... .......... .......... .......... .......... 56% 37.3M 4s
    ##  75550K .......... .......... .......... .......... .......... 56% 31.6M 4s
    ##  75600K .......... .......... .......... .......... .......... 56% 36.1M 4s
    ##  75650K .......... .......... .......... .......... .......... 56% 36.9M 4s
    ##  75700K .......... .......... .......... .......... .......... 56% 35.4M 4s
    ##  75750K .......... .......... .......... .......... .......... 56% 31.1M 4s
    ##  75800K .......... .......... .......... .......... .......... 56% 35.7M 4s
    ##  75850K .......... .......... .......... .......... .......... 56% 35.6M 4s
    ##  75900K .......... .......... .......... .......... .......... 56% 36.3M 4s
    ##  75950K .......... .......... .......... .......... .......... 56% 29.8M 4s
    ##  76000K .......... .......... .......... .......... .......... 56% 35.2M 4s
    ##  76050K .......... .......... .......... .......... .......... 56% 35.9M 4s
    ##  76100K .......... .......... .......... .......... .......... 56% 37.0M 4s
    ##  76150K .......... .......... .......... .......... .......... 56% 31.9M 4s
    ##  76200K .......... .......... .......... .......... .......... 56% 36.8M 4s
    ##  76250K .......... .......... .......... .......... .......... 56% 37.3M 4s
    ##  76300K .......... .......... .......... .......... .......... 56% 36.7M 4s
    ##  76350K .......... .......... .......... .......... .......... 56% 32.2M 4s
    ##  76400K .......... .......... .......... .......... .......... 57% 35.6M 4s
    ##  76450K .......... .......... .......... .......... .......... 57% 35.5M 4s
    ##  76500K .......... .......... .......... .......... .......... 57% 37.3M 4s
    ##  76550K .......... .......... .......... .......... .......... 57% 31.1M 4s
    ##  76600K .......... .......... .......... .......... .......... 57% 36.8M 4s
    ##  76650K .......... .......... .......... .......... .......... 57% 37.3M 4s
    ##  76700K .......... .......... .......... .......... .......... 57% 35.4M 4s
    ##  76750K .......... .......... .......... .......... .......... 57% 31.9M 4s
    ##  76800K .......... .......... .......... .......... .......... 57% 37.4M 4s
    ##  76850K .......... .......... .......... .......... .......... 57% 36.5M 4s
    ##  76900K .......... .......... .......... .......... .......... 57% 36.0M 4s
    ##  76950K .......... .......... .......... .......... .......... 57% 32.2M 4s
    ##  77000K .......... .......... .......... .......... .......... 57% 36.7M 4s
    ##  77050K .......... .......... .......... .......... .......... 57% 37.2M 4s
    ##  77100K .......... .......... .......... .......... .......... 57% 35.5M 4s
    ##  77150K .......... .......... .......... .......... .......... 57% 31.6M 4s
    ##  77200K .......... .......... .......... .......... .......... 57% 37.2M 4s
    ##  77250K .......... .......... .......... .......... .......... 57% 35.6M 4s
    ##  77300K .......... .......... .......... .......... .......... 57% 35.9M 4s
    ##  77350K .......... .......... .......... .......... .......... 57% 32.8M 4s
    ##  77400K .......... .......... .......... .......... .......... 57% 36.7M 4s
    ##  77450K .......... .......... .......... .......... .......... 57% 33.6M 4s
    ##  77500K .......... .......... .......... .......... .......... 57% 37.1M 4s
    ##  77550K .......... .......... .......... .......... .......... 57% 32.0M 4s
    ##  77600K .......... .......... .......... .......... .......... 57% 34.7M 4s
    ##  77650K .......... .......... .......... .......... .......... 57% 37.5M 4s
    ##  77700K .......... .......... .......... .......... .......... 57% 36.4M 4s
    ##  77750K .......... .......... .......... .......... .......... 58% 31.1M 4s
    ##  77800K .......... .......... .......... .......... .......... 58% 37.0M 4s
    ##  77850K .......... .......... .......... .......... .......... 58% 36.6M 4s
    ##  77900K .......... .......... .......... .......... .......... 58% 3.92M 4s
    ##  77950K .......... .......... .......... .......... .......... 58% 32.2M 4s
    ##  78000K .......... .......... .......... .......... .......... 58% 36.5M 4s
    ##  78050K .......... .......... .......... .......... .......... 58% 35.9M 4s
    ##  78100K .......... .......... .......... .......... .......... 58% 37.7M 4s
    ##  78150K .......... .......... .......... .......... .......... 58% 32.1M 4s
    ##  78200K .......... .......... .......... .......... .......... 58% 37.0M 4s
    ##  78250K .......... .......... .......... .......... .......... 58% 35.6M 4s
    ##  78300K .......... .......... .......... .......... .......... 58% 33.5M 4s
    ##  78350K .......... .......... .......... .......... .......... 58% 31.9M 4s
    ##  78400K .......... .......... .......... .......... .......... 58% 37.3M 4s
    ##  78450K .......... .......... .......... .......... .......... 58% 36.5M 4s
    ##  78500K .......... .......... .......... .......... .......... 58% 34.8M 4s
    ##  78550K .......... .......... .......... .......... .......... 58% 32.6M 4s
    ##  78600K .......... .......... .......... .......... .......... 58% 36.2M 4s
    ##  78650K .......... .......... .......... .......... .......... 58% 36.9M 4s
    ##  78700K .......... .......... .......... .......... .......... 58% 37.5M 4s
    ##  78750K .......... .......... .......... .......... .......... 58% 29.9M 4s
    ##  78800K .......... .......... .......... .......... .......... 58% 37.3M 4s
    ##  78850K .......... .......... .......... .......... .......... 58% 37.4M 4s
    ##  78900K .......... .......... .......... .......... .......... 58% 36.2M 4s
    ##  78950K .......... .......... .......... .......... .......... 58% 30.2M 4s
    ##  79000K .......... .......... .......... .......... .......... 58% 35.8M 4s
    ##  79050K .......... .......... .......... .......... .......... 59% 36.8M 4s
    ##  79100K .......... .......... .......... .......... .......... 59% 37.1M 4s
    ##  79150K .......... .......... .......... .......... .......... 59% 32.2M 4s
    ##  79200K .......... .......... .......... .......... .......... 59% 35.6M 4s
    ##  79250K .......... .......... .......... .......... .......... 59% 34.1M 4s
    ##  79300K .......... .......... .......... .......... .......... 59% 36.0M 4s
    ##  79350K .......... .......... .......... .......... .......... 59% 32.7M 4s
    ##  79400K .......... .......... .......... .......... .......... 59% 34.9M 4s
    ##  79450K .......... .......... .......... .......... .......... 59% 34.9M 4s
    ##  79500K .......... .......... .......... .......... .......... 59% 36.0M 4s
    ##  79550K .......... .......... .......... .......... .......... 59% 31.2M 4s
    ##  79600K .......... .......... .......... .......... .......... 59% 36.3M 4s
    ##  79650K .......... .......... .......... .......... .......... 59% 37.7M 4s
    ##  79700K .......... .......... .......... .......... .......... 59% 36.1M 4s
    ##  79750K .......... .......... .......... .......... .......... 59% 30.1M 4s
    ##  79800K .......... .......... .......... .......... .......... 59% 36.4M 4s
    ##  79850K .......... .......... .......... .......... .......... 59% 35.8M 4s
    ##  79900K .......... .......... .......... .......... .......... 59% 34.3M 4s
    ##  79950K .......... .......... .......... .......... .......... 59% 31.3M 4s
    ##  80000K .......... .......... .......... .......... .......... 59% 35.3M 4s
    ##  80050K .......... .......... .......... .......... .......... 59% 32.1M 4s
    ##  80100K .......... .......... .......... .......... .......... 59% 33.8M 4s
    ##  80150K .......... .......... .......... .......... .......... 59% 31.8M 4s
    ##  80200K .......... .......... .......... .......... .......... 59% 33.1M 4s
    ##  80250K .......... .......... .......... .......... .......... 59% 32.6M 4s
    ##  80300K .......... .......... .......... .......... .......... 59% 35.1M 4s
    ##  80350K .......... .......... .......... .......... .......... 59% 32.5M 4s
    ##  80400K .......... .......... .......... .......... .......... 60% 36.7M 4s
    ##  80450K .......... .......... .......... .......... .......... 60% 37.0M 4s
    ##  80500K .......... .......... .......... .......... .......... 60% 37.0M 4s
    ##  80550K .......... .......... .......... .......... .......... 60% 32.5M 4s
    ##  80600K .......... .......... .......... .......... .......... 60% 35.2M 4s
    ##  80650K .......... .......... .......... .......... .......... 60% 35.8M 4s
    ##  80700K .......... .......... .......... .......... .......... 60% 37.4M 4s
    ##  80750K .......... .......... .......... .......... .......... 60% 31.0M 4s
    ##  80800K .......... .......... .......... .......... .......... 60% 36.6M 4s
    ##  80850K .......... .......... .......... .......... .......... 60% 36.0M 4s
    ##  80900K .......... .......... .......... .......... .......... 60% 37.4M 4s
    ##  80950K .......... .......... .......... .......... .......... 60% 32.3M 4s
    ##  81000K .......... .......... .......... .......... .......... 60% 36.7M 4s
    ##  81050K .......... .......... .......... .......... .......... 60% 36.8M 4s
    ##  81100K .......... .......... .......... .......... .......... 60% 35.4M 4s
    ##  81150K .......... .......... .......... .......... .......... 60% 30.7M 4s
    ##  81200K .......... .......... .......... .......... .......... 60% 37.1M 4s
    ##  81250K .......... .......... .......... .......... .......... 60% 37.2M 4s
    ##  81300K .......... .......... .......... .......... .......... 60% 35.9M 4s
    ##  81350K .......... .......... .......... .......... .......... 60% 32.8M 4s
    ##  81400K .......... .......... .......... .......... .......... 60% 37.5M 4s
    ##  81450K .......... .......... .......... .......... .......... 60% 34.2M 4s
    ##  81500K .......... .......... .......... .......... .......... 60% 37.4M 4s
    ##  81550K .......... .......... .......... .......... .......... 60% 31.9M 4s
    ##  81600K .......... .......... .......... .......... .......... 60% 34.6M 4s
    ##  81650K .......... .......... .......... .......... .......... 60% 37.8M 4s
    ##  81700K .......... .......... .......... .......... .......... 60% 36.6M 4s
    ##  81750K .......... .......... .......... .......... .......... 61% 30.6M 4s
    ##  81800K .......... .......... .......... .......... .......... 61% 36.9M 4s
    ##  81850K .......... .......... .......... .......... .......... 61% 34.4M 4s
    ##  81900K .......... .......... .......... .......... .......... 61% 34.2M 4s
    ##  81950K .......... .......... .......... .......... .......... 61% 31.2M 4s
    ##  82000K .......... .......... .......... .......... .......... 61% 35.1M 4s
    ##  82050K .......... .......... .......... .......... .......... 61% 36.8M 4s
    ##  82100K .......... .......... .......... .......... .......... 61% 36.8M 4s
    ##  82150K .......... .......... .......... .......... .......... 61% 32.0M 4s
    ##  82200K .......... .......... .......... .......... .......... 61% 36.9M 4s
    ##  82250K .......... .......... .......... .......... .......... 61% 37.5M 4s
    ##  82300K .......... .......... .......... .......... .......... 61% 35.5M 4s
    ##  82350K .......... .......... .......... .......... .......... 61% 31.6M 4s
    ##  82400K .......... .......... .......... .......... .......... 61% 37.6M 4s
    ##  82450K .......... .......... .......... .......... .......... 61% 33.1M 4s
    ##  82500K .......... .......... .......... .......... .......... 61% 23.1M 4s
    ##  82550K .......... .......... .......... .......... .......... 61% 23.1M 4s
    ##  82600K .......... .......... .......... .......... .......... 61% 37.5M 4s
    ##  82650K .......... .......... .......... .......... .......... 61% 19.9M 4s
    ##  82700K .......... .......... .......... .......... .......... 61% 25.4M 4s
    ##  82750K .......... .......... .......... .......... .......... 61% 17.9M 4s
    ##  82800K .......... .......... .......... .......... .......... 61% 29.0M 4s
    ##  82850K .......... .......... .......... .......... .......... 61% 24.2M 4s
    ##  82900K .......... .......... .......... .......... .......... 61% 36.8M 4s
    ##  82950K .......... .......... .......... .......... .......... 61% 20.4M 4s
    ##  83000K .......... .......... .......... .......... .......... 61% 37.0M 4s
    ##  83050K .......... .......... .......... .......... .......... 61% 18.9M 4s
    ##  83100K .......... .......... .......... .......... .......... 62% 25.4M 4s
    ##  83150K .......... .......... .......... .......... .......... 62% 20.2M 4s
    ##  83200K .......... .......... .......... .......... .......... 62% 25.9M 4s
    ##  83250K .......... .......... .......... .......... .......... 62% 25.3M 4s
    ##  83300K .......... .......... .......... .......... .......... 62% 28.9M 4s
    ##  83350K .......... .......... .......... .......... .......... 62% 23.8M 4s
    ##  83400K .......... .......... .......... .......... .......... 62% 25.4M 4s
    ##  83450K .......... .......... .......... .......... .......... 62% 29.9M 4s
    ##  83500K .......... .......... .......... .......... .......... 62% 25.4M 4s
    ##  83550K .......... .......... .......... .......... .......... 62% 21.5M 4s
    ##  83600K .......... .......... .......... .......... .......... 62% 23.5M 4s
    ##  83650K .......... .......... .......... .......... .......... 62% 29.5M 4s
    ##  83700K .......... .......... .......... .......... .......... 62% 28.9M 4s
    ##  83750K .......... .......... .......... .......... .......... 62% 22.1M 4s
    ##  83800K .......... .......... .......... .......... .......... 62% 25.0M 4s
    ##  83850K .......... .......... .......... .......... .......... 62% 28.0M 4s
    ##  83900K .......... .......... .......... .......... .......... 62% 28.8M 4s
    ##  83950K .......... .......... .......... .......... .......... 62% 22.4M 4s
    ##  84000K .......... .......... .......... .......... .......... 62% 27.9M 4s
    ##  84050K .......... .......... .......... .......... .......... 62% 24.3M 4s
    ##  84100K .......... .......... .......... .......... .......... 62% 28.8M 4s
    ##  84150K .......... .......... .......... .......... .......... 62% 25.9M 4s
    ##  84200K .......... .......... .......... .......... .......... 62% 20.9M 3s
    ##  84250K .......... .......... .......... .......... .......... 62% 27.7M 3s
    ##  84300K .......... .......... .......... .......... .......... 62% 29.1M 3s
    ##  84350K .......... .......... .......... .......... .......... 62% 22.6M 3s
    ##  84400K .......... .......... .......... .......... .......... 62% 25.1M 3s
    ##  84450K .......... .......... .......... .......... .......... 63% 36.5M 3s
    ##  84500K .......... .......... .......... .......... .......... 63% 22.0M 3s
    ##  84550K .......... .......... .......... .......... .......... 63% 32.1M 3s
    ##  84600K .......... .......... .......... .......... .......... 63% 21.9M 3s
    ##  84650K .......... .......... .......... .......... .......... 63% 26.1M 3s
    ##  84700K .......... .......... .......... .......... .......... 63% 25.0M 3s
    ##  84750K .......... .......... .......... .......... .......... 63% 26.9M 3s
    ##  84800K .......... .......... .......... .......... .......... 63% 21.2M 3s
    ##  84850K .......... .......... .......... .......... .......... 63% 36.8M 3s
    ##  84900K .......... .......... .......... .......... .......... 63% 23.3M 3s
    ##  84950K .......... .......... .......... .......... .......... 63% 23.4M 3s
    ##  85000K .......... .......... .......... .......... .......... 63% 28.0M 3s
    ##  85050K .......... .......... .......... .......... .......... 63% 25.1M 3s
    ##  85100K .......... .......... .......... .......... .......... 63% 22.8M 3s
    ##  85150K .......... .......... .......... .......... .......... 63% 23.5M 3s
    ##  85200K .......... .......... .......... .......... .......... 63% 32.0M 3s
    ##  85250K .......... .......... .......... .......... .......... 63% 24.9M 3s
    ##  85300K .......... .......... .......... .......... .......... 63% 37.5M 3s
    ##  85350K .......... .......... .......... .......... .......... 63% 19.6M 3s
    ##  85400K .......... .......... .......... .......... .......... 63% 37.4M 3s
    ##  85450K .......... .......... .......... .......... .......... 63% 21.4M 3s
    ##  85500K .......... .......... .......... .......... .......... 63% 24.2M 3s
    ##  85550K .......... .......... .......... .......... .......... 63% 32.0M 3s
    ##  85600K .......... .......... .......... .......... .......... 63% 21.5M 3s
    ##  85650K .......... .......... .......... .......... .......... 63% 36.7M 3s
    ##  85700K .......... .......... .......... .......... .......... 63% 29.7M 3s
    ##  85750K .......... .......... .......... .......... .......... 63% 23.2M 3s
    ##  85800K .......... .......... .......... .......... .......... 64% 19.6M 3s
    ##  85850K .......... .......... .......... .......... .......... 64% 25.7M 3s
    ##  85900K .......... .......... .......... .......... .......... 64% 23.2M 3s
    ##  85950K .......... .......... .......... .......... .......... 64% 20.9M 3s
    ##  86000K .......... .......... .......... .......... .......... 64% 29.4M 3s
    ##  86050K .......... .......... .......... .......... .......... 64% 13.8M 3s
    ##  86100K .......... .......... .......... .......... .......... 64% 36.2M 3s
    ##  86150K .......... .......... .......... .......... .......... 64% 21.6M 3s
    ##  86200K .......... .......... .......... .......... .......... 64% 24.7M 3s
    ##  86250K .......... .......... .......... .......... .......... 64% 28.8M 3s
    ##  86300K .......... .......... .......... .......... .......... 64% 24.2M 3s
    ##  86350K .......... .......... .......... .......... .......... 64% 30.6M 3s
    ##  86400K .......... .......... .......... .......... .......... 64% 21.3M 3s
    ##  86450K .......... .......... .......... .......... .......... 64% 36.6M 3s
    ##  86500K .......... .......... .......... .......... .......... 64% 21.9M 3s
    ##  86550K .......... .......... .......... .......... .......... 64% 32.3M 3s
    ##  86600K .......... .......... .......... .......... .......... 64% 35.5M 3s
    ##  86650K .......... .......... .......... .......... .......... 64% 37.0M 3s
    ##  86700K .......... .......... .......... .......... .......... 64% 20.9M 3s
    ##  86750K .......... .......... .......... .......... .......... 64% 32.1M 3s
    ##  86800K .......... .......... .......... .......... .......... 64% 37.3M 3s
    ##  86850K .......... .......... .......... .......... .......... 64% 18.6M 3s
    ##  86900K .......... .......... .......... .......... .......... 64% 23.9M 3s
    ##  86950K .......... .......... .......... .......... .......... 64% 20.2M 3s
    ##  87000K .......... .......... .......... .......... .......... 64% 36.4M 3s
    ##  87050K .......... .......... .......... .......... .......... 64% 27.3M 3s
    ##  87100K .......... .......... .......... .......... .......... 65% 25.1M 3s
    ##  87150K .......... .......... .......... .......... .......... 65% 17.9M 3s
    ##  87200K .......... .......... .......... .......... .......... 65% 25.6M 3s
    ##  87250K .......... .......... .......... .......... .......... 65% 23.1M 3s
    ##  87300K .......... .......... .......... .......... .......... 65% 23.9M 3s
    ##  87350K .......... .......... .......... .......... .......... 65% 23.5M 3s
    ##  87400K .......... .......... .......... .......... .......... 65% 30.9M 3s
    ##  87450K .......... .......... .......... .......... .......... 65% 28.1M 3s
    ##  87500K .......... .......... .......... .......... .......... 65% 30.8M 3s
    ##  87550K .......... .......... .......... .......... .......... 65% 21.6M 3s
    ##  87600K .......... .......... .......... .......... .......... 65% 25.5M 3s
    ##  87650K .......... .......... .......... .......... .......... 65% 29.6M 3s
    ##  87700K .......... .......... .......... .......... .......... 65% 21.5M 3s
    ##  87750K .......... .......... .......... .......... .......... 65% 22.2M 3s
    ##  87800K .......... .......... .......... .......... .......... 65% 37.2M 3s
    ##  87850K .......... .......... .......... .......... .......... 65% 23.4M 3s
    ##  87900K .......... .......... .......... .......... .......... 65% 37.5M 3s
    ##  87950K .......... .......... .......... .......... .......... 65% 19.4M 3s
    ##  88000K .......... .......... .......... .......... .......... 65% 21.8M 3s
    ##  88050K .......... .......... .......... .......... .......... 65% 30.1M 3s
    ##  88100K .......... .......... .......... .......... .......... 65% 27.6M 3s
    ##  88150K .......... .......... .......... .......... .......... 65% 22.8M 3s
    ##  88200K .......... .......... .......... .......... .......... 65% 27.7M 3s
    ##  88250K .......... .......... .......... .......... .......... 65% 21.7M 3s
    ##  88300K .......... .......... .......... .......... .......... 65% 35.8M 3s
    ##  88350K .......... .......... .......... .......... .......... 65% 24.7M 3s
    ##  88400K .......... .......... .......... .......... .......... 65% 37.6M 3s
    ##  88450K .......... .......... .......... .......... .......... 66% 29.5M 3s
    ##  88500K .......... .......... .......... .......... .......... 66% 24.9M 3s
    ##  88550K .......... .......... .......... .......... .......... 66% 23.3M 3s
    ##  88600K .......... .......... .......... .......... .......... 66% 22.0M 3s
    ##  88650K .......... .......... .......... .......... .......... 66%  135M 3s
    ##  88700K .......... .......... .......... .......... .......... 66%  140M 3s
    ##  88750K .......... .......... .......... .......... .......... 66%  119M 3s
    ##  88800K .......... .......... .......... .......... .......... 66%  120M 3s
    ##  88850K .......... .......... .......... .......... .......... 66% 35.0M 3s
    ##  88900K .......... .......... .......... .......... .......... 66%  116M 3s
    ##  88950K .......... .......... .......... .......... .......... 66%  119M 3s
    ##  89000K .......... .......... .......... .......... .......... 66% 43.8M 3s
    ##  89050K .......... .......... .......... .......... .......... 66%  142M 3s
    ##  89100K .......... .......... .......... .......... .......... 66%  128M 3s
    ##  89150K .......... .......... .......... .......... .......... 66% 47.1M 3s
    ##  89200K .......... .......... .......... .......... .......... 66%  132M 3s
    ##  89250K .......... .......... .......... .......... .......... 66% 50.4M 3s
    ##  89300K .......... .......... .......... .......... .......... 66% 51.9M 3s
    ##  89350K .......... .......... .......... .......... .......... 66% 22.7M 3s
    ##  89400K .......... .......... .......... .......... .......... 66% 49.6M 3s
    ##  89450K .......... .......... .......... .......... .......... 66%  131M 3s
    ##  89500K .......... .......... .......... .......... .......... 66%  143M 3s
    ##  89550K .......... .......... .......... .......... .......... 66% 31.3M 3s
    ##  89600K .......... .......... .......... .......... .......... 66% 50.9M 3s
    ##  89650K .......... .......... .......... .......... .......... 66%  141M 3s
    ##  89700K .......... .......... .......... .......... .......... 66%  137M 3s
    ##  89750K .......... .......... .......... .......... .......... 66% 30.5M 3s
    ##  89800K .......... .......... .......... .......... .......... 67% 50.8M 3s
    ##  89850K .......... .......... .......... .......... .......... 67%  130M 3s
    ##  89900K .......... .......... .......... .......... .......... 67% 50.4M 3s
    ##  89950K .......... .......... .......... .......... .......... 67% 47.3M 3s
    ##  90000K .......... .......... .......... .......... .......... 67%  134M 3s
    ##  90050K .......... .......... .......... .......... .......... 67%  142M 3s
    ##  90100K .......... .......... .......... .......... .......... 67% 30.7M 3s
    ##  90150K .......... .......... .......... .......... .......... 67% 47.3M 3s
    ##  90200K .......... .......... .......... .......... .......... 67%  139M 3s
    ##  90250K .......... .......... .......... .......... .......... 67% 33.1M 3s
    ##  90300K .......... .......... .......... .......... .......... 67% 52.9M 3s
    ##  90350K .......... .......... .......... .......... .......... 67%  101M 3s
    ##  90400K .......... .......... .......... .......... .......... 67%  141M 3s
    ##  90450K .......... .......... .......... .......... .......... 67%  144M 3s
    ##  90500K .......... .......... .......... .......... .......... 67%  134M 3s
    ##  90550K .......... .......... .......... .......... .......... 67%  123M 3s
    ##  90600K .......... .......... .......... .......... .......... 67%  139M 3s
    ##  90650K .......... .......... .......... .......... .......... 67%  135M 3s
    ##  90700K .......... .......... .......... .......... .......... 67%  128M 3s
    ##  90750K .......... .......... .......... .......... .......... 67%  118M 3s
    ##  90800K .......... .......... .......... .......... .......... 67%  134M 3s
    ##  90850K .......... .......... .......... .......... .......... 67%  144M 3s
    ##  90900K .......... .......... .......... .......... .......... 67%  137M 3s
    ##  90950K .......... .......... .......... .......... .......... 67% 15.9M 3s
    ##  91000K .......... .......... .......... .......... .......... 67% 73.5M 3s
    ##  91050K .......... .......... .......... .......... .......... 67%  137M 3s
    ##  91100K .......... .......... .......... .......... .......... 67%  140M 3s
    ##  91150K .......... .......... .......... .......... .......... 68%  117M 3s
    ##  91200K .......... .......... .......... .......... .......... 68%  146M 3s
    ##  91250K .......... .......... .......... .......... .......... 68%  141M 3s
    ##  91300K .......... .......... .......... .......... .......... 68% 16.0M 3s
    ##  91350K .......... .......... .......... .......... .......... 68% 68.5M 3s
    ##  91400K .......... .......... .......... .......... .......... 68% 76.7M 3s
    ##  91450K .......... .......... .......... .......... .......... 68% 81.0M 3s
    ##  91500K .......... .......... .......... .......... .......... 68% 79.7M 3s
    ##  91550K .......... .......... .......... .......... .......... 68% 57.2M 3s
    ##  91600K .......... .......... .......... .......... .......... 68% 73.0M 3s
    ##  91650K .......... .......... .......... .......... .......... 68% 61.6M 3s
    ##  91700K .......... .......... .......... .......... .......... 68% 65.9M 3s
    ##  91750K .......... .......... .......... .......... .......... 68% 68.0M 3s
    ##  91800K .......... .......... .......... .......... .......... 68% 74.2M 3s
    ##  91850K .......... .......... .......... .......... .......... 68% 66.5M 3s
    ##  91900K .......... .......... .......... .......... .......... 68% 71.2M 3s
    ##  91950K .......... .......... .......... .......... .......... 68% 54.9M 3s
    ##  92000K .......... .......... .......... .......... .......... 68% 77.7M 3s
    ##  92050K .......... .......... .......... .......... .......... 68% 62.8M 3s
    ##  92100K .......... .......... .......... .......... .......... 68% 76.0M 3s
    ##  92150K .......... .......... .......... .......... .......... 68% 75.7M 3s
    ##  92200K .......... .......... .......... .......... .......... 68% 82.8M 3s
    ##  92250K .......... .......... .......... .......... .......... 68% 84.8M 3s
    ##  92300K .......... .......... .......... .......... .......... 68% 84.1M 3s
    ##  92350K .......... .......... .......... .......... .......... 68% 59.9M 3s
    ##  92400K .......... .......... .......... .......... .......... 68% 64.2M 3s
    ##  92450K .......... .......... .......... .......... .......... 68% 85.3M 3s
    ##  92500K .......... .......... .......... .......... .......... 69% 77.2M 3s
    ##  92550K .......... .......... .......... .......... .......... 69% 74.1M 3s
    ##  92600K .......... .......... .......... .......... .......... 69% 86.1M 3s
    ##  92650K .......... .......... .......... .......... .......... 69% 84.2M 3s
    ##  92700K .......... .......... .......... .......... .......... 69% 75.2M 3s
    ##  92750K .......... .......... .......... .......... .......... 69% 64.9M 3s
    ##  92800K .......... .......... .......... .......... .......... 69% 83.6M 3s
    ##  92850K .......... .......... .......... .......... .......... 69% 85.6M 3s
    ##  92900K .......... .......... .......... .......... .......... 69% 85.9M 3s
    ##  92950K .......... .......... .......... .......... .......... 69% 82.4M 3s
    ##  93000K .......... .......... .......... .......... .......... 69% 88.4M 3s
    ##  93050K .......... .......... .......... .......... .......... 69% 80.2M 3s
    ##  93100K .......... .......... .......... .......... .......... 69% 95.9M 3s
    ##  93150K .......... .......... .......... .......... .......... 69% 74.9M 3s
    ##  93200K .......... .......... .......... .......... .......... 69% 93.7M 3s
    ##  93250K .......... .......... .......... .......... .......... 69% 91.0M 3s
    ##  93300K .......... .......... .......... .......... .......... 69% 93.7M 3s
    ##  93350K .......... .......... .......... .......... .......... 69% 81.9M 3s
    ##  93400K .......... .......... .......... .......... .......... 69% 88.0M 3s
    ##  93450K .......... .......... .......... .......... .......... 69% 89.7M 3s
    ##  93500K .......... .......... .......... .......... .......... 69% 90.5M 3s
    ##  93550K .......... .......... .......... .......... .......... 69% 74.3M 3s
    ##  93600K .......... .......... .......... .......... .......... 69% 81.4M 3s
    ##  93650K .......... .......... .......... .......... .......... 69% 91.7M 3s
    ##  93700K .......... .......... .......... .......... .......... 69% 91.3M 3s
    ##  93750K .......... .......... .......... .......... .......... 69% 73.5M 3s
    ##  93800K .......... .......... .......... .......... .......... 70% 94.6M 3s
    ##  93850K .......... .......... .......... .......... .......... 70% 97.4M 3s
    ##  93900K .......... .......... .......... .......... .......... 70% 96.5M 3s
    ##  93950K .......... .......... .......... .......... .......... 70% 81.8M 3s
    ##  94000K .......... .......... .......... .......... .......... 70% 98.8M 3s
    ##  94050K .......... .......... .......... .......... .......... 70% 96.5M 3s
    ##  94100K .......... .......... .......... .......... .......... 70% 88.3M 3s
    ##  94150K .......... .......... .......... .......... .......... 70% 69.4M 3s
    ##  94200K .......... .......... .......... .......... .......... 70% 79.0M 3s
    ##  94250K .......... .......... .......... .......... .......... 70% 94.8M 3s
    ##  94300K .......... .......... .......... .......... .......... 70% 96.5M 3s
    ##  94350K .......... .......... .......... .......... .......... 70% 84.3M 3s
    ##  94400K .......... .......... .......... .......... .......... 70% 95.3M 3s
    ##  94450K .......... .......... .......... .......... .......... 70% 96.3M 3s
    ##  94500K .......... .......... .......... .......... .......... 70% 96.0M 3s
    ##  94550K .......... .......... .......... .......... .......... 70% 86.6M 3s
    ##  94600K .......... .......... .......... .......... .......... 70% 96.7M 3s
    ##  94650K .......... .......... .......... .......... .......... 70% 88.3M 3s
    ##  94700K .......... .......... .......... .......... .......... 70% 92.0M 3s
    ##  94750K .......... .......... .......... .......... .......... 70% 87.1M 3s
    ##  94800K .......... .......... .......... .......... .......... 70%  103M 3s
    ##  94850K .......... .......... .......... .......... .......... 70% 91.7M 3s
    ##  94900K .......... .......... .......... .......... .......... 70% 84.3M 3s
    ##  94950K .......... .......... .......... .......... .......... 70% 78.0M 3s
    ##  95000K .......... .......... .......... .......... .......... 70%  101M 3s
    ##  95050K .......... .......... .......... .......... .......... 70% 91.8M 3s
    ##  95100K .......... .......... .......... .......... .......... 70% 80.7M 3s
    ##  95150K .......... .......... .......... .......... .......... 71% 83.8M 3s
    ##  95200K .......... .......... .......... .......... .......... 71%  101M 3s
    ##  95250K .......... .......... .......... .......... .......... 71%  103M 3s
    ##  95300K .......... .......... .......... .......... .......... 71%  104M 3s
    ##  95350K .......... .......... .......... .......... .......... 71% 80.4M 3s
    ##  95400K .......... .......... .......... .......... .......... 71% 86.4M 2s
    ##  95450K .......... .......... .......... .......... .......... 71% 79.3M 2s
    ##  95500K .......... .......... .......... .......... .......... 71% 93.1M 2s
    ##  95550K .......... .......... .......... .......... .......... 71% 7.11M 2s
    ##  95600K .......... .......... .......... .......... .......... 71% 82.7M 2s
    ##  95650K .......... .......... .......... .......... .......... 71%  102M 2s
    ##  95700K .......... .......... .......... .......... .......... 71% 86.1M 2s
    ##  95750K .......... .......... .......... .......... .......... 71% 74.1M 2s
    ##  95800K .......... .......... .......... .......... .......... 71% 96.1M 2s
    ##  95850K .......... .......... .......... .......... .......... 71% 99.3M 2s
    ##  95900K .......... .......... .......... .......... .......... 71% 34.3M 2s
    ##  95950K .......... .......... .......... .......... .......... 71% 72.7M 2s
    ##  96000K .......... .......... .......... .......... .......... 71% 99.7M 2s
    ##  96050K .......... .......... .......... .......... .......... 71% 54.5M 2s
    ##  96100K .......... .......... .......... .......... .......... 71% 65.8M 2s
    ##  96150K .......... .......... .......... .......... .......... 71% 73.2M 2s
    ##  96200K .......... .......... .......... .......... .......... 71% 20.3M 2s
    ##  96250K .......... .......... .......... .......... .......... 71% 64.7M 2s
    ##  96300K .......... .......... .......... .......... .......... 71% 86.8M 2s
    ##  96350K .......... .......... .......... .......... .......... 71% 65.7M 2s
    ##  96400K .......... .......... .......... .......... .......... 71% 83.2M 2s
    ##  96450K .......... .......... .......... .......... .......... 71% 75.0M 2s
    ##  96500K .......... .......... .......... .......... .......... 72% 83.6M 2s
    ##  96550K .......... .......... .......... .......... .......... 72% 68.8M 2s
    ##  96600K .......... .......... .......... .......... .......... 72% 86.1M 2s
    ##  96650K .......... .......... .......... .......... .......... 72% 78.9M 2s
    ##  96700K .......... .......... .......... .......... .......... 72% 88.3M 2s
    ##  96750K .......... .......... .......... .......... .......... 72% 68.1M 2s
    ##  96800K .......... .......... .......... .......... .......... 72% 78.6M 2s
    ##  96850K .......... .......... .......... .......... .......... 72% 85.3M 2s
    ##  96900K .......... .......... .......... .......... .......... 72% 80.1M 2s
    ##  96950K .......... .......... .......... .......... .......... 72% 69.4M 2s
    ##  97000K .......... .......... .......... .......... .......... 72% 91.1M 2s
    ##  97050K .......... .......... .......... .......... .......... 72% 89.8M 2s
    ##  97100K .......... .......... .......... .......... .......... 72% 86.5M 2s
    ##  97150K .......... .......... .......... .......... .......... 72% 67.1M 2s
    ##  97200K .......... .......... .......... .......... .......... 72% 90.1M 2s
    ##  97250K .......... .......... .......... .......... .......... 72% 76.4M 2s
    ##  97300K .......... .......... .......... .......... .......... 72% 88.0M 2s
    ##  97350K .......... .......... .......... .......... .......... 72% 78.8M 2s
    ##  97400K .......... .......... .......... .......... .......... 72% 78.4M 2s
    ##  97450K .......... .......... .......... .......... .......... 72% 67.8M 2s
    ##  97500K .......... .......... .......... .......... .......... 72% 81.7M 2s
    ##  97550K .......... .......... .......... .......... .......... 72% 77.1M 2s
    ##  97600K .......... .......... .......... .......... .......... 72% 92.1M 2s
    ##  97650K .......... .......... .......... .......... .......... 72% 9.74M 2s
    ##  97700K .......... .......... .......... .......... .......... 72% 53.1M 2s
    ##  97750K .......... .......... .......... .......... .......... 72% 67.0M 2s
    ##  97800K .......... .......... .......... .......... .......... 72% 65.3M 2s
    ##  97850K .......... .......... .......... .......... .......... 73% 76.7M 2s
    ##  97900K .......... .......... .......... .......... .......... 73% 96.7M 2s
    ##  97950K .......... .......... .......... .......... .......... 73% 82.1M 2s
    ##  98000K .......... .......... .......... .......... .......... 73% 91.2M 2s
    ##  98050K .......... .......... .......... .......... .......... 73% 89.7M 2s
    ##  98100K .......... .......... .......... .......... .......... 73% 32.8M 2s
    ##  98150K .......... .......... .......... .......... .......... 73% 48.3M 2s
    ##  98200K .......... .......... .......... .......... .......... 73% 58.9M 2s
    ##  98250K .......... .......... .......... .......... .......... 73% 64.1M 2s
    ##  98300K .......... .......... .......... .......... .......... 73% 77.0M 2s
    ##  98350K .......... .......... .......... .......... .......... 73% 53.5M 2s
    ##  98400K .......... .......... .......... .......... .......... 73% 76.4M 2s
    ##  98450K .......... .......... .......... .......... .......... 73% 78.5M 2s
    ##  98500K .......... .......... .......... .......... .......... 73% 70.8M 2s
    ##  98550K .......... .......... .......... .......... .......... 73% 8.79M 2s
    ##  98600K .......... .......... .......... .......... .......... 73% 49.4M 2s
    ##  98650K .......... .......... .......... .......... .......... 73% 63.4M 2s
    ##  98700K .......... .......... .......... .......... .......... 73% 69.9M 2s
    ##  98750K .......... .......... .......... .......... .......... 73% 62.5M 2s
    ##  98800K .......... .......... .......... .......... .......... 73% 64.4M 2s
    ##  98850K .......... .......... .......... .......... .......... 73% 66.0M 2s
    ##  98900K .......... .......... .......... .......... .......... 73% 62.8M 2s
    ##  98950K .......... .......... .......... .......... .......... 73% 62.4M 2s
    ##  99000K .......... .......... .......... .......... .......... 73% 81.3M 2s
    ##  99050K .......... .......... .......... .......... .......... 73% 82.1M 2s
    ##  99100K .......... .......... .......... .......... .......... 73% 75.7M 2s
    ##  99150K .......... .......... .......... .......... .......... 73% 60.4M 2s
    ##  99200K .......... .......... .......... .......... .......... 74% 83.3M 2s
    ##  99250K .......... .......... .......... .......... .......... 74% 90.6M 2s
    ##  99300K .......... .......... .......... .......... .......... 74% 6.01M 2s
    ##  99350K .......... .......... .......... .......... .......... 74% 55.7M 2s
    ##  99400K .......... .......... .......... .......... .......... 74% 71.1M 2s
    ##  99450K .......... .......... .......... .......... .......... 74% 66.0M 2s
    ##  99500K .......... .......... .......... .......... .......... 74% 71.9M 2s
    ##  99550K .......... .......... .......... .......... .......... 74% 70.7M 2s
    ##  99600K .......... .......... .......... .......... .......... 74% 72.3M 2s
    ##  99650K .......... .......... .......... .......... .......... 74% 73.0M 2s
    ##  99700K .......... .......... .......... .......... .......... 74% 84.9M 2s
    ##  99750K .......... .......... .......... .......... .......... 74% 26.6M 2s
    ##  99800K .......... .......... .......... .......... .......... 74% 85.3M 2s
    ##  99850K .......... .......... .......... .......... .......... 74% 81.1M 2s
    ##  99900K .......... .......... .......... .......... .......... 74% 80.6M 2s
    ##  99950K .......... .......... .......... .......... .......... 74% 57.4M 2s
    ## 100000K .......... .......... .......... .......... .......... 74% 67.6M 2s
    ## 100050K .......... .......... .......... .......... .......... 74% 70.9M 2s
    ## 100100K .......... .......... .......... .......... .......... 74% 65.8M 2s
    ## 100150K .......... .......... .......... .......... .......... 74% 68.8M 2s
    ## 100200K .......... .......... .......... .......... .......... 74% 84.9M 2s
    ## 100250K .......... .......... .......... .......... .......... 74% 71.9M 2s
    ## 100300K .......... .......... .......... .......... .......... 74% 73.2M 2s
    ## 100350K .......... .......... .......... .......... .......... 74% 65.6M 2s
    ## 100400K .......... .......... .......... .......... .......... 74% 86.4M 2s
    ## 100450K .......... .......... .......... .......... .......... 74% 69.5M 2s
    ## 100500K .......... .......... .......... .......... .......... 75% 76.2M 2s
    ## 100550K .......... .......... .......... .......... .......... 75% 71.8M 2s
    ## 100600K .......... .......... .......... .......... .......... 75% 84.2M 2s
    ## 100650K .......... .......... .......... .......... .......... 75% 79.7M 2s
    ## 100700K .......... .......... .......... .......... .......... 75% 77.4M 2s
    ## 100750K .......... .......... .......... .......... .......... 75% 62.8M 2s
    ## 100800K .......... .......... .......... .......... .......... 75% 51.9M 2s
    ## 100850K .......... .......... .......... .......... .......... 75% 61.0M 2s
    ## 100900K .......... .......... .......... .......... .......... 75% 84.0M 2s
    ## 100950K .......... .......... .......... .......... .......... 75% 72.3M 2s
    ## 101000K .......... .......... .......... .......... .......... 75% 92.2M 2s
    ## 101050K .......... .......... .......... .......... .......... 75% 84.4M 2s
    ## 101100K .......... .......... .......... .......... .......... 75% 85.0M 2s
    ## 101150K .......... .......... .......... .......... .......... 75% 77.8M 2s
    ## 101200K .......... .......... .......... .......... .......... 75% 96.7M 2s
    ## 101250K .......... .......... .......... .......... .......... 75% 94.4M 2s
    ## 101300K .......... .......... .......... .......... .......... 75% 92.1M 2s
    ## 101350K .......... .......... .......... .......... .......... 75% 85.4M 2s
    ## 101400K .......... .......... .......... .......... .......... 75% 97.7M 2s
    ## 101450K .......... .......... .......... .......... .......... 75% 93.8M 2s
    ## 101500K .......... .......... .......... .......... .......... 75% 93.1M 2s
    ## 101550K .......... .......... .......... .......... .......... 75% 84.6M 2s
    ## 101600K .......... .......... .......... .......... .......... 75% 24.5M 2s
    ## 101650K .......... .......... .......... .......... .......... 75% 89.7M 2s
    ## 101700K .......... .......... .......... .......... .......... 75% 80.1M 2s
    ## 101750K .......... .......... .......... .......... .......... 75% 83.2M 2s
    ## 101800K .......... .......... .......... .......... .......... 75%  110M 2s
    ## 101850K .......... .......... .......... .......... .......... 76%  103M 2s
    ## 101900K .......... .......... .......... .......... .......... 76% 38.7M 2s
    ## 101950K .......... .......... .......... .......... .......... 76% 47.8M 2s
    ## 102000K .......... .......... .......... .......... .......... 76% 63.3M 2s
    ## 102050K .......... .......... .......... .......... .......... 76% 64.8M 2s
    ## 102100K .......... .......... .......... .......... .......... 76%  104M 2s
    ## 102150K .......... .......... .......... .......... .......... 76% 75.8M 2s
    ## 102200K .......... .......... .......... .......... .......... 76% 31.2M 2s
    ## 102250K .......... .......... .......... .......... .......... 76% 77.5M 2s
    ## 102300K .......... .......... .......... .......... .......... 76% 55.2M 2s
    ## 102350K .......... .......... .......... .......... .......... 76% 60.0M 2s
    ## 102400K .......... .......... .......... .......... .......... 76% 73.4M 2s
    ## 102450K .......... .......... .......... .......... .......... 76% 36.1M 2s
    ## 102500K .......... .......... .......... .......... .......... 76%  104M 2s
    ## 102550K .......... .......... .......... .......... .......... 76% 78.6M 2s
    ## 102600K .......... .......... .......... .......... .......... 76% 91.4M 2s
    ## 102650K .......... .......... .......... .......... .......... 76% 98.1M 2s
    ## 102700K .......... .......... .......... .......... .......... 76% 58.0M 2s
    ## 102750K .......... .......... .......... .......... .......... 76% 76.8M 2s
    ## 102800K .......... .......... .......... .......... .......... 76% 86.2M 2s
    ## 102850K .......... .......... .......... .......... .......... 76%  101M 2s
    ## 102900K .......... .......... .......... .......... .......... 76% 63.1M 2s
    ## 102950K .......... .......... .......... .......... .......... 76% 71.3M 2s
    ## 103000K .......... .......... .......... .......... .......... 76% 99.6M 2s
    ## 103050K .......... .......... .......... .......... .......... 76% 84.3M 2s
    ## 103100K .......... .......... .......... .......... .......... 76% 76.5M 2s
    ## 103150K .......... .......... .......... .......... .......... 76% 81.3M 2s
    ## 103200K .......... .......... .......... .......... .......... 77% 89.0M 2s
    ## 103250K .......... .......... .......... .......... .......... 77% 50.5M 2s
    ## 103300K .......... .......... .......... .......... .......... 77% 73.2M 2s
    ## 103350K .......... .......... .......... .......... .......... 77% 85.5M 2s
    ## 103400K .......... .......... .......... .......... .......... 77%  103M 2s
    ## 103450K .......... .......... .......... .......... .......... 77% 54.6M 2s
    ## 103500K .......... .......... .......... .......... .......... 77% 85.9M 2s
    ## 103550K .......... .......... .......... .......... .......... 77% 64.5M 2s
    ## 103600K .......... .......... .......... .......... .......... 77% 80.2M 2s
    ## 103650K .......... .......... .......... .......... .......... 77% 72.6M 2s
    ## 103700K .......... .......... .......... .......... .......... 77% 91.9M 2s
    ## 103750K .......... .......... .......... .......... .......... 77% 34.2M 2s
    ## 103800K .......... .......... .......... .......... .......... 77% 86.5M 2s
    ## 103850K .......... .......... .......... .......... .......... 77%  113M 2s
    ## 103900K .......... .......... .......... .......... .......... 77% 89.6M 2s
    ## 103950K .......... .......... .......... .......... .......... 77% 73.5M 2s
    ## 104000K .......... .......... .......... .......... .......... 77% 73.2M 2s
    ## 104050K .......... .......... .......... .......... .......... 77% 41.8M 2s
    ## 104100K .......... .......... .......... .......... .......... 77%  106M 2s
    ## 104150K .......... .......... .......... .......... .......... 77%  101M 2s
    ## 104200K .......... .......... .......... .......... .......... 77% 33.9M 2s
    ## 104250K .......... .......... .......... .......... .......... 77%  105M 2s
    ## 104300K .......... .......... .......... .......... .......... 77% 90.3M 2s
    ## 104350K .......... .......... .......... .......... .......... 77% 98.8M 2s
    ## 104400K .......... .......... .......... .......... .......... 77% 79.1M 2s
    ## 104450K .......... .......... .......... .......... .......... 77%  117M 2s
    ## 104500K .......... .......... .......... .......... .......... 77% 21.9M 2s
    ## 104550K .......... .......... .......... .......... .......... 78% 78.1M 2s
    ## 104600K .......... .......... .......... .......... .......... 78%  128M 2s
    ## 104650K .......... .......... .......... .......... .......... 78%  115M 2s
    ## 104700K .......... .......... .......... .......... .......... 78% 91.4M 2s
    ## 104750K .......... .......... .......... .......... .......... 78% 88.9M 2s
    ## 104800K .......... .......... .......... .......... .......... 78% 39.2M 2s
    ## 104850K .......... .......... .......... .......... .......... 78% 97.2M 2s
    ## 104900K .......... .......... .......... .......... .......... 78% 76.6M 2s
    ## 104950K .......... .......... .......... .......... .......... 78% 82.8M 2s
    ## 105000K .......... .......... .......... .......... .......... 78% 60.0M 2s
    ## 105050K .......... .......... .......... .......... .......... 78% 95.3M 2s
    ## 105100K .......... .......... .......... .......... .......... 78% 57.6M 2s
    ## 105150K .......... .......... .......... .......... .......... 78% 48.7M 2s
    ## 105200K .......... .......... .......... .......... .......... 78% 60.3M 2s
    ## 105250K .......... .......... .......... .......... .......... 78% 64.3M 2s
    ## 105300K .......... .......... .......... .......... .......... 78% 70.1M 2s
    ## 105350K .......... .......... .......... .......... .......... 78% 89.9M 2s
    ## 105400K .......... .......... .......... .......... .......... 78% 64.9M 2s
    ## 105450K .......... .......... .......... .......... .......... 78% 58.4M 2s
    ## 105500K .......... .......... .......... .......... .......... 78% 70.0M 2s
    ## 105550K .......... .......... .......... .......... .......... 78% 87.7M 2s
    ## 105600K .......... .......... .......... .......... .......... 78% 45.2M 2s
    ## 105650K .......... .......... .......... .......... .......... 78%  100M 2s
    ## 105700K .......... .......... .......... .......... .......... 78%  101M 2s
    ## 105750K .......... .......... .......... .......... .......... 78% 80.9M 2s
    ## 105800K .......... .......... .......... .......... .......... 78% 61.4M 2s
    ## 105850K .......... .......... .......... .......... .......... 78% 87.7M 2s
    ## 105900K .......... .......... .......... .......... .......... 79% 54.4M 2s
    ## 105950K .......... .......... .......... .......... .......... 79% 65.0M 2s
    ## 106000K .......... .......... .......... .......... .......... 79%  117M 2s
    ## 106050K .......... .......... .......... .......... .......... 79% 79.9M 2s
    ## 106100K .......... .......... .......... .......... .......... 79% 58.8M 2s
    ## 106150K .......... .......... .......... .......... .......... 79% 62.2M 2s
    ## 106200K .......... .......... .......... .......... .......... 79%  106M 2s
    ## 106250K .......... .......... .......... .......... .......... 79% 98.7M 2s
    ## 106300K .......... .......... .......... .......... .......... 79% 12.9M 2s
    ## 106350K .......... .......... .......... .......... .......... 79% 50.1M 2s
    ## 106400K .......... .......... .......... .......... .......... 79% 53.2M 2s
    ## 106450K .......... .......... .......... .......... .......... 79% 97.5M 2s
    ## 106500K .......... .......... .......... .......... .......... 79% 96.7M 2s
    ## 106550K .......... .......... .......... .......... .......... 79% 87.7M 2s
    ## 106600K .......... .......... .......... .......... .......... 79%  137M 2s
    ## 106650K .......... .......... .......... .......... .......... 79% 68.5M 2s
    ## 106700K .......... .......... .......... .......... .......... 79% 54.1M 2s
    ## 106750K .......... .......... .......... .......... .......... 79% 23.7M 2s
    ## 106800K .......... .......... .......... .......... .......... 79%  118M 2s
    ## 106850K .......... .......... .......... .......... .......... 79%  131M 2s
    ## 106900K .......... .......... .......... .......... .......... 79% 15.2M 2s
    ## 106950K .......... .......... .......... .......... .......... 79% 61.3M 2s
    ## 107000K .......... .......... .......... .......... .......... 79% 89.9M 2s
    ## 107050K .......... .......... .......... .......... .......... 79% 44.2M 2s
    ## 107100K .......... .......... .......... .......... .......... 79% 95.0M 2s
    ## 107150K .......... .......... .......... .......... .......... 79% 95.0M 2s
    ## 107200K .......... .......... .......... .......... .......... 79%  102M 2s
    ## 107250K .......... .......... .......... .......... .......... 80%  135M 2s
    ## 107300K .......... .......... .......... .......... .......... 80% 77.6M 2s
    ## 107350K .......... .......... .......... .......... .......... 80%  118M 2s
    ## 107400K .......... .......... .......... .......... .......... 80% 30.9M 2s
    ## 107450K .......... .......... .......... .......... .......... 80%  112M 2s
    ## 107500K .......... .......... .......... .......... .......... 80%  120M 2s
    ## 107550K .......... .......... .......... .......... .......... 80%  104M 2s
    ## 107600K .......... .......... .......... .......... .......... 80%  103M 2s
    ## 107650K .......... .......... .......... .......... .......... 80%  119M 2s
    ## 107700K .......... .......... .......... .......... .......... 80% 26.7M 2s
    ## 107750K .......... .......... .......... .......... .......... 80% 85.0M 2s
    ## 107800K .......... .......... .......... .......... .......... 80% 37.7M 2s
    ## 107850K .......... .......... .......... .......... .......... 80%  131M 2s
    ## 107900K .......... .......... .......... .......... .......... 80% 35.0M 2s
    ## 107950K .......... .......... .......... .......... .......... 80%  110M 2s
    ## 108000K .......... .......... .......... .......... .......... 80%  106M 2s
    ## 108050K .......... .......... .......... .......... .......... 80%  128M 2s
    ## 108100K .......... .......... .......... .......... .......... 80% 41.9M 2s
    ## 108150K .......... .......... .......... .......... .......... 80%  105M 2s
    ## 108200K .......... .......... .......... .......... .......... 80%  112M 2s
    ## 108250K .......... .......... .......... .......... .......... 80%  102M 2s
    ## 108300K .......... .......... .......... .......... .......... 80% 47.6M 2s
    ## 108350K .......... .......... .......... .......... .......... 80% 86.1M 2s
    ## 108400K .......... .......... .......... .......... .......... 80% 36.5M 2s
    ## 108450K .......... .......... .......... .......... .......... 80%  115M 2s
    ## 108500K .......... .......... .......... .......... .......... 80%  138M 2s
    ## 108550K .......... .......... .......... .......... .......... 81% 22.3M 1s
    ## 108600K .......... .......... .......... .......... .......... 81%  128M 1s
    ## 108650K .......... .......... .......... .......... .......... 81%  120M 1s
    ## 108700K .......... .......... .......... .......... .......... 81%  139M 1s
    ## 108750K .......... .......... .......... .......... .......... 81%  112M 1s
    ## 108800K .......... .......... .......... .......... .......... 81% 27.7M 1s
    ## 108850K .......... .......... .......... .......... .......... 81% 35.7M 1s
    ## 108900K .......... .......... .......... .......... .......... 81%  116M 1s
    ## 108950K .......... .......... .......... .......... .......... 81% 48.4M 1s
    ## 109000K .......... .......... .......... .......... .......... 81% 96.4M 1s
    ## 109050K .......... .......... .......... .......... .......... 81%  102M 1s
    ## 109100K .......... .......... .......... .......... .......... 81%  116M 1s
    ## 109150K .......... .......... .......... .......... .......... 81% 34.8M 1s
    ## 109200K .......... .......... .......... .......... .......... 81%  125M 1s
    ## 109250K .......... .......... .......... .......... .......... 81% 43.9M 1s
    ## 109300K .......... .......... .......... .......... .......... 81%  112M 1s
    ## 109350K .......... .......... .......... .......... .......... 81%  128M 1s
    ## 109400K .......... .......... .......... .......... .......... 81%  145M 1s
    ## 109450K .......... .......... .......... .......... .......... 81% 61.1M 1s
    ## 109500K .......... .......... .......... .......... .......... 81%  128M 1s
    ## 109550K .......... .......... .......... .......... .......... 81% 23.0M 1s
    ## 109600K .......... .......... .......... .......... .......... 81%  145M 1s
    ## 109650K .......... .......... .......... .......... .......... 81%  138M 1s
    ## 109700K .......... .......... .......... .......... .......... 81%  154M 1s
    ## 109750K .......... .......... .......... .......... .......... 81%  112M 1s
    ## 109800K .......... .......... .......... .......... .......... 81%  137M 1s
    ## 109850K .......... .......... .......... .......... .......... 81%  134M 1s
    ## 109900K .......... .......... .......... .......... .......... 82% 21.3M 1s
    ## 109950K .......... .......... .......... .......... .......... 82%  114M 1s
    ## 110000K .......... .......... .......... .......... .......... 82%  124M 1s
    ## 110050K .......... .......... .......... .......... .......... 82%  136M 1s
    ## 110100K .......... .......... .......... .......... .......... 82%  142M 1s
    ## 110150K .......... .......... .......... .......... .......... 82%  114M 1s
    ## 110200K .......... .......... .......... .......... .......... 82% 63.4M 1s
    ## 110250K .......... .......... .......... .......... .......... 82% 31.6M 1s
    ## 110300K .......... .......... .......... .......... .......... 82%  114M 1s
    ## 110350K .......... .......... .......... .......... .......... 82% 85.4M 1s
    ## 110400K .......... .......... .......... .......... .......... 82%  122M 1s
    ## 110450K .......... .......... .......... .......... .......... 82%  109M 1s
    ## 110500K .......... .......... .......... .......... .......... 82%  139M 1s
    ## 110550K .......... .......... .......... .......... .......... 82%  124M 1s
    ## 110600K .......... .......... .......... .......... .......... 82% 36.1M 1s
    ## 110650K .......... .......... .......... .......... .......... 82%  115M 1s
    ## 110700K .......... .......... .......... .......... .......... 82%  112M 1s
    ## 110750K .......... .......... .......... .......... .......... 82% 95.6M 1s
    ## 110800K .......... .......... .......... .......... .......... 82%  122M 1s
    ## 110850K .......... .......... .......... .......... .......... 82%  139M 1s
    ## 110900K .......... .......... .......... .......... .......... 82% 23.7M 1s
    ## 110950K .......... .......... .......... .......... .......... 82%  119M 1s
    ## 111000K .......... .......... .......... .......... .......... 82% 88.5M 1s
    ## 111050K .......... .......... .......... .......... .......... 82% 34.5M 1s
    ## 111100K .......... .......... .......... .......... .......... 82%  115M 1s
    ## 111150K .......... .......... .......... .......... .......... 82%  111M 1s
    ## 111200K .......... .......... .......... .......... .......... 82%  134M 1s
    ## 111250K .......... .......... .......... .......... .......... 83%  120M 1s
    ## 111300K .......... .......... .......... .......... .......... 83%  113M 1s
    ## 111350K .......... .......... .......... .......... .......... 83% 95.8M 1s
    ## 111400K .......... .......... .......... .......... .......... 83%  142M 1s
    ## 111450K .......... .......... .......... .......... .......... 83%  135M 1s
    ## 111500K .......... .......... .......... .......... .......... 83% 38.5M 1s
    ## 111550K .......... .......... .......... .......... .......... 83% 94.1M 1s
    ## 111600K .......... .......... .......... .......... .......... 83%  131M 1s
    ## 111650K .......... .......... .......... .......... .......... 83% 37.7M 1s
    ## 111700K .......... .......... .......... .......... .......... 83% 94.8M 1s
    ## 111750K .......... .......... .......... .......... .......... 83% 95.7M 1s
    ## 111800K .......... .......... .......... .......... .......... 83%  139M 1s
    ## 111850K .......... .......... .......... .......... .......... 83% 95.3M 1s
    ## 111900K .......... .......... .......... .......... .......... 83%  115M 1s
    ## 111950K .......... .......... .......... .......... .......... 83% 19.1M 1s
    ## 112000K .......... .......... .......... .......... .......... 83% 75.0M 1s
    ## 112050K .......... .......... .......... .......... .......... 83%  102M 1s
    ## 112100K .......... .......... .......... .......... .......... 83%  113M 1s
    ## 112150K .......... .......... .......... .......... .......... 83%  117M 1s
    ## 112200K .......... .......... .......... .......... .......... 83%  139M 1s
    ## 112250K .......... .......... .......... .......... .......... 83%  134M 1s
    ## 112300K .......... .......... .......... .......... .......... 83% 28.1M 1s
    ## 112350K .......... .......... .......... .......... .......... 83% 53.8M 1s
    ## 112400K .......... .......... .......... .......... .......... 83% 52.7M 1s
    ## 112450K .......... .......... .......... .......... .......... 83%  116M 1s
    ## 112500K .......... .......... .......... .......... .......... 83%  114M 1s
    ## 112550K .......... .......... .......... .......... .......... 83%  125M 1s
    ## 112600K .......... .......... .......... .......... .......... 84% 95.0M 1s
    ## 112650K .......... .......... .......... .......... .......... 84% 66.7M 1s
    ## 112700K .......... .......... .......... .......... .......... 84%  129M 1s
    ## 112750K .......... .......... .......... .......... .......... 84% 16.9M 1s
    ## 112800K .......... .......... .......... .......... .......... 84% 41.2M 1s
    ## 112850K .......... .......... .......... .......... .......... 84% 94.8M 1s
    ## 112900K .......... .......... .......... .......... .......... 84% 99.4M 1s
    ## 112950K .......... .......... .......... .......... .......... 84% 81.1M 1s
    ## 113000K .......... .......... .......... .......... .......... 84%  109M 1s
    ## 113050K .......... .......... .......... .......... .......... 84% 97.8M 1s
    ## 113100K .......... .......... .......... .......... .......... 84% 87.6M 1s
    ## 113150K .......... .......... .......... .......... .......... 84% 46.7M 1s
    ## 113200K .......... .......... .......... .......... .......... 84%  109M 1s
    ## 113250K .......... .......... .......... .......... .......... 84%  120M 1s
    ## 113300K .......... .......... .......... .......... .......... 84%  119M 1s
    ## 113350K .......... .......... .......... .......... .......... 84%  113M 1s
    ## 113400K .......... .......... .......... .......... .......... 84%  102M 1s
    ## 113450K .......... .......... .......... .......... .......... 84%  112M 1s
    ## 113500K .......... .......... .......... .......... .......... 84% 75.6M 1s
    ## 113550K .......... .......... .......... .......... .......... 84% 98.5M 1s
    ## 113600K .......... .......... .......... .......... .......... 84%  109M 1s
    ## 113650K .......... .......... .......... .......... .......... 84% 27.5M 1s
    ## 113700K .......... .......... .......... .......... .......... 84%  129M 1s
    ## 113750K .......... .......... .......... .......... .......... 84% 97.5M 1s
    ## 113800K .......... .......... .......... .......... .......... 84%  120M 1s
    ## 113850K .......... .......... .......... .......... .......... 84%  107M 1s
    ## 113900K .......... .......... .......... .......... .......... 84%  127M 1s
    ## 113950K .......... .......... .......... .......... .......... 85%  103M 1s
    ## 114000K .......... .......... .......... .......... .......... 85% 52.2M 1s
    ## 114050K .......... .......... .......... .......... .......... 85%  103M 1s
    ## 114100K .......... .......... .......... .......... .......... 85%  119M 1s
    ## 114150K .......... .......... .......... .......... .......... 85% 65.0M 1s
    ## 114200K .......... .......... .......... .......... .......... 85%  102M 1s
    ## 114250K .......... .......... .......... .......... .......... 85%  102M 1s
    ## 114300K .......... .......... .......... .......... .......... 85%  130M 1s
    ## 114350K .......... .......... .......... .......... .......... 85% 60.0M 1s
    ## 114400K .......... .......... .......... .......... .......... 85%  120M 1s
    ## 114450K .......... .......... .......... .......... .......... 85%  122M 1s
    ## 114500K .......... .......... .......... .......... .......... 85% 69.9M 1s
    ## 114550K .......... .......... .......... .......... .......... 85% 5.17M 1s
    ## 114600K .......... .......... .......... .......... .......... 85% 86.9M 1s
    ## 114650K .......... .......... .......... .......... .......... 85% 59.6M 1s
    ## 114700K .......... .......... .......... .......... .......... 85%  103M 1s
    ## 114750K .......... .......... .......... .......... .......... 85% 99.0M 1s
    ## 114800K .......... .......... .......... .......... .......... 85%  130M 1s
    ## 114850K .......... .......... .......... .......... .......... 85%  128M 1s
    ## 114900K .......... .......... .......... .......... .......... 85% 21.0M 1s
    ## 114950K .......... .......... .......... .......... .......... 85% 97.3M 1s
    ## 115000K .......... .......... .......... .......... .......... 85%  120M 1s
    ## 115050K .......... .......... .......... .......... .......... 85%  115M 1s
    ## 115100K .......... .......... .......... .......... .......... 85%  106M 1s
    ## 115150K .......... .......... .......... .......... .......... 85%  104M 1s
    ## 115200K .......... .......... .......... .......... .......... 85% 41.0M 1s
    ## 115250K .......... .......... .......... .......... .......... 86%  106M 1s
    ## 115300K .......... .......... .......... .......... .......... 86% 26.8M 1s
    ## 115350K .......... .......... .......... .......... .......... 86% 38.3M 1s
    ## 115400K .......... .......... .......... .......... .......... 86% 58.9M 1s
    ## 115450K .......... .......... .......... .......... .......... 86%  105M 1s
    ## 115500K .......... .......... .......... .......... .......... 86%  109M 1s
    ## 115550K .......... .......... .......... .......... .......... 86%  101M 1s
    ## 115600K .......... .......... .......... .......... .......... 86% 26.9M 1s
    ## 115650K .......... .......... .......... .......... .......... 86%  105M 1s
    ## 115700K .......... .......... .......... .......... .......... 86%  105M 1s
    ## 115750K .......... .......... .......... .......... .......... 86% 38.8M 1s
    ## 115800K .......... .......... .......... .......... .......... 86%  130M 1s
    ## 115850K .......... .......... .......... .......... .......... 86%  128M 1s
    ## 115900K .......... .......... .......... .......... .......... 86% 76.5M 1s
    ## 115950K .......... .......... .......... .......... .......... 86% 97.7M 1s
    ## 116000K .......... .......... .......... .......... .......... 86% 80.9M 1s
    ## 116050K .......... .......... .......... .......... .......... 86%  104M 1s
    ## 116100K .......... .......... .......... .......... .......... 86% 12.5M 1s
    ## 116150K .......... .......... .......... .......... .......... 86%  101M 1s
    ## 116200K .......... .......... .......... .......... .......... 86%  106M 1s
    ## 116250K .......... .......... .......... .......... .......... 86%  109M 1s
    ## 116300K .......... .......... .......... .......... .......... 86%  125M 1s
    ## 116350K .......... .......... .......... .......... .......... 86% 81.6M 1s
    ## 116400K .......... .......... .......... .......... .......... 86%  123M 1s
    ## 116450K .......... .......... .......... .......... .......... 86%  134M 1s
    ## 116500K .......... .......... .......... .......... .......... 86% 59.6M 1s
    ## 116550K .......... .......... .......... .......... .......... 86% 97.8M 1s
    ## 116600K .......... .......... .......... .......... .......... 87% 27.1M 1s
    ## 116650K .......... .......... .......... .......... .......... 87%  124M 1s
    ## 116700K .......... .......... .......... .......... .......... 87%  131M 1s
    ## 116750K .......... .......... .......... .......... .......... 87% 34.1M 1s
    ## 116800K .......... .......... .......... .......... .......... 87%  103M 1s
    ## 116850K .......... .......... .......... .......... .......... 87%  106M 1s
    ## 116900K .......... .......... .......... .......... .......... 87%  132M 1s
    ## 116950K .......... .......... .......... .......... .......... 87%  102M 1s
    ## 117000K .......... .......... .......... .......... .......... 87%  126M 1s
    ## 117050K .......... .......... .......... .......... .......... 87%  135M 1s
    ## 117100K .......... .......... .......... .......... .......... 87% 48.8M 1s
    ## 117150K .......... .......... .......... .......... .......... 87% 48.9M 1s
    ## 117200K .......... .......... .......... .......... .......... 87% 45.2M 1s
    ## 117250K .......... .......... .......... .......... .......... 87%  102M 1s
    ## 117300K .......... .......... .......... .......... .......... 87%  105M 1s
    ## 117350K .......... .......... .......... .......... .......... 87% 83.8M 1s
    ## 117400K .......... .......... .......... .......... .......... 87% 35.7M 1s
    ## 117450K .......... .......... .......... .......... .......... 87% 98.7M 1s
    ## 117500K .......... .......... .......... .......... .......... 87%  133M 1s
    ## 117550K .......... .......... .......... .......... .......... 87% 41.0M 1s
    ## 117600K .......... .......... .......... .......... .......... 87% 87.1M 1s
    ## 117650K .......... .......... .......... .......... .......... 87% 49.3M 1s
    ## 117700K .......... .......... .......... .......... .......... 87%  108M 1s
    ## 117750K .......... .......... .......... .......... .......... 87% 28.3M 1s
    ## 117800K .......... .......... .......... .......... .......... 87% 57.0M 1s
    ## 117850K .......... .......... .......... .......... .......... 87% 77.0M 1s
    ## 117900K .......... .......... .......... .......... .......... 87% 96.9M 1s
    ## 117950K .......... .......... .......... .......... .......... 88% 56.1M 1s
    ## 118000K .......... .......... .......... .......... .......... 88% 58.8M 1s
    ## 118050K .......... .......... .......... .......... .......... 88% 57.5M 1s
    ## 118100K .......... .......... .......... .......... .......... 88% 95.4M 1s
    ## 118150K .......... .......... .......... .......... .......... 88% 47.4M 1s
    ## 118200K .......... .......... .......... .......... .......... 88% 51.4M 1s
    ## 118250K .......... .......... .......... .......... .......... 88% 94.9M 1s
    ## 118300K .......... .......... .......... .......... .......... 88%  114M 1s
    ## 118350K .......... .......... .......... .......... .......... 88% 51.8M 1s
    ## 118400K .......... .......... .......... .......... .......... 88%  120M 1s
    ## 118450K .......... .......... .......... .......... .......... 88% 51.0M 1s
    ## 118500K .......... .......... .......... .......... .......... 88%  115M 1s
    ## 118550K .......... .......... .......... .......... .......... 88% 79.7M 1s
    ## 118600K .......... .......... .......... .......... .......... 88%  107M 1s
    ## 118650K .......... .......... .......... .......... .......... 88%  121M 1s
    ## 118700K .......... .......... .......... .......... .......... 88%  108M 1s
    ## 118750K .......... .......... .......... .......... .......... 88% 31.7M 1s
    ## 118800K .......... .......... .......... .......... .......... 88%  102M 1s
    ## 118850K .......... .......... .......... .......... .......... 88% 69.5M 1s
    ## 118900K .......... .......... .......... .......... .......... 88%  109M 1s
    ## 118950K .......... .......... .......... .......... .......... 88% 44.7M 1s
    ## 119000K .......... .......... .......... .......... .......... 88%  130M 1s
    ## 119050K .......... .......... .......... .......... .......... 88%  126M 1s
    ## 119100K .......... .......... .......... .......... .......... 88%  116M 1s
    ## 119150K .......... .......... .......... .......... .......... 88% 33.5M 1s
    ## 119200K .......... .......... .......... .......... .......... 88% 80.8M 1s
    ## 119250K .......... .......... .......... .......... .......... 88%  107M 1s
    ## 119300K .......... .......... .......... .......... .......... 89% 24.9M 1s
    ## 119350K .......... .......... .......... .......... .......... 89% 96.0M 1s
    ## 119400K .......... .......... .......... .......... .......... 89%  125M 1s
    ## 119450K .......... .......... .......... .......... .......... 89%  129M 1s
    ## 119500K .......... .......... .......... .......... .......... 89%  114M 1s
    ## 119550K .......... .......... .......... .......... .......... 89%  110M 1s
    ## 119600K .......... .......... .......... .......... .......... 89%  131M 1s
    ## 119650K .......... .......... .......... .......... .......... 89% 23.4M 1s
    ## 119700K .......... .......... .......... .......... .......... 89%  127M 1s
    ## 119750K .......... .......... .......... .......... .......... 89% 92.1M 1s
    ## 119800K .......... .......... .......... .......... .......... 89% 66.7M 1s
    ## 119850K .......... .......... .......... .......... .......... 89% 73.5M 1s
    ## 119900K .......... .......... .......... .......... .......... 89% 97.7M 1s
    ## 119950K .......... .......... .......... .......... .......... 89% 82.8M 1s
    ## 120000K .......... .......... .......... .......... .......... 89%  132M 1s
    ## 120050K .......... .......... .......... .......... .......... 89%  110M 1s
    ## 120100K .......... .......... .......... .......... .......... 89%  122M 1s
    ## 120150K .......... .......... .......... .......... .......... 89% 33.3M 1s
    ## 120200K .......... .......... .......... .......... .......... 89% 79.6M 1s
    ## 120250K .......... .......... .......... .......... .......... 89% 52.6M 1s
    ## 120300K .......... .......... .......... .......... .......... 89% 63.9M 1s
    ## 120350K .......... .......... .......... .......... .......... 89% 84.9M 1s
    ## 120400K .......... .......... .......... .......... .......... 89% 29.2M 1s
    ## 120450K .......... .......... .......... .......... .......... 89% 95.0M 1s
    ## 120500K .......... .......... .......... .......... .......... 89%  127M 1s
    ## 120550K .......... .......... .......... .......... .......... 89%  111M 1s
    ## 120600K .......... .......... .......... .......... .......... 89% 55.6M 1s
    ## 120650K .......... .......... .......... .......... .......... 90%  109M 1s
    ## 120700K .......... .......... .......... .......... .......... 90%  111M 1s
    ## 120750K .......... .......... .......... .......... .......... 90% 95.6M 1s
    ## 120800K .......... .......... .......... .......... .......... 90%  120M 1s
    ## 120850K .......... .......... .......... .......... .......... 90%  110M 1s
    ## 120900K .......... .......... .......... .......... .......... 90%  104M 1s
    ## 120950K .......... .......... .......... .......... .......... 90% 96.9M 1s
    ## 121000K .......... .......... .......... .......... .......... 90% 62.0M 1s
    ## 121050K .......... .......... .......... .......... .......... 90% 57.8M 1s
    ## 121100K .......... .......... .......... .......... .......... 90%  101M 1s
    ## 121150K .......... .......... .......... .......... .......... 90% 99.2M 1s
    ## 121200K .......... .......... .......... .......... .......... 90%  104M 1s
    ## 121250K .......... .......... .......... .......... .......... 90% 76.6M 1s
    ## 121300K .......... .......... .......... .......... .......... 90% 88.3M 1s
    ## 121350K .......... .......... .......... .......... .......... 90% 79.7M 1s
    ## 121400K .......... .......... .......... .......... .......... 90%  108M 1s
    ## 121450K .......... .......... .......... .......... .......... 90% 68.9M 1s
    ## 121500K .......... .......... .......... .......... .......... 90% 62.7M 1s
    ## 121550K .......... .......... .......... .......... .......... 90% 63.9M 1s
    ## 121600K .......... .......... .......... .......... .......... 90%  103M 1s
    ## 121650K .......... .......... .......... .......... .......... 90%  103M 1s
    ## 121700K .......... .......... .......... .......... .......... 90% 55.7M 1s
    ## 121750K .......... .......... .......... .......... .......... 90% 87.2M 1s
    ## 121800K .......... .......... .......... .......... .......... 90%  117M 1s
    ## 121850K .......... .......... .......... .......... .......... 90%  105M 1s
    ## 121900K .......... .......... .......... .......... .......... 90% 65.4M 1s
    ## 121950K .......... .......... .......... .......... .......... 91% 86.0M 1s
    ## 122000K .......... .......... .......... .......... .......... 91% 89.6M 1s
    ## 122050K .......... .......... .......... .......... .......... 91%  107M 1s
    ## 122100K .......... .......... .......... .......... .......... 91%  105M 1s
    ## 122150K .......... .......... .......... .......... .......... 91% 89.4M 1s
    ## 122200K .......... .......... .......... .......... .......... 91%  104M 1s
    ## 122250K .......... .......... .......... .......... .......... 91%  118M 1s
    ## 122300K .......... .......... .......... .......... .......... 91%  103M 1s
    ## 122350K .......... .......... .......... .......... .......... 91% 90.8M 1s
    ## 122400K .......... .......... .......... .......... .......... 91% 64.3M 1s
    ## 122450K .......... .......... .......... .......... .......... 91% 60.8M 1s
    ## 122500K .......... .......... .......... .......... .......... 91%  107M 1s
    ## 122550K .......... .......... .......... .......... .......... 91% 97.5M 1s
    ## 122600K .......... .......... .......... .......... .......... 91%  108M 1s
    ## 122650K .......... .......... .......... .......... .......... 91% 55.0M 1s
    ## 122700K .......... .......... .......... .......... .......... 91% 81.2M 1s
    ## 122750K .......... .......... .......... .......... .......... 91% 33.6M 1s
    ## 122800K .......... .......... .......... .......... .......... 91% 77.5M 1s
    ## 122850K .......... .......... .......... .......... .......... 91% 89.0M 1s
    ## 122900K .......... .......... .......... .......... .......... 91% 89.6M 1s
    ## 122950K .......... .......... .......... .......... .......... 91% 65.4M 1s
    ## 123000K .......... .......... .......... .......... .......... 91% 77.5M 1s
    ## 123050K .......... .......... .......... .......... .......... 91% 88.6M 1s
    ## 123100K .......... .......... .......... .......... .......... 91% 83.0M 1s
    ## 123150K .......... .......... .......... .......... .......... 91% 46.7M 1s
    ## 123200K .......... .......... .......... .......... .......... 91% 93.8M 1s
    ## 123250K .......... .......... .......... .......... .......... 91% 92.6M 1s
    ## 123300K .......... .......... .......... .......... .......... 92% 95.0M 1s
    ## 123350K .......... .......... .......... .......... .......... 92% 17.0M 1s
    ## 123400K .......... .......... .......... .......... .......... 92% 94.5M 1s
    ## 123450K .......... .......... .......... .......... .......... 92%  103M 1s
    ## 123500K .......... .......... .......... .......... .......... 92% 88.4M 1s
    ## 123550K .......... .......... .......... .......... .......... 92% 76.4M 1s
    ## 123600K .......... .......... .......... .......... .......... 92% 93.6M 1s
    ## 123650K .......... .......... .......... .......... .......... 92%  104M 1s
    ## 123700K .......... .......... .......... .......... .......... 92% 11.1M 1s
    ## 123750K .......... .......... .......... .......... .......... 92% 71.2M 1s
    ## 123800K .......... .......... .......... .......... .......... 92% 41.3M 1s
    ## 123850K .......... .......... .......... .......... .......... 92% 80.2M 1s
    ## 123900K .......... .......... .......... .......... .......... 92% 80.7M 1s
    ## 123950K .......... .......... .......... .......... .......... 92% 57.4M 1s
    ## 124000K .......... .......... .......... .......... .......... 92% 99.3M 1s
    ## 124050K .......... .......... .......... .......... .......... 92% 20.6M 1s
    ## 124100K .......... .......... .......... .......... .......... 92%  103M 1s
    ## 124150K .......... .......... .......... .......... .......... 92% 66.7M 1s
    ## 124200K .......... .......... .......... .......... .......... 92% 46.8M 1s
    ## 124250K .......... .......... .......... .......... .......... 92% 67.7M 1s
    ## 124300K .......... .......... .......... .......... .......... 92%  105M 1s
    ## 124350K .......... .......... .......... .......... .......... 92% 90.8M 1s
    ## 124400K .......... .......... .......... .......... .......... 92%  107M 1s
    ## 124450K .......... .......... .......... .......... .......... 92% 10.9M 1s
    ## 124500K .......... .......... .......... .......... .......... 92%  101M 1s
    ## 124550K .......... .......... .......... .......... .......... 92% 75.9M 1s
    ## 124600K .......... .......... .......... .......... .......... 92% 87.3M 1s
    ## 124650K .......... .......... .......... .......... .......... 93% 91.7M 0s
    ## 124700K .......... .......... .......... .......... .......... 93% 83.4M 0s
    ## 124750K .......... .......... .......... .......... .......... 93% 81.9M 0s
    ## 124800K .......... .......... .......... .......... .......... 93%  114M 0s
    ## 124850K .......... .......... .......... .......... .......... 93% 37.8M 0s
    ## 124900K .......... .......... .......... .......... .......... 93% 48.7M 0s
    ## 124950K .......... .......... .......... .......... .......... 93% 76.2M 0s
    ## 125000K .......... .......... .......... .......... .......... 93% 60.9M 0s
    ## 125050K .......... .......... .......... .......... .......... 93% 65.6M 0s
    ## 125100K .......... .......... .......... .......... .......... 93% 63.0M 0s
    ## 125150K .......... .......... .......... .......... .......... 93% 88.6M 0s
    ## 125200K .......... .......... .......... .......... .......... 93%  118M 0s
    ## 125250K .......... .......... .......... .......... .......... 93%  110M 0s
    ## 125300K .......... .......... .......... .......... .......... 93%  100M 0s
    ## 125350K .......... .......... .......... .......... .......... 93% 78.0M 0s
    ## 125400K .......... .......... .......... .......... .......... 93% 41.7M 0s
    ## 125450K .......... .......... .......... .......... .......... 93% 69.4M 0s
    ## 125500K .......... .......... .......... .......... .......... 93% 99.0M 0s
    ## 125550K .......... .......... .......... .......... .......... 93% 73.9M 0s
    ## 125600K .......... .......... .......... .......... .......... 93% 70.9M 0s
    ## 125650K .......... .......... .......... .......... .......... 93% 53.9M 0s
    ## 125700K .......... .......... .......... .......... .......... 93% 71.9M 0s
    ## 125750K .......... .......... .......... .......... .......... 93% 32.5M 0s
    ## 125800K .......... .......... .......... .......... .......... 93%  106M 0s
    ## 125850K .......... .......... .......... .......... .......... 93%  122M 0s
    ## 125900K .......... .......... .......... .......... .......... 93% 96.5M 0s
    ## 125950K .......... .......... .......... .......... .......... 93% 90.7M 0s
    ## 126000K .......... .......... .......... .......... .......... 94%  115M 0s
    ## 126050K .......... .......... .......... .......... .......... 94%  112M 0s
    ## 126100K .......... .......... .......... .......... .......... 94% 97.9M 0s
    ## 126150K .......... .......... .......... .......... .......... 94% 92.9M 0s
    ## 126200K .......... .......... .......... .......... .......... 94% 50.7M 0s
    ## 126250K .......... .......... .......... .......... .......... 94% 18.0M 0s
    ## 126300K .......... .......... .......... .......... .......... 94% 31.2M 0s
    ## 126350K .......... .......... .......... .......... .......... 94% 99.0M 0s
    ## 126400K .......... .......... .......... .......... .......... 94%  121M 0s
    ## 126450K .......... .......... .......... .......... .......... 94%  113M 0s
    ## 126500K .......... .......... .......... .......... .......... 94%  117M 0s
    ## 126550K .......... .......... .......... .......... .......... 94%  110M 0s
    ## 126600K .......... .......... .......... .......... .......... 94% 95.9M 0s
    ## 126650K .......... .......... .......... .......... .......... 94%  122M 0s
    ## 126700K .......... .......... .......... .......... .......... 94%  126M 0s
    ## 126750K .......... .......... .......... .......... .......... 94% 60.8M 0s
    ## 126800K .......... .......... .......... .......... .......... 94% 54.9M 0s
    ## 126850K .......... .......... .......... .......... .......... 94% 75.8M 0s
    ## 126900K .......... .......... .......... .......... .......... 94% 88.4M 0s
    ## 126950K .......... .......... .......... .......... .......... 94% 67.4M 0s
    ## 127000K .......... .......... .......... .......... .......... 94% 49.8M 0s
    ## 127050K .......... .......... .......... .......... .......... 94% 37.4M 0s
    ## 127100K .......... .......... .......... .......... .......... 94%  122M 0s
    ## 127150K .......... .......... .......... .......... .......... 94% 79.3M 0s
    ## 127200K .......... .......... .......... .......... .......... 94% 6.35M 0s
    ## 127250K .......... .......... .......... .......... .......... 94% 61.8M 0s
    ## 127300K .......... .......... .......... .......... .......... 94% 50.5M 0s
    ## 127350K .......... .......... .......... .......... .......... 95% 52.1M 0s
    ## 127400K .......... .......... .......... .......... .......... 95% 75.5M 0s
    ## 127450K .......... .......... .......... .......... .......... 95% 68.5M 0s
    ## 127500K .......... .......... .......... .......... .......... 95% 90.4M 0s
    ## 127550K .......... .......... .......... .......... .......... 95% 87.6M 0s
    ## 127600K .......... .......... .......... .......... .......... 95%  106M 0s
    ## 127650K .......... .......... .......... .......... .......... 95%  105M 0s
    ## 127700K .......... .......... .......... .......... .......... 95% 97.1M 0s
    ## 127750K .......... .......... .......... .......... .......... 95%  105M 0s
    ## 127800K .......... .......... .......... .......... .......... 95% 89.4M 0s
    ## 127850K .......... .......... .......... .......... .......... 95%  103M 0s
    ## 127900K .......... .......... .......... .......... .......... 95% 89.5M 0s
    ## 127950K .......... .......... .......... .......... .......... 95% 60.8M 0s
    ## 128000K .......... .......... .......... .......... .......... 95% 88.3M 0s
    ## 128050K .......... .......... .......... .......... .......... 95% 59.1M 0s
    ## 128100K .......... .......... .......... .......... .......... 95% 93.2M 0s
    ## 128150K .......... .......... .......... .......... .......... 95% 95.7M 0s
    ## 128200K .......... .......... .......... .......... .......... 95% 75.6M 0s
    ## 128250K .......... .......... .......... .......... .......... 95% 68.4M 0s
    ## 128300K .......... .......... .......... .......... .......... 95% 79.1M 0s
    ## 128350K .......... .......... .......... .......... .......... 95% 67.0M 0s
    ## 128400K .......... .......... .......... .......... .......... 95% 91.4M 0s
    ## 128450K .......... .......... .......... .......... .......... 95% 97.0M 0s
    ## 128500K .......... .......... .......... .......... .......... 95%  110M 0s
    ## 128550K .......... .......... .......... .......... .......... 95%  102M 0s
    ## 128600K .......... .......... .......... .......... .......... 95% 69.3M 0s
    ## 128650K .......... .......... .......... .......... .......... 95%  120M 0s
    ## 128700K .......... .......... .......... .......... .......... 96% 85.8M 0s
    ## 128750K .......... .......... .......... .......... .......... 96% 70.3M 0s
    ## 128800K .......... .......... .......... .......... .......... 96% 85.7M 0s
    ## 128850K .......... .......... .......... .......... .......... 96% 98.2M 0s
    ## 128900K .......... .......... .......... .......... .......... 96% 87.4M 0s
    ## 128950K .......... .......... .......... .......... .......... 96% 74.5M 0s
    ## 129000K .......... .......... .......... .......... .......... 96% 83.7M 0s
    ## 129050K .......... .......... .......... .......... .......... 96% 76.5M 0s
    ## 129100K .......... .......... .......... .......... .......... 96% 95.4M 0s
    ## 129150K .......... .......... .......... .......... .......... 96% 32.1M 0s
    ## 129200K .......... .......... .......... .......... .......... 96% 57.2M 0s
    ## 129250K .......... .......... .......... .......... .......... 96% 86.0M 0s
    ## 129300K .......... .......... .......... .......... .......... 96% 72.2M 0s
    ## 129350K .......... .......... .......... .......... .......... 96% 66.1M 0s
    ## 129400K .......... .......... .......... .......... .......... 96% 75.5M 0s
    ## 129450K .......... .......... .......... .......... .......... 96% 81.4M 0s
    ## 129500K .......... .......... .......... .......... .......... 96% 89.9M 0s
    ## 129550K .......... .......... .......... .......... .......... 96% 6.75M 0s
    ## 129600K .......... .......... .......... .......... .......... 96% 83.2M 0s
    ## 129650K .......... .......... .......... .......... .......... 96% 38.6M 0s
    ## 129700K .......... .......... .......... .......... .......... 96% 57.4M 0s
    ## 129750K .......... .......... .......... .......... .......... 96% 49.2M 0s
    ## 129800K .......... .......... .......... .......... .......... 96% 79.5M 0s
    ## 129850K .......... .......... .......... .......... .......... 96%  122M 0s
    ## 129900K .......... .......... .......... .......... .......... 96%  107M 0s
    ## 129950K .......... .......... .......... .......... .......... 96% 77.7M 0s
    ## 130000K .......... .......... .......... .......... .......... 97%  133M 0s
    ## 130050K .......... .......... .......... .......... .......... 97% 26.8M 0s
    ## 130100K .......... .......... .......... .......... .......... 97% 40.9M 0s
    ## 130150K .......... .......... .......... .......... .......... 97% 40.0M 0s
    ## 130200K .......... .......... .......... .......... .......... 97%  103M 0s
    ## 130250K .......... .......... .......... .......... .......... 97% 67.3M 0s
    ## 130300K .......... .......... .......... .......... .......... 97% 82.9M 0s
    ## 130350K .......... .......... .......... .......... .......... 97% 82.3M 0s
    ## 130400K .......... .......... .......... .......... .......... 97% 95.4M 0s
    ## 130450K .......... .......... .......... .......... .......... 97% 9.23M 0s
    ## 130500K .......... .......... .......... .......... .......... 97% 78.1M 0s
    ## 130550K .......... .......... .......... .......... .......... 97% 49.7M 0s
    ## 130600K .......... .......... .......... .......... .......... 97% 73.4M 0s
    ## 130650K .......... .......... .......... .......... .......... 97% 74.8M 0s
    ## 130700K .......... .......... .......... .......... .......... 97% 80.2M 0s
    ## 130750K .......... .......... .......... .......... .......... 97% 75.0M 0s
    ## 130800K .......... .......... .......... .......... .......... 97% 38.1M 0s
    ## 130850K .......... .......... .......... .......... .......... 97%  128M 0s
    ## 130900K .......... .......... .......... .......... .......... 97%  121M 0s
    ## 130950K .......... .......... .......... .......... .......... 97%  119M 0s
    ## 131000K .......... .......... .......... .......... .......... 97% 27.9M 0s
    ## 131050K .......... .......... .......... .......... .......... 97% 82.0M 0s
    ## 131100K .......... .......... .......... .......... .......... 97% 85.3M 0s
    ## 131150K .......... .......... .......... .......... .......... 97% 67.9M 0s
    ## 131200K .......... .......... .......... .......... .......... 97% 24.5M 0s
    ## 131250K .......... .......... .......... .......... .......... 97% 17.6M 0s
    ## 131300K .......... .......... .......... .......... .......... 97%  103M 0s
    ## 131350K .......... .......... .......... .......... .......... 98%  107M 0s
    ## 131400K .......... .......... .......... .......... .......... 98%  133M 0s
    ## 131450K .......... .......... .......... .......... .......... 98%  123M 0s
    ## 131500K .......... .......... .......... .......... .......... 98%  118M 0s
    ## 131550K .......... .......... .......... .......... .......... 98%  113M 0s
    ## 131600K .......... .......... .......... .......... .......... 98% 10.7M 0s
    ## 131650K .......... .......... .......... .......... .......... 98% 44.2M 0s
    ## 131700K .......... .......... .......... .......... .......... 98% 57.1M 0s
    ## 131750K .......... .......... .......... .......... .......... 98% 64.0M 0s
    ## 131800K .......... .......... .......... .......... .......... 98% 21.3M 0s
    ## 131850K .......... .......... .......... .......... .......... 98%  127M 0s
    ## 131900K .......... .......... .......... .......... .......... 98%  129M 0s
    ## 131950K .......... .......... .......... .......... .......... 98%  107M 0s
    ## 132000K .......... .......... .......... .......... .......... 98% 16.7M 0s
    ## 132050K .......... .......... .......... .......... .......... 98%  118M 0s
    ## 132100K .......... .......... .......... .......... .......... 98%  132M 0s
    ## 132150K .......... .......... .......... .......... .......... 98% 83.2M 0s
    ## 132200K .......... .......... .......... .......... .......... 98% 90.3M 0s
    ## 132250K .......... .......... .......... .......... .......... 98% 96.9M 0s
    ## 132300K .......... .......... .......... .......... .......... 98% 88.8M 0s
    ## 132350K .......... .......... .......... .......... .......... 98% 25.0M 0s
    ## 132400K .......... .......... .......... .......... .......... 98% 38.6M 0s
    ## 132450K .......... .......... .......... .......... .......... 98% 39.4M 0s
    ## 132500K .......... .......... .......... .......... .......... 98% 80.3M 0s
    ## 132550K .......... .......... .......... .......... .......... 98% 79.2M 0s
    ## 132600K .......... .......... .......... .......... .......... 98% 99.4M 0s
    ## 132650K .......... .......... .......... .......... .......... 98% 94.9M 0s
    ## 132700K .......... .......... .......... .......... .......... 99%  107M 0s
    ## 132750K .......... .......... .......... .......... .......... 99%  109M 0s
    ## 132800K .......... .......... .......... .......... .......... 99% 46.7M 0s
    ## 132850K .......... .......... .......... .......... .......... 99% 79.2M 0s
    ## 132900K .......... .......... .......... .......... .......... 99% 44.9M 0s
    ## 132950K .......... .......... .......... .......... .......... 99% 40.8M 0s
    ## 133000K .......... .......... .......... .......... .......... 99%  104M 0s
    ## 133050K .......... .......... .......... .......... .......... 99% 68.0M 0s
    ## 133100K .......... .......... .......... .......... .......... 99% 74.5M 0s
    ## 133150K .......... .......... .......... .......... .......... 99% 65.4M 0s
    ## 133200K .......... .......... .......... .......... .......... 99% 88.0M 0s
    ## 133250K .......... .......... .......... .......... .......... 99% 6.05M 0s
    ## 133300K .......... .......... .......... .......... .......... 99% 39.1M 0s
    ## 133350K .......... .......... .......... .......... .......... 99% 59.9M 0s
    ## 133400K .......... .......... .......... .......... .......... 99% 84.1M 0s
    ## 133450K .......... .......... .......... .......... .......... 99% 95.5M 0s
    ## 133500K .......... .......... .......... .......... .......... 99% 73.1M 0s
    ## 133550K .......... .......... .......... .......... .......... 99% 66.4M 0s
    ## 133600K .......... .......... .......... .......... .......... 99% 88.7M 0s
    ## 133650K .......... .......... .......... .......... .......... 99% 16.9M 0s
    ## 133700K .......... .......... .......... .......... .......... 99% 48.9M 0s
    ## 133750K .......... .......... .......... .......... .......... 99% 62.6M 0s
    ## 133800K .......... .......... .......... .......... .......... 99% 73.4M 0s
    ## 133850K .......... .......... .......... .......... .......... 99% 68.5M 0s
    ## 133900K .......... .......... .......... .......... .......... 99% 80.2M 0s
    ## 133950K .......... .......... .......... .......... .......... 99% 76.5M 0s
    ## 134000K .......... .......... .......... .......... .......... 99% 86.0M 0s
    ## 134050K .......... .....                                      100% 88.1M=6.8s
    ## 
    ## 2021-12-15 08:50:07 (19.3 MB/s) - ‘silva_nr99_v138.1_train_set.fa.gz.1’ saved [137283333/137283333]

``` r
fastaRef <- "/home/rstudio/silva_nr99_v138.1_train_set.fa.gz"
taxTab <- assignTaxonomy(seqtabNoC, refFasta=fastaRef, multithread=TRUE)
unname(head(taxTab))
```

    ##      [,1]       [,2]           [,3]          [,4]            [,5]            
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      [,6]         
    ## [1,] NA           
    ## [2,] NA           
    ## [3,] NA           
    ## [4,] NA           
    ## [5,] "Bacteroides"
    ## [6,] NA

``` r
seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
```

``` r
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
```

``` r
samdf <- read.csv("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/MIMARKS_Data_combined.csv",header=TRUE)
samdf$SampleID <- paste0(gsub("00", "", samdf$host_subject_id), "D", samdf$age-21)
samdf <- samdf[!duplicated(samdf$SampleID),] # Remove dupicate entries for reverse reads
rownames(seqtabAll) <- gsub("124", "125", rownames(seqtabAll)) # Fix discrepancy
all(rownames(seqtabAll) %in% samdf$SampleID) # TRUE
```

    ## [1] TRUE

``` r
rownames(samdf) <- samdf$SampleID
keep.cols <- c("collection_date", "biome", "target_gene", "target_subfragment",
"host_common_name", "host_subject_id", "age", "sex", "body_product", "tot_mass",
"diet", "family_relationship", "genotype", "SampleID") 
samdf <- samdf[rownames(seqtabAll), keep.cols]
```

``` r
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxTab),phy_tree(fitGTR$tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 218 taxa and 19 samples ]
    ## sample_data() Sample Data:       [ 19 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 218 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 218 tips and 216 internal nodes ]

``` r
ps_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/ps.rds")
ps = readRDS(ps_connect)
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 389 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 389 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 389 tips and 387 internal nodes ]

``` r
# Show available ranks in the dataset
rank_names(ps)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"

``` r
# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)
```

    ## 
    ##              Actinobacteria               Bacteroidetes 
    ##                          13                          23 
    ## Candidatus_Saccharibacteria   Cyanobacteria/Chloroplast 
    ##                           1                           4 
    ##         Deinococcus-Thermus                  Firmicutes 
    ##                           1                         327 
    ##                Fusobacteria              Proteobacteria 
    ##                           1                          11 
    ##                 Tenericutes             Verrucomicrobia 
    ##                           1                           1 
    ##                        <NA> 
    ##                           6

``` r
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
```

``` r
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
```

``` r
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

    ##                         Phylum         1     2
    ## 1               Actinobacteria 120.15385  1562
    ## 2                Bacteroidetes 265.52174  6107
    ## 3  Candidatus_Saccharibacteria 280.00000   280
    ## 4    Cyanobacteria/Chloroplast  64.25000   257
    ## 5          Deinococcus-Thermus  52.00000    52
    ## 6                   Firmicutes 179.24771 58614
    ## 7                 Fusobacteria   2.00000     2
    ## 8               Proteobacteria  59.09091   650
    ## 9                  Tenericutes 234.00000   234
    ## 10             Verrucomicrobia 104.00000   104

``` r
# Define phyla to filter
filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)
ps1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 381 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 381 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 381 tips and 379 internal nodes ]

``` r
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](cc1_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

``` r
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
```

    ## [1] 18

``` r
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
```

``` r
# Show available ranks in the dataset
rank_names(ps)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"

``` r
# Create table, number of features for each phyla
```

``` r
table(tax_table(ps)[, "Phylum"], exclude = NULL)
```

    ## 
    ##              Actinobacteria               Bacteroidetes 
    ##                          13                          23 
    ## Candidatus_Saccharibacteria   Cyanobacteria/Chloroplast 
    ##                           1                           4 
    ##         Deinococcus-Thermus                  Firmicutes 
    ##                           1                         327 
    ##                Fusobacteria              Proteobacteria 
    ##                           1                          11 
    ##                 Tenericutes             Verrucomicrobia 
    ##                           1                           1

``` r
## 
##              Actinobacteria               Bacteroidetes Candidatus_Saccharibacteria 
##                          13                          23                           1 
##   Cyanobacteria/Chloroplast         Deinococcus-Thermus                  Firmicutes 
##                           4                           1                         327 
##                Fusobacteria              Proteobacteria                 Tenericutes 
##                           1                          11                           1 
##             Verrucomicrobia                        <NA> 
##                           1                           6
```

``` r
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
```

``` r
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
```

``` r
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

    ##                         Phylum         1     2
    ## 1               Actinobacteria 120.15385  1562
    ## 2                Bacteroidetes 265.52174  6107
    ## 3  Candidatus_Saccharibacteria 280.00000   280
    ## 4    Cyanobacteria/Chloroplast  64.25000   257
    ## 5          Deinococcus-Thermus  52.00000    52
    ## 6                   Firmicutes 179.24771 58614
    ## 7                 Fusobacteria   2.00000     2
    ## 8               Proteobacteria  59.09091   650
    ## 9                  Tenericutes 234.00000   234
    ## 10             Verrucomicrobia 104.00000   104

``` r
# Define phyla to filter
filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)
ps1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 381 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 381 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 381 tips and 379 internal nodes ]

``` r
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](cc1_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

``` r
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
```

    ## [1] 18

``` r
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
```

``` r
# How many genera would be present after filtering?
```

``` r
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
```

    ## [1] 49

``` r
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
```

``` r
h1 = 0.4
ps4 = tip_glom(ps2, h = h1)
```

``` r
multiPlotTitleTextSize = 15
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
```

``` r
# group plots together
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)
```

![](cc1_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

``` r
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
```

``` r
# Transform to relative abundance. Save as new object.
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})
```

``` r
plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2,  plotBefore, plotAfter)
```

![](cc1_files/figure-gfm/unnamed-chunk-55-1.png)<!-- -->

``` r
psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)
```

![](cc1_files/figure-gfm/unnamed-chunk-56-1.png)<!-- -->

``` r
.cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools",
                  "reshape2", "PMA", "structSSI", "ade4",
                  "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("genefilter", "impute")
```

``` r
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)){
  install.packages(.cran_packages[!.inst],repos = "http://cran.rstudio.com/")
}
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

    ## Warning: package 'structSSI' is not available for this version of R
    ## 
    ## A version of this package for your version of R might be available elsewhere,
    ## see the ideas at
    ## https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages

``` r
.inst <- .github_packages %in% installed.packages()
if (any(!.inst)){
  devtools::install_github(.github_packages[!.inst])
}
```

    ## Skipping install of 'phyloseqGraphTest' from a github remote, the SHA1 (3fb6c274) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)){BiocManager::install(.bioc_packages[!.inst])
}
```

``` r
qplot(sample_data(ps)$age, geom = "histogram",binwidth=20) + xlab("age")
```

![](cc1_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->

``` r
qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
```

![](cc1_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->

``` r
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship=gsub(" ","",sample_data(ps)$family_relationship)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTATCCGGAATGACTGGGCGTAAAGGGTGCGTAGGTGGTTTGGCAAGTTGGTAGCGTAATTCCGGGGCTCAACCTCGGCGCTACTACCAAAACTGCTGGACTTGAGTGCAGGAGGGGTGAATGGAATTCCTAGTGTAGCGGTGGAATGCGTAGATATTAGGAAGAACACCAGCGGCGAAGGCGATTCACTGGACTGTAACTGACACTGAGGCACGAAAGCGTGGGGAG
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](cc1_files/figure-gfm/unnamed-chunk-63-1.png)<!-- -->

``` r
rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram",binwidth=0.05) +
  xlab("Relative abundance")
```

![](cc1_files/figure-gfm/unnamed-chunk-64-1.png)<!-- -->

``` r
outliers <- c("F5D165", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
ps <- prune_samples(!(sample_names(ps) %in% outliers), ps)
```

``` r
which(!rowSums(otu_table(ps)) > 1000)
```

    ## F5D145 M1D149   M1D9 M2D125  M2D19 M3D148 M3D149   M3D3   M3D5   M3D8 
    ##     69    185    200    204    218    243    244    252    256    260

``` r
ps <- prune_samples(rowSums(otu_table(ps)) > 1000, ps)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
```

``` r
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "age_binned",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](cc1_files/figure-gfm/unnamed-chunk-68-1.png)<!-- -->

``` r
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "age_binned", label= "SampleID",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](cc1_files/figure-gfm/unnamed-chunk-69-1.png)<!-- -->

``` r
plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](cc1_files/figure-gfm/unnamed-chunk-70-1.png)<!-- -->

``` r
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCTTGATAAGTCTGAAGTGAAAGGCCAAGGCTTAACCATGGAACTGCTTTGGAAACTATGAGGCTAGAGTGCTGGAGAGGTAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGCTTACTGGACAGAAACTGACGTTGAGGCTCGAAAGCGTGGGGAG
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned",
                  shape = "family_relationship") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age", shape = "Litter")
```

![](cc1_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->

``` r
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))
```

``` r
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1
```

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(reshape2)
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")
```

![](cc1_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->

``` r
library(ade4)
```

    ## 
    ## Attaching package: 'ade4'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     score

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     score

``` r
ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         SampleID = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))
tax <- tax_table(ps) %>%
  data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))
row_scores <- row_scores %>%
  left_join(sample_data(pslog))
```

    ## Joining, by = "SampleID"

``` r
col_scores <- col_scores %>%
  left_join(tax)
```

    ## Joining, by = "seq"

``` r
evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
ggplot() +
  geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2), shape = 2) +
  geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
             size = .3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  facet_grid(~ age_binned) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](cc1_files/figure-gfm/unnamed-chunk-76-1.png)<!-- -->

``` r
ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ age_binned + family_relationship)
```

``` r
library(ggrepel)
ps_scores <- vegan::scores(ps_ccpna)
sites <- data.frame(ps_scores$sites)
sites$SampleID <- rownames(sites)
sites <- sites %>%
  left_join(sample_data(ps))
```

    ## Joining, by = "SampleID"

``` r
species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(ps)))
species <- species %>%
  left_join(tax)
```

    ## Joining, by = "otu_id"

``` r
evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                    aes(x = CCA1, y = CCA2, label = otu_id),
            size = 1.5, segment.size = 0.1) +
  facet_grid(. ~ family_relationship) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
        y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45   ) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

    ## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](cc1_files/figure-gfm/unnamed-chunk-78-1.png)<!-- -->

``` r
library(caret)
```

    ## Loading required package: lattice

``` r
sample_data(pslog)$age2 <- cut(sample_data(pslog)$age, c(0, 100, 400))
dataMatrix <- data.frame(age = sample_data(pslog)$age2, otu_table(pslog))
# take 8 mice at random to be the training set, and the remaining 4 the test set
trainingMice <- sample(unique(sample_data(pslog)$host_subject_id), size = 8)
inTrain <- which(sample_data(pslog)$host_subject_id %in% trainingMice)
training <- dataMatrix[inTrain,]
testing <- dataMatrix[-inTrain,]
plsFit <- train(age ~ ., data = training,
                method = "pls", preProc = "center")
```

``` r
plsClasses <- predict(plsFit, newdata = testing)
table(plsClasses, testing$age)
```

    ##            
    ## plsClasses  (0,100] (100,400]
    ##   (0,100]        73         2
    ##   (100,400]       3        45

``` r
library(randomForest)
```

    ## randomForest 4.6-14

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

``` r
rfFit <- train(age ~ ., data = training, method = "rf",
               preProc = "center", proximity = TRUE)
rfClasses <- predict(rfFit, newdata = testing)
table(rfClasses, testing$age)
```

    ##            
    ## rfClasses   (0,100] (100,400]
    ##   (0,100]        75         1
    ##   (100,400]       1        46

``` r
library(vegan)
```

    ## Loading required package: permute

    ## This is vegan 2.5-7

    ## 
    ## Attaching package: 'vegan'

    ## The following object is masked from 'package:caret':
    ## 
    ##     tolerance

``` r
pls_biplot <- list("loadings" = loadings(plsFit$finalModel),
                   "scores" = scores(plsFit$finalModel))
class(pls_biplot$scores) <- "matrix"

pls_biplot$scores <- data.frame(sample_data(pslog)[inTrain, ],pls_biplot$scores)

tax <- tax_table(ps)@.Data %>%
  data.frame(stringsAsFactors = FALSE)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
class(pls_biplot$loadings) <- "matrix"
pls_biplot$loadings <- data.frame(tax, pls_biplot$loadings)
```

``` r
ggplot() +
  geom_point(data = pls_biplot$scores,
             aes(x = Comp.1, y = Comp.2), shape = 2) +
  geom_point(data = pls_biplot$loadings,
             aes(x = 25 * Comp.1, y = 25 * Comp.2, col = Order),
             size = 0.3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Axis1", y = "Axis2", col = "Binned Age") +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  facet_grid( ~ age2) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](cc1_files/figure-gfm/unnamed-chunk-83-1.png)<!-- -->

``` r
rf_prox <- cmdscale(1 - rfFit$finalModel$proximity) %>%
  data.frame(sample_data(pslog)[inTrain, ])

ggplot(rf_prox) +
  geom_point(aes(x = X1, y = X2, col = age_binned),
             size = 1, alpha = 0.7) +
  scale_color_manual(values = c("#A66EB8", "#238DB5", "#748B4F")) +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  labs(col = "Binned Age", x = "Axis1", y = "Axis2")
```

![](cc1_files/figure-gfm/unnamed-chunk-84-1.png)<!-- -->

``` r
as.vector(tax_table(ps)[which.max(importance(rfFit$finalModel)), c("Family", "Genus")])
```

    ## [1] "Lachnospiraceae" "Roseburia"

``` r
impOtu <- as.vector(otu_table(pslog)[,which.max(importance(rfFit$finalModel))])
maxImpDF <- data.frame(sample_data(pslog), abund = impOtu)
ggplot(maxImpDF) +   geom_histogram(aes(x = abund)) +
  facet_grid(age2 ~ .) +
  labs(x = "Abundance of discriminative bacteria", y = "Number of samples")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](cc1_files/figure-gfm/unnamed-chunk-86-1.png)<!-- -->
