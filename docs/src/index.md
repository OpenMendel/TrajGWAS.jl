# TrajGWAS.jl

TrajGWAS.jl is a Julia package for performing genome-wide association studies (GWAS) for continuous longitudinal phenotypes using a modified linear mixed effects model. It builds upon the [within-subject variance estimation by robust regression (WiSER) method](https://github.com/OpenMendel/WiSER.jl) and can be used to identify variants associated with changes in the mean and within-subject variability of the longitduinal trait. The estimation procedure is robust to distributional misspecifications of both the random effects and the response. A saddlepoint approximation (SPA) option is implemented to provide improved power and calibrated type I error for rare variants. 


TrajGWAS.jl currently supports [PLINK](https://zzz.bwh.harvard.edu/plink/), [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) (both dosage and genotype data), and [BGEN](https://www.well.ox.ac.uk/~gav/bgen_format/) file formats. We plan to add [PGEN](https://www.cog-genomics.org/plink/2.0/formats#pgen) support in the future. 

## Installation

This package requires Julia v1.6 or later. The package is registered and can be installed using the command:

```julia
using Pkg
pkg"add TrajGWAS"
```

To run the code in this document, the packages installed by the following command are also necessary:
```julia
pkg"add BenchmarkTools CSV Glob PrettyTables"
```


```julia
# machine information for this tutorial
versioninfo()
```

    Julia Version 1.7.1
    Commit ac5cc99908 (2021-12-22 19:35 UTC)
    Platform Info:
      OS: macOS (x86_64-apple-darwin19.5.0)
      CPU: Intel(R) Core(TM) i7-7820HQ CPU @ 2.90GHz
      WORD_SIZE: 64
      LIBM: libopenlibm
      LLVM: libLLVM-12.0.1 (ORCJIT, skylake)



```julia
# for use in this tutorial
ENV["COLUMNS"] = 400
using BenchmarkTools, CSV, Glob, SnpArrays, TrajGWAS, PrettyTables
```

## Example data sets

The `data` folder of the package contains the example data sets for use with PLINK and VCF Files. In general, the user can locate this folder by command:


```julia
pvalpath = "trajgwas.pval.txt"
nullpath = "trajgwas.null.txt"
const datadir = normpath(joinpath(dirname(pathof(TrajGWAS)), "../data/"))
```




    "/Users/xyz/.julia/dev/TrajGWAS/data/"




```julia
# content of the data folder
readdir(glob"*.*", datadir)
```




    14-element Vector{String}:
     "/Users/xyz/.julia/dev/TrajGWAS/data/bgen_snpsetfile.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/covariate.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/example.8bits.bgen"
     "/Users/xyz/.julia/dev/TrajGWAS/data/example.8bits.bgen.bgi"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.bed"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.bim"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.fam"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap_snpsetfile.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/sim_data.jl"
     "/Users/xyz/.julia/dev/TrajGWAS/data/snpsetfile_vcf.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/test_vcf.vcf.gz"
     "/Users/xyz/.julia/dev/TrajGWAS/data/trajgwas_bgen_ex.csv"
     "/Users/xyz/.julia/dev/TrajGWAS/data/trajgwas_plinkex.csv"
     "/Users/xyz/.julia/dev/TrajGWAS/data/trajgwas_vcfex.csv"



The `hapmap3` files `trajgwas_plinkex.csv` file correspond to data examples using PLINK formatted files (.bed, .bim, .fam). 

The `test_vcf.vcf.gz` and `trajgwas_vcfex.csv` files are for an example analysis using VCF formatted files. 

The `example.8bits.bgen` and `trajgwas_bgen_ex.csv` files are for an example analysis using VCF formatted files. 


## Basic usage

The following command performs GWAS using for the hapmap3 PLINK files. The output is the fitted null model.

The default type of GWAS performed is a single-snp significance genome-wide scan, this can be changed by the keyword `analysistype` (default is "singlesnp"). Other types of analyses are covered later. It outputs the null model, runtime of fitting the null model, and convergence metrics. 


```julia
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = pvalpath,
        nullfile = nullpath)
```

    
    ******************************************************************************
    This program contains Ipopt, a library for large-scale nonlinear optimization.
     Ipopt is released as open source code under the Eclipse Public License (EPL).
             For more information visit https://github.com/coin-or/Ipopt
    ******************************************************************************
    
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.602286
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.040254





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 324
    Total observations: 3240
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  13.2282     0.146459    90.32    <1e-99
    β2: sex          -3.29295    0.2101     -15.67    <1e-54
    β3: onMeds        0.459585   0.0596002    7.71    <1e-13
    τ1: (Intercept)   0.792508   0.0850728    9.32    <1e-19
    τ2: sex          -0.2865     0.0970732   -2.95    0.0032
    τ3: onMeds        0.422303   0.063825     6.62    <1e-10
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  3.32057
    




For documentation of the `trajgwas` function, type `?trajgwas` in Julia REPL.
```@docs
trajgwas
```

### Formula for null model

The first three arguments specify the null model without SNP effects. The fourth argument specifies the grouping variable (subject ID, cluster, etc.).

- The first term is the mean fixed effects formula e.g., `@formula(y ~ 1 + sex + onMeds)`.
- The second argument is the random (location) effects -- `@formula(y ~ 1)` specifies a random intercept, whereas `@formula(y ~ 1 + time)` would specify a random intercept and slope (time).
- The third argument specifies the within-subject variance (WSV) formula `@formula(y ~ 1 + sex + onMeds)` models the intra-individual variability of `y` with an intercept, sex, and medication. Note: The WSV formula must have an intercept. 
- The fourth argument specifies the grouping variable of the longitudinal data. In the example data, the `id` variable specifies the grouping variable the repeated measures are collected on. 



### Input files

`trajgwas` expects two input files: one for responses plus covariates (fifth argument) in long format, the other the genetic file(s) for dosages/genotypes (sixth argument).

#### Covariate and trait file

Covariates and phenotype are provided in a csv file, e.g., `trajgwas_plinkex.csv`, which has one header line for variable names. In this example, variable `y` is continuous variable with repeated measures on individuals specified by `id`. We want to include variable `sex` and `onMeds` as covariates for GWAS.


```julia
run(`head $(datadir)trajgwas_plinkex.csv`);
```

    sex,onMeds,snp1,snp2,snp3,snp4,y,id
    0.0,1.0,0.0,1.0,2.0,0.0,12.26667411332518,A1
    0.0,0.0,0.0,1.0,2.0,0.0,10.268123812744903,A1
    0.0,0.0,0.0,1.0,2.0,0.0,12.165997408822557,A1
    0.0,0.0,0.0,1.0,2.0,0.0,11.879709602937222,A1
    0.0,0.0,0.0,1.0,2.0,0.0,12.812705165990007,A1
    0.0,0.0,0.0,1.0,2.0,0.0,9.987659617201372,A1
    0.0,0.0,0.0,1.0,2.0,0.0,12.140779426464974,A1
    0.0,1.0,0.0,1.0,2.0,0.0,13.205778146177705,A1
    0.0,0.0,0.0,1.0,2.0,0.0,11.363145919060207,A1


#### Genetic file(s)

TrajGWAS supports PLINK files, VCF files, and BGEN files.

- Genotype data is available as binary PLINK files.

- TrajGWAS can use dosage or genotype data from VCF files. 

- TrajGWAS uses dosage data from BGEN files.

!!! note

    By default, TrajGWAS assumes a set of PLINK files will be used. When using a VCF File, VCF file and type of data (dosage, genotype) must be specified by the `geneticformat` and `vcftype` options (as shown later). Similarly, BGEN must be specified as the `geneticformat` if using a BGEN file. 


```julia
readdir(glob"hapmap3.*", datadir)
```




    3-element Vector{String}:
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.bed"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.bim"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.fam"



In this example, there are 324 samples at 13,928 SNPs.


```julia
size(SnpArray(datadir * "hapmap3.bed"))
```




    (324, 13928)



Compressed PLINK and VCF files are supported. For example, if Plink files are `hapmap3.bed.gz`, `hapmap3.bim.gz` and `hapmap3.fam.gz`, the same command
```julia
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = pvalpath,
        nullfile = nullpath)
```
still works. Check all supported compression format by


```julia
SnpArrays.ALLOWED_FORMAT
```




    6-element Vector{String}:
     "gz"
     "zlib"
     "zz"
     "xz"
     "zst"
     "bz2"



### Output files

`trajgwas` outputs two files: `trajgwas.null.txt` and `trajgwas.pval.txt`. 

* `trajgwas.null.txt` lists the estimated null model (without SNPs). 


```julia
run(`cat trajgwas.null.txt`);
```

    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 324
    Total observations: 3240
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  13.2282     0.146459    90.32    <1e-99
    β2: sex          -3.29295    0.2101     -15.67    <1e-54
    β3: onMeds        0.459585   0.0596002    7.71    <1e-13
    τ1: (Intercept)   0.792508   0.0850728    9.32    <1e-19
    τ2: sex          -0.2865     0.0970732   -2.95    0.0032
    τ3: onMeds        0.422303   0.063825     6.62    <1e-10
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  3.32057
    


* `trajgwas.pval.txt` tallies the SNPs, their pvalues, and relevant information on each SNP.
    - `betapval` represents the p-value of the SNP's effect on the mean of the trait. If `spa=true` (default), then this is the SPA p-value. If `spa=false`, then this is the score test beta p-value. 
    - `taupval` represents the p-value of the SNP's effect on the within-subject variability of the trait. If `spa=true` (default), then this is the SPA p-value. If `spa=false`, then this is the score test tau p-value. 
    - `jointpval` represents a joint p-value of the SNP's effect on both the mean and variance. By default `spa=true` this is the harmonic mean of the saddlepoint approximated p-values for beta and tau. If `spa=false`, this is the joint score test p-value.


```julia
pretty_table(first(CSV.read("trajgwas.pval.txt", DataFrame), 8), tf = tf_markdown)
```

    |[1m   chr [0m|[1m     pos [0m|[1m      snpid [0m|[1m        maf [0m|[1m     hwepval [0m|[1m    betapval [0m|[1m betadir [0m|[1m     taupval [0m|[1m taudir [0m|[1m   jointpval [0m|
    |[90m Int64 [0m|[90m   Int64 [0m|[90m   String31 [0m|[90m    Float64 [0m|[90m     Float64 [0m|[90m     Float64 [0m|[90m   Int64 [0m|[90m     Float64 [0m|[90m  Int64 [0m|[90m     Float64 [0m|
    |-------|---------|------------|------------|-------------|-------------|---------|-------------|--------|-------------|
    |     1 |  554484 | rs10458597 |        0.0 |         1.0 |         1.0 |       0 |         1.0 |      0 |         1.0 |
    |     1 |  758311 | rs12562034 |  0.0776398 |    0.409876 |     0.38939 |      -1 |    0.922822 |     -1 |    0.547682 |
    |     1 |  967643 |  rs2710875 |   0.324074 |  4.07625e-7 |  5.83856e-6 |       1 |  3.35408e-6 |     -1 |  4.93598e-7 |
    |     1 | 1168108 | rs11260566 |   0.191589 |    0.128568 | 0.000850889 |       1 |   0.0559055 |      1 |  0.00101827 |
    |     1 | 1375074 |  rs1312568 |   0.441358 |  2.5376e-19 | 0.000171683 |      -1 | 1.84811e-13 |      1 | 2.07091e-15 |
    |     1 | 1588771 | rs35154105 |        0.0 |         1.0 |         1.0 |       0 |         1.0 |      0 |         1.0 |
    |     1 | 1789051 | rs16824508 | 0.00462963 |    0.933278 |    0.295035 |       1 |    0.304109 |      1 |    0.299504 |
    |     1 | 1990452 |  rs2678939 |   0.453704 | 5.07696e-11 |  1.27871e-7 |       1 |  7.73809e-9 |     -1 |   3.3986e-9 |


Output file names can be changed by the `nullfile` and `pvalfile` keywords respectively. For example, 
```julia
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = "trajgwas.pval.txt.gz",
        nullfile = nullpath)

```
will output the p-value file in compressed gz format.

### Subsamples

Use the keyword `covrowinds` to specify selected samples in the covarite file. Use the keyword `geneticrowinds` to specify selected samples in the Plink (.bed), VCF, or BGEN file. For example, to use the first 300 samples in both covariate and bed file:
```julia
trajgwas(@formula(trait ~ 1+ sex), 
    @formula(trait ~ 1), 
    @formula(trait ~ 1 + sex), 
    :id,
    covfile, geneticfile, 
    covrowinds=1:300, geneticrowinds=1:300)
```
!!! note

    Users should make sure that the selected samples in covariate file match exactly those in bed file, otherwise an error will be displayed. 

### Input non-genetic data as DataFrame

Internally `trajgwas` parses the covariate file as a DataFrame by `df = CSV.read(covfile, DataFrame)`. For covariate file of other formats, or for specifying default levels of categorical variables, users can parse their datafile as a DataFrame and then input the DataFrame to `trajgwas` directly.
```julia
trajgwas(@formula(trait ~ 1+ sex), 
    @formula(trait ~ 1), 
    @formula(trait ~ 1 + sex),
    :id,
    df, geneticfile)
```
!!! note

    Users should always make sure that individuals in covariate file or DataFrame match those in Plink fam file/VCF File/BGEN file. 

For example, following code checks that the subject ID of the `covariate.txt` file match that of the `hapmap3.fam` file exactly.


```julia
covdf = CSV.read(datadir * "covariate.txt", DataFrame)
plkfam = CSV.read(datadir * "hapmap3.fam", DataFrame, header=0, delim=' ')
all(covdf[!, 2] .== plkfam[!, 2])
```




    true




```julia
pretty_table(first(covdf, 8), tf = tf_markdown)
```

    |[1m   famid [0m|[1m   perid [0m|[1m  faid [0m|[1m  moid [0m|[1m   sex [0m|[1m trait [0m|
    |[90m String7 [0m|[90m String7 [0m|[90m Int64 [0m|[90m Int64 [0m|[90m Int64 [0m|[90m Int64 [0m|
    |---------|---------|-------|-------|-------|-------|
    |    2431 | NA19916 |     0 |     0 |     1 |     4 |
    |    2424 | NA19835 |     0 |     0 |     2 |     4 |
    |    2469 | NA20282 |     0 |     0 |     2 |     4 |
    |    2368 | NA19703 |     0 |     0 |     1 |     3 |
    |    2425 | NA19901 |     0 |     0 |     2 |     3 |
    |    2427 | NA19908 |     0 |     0 |     1 |     4 |
    |    2430 | NA19914 |     0 |     0 |     2 |     4 |
    |    2470 | NA20287 |     0 |     0 |     2 |     1 |



```julia
pretty_table(first(plkfam, 8), tf = tf_markdown)
```

    |[1m Column1 [0m|[1m Column2 [0m|[1m Column3 [0m|[1m Column4 [0m|[1m Column5 [0m|[1m Column6 [0m|
    |[90m String7 [0m|[90m String7 [0m|[90m   Int64 [0m|[90m   Int64 [0m|[90m   Int64 [0m|[90m   Int64 [0m|
    |---------|---------|---------|---------|---------|---------|
    |      A1 | NA19916 |       0 |       0 |       1 |      -9 |
    |       2 | NA19835 |       0 |       0 |       2 |      -9 |
    |       3 | NA20282 |       0 |       0 |       2 |      -9 |
    |       4 | NA19703 |       0 |       0 |       1 |      -9 |
    |       5 | NA19901 |       0 |       0 |       2 |      -9 |
    |       6 | NA19908 |       0 |       0 |       1 |      -9 |
    |       7 | NA19914 |       0 |       0 |       2 |      -9 |
    |       8 | NA20287 |       0 |       0 |       2 |      -9 |


### Timing

For this moderate-sized data set, `trajgwas` takes around 1 second without applying SPA. With SPA, it takes a little longer. 


```julia
@btime(trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = pvalpath,
        nullfile = nullpath, 
        usespa = false));
```

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044445
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.042232
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043275
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039053
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044923
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039103
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044385
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039346
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043468
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039741
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043577
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.040535
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043344
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039667
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043612
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039483
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.045086
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039441
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044173
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.040099
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044949
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.040097
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044071
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.041057
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044740
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.040537
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044379
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.042448
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.045716
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.042732
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043419
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.040273
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044672
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.040251
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.045095
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.041001
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043335
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039804
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044344
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039873
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043518
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039340
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044033
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039883
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043352
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.044078
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043401
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039271
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044151
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039558
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043364
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039549
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044107
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.040570
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043527
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.040298
      459.587 ms (1433770 allocations: 127.29 MiB)



```julia
# clean up
rm("trajgwas.null.txt", force=true)
rm("trajgwas.pval.txt", force=true)
```

## VCF Formatted Files

By default, TrajGWAS.jl will assume you are using PLINK files. It also supports VCF (and BGEN) Files. `vcf.gz` files are supported as well. To use vcf files in any of the analysis options detailed in this documentation, you simply need to add two keyword options to the `trajgwas` function:
* `geneticformat`: Choices are "VCF" or "PLINK". If you are using a VCF file, use `geneticformat = "VCF"`.
* `vcftype`: Choices are :GT (for genotypes) or :DS (for dosages). This tells TrajGWAS which type of data to extract from the VCF file.

Using a VCF File does not output minor allele frequency or hardy weinberg equillibrium p-values for each SNP tested since they may be dosages. 

The following shows how to run an analysis with a VCF file using the dosage information. 


```julia
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_vcfex.csv",
        datadir * "test_vcf",
        pvalfile = pvalpath,
        geneticformat = "VCF",
        vcftype = :DS)
```

    run = 1, ‖Δβ‖ = 0.003459, ‖Δτ‖ = 0.421747, ‖ΔL‖ = 0.002492, status = LOCALLY_SOLVED, time(s) = 0.020257
    run = 2, ‖Δβ‖ = 0.001369, ‖Δτ‖ = 0.015998, ‖ΔL‖ = 0.003244, status = LOCALLY_SOLVED, time(s) = 0.017106





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 191
    Total observations: 1910
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  12.1708     0.18351     66.32    <1e-99
    β2: sex          -3.41012    0.253155   -13.47    <1e-40
    β3: onMeds        0.485447   0.0701423    6.92    <1e-11
    τ1: (Intercept)   0.774299   0.0849867    9.11    <1e-19
    τ2: sex          -0.524117   0.106719    -4.91    <1e-06
    τ3: onMeds        0.633337   0.0814974    7.77    <1e-14
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  2.86452
    





```julia
pretty_table(first(CSV.read("trajgwas.pval.txt", DataFrame), 8), tf = tf_markdown)
```

    |[1m   chr [0m|[1m      pos [0m|[1m       snpid [0m|[1m  betapval [0m|[1m betadir [0m|[1m    taupval [0m|[1m taudir [0m|[1m  jointpval [0m|
    |[90m Int64 [0m|[90m    Int64 [0m|[90m    String15 [0m|[90m   Float64 [0m|[90m   Int64 [0m|[90m    Float64 [0m|[90m  Int64 [0m|[90m    Float64 [0m|
    |-------|----------|-------------|-----------|---------|------------|--------|------------|
    |    22 | 20000086 | rs138720731 |  0.322684 |      -1 |    0.61956 |     -1 |   0.424353 |
    |    22 | 20000146 |  rs73387790 |       1.0 |       0 |        1.0 |      0 |        1.0 |
    |    22 | 20000199 | rs183293480 |  0.214207 |       1 |   0.190785 |      1 |   0.201819 |
    |    22 | 20000291 | rs185807825 |  0.292231 |       1 |   0.200029 |      1 |   0.237495 |
    |    22 | 20000428 |  rs55902548 | 0.0140114 |       1 | 0.00369005 |      1 | 0.00576184 |
    |    22 | 20000683 | rs142720028 |       1.0 |       0 |        1.0 |      0 |        1.0 |
    |    22 | 20000771 | rs114690707 |  0.536406 |      -1 |   0.145113 |      1 |    0.22843 |
    |    22 | 20000793 | rs189842693 |  0.161414 |      -1 |   0.757017 |      1 |   0.266091 |


## BGEN Formatted Files

By default, TrajGWAS.jl will assume you are using PLINK files. It also supports BGEN Files. To use BGEN files in any of the analysis options detailed in this documentation, you simply need to add the following keyword option to the `trajgwas` function:
* `geneticformat`: Choices are "VCF" or "PLINK" or "BGEN". If you are using a BGEN file, use `geneticformat = "BGEN"`.

Using a BGEN File does not output minor allele frequency or hardy weinberg equillibrium p-values for each SNP tested.

Some features, such as SNP-set analyses, are only available when there's an indexing `.bgi` file available. Additionally, if the BGEN file does not contain sample information (i.e. sample IDs), then a sample file is necessary and its path can be specified with the `samplepath` keyword. 

!!! note

    BGEN files can contain an optional index file (`.bgi` file) that allows the variants to be specified in order of position. TrajGWAS will automatically look for a file in the same directory as the `BGENFILENAME` with the name `BGENFILENAME.bgi`. The BGEN file is read either in the `.bgi` file order if `BGENFILENAME.bgi` is supplied in the same directory as the BGEN file, otherwise it will use the order in the BGEN file. This is important in analyses specifying `snpinds` as well as annotation groupings. You must make sure this matches the way the BGEN file will be read in. 

The following shows how to run an analysis with a BGEN file using the dosage information. 


```julia
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id, 
        datadir * "trajgwas_bgen_ex.csv", 
        datadir * "example.8bits", 
        geneticformat = "BGEN", 
        pvalfile = pvalpath)
```

    run = 1, ‖Δβ‖ = 0.096578, ‖Δτ‖ = 0.162563, ‖ΔL‖ = 0.007625, status = LOCALLY_SOLVED, time(s) = 0.055072
    run = 2, ‖Δβ‖ = 0.003173, ‖Δτ‖ = 0.005584, ‖ΔL‖ = 0.001468, status = LOCALLY_SOLVED, time(s) = 0.040037





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 500
    Total observations: 5000
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  10.5275     0.0992675  106.05    <1e-99
    β2: sex          -3.38067    0.146416   -23.09    <1e-99
    β3: onMeds        0.522162   0.0358727   14.56    <1e-47
    τ1: (Intercept)   0.294364   0.0443191    6.64    <1e-10
    τ2: sex          -0.365009   0.0503589   -7.25    <1e-12
    τ3: onMeds        0.559769   0.0467046   11.99    <1e-32
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  2.51561
    





```julia
pretty_table(first(CSV.read("trajgwas.pval.txt", DataFrame), 8), tf = tf_markdown)
```

    |[1m   chr [0m|[1m   pos [0m|[1m    snpid [0m|[1m     varid [0m|[1m  hwepval [0m|[1m      maf [0m|[1m infoscore [0m|[1m    betapval [0m|[1m betadir [0m|[1m     taupval [0m|[1m taudir [0m|[1m   jointpval [0m|
    |[90m Int64 [0m|[90m Int64 [0m|[90m String15 [0m|[90m  String15 [0m|[90m  Float64 [0m|[90m  Float64 [0m|[90m   Float64 [0m|[90m     Float64 [0m|[90m   Int64 [0m|[90m     Float64 [0m|[90m  Int64 [0m|[90m     Float64 [0m|
    |-------|-------|----------|-----------|----------|----------|-----------|-------------|---------|-------------|--------|-------------|
    |     1 |  1001 | RSID_101 | SNPID_101 |  0.86274 | 0.416977 |  0.984639 |    0.510629 |      -1 |    0.425112 |     -1 |    0.463963 |
    |     1 |  2000 |   RSID_2 |   SNPID_2 | 0.192181 |  0.19751 |       9.0 |  4.5062e-10 |       1 | 1.40173e-21 |     -1 | 6.67473e-21 |
    |     1 |  2001 | RSID_102 | SNPID_102 |   0.1844 | 0.197667 |  0.727308 | 4.47847e-10 |      -1 | 1.79121e-21 |      1 | 8.45783e-21 |
    |     1 |  3000 |   RSID_3 |   SNPID_3 | 0.965354 | 0.483396 |  0.955355 |   0.0721822 |      -1 |    0.600949 |     -1 |    0.128884 |
    |     1 |  3001 | RSID_103 | SNPID_103 | 0.965354 | 0.483396 |  0.955355 |   0.0721821 |       1 |    0.600949 |      1 |    0.128884 |
    |     1 |  4000 |   RSID_4 |   SNPID_4 | 0.371927 |  0.21671 |  0.991768 |    0.112472 |      -1 |    0.133752 |     -1 |    0.122192 |
    |     1 |  4001 | RSID_104 | SNPID_104 | 0.371928 |  0.21671 |  0.991768 |    0.112471 |       1 |    0.133752 |      1 |    0.122192 |
    |     1 |  5000 |   RSID_5 |   SNPID_5 | 0.587013 | 0.388082 |  0.968258 |    0.526325 |      -1 |    0.542332 |     -1 |    0.534209 |


## SNP models

Genotypes are translated into numeric values according to different genetic model, which is specified by the `snpmodel` keyword. Default is `ADDITIVE_MODEL`.

| Genotype | `SnpArray` | `ADDITIVE_MODEL` | `DOMINANT_MODEL` | `RECESSIVE_MODEL` |    
|:---:|:---:|:---:|:---:|:---:|  
| A1,A1 | 0x00 | 0 | 0 | 0 |  
| missing | 0x01 | NaN | NaN | NaN |
| A1,A2 | 0x02 | 1 | 1 | 0 |  
| A2,A2 | 0x03 | 2 | 1 | 1 |  

!!! note

    `trajgwas` imputes missing genotypes according to minor allele frequencies. 
    
Users are advised to impute genotypes using more sophiscated methods before GWAS.

## SNP and/or sample masks

In practice, we often perform GWAS on selected SNPs and/or selected samples. They can be specified by the `snpinds`, `covrowinds` and `geneticrowinds` keywords of `trajgwas` function. 

For example, to perform GWAS on SNPs with minor allele frequency (MAF) above 0.05


```julia
# create SNP mask
snpinds = maf(SnpArray("../data/hapmap3.bed")) .≥ 0.05
# GWAS on selected SNPs
@time trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = "commonvariant.pval.txt",
        snpinds = snpinds)
```

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043349
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039752
      2.733297 seconds (8.07 M allocations: 492.816 MiB, 5.49% gc time, 65.58% compilation time)





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 324
    Total observations: 3240
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  13.2282     0.146459    90.32    <1e-99
    β2: sex          -3.29295    0.2101     -15.67    <1e-54
    β3: onMeds        0.459585   0.0596002    7.71    <1e-13
    τ1: (Intercept)   0.792508   0.0850728    9.32    <1e-19
    τ2: sex          -0.2865     0.0970732   -2.95    0.0032
    τ3: onMeds        0.422303   0.063825     6.62    <1e-10
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  3.32057
    





```julia
run(`head commonvariant.pval.txt`);
```

    chr	pos	snpid	maf	hwepval	betapval	betadir	taupval	taudir	jointpval
    1	758311	rs12562034	0.07763975155279501	0.4098763332666681	0.3893898715952268	-1	0.9228220976232073	-1	0.5476822137398478
    1	967643	rs2710875	0.32407407407407407	4.076249100705747e-7	5.838562452864651e-6	1	3.354075049151666e-6	-1	4.935983727867445e-7
    1	1168108	rs11260566	0.19158878504672894	0.1285682279446898	0.0008508894259509042	1	0.05590549977468929	1	0.0010182674469333282
    1	1375074	rs1312568	0.441358024691358	2.5376019650614977e-19	0.00017168324418138014	-1	1.8481146312336014e-13	1	2.0709090755361924e-15
    1	1990452	rs2678939	0.4537037037037037	5.07695957708431e-11	1.278709529589261e-7	1	7.7380891475822e-9	-1	3.3986012947082272e-9
    1	2194615	rs7553178	0.22685185185185186	0.17056143157457776	0.07404001756940115	-1	0.00023109048716816823	1	0.00010021413726889922
    1	2396747	rs13376356	0.1448598130841121	0.9053079215078139	0.48972088680069503	-1	0.6595593862223496	1	0.5620909278620404
    1	2823603	rs1563468	0.4830246913580247	4.23065537243926e-9	2.581476982676523e-7	-1	1.1700736777037956e-8	1	6.36797979359385e-8
    1	3025087	rs6690373	0.2538699690402477	9.238641887192776e-8	3.111965106232366e-5	1	2.0479935907554532e-7	-1	2.123094937745341e-5



```julia
# extra header line in commonvariant.pval.txt
countlines("commonvariant.pval.txt"), count(snpinds)
```




    (12086, 12085)




```julia
# clean up
rm("trajgwas.null.txt", force=true)
rm("commonvariant.pval.txt", force=true)
```

`covrowinds` specify the samples in the covariate file and `geneticrowinds` for PLINK or VCF File. User should be particularly careful when using these two keyword. Selected rows in SnpArray should exactly match the samples in the null model. Otherwise the results are meaningless.

## Estimating Effect Sizes (Wald)

By default, `trajgwas` calculates p-value for each SNP using SPA/score test. Score test is fast because it doesn't require fitting alternative model for each SNP. User can request the Wald p-values and the estimated effect size of each SNP using keyword `test=:wald`. Wald is much slower but will give you estimated effect sizes from the WiSER model.


```julia
@time trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = "wald.pval.txt",
        snpinds = 1:5,
        test = :wald)
```

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.051852
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.042647
    run = 1, ‖Δβ‖ = 0.033243, ‖Δτ‖ = 0.146511, ‖ΔL‖ = 0.005476, status = LOCALLY_SOLVED, time(s) = 0.052239
    run = 2, ‖Δβ‖ = 0.005774, ‖Δτ‖ = 0.042246, ‖ΔL‖ = 0.001784, status = LOCALLY_SOLVED, time(s) = 0.048045
    run = 1, ‖Δβ‖ = 0.013090, ‖Δτ‖ = 0.130781, ‖ΔL‖ = 0.005011, status = LOCALLY_SOLVED, time(s) = 0.057742
    run = 2, ‖Δβ‖ = 0.003913, ‖Δτ‖ = 0.037309, ‖ΔL‖ = 0.001516, status = LOCALLY_SOLVED, time(s) = 0.052004
    run = 1, ‖Δβ‖ = 0.022159, ‖Δτ‖ = 0.141135, ‖ΔL‖ = 0.005554, status = LOCALLY_SOLVED, time(s) = 0.048711
    run = 2, ‖Δβ‖ = 0.001482, ‖Δτ‖ = 0.021700, ‖ΔL‖ = 0.001435, status = LOCALLY_SOLVED, time(s) = 0.044404
    run = 1, ‖Δβ‖ = 0.026764, ‖Δτ‖ = 0.368620, ‖ΔL‖ = 0.000317, status = LOCALLY_SOLVED, time(s) = 0.047687
    run = 2, ‖Δβ‖ = 0.003023, ‖Δτ‖ = 0.030938, ‖ΔL‖ = 0.003568, status = LOCALLY_SOLVED, time(s) = 0.043295
      4.293410 seconds (7.17 M allocations: 431.003 MiB, 2.57% gc time, 86.00% compilation time)





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 324
    Total observations: 3240
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  13.2282     0.146459    90.32    <1e-99
    β2: sex          -3.29295    0.2101     -15.67    <1e-54
    β3: onMeds        0.459585   0.0596002    7.71    <1e-13
    τ1: (Intercept)   0.792508   0.0850728    9.32    <1e-19
    τ2: sex          -0.2865     0.0970732   -2.95    0.0032
    τ3: onMeds        0.422303   0.063825     6.62    <1e-10
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  3.32057
    




Note the extra `effect` column in pvalfile, which is the effect size (regression coefficient) for each SNP. 


```julia
run(`head wald.pval.txt`);
```

    chr	pos	snpid	maf	hwepval	betaeffect	betapval	taueffect	taupval
    1	554484	rs10458597	0.0	1.0	0.0	1.0	0.0	1.0
    1	758311	rs12562034	0.07763975155279501	0.4098763332666681	-0.231072479332816	0.3929844692486013	0.01340575997086949	0.9216488777690096
    1	967643	rs2710875	0.32407407407407407	4.076249100705747e-7	0.6406954913590558	8.824209286611495e-7	-0.34763283794106525	1.0836367479719384e-11
    1	1168108	rs11260566	0.19158878504672894	0.1285682279446898	0.6188250819875118	0.0002537468407208455	-0.14748437551458396	0.0502802758129169
    1	1375074	rs1312568	0.441358024691358	2.5376019650614977e-19	-0.4560323801786057	0.00016592110001810087	0.4964521486169327	1.177346561221082e-36



```julia
# clean up
rm("wald.pval.txt", force=true)
rm("trajgwas.null.txt", force=true)
```

## Score test for screening, Wald for effect size estimates 

For large data sets, a practical solution is to perform the score test first across all SNPs, then re-do Wald for the most promising SNPs according to score test p-values in order to get estimated effect sizes.

**Step 1**: Perform score test GWAS, results in `score.pval.txt`.


```julia
@time trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = "score.pval.txt")
```

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044256
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.041231
      1.121874 seconds (4.00 M allocations: 259.048 MiB, 13.18% gc time)





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 324
    Total observations: 3240
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  13.2282     0.146459    90.32    <1e-99
    β2: sex          -3.29295    0.2101     -15.67    <1e-54
    β3: onMeds        0.459585   0.0596002    7.71    <1e-13
    τ1: (Intercept)   0.792508   0.0850728    9.32    <1e-19
    τ2: sex          -0.2865     0.0970732   -2.95    0.0032
    τ3: onMeds        0.422303   0.063825     6.62    <1e-10
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  3.32057
    





```julia
pretty_table(first(CSV.read("score.pval.txt", DataFrame), 8), tf = tf_markdown)
```

    |[1m   chr [0m|[1m     pos [0m|[1m      snpid [0m|[1m        maf [0m|[1m     hwepval [0m|[1m    betapval [0m|[1m betadir [0m|[1m     taupval [0m|[1m taudir [0m|[1m   jointpval [0m|
    |[90m Int64 [0m|[90m   Int64 [0m|[90m   String31 [0m|[90m    Float64 [0m|[90m     Float64 [0m|[90m     Float64 [0m|[90m   Int64 [0m|[90m     Float64 [0m|[90m  Int64 [0m|[90m     Float64 [0m|
    |-------|---------|------------|------------|-------------|-------------|---------|-------------|--------|-------------|
    |     1 |  554484 | rs10458597 |        0.0 |         1.0 |         1.0 |       0 |         1.0 |      0 |         1.0 |
    |     1 |  758311 | rs12562034 |  0.0776398 |    0.409876 |     0.38939 |      -1 |    0.922822 |     -1 |    0.547682 |
    |     1 |  967643 |  rs2710875 |   0.324074 |  4.07625e-7 |  5.83856e-6 |       1 |  3.35408e-6 |     -1 |  4.93598e-7 |
    |     1 | 1168108 | rs11260566 |   0.191589 |    0.128568 | 0.000850889 |       1 |   0.0559055 |      1 |  0.00101827 |
    |     1 | 1375074 |  rs1312568 |   0.441358 |  2.5376e-19 | 0.000171683 |      -1 | 1.84811e-13 |      1 | 2.07091e-15 |
    |     1 | 1588771 | rs35154105 |        0.0 |         1.0 |         1.0 |       0 |         1.0 |      0 |         1.0 |
    |     1 | 1789051 | rs16824508 | 0.00462963 |    0.933278 |    0.295035 |       1 |    0.304109 |      1 |    0.299504 |
    |     1 | 1990452 |  rs2678939 |   0.453704 | 5.07696e-11 |  1.27871e-7 |       1 |  7.73809e-9 |     -1 |   3.3986e-9 |


**Step 2**: Sort score test p-values and find top 10 SNPs.


```julia
scorepvals = CSV.read("score.pval.txt", DataFrame)[!, :taupval] # tau p-values 
tophits = sortperm(scorepvals)[1:10] # indices of 10 SNPs with smallest p-values
scorepvals[tophits] # smallest 10 p-values
```




    10-element Vector{Float64}:
     7.872202557706978e-17
     1.8481146312336014e-13
     1.497868141369434e-12
     9.730062616490452e-12
     1.2892059340654512e-11
     1.4844638660225815e-11
     1.5167114834215392e-11
     1.8761256218261456e-11
     2.3491352960094672e-11
     2.6075629334507108e-11



**Step 3**: Re-do LRT on top hits.


```julia
@time trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = "wald.pval.txt",
        snpinds = tophits,
        test = :wald)
```

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043852
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.038897
    run = 1, ‖Δβ‖ = 0.026764, ‖Δτ‖ = 0.368620, ‖ΔL‖ = 0.000317, status = LOCALLY_SOLVED, time(s) = 0.051141
    run = 2, ‖Δβ‖ = 0.003023, ‖Δτ‖ = 0.030938, ‖ΔL‖ = 0.003568, status = LOCALLY_SOLVED, time(s) = 0.045226
    run = 1, ‖Δβ‖ = 0.039946, ‖Δτ‖ = 0.185049, ‖ΔL‖ = 0.005690, status = LOCALLY_SOLVED, time(s) = 0.047814
    run = 2, ‖Δβ‖ = 0.001640, ‖Δτ‖ = 0.029639, ‖ΔL‖ = 0.002028, status = LOCALLY_SOLVED, time(s) = 0.047543
    run = 1, ‖Δβ‖ = 0.040207, ‖Δτ‖ = 0.407863, ‖ΔL‖ = 0.001924, status = LOCALLY_SOLVED, time(s) = 0.049603
    run = 2, ‖Δβ‖ = 0.002632, ‖Δτ‖ = 0.039609, ‖ΔL‖ = 0.005131, status = LOCALLY_SOLVED, time(s) = 0.049566
    run = 1, ‖Δβ‖ = 0.052372, ‖Δτ‖ = 0.828946, ‖ΔL‖ = 0.001277, status = LOCALLY_SOLVED, time(s) = 0.061910
    run = 2, ‖Δβ‖ = 0.013694, ‖Δτ‖ = 0.185493, ‖ΔL‖ = 0.007118, status = LOCALLY_SOLVED, time(s) = 0.052748
    run = 1, ‖Δβ‖ = 0.032435, ‖Δτ‖ = 0.199659, ‖ΔL‖ = 0.004973, status = LOCALLY_SOLVED, time(s) = 0.049459
    run = 2, ‖Δβ‖ = 0.002847, ‖Δτ‖ = 0.032438, ‖ΔL‖ = 0.002941, status = LOCALLY_SOLVED, time(s) = 0.048901
    run = 1, ‖Δβ‖ = 0.009362, ‖Δτ‖ = 0.625969, ‖ΔL‖ = 0.002237, status = LOCALLY_SOLVED, time(s) = 0.061187
    run = 2, ‖Δβ‖ = 0.005893, ‖Δτ‖ = 0.176417, ‖ΔL‖ = 0.005017, status = LOCALLY_SOLVED, time(s) = 0.052250
    run = 1, ‖Δβ‖ = 0.012805, ‖Δτ‖ = 0.243313, ‖ΔL‖ = 0.006289, status = LOCALLY_SOLVED, time(s) = 0.051522
    run = 2, ‖Δβ‖ = 0.002636, ‖Δτ‖ = 0.044210, ‖ΔL‖ = 0.003312, status = LOCALLY_SOLVED, time(s) = 0.053374
    run = 1, ‖Δβ‖ = 0.014062, ‖Δτ‖ = 0.225841, ‖ΔL‖ = 0.003899, status = LOCALLY_SOLVED, time(s) = 0.054011
    run = 2, ‖Δβ‖ = 0.002112, ‖Δτ‖ = 0.030955, ‖ΔL‖ = 0.003860, status = LOCALLY_SOLVED, time(s) = 0.053688
    run = 1, ‖Δβ‖ = 0.026425, ‖Δτ‖ = 0.258805, ‖ΔL‖ = 0.004204, status = LOCALLY_SOLVED, time(s) = 0.051707
    run = 2, ‖Δβ‖ = 0.001090, ‖Δτ‖ = 0.051281, ‖ΔL‖ = 0.003822, status = LOCALLY_SOLVED, time(s) = 0.054497
    run = 1, ‖Δβ‖ = 0.035737, ‖Δτ‖ = 0.174702, ‖ΔL‖ = 0.001439, status = LOCALLY_SOLVED, time(s) = 0.053476
    run = 2, ‖Δβ‖ = 0.004776, ‖Δτ‖ = 0.019790, ‖ΔL‖ = 0.001867, status = LOCALLY_SOLVED, time(s) = 0.051735
      3.154140 seconds (5.48 M allocations: 365.218 MiB, 3.34% gc time, 58.54% compilation time)





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 324
    Total observations: 3240
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  13.2282     0.146459    90.32    <1e-99
    β2: sex          -3.29295    0.2101     -15.67    <1e-54
    β3: onMeds        0.459585   0.0596002    7.71    <1e-13
    τ1: (Intercept)   0.792508   0.0850728    9.32    <1e-19
    τ2: sex          -0.2865     0.0970732   -2.95    0.0032
    τ3: onMeds        0.422303   0.063825     6.62    <1e-10
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  3.32057
    





```julia
pretty_table(CSV.read("wald.pval.txt", DataFrame), tf = tf_markdown)
```

    |[1m   chr [0m|[1m       pos [0m|[1m     snpid [0m|[1m      maf [0m|[1m     hwepval [0m|[1m betaeffect [0m|[1m    betapval [0m|[1m taueffect [0m|[1m     taupval [0m|
    |[90m Int64 [0m|[90m     Int64 [0m|[90m  String15 [0m|[90m  Float64 [0m|[90m     Float64 [0m|[90m    Float64 [0m|[90m     Float64 [0m|[90m   Float64 [0m|[90m     Float64 [0m|
    |-------|-----------|-----------|----------|-------------|------------|-------------|-----------|-------------|
    |     1 |   1375074 | rs1312568 | 0.441358 |  2.5376e-19 |  -0.456032 | 0.000165921 |  0.496452 | 1.17735e-36 |
    |     1 |  11552817 | rs2745282 | 0.421053 |  1.3537e-14 |   -0.67282 |   7.6492e-9 |  0.425653 | 1.45395e-19 |
    |     1 | 120276030 | rs6688004 | 0.469136 | 1.75473e-27 |  -0.637817 |  5.13062e-9 |  0.430562 | 6.07457e-21 |
    |     2 | 135623558 | rs6730157 | 0.313272 | 7.76753e-16 |  -0.556166 |  1.10163e-5 |  0.462346 | 2.88628e-20 |
    |     9 | 126307510 | rs3814134 | 0.427469 | 4.28276e-33 |   0.723472 | 1.41601e-11 | -0.518216 | 4.75338e-39 |
    |    15 |  40603025 | rs2617236 | 0.390966 | 1.16466e-14 |  -0.644375 |  1.44948e-7 |  0.445647 | 9.43451e-19 |
    |    15 |  40803767 | rs3742988 |  0.42284 | 4.26181e-30 |  -0.796743 | 2.24785e-13 |  0.389234 | 7.09349e-17 |
    |    17 |  56509992 | rs8064681 | 0.481481 | 1.37485e-21 |   -0.58961 |   1.9256e-7 |  0.425417 | 5.60171e-25 |
    |    17 |  71293786 | rs2125345 | 0.339009 | 1.75527e-14 |  -0.635669 |  1.82739e-7 |  0.449612 |  1.7146e-21 |
    |    23 |  64815688 | rs5964999 | 0.475078 | 3.63212e-56 |    0.75886 | 5.04109e-14 |  -0.38849 | 5.39098e-23 |



```julia
# clean up
rm("trajgwas.null.txt", force=true)
rm("trajgwas.pval.txt", force=true)
rm("score.pval.txt", force=true)
rm("wald.pval.txt", force=true)
```

## GxE or other interactions

### Testing jointly G + GxE 

In many applications, we want to test SNP effect and/or its interaction with other terms. `testformula` keyword specifies the test unit **besides** the covariates in `nullformula`. 

In following example, keyword `testformula=@formula(trait ~ snp + snp & sex)` instructs `trajgwas` to test joint effect of `snp` and `snp & sex` interaction.


```julia
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = "GxE.pval.txt",
        testformula=@formula(trait ~ snp + snp & sex))
```

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.046905
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.043433





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 324
    Total observations: 3240
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  13.2282     0.146459    90.32    <1e-99
    β2: sex          -3.29295    0.2101     -15.67    <1e-54
    β3: onMeds        0.459585   0.0596002    7.71    <1e-13
    τ1: (Intercept)   0.792508   0.0850728    9.32    <1e-19
    τ2: sex          -0.2865     0.0970732   -2.95    0.0032
    τ3: onMeds        0.422303   0.063825     6.62    <1e-10
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  3.32057
    





```julia
pretty_table(first(CSV.read("GxE.pval.txt", DataFrame), 5), tf = tf_markdown)
```

    |[1m   chr [0m|[1m     pos [0m|[1m      snpid [0m|[1m       maf [0m|[1m    hwepval [0m|[1m    betapval [0m|[1m     taupval [0m|[1m   jointpval [0m|
    |[90m Int64 [0m|[90m   Int64 [0m|[90m   String31 [0m|[90m   Float64 [0m|[90m    Float64 [0m|[90m     Float64 [0m|[90m     Float64 [0m|[90m     Float64 [0m|
    |-------|---------|------------|-----------|------------|-------------|-------------|-------------|
    |     1 |  554484 | rs10458597 |       0.0 |        1.0 |         1.0 |         1.0 |         1.0 |
    |     1 |  758311 | rs12562034 | 0.0776398 |   0.409876 |       0.383 |     0.33511 |    0.357458 |
    |     1 |  967643 |  rs2710875 |  0.324074 | 4.07625e-7 |  9.98217e-6 |  5.73202e-7 |  1.08415e-6 |
    |     1 | 1168108 | rs11260566 |  0.191589 |   0.128568 |  0.00225789 |    0.105482 |  0.00442114 |
    |     1 | 1375074 |  rs1312568 |  0.441358 | 2.5376e-19 | 0.000928728 | 3.48807e-15 | 6.97615e-15 |



```julia
# clean up
rm("trajgwas.null.txt", force=true)
rm("GxE.pval.txt",  force=true)
```

### Testing only GxE interaction term

For some applications, the user may want to simply test the GxE interaction effect. This requires fitting the SNP in the null model and is much slower, but the command `trajgwas()` with keyword `analysistype = "gxe"` can be used test the interaction effect.
The environmental variable must be specified in the command using the keyword argument `e`, either as a symbol, such as `:age` or as a string `"age"`. 


```julia
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = "gxe_snp.pval.txt", 
        analysistype = "gxe",
        e = :sex, 
        snpinds=1:5)

```

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.045341
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.040165
    run = 1, ‖Δβ‖ = 0.480377, ‖Δτ‖ = 0.023510, ‖ΔL‖ = 0.001822, status = LOCALLY_SOLVED, time(s) = 0.045323
    run = 2, ‖Δβ‖ = 0.000323, ‖Δτ‖ = 0.002136, ‖ΔL‖ = 0.000025, status = LOCALLY_SOLVED, time(s) = 0.044147
    run = 1, ‖Δβ‖ = 1.059989, ‖Δτ‖ = 0.458239, ‖ΔL‖ = 0.063093, status = LOCALLY_SOLVED, time(s) = 0.060090
    run = 2, ‖Δβ‖ = 0.009118, ‖Δτ‖ = 0.105352, ‖ΔL‖ = 0.000056, status = LOCALLY_SOLVED, time(s) = 0.056690
    run = 1, ‖Δβ‖ = 1.144781, ‖Δτ‖ = 0.251459, ‖ΔL‖ = 0.034318, status = LOCALLY_SOLVED, time(s) = 0.047379
    run = 2, ‖Δβ‖ = 0.000830, ‖Δτ‖ = 0.025631, ‖ΔL‖ = 0.000042, status = LOCALLY_SOLVED, time(s) = 0.043087
    run = 1, ‖Δβ‖ = 0.683108, ‖Δτ‖ = 0.773919, ‖ΔL‖ = 0.042397, status = LOCALLY_SOLVED, time(s) = 0.044666
    run = 2, ‖Δβ‖ = 0.008159, ‖Δτ‖ = 0.069178, ‖ΔL‖ = 0.003574, status = LOCALLY_SOLVED, time(s) = 0.047489





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 324
    Total observations: 3240
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  13.2282     0.146459    90.32    <1e-99
    β2: sex          -3.29295    0.2101     -15.67    <1e-54
    β3: onMeds        0.459585   0.0596002    7.71    <1e-13
    τ1: (Intercept)   0.792508   0.0850728    9.32    <1e-19
    τ2: sex          -0.2865     0.0970732   -2.95    0.0032
    τ3: onMeds        0.422303   0.063825     6.62    <1e-10
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  3.32057
    





```julia
pretty_table(CSV.read("gxe_snp.pval.txt", DataFrame), tf = tf_markdown)
```

    |[1m   chr [0m|[1m     pos [0m|[1m      snpid [0m|[1m       maf [0m|[1m    hwepval [0m|[1m snpeffectnullbeta [0m|[1m snpeffectnulltau [0m|[1m betapval [0m|[1m  taupval [0m|
    |[90m Int64 [0m|[90m   Int64 [0m|[90m   String15 [0m|[90m   Float64 [0m|[90m    Float64 [0m|[90m           Float64 [0m|[90m          Float64 [0m|[90m  Float64 [0m|[90m  Float64 [0m|
    |-------|---------|------------|-----------|------------|-------------------|------------------|----------|----------|
    |     1 |  554484 | rs10458597 |       0.0 |        1.0 |               0.0 |              0.0 |      1.0 |      1.0 |
    |     1 |  758311 | rs12562034 | 0.0776398 |   0.409876 |         -0.230593 |        0.0102355 | 0.448376 | 0.149797 |
    |     1 |  967643 |  rs2710875 |  0.324074 | 4.07625e-7 |          0.639871 |        -0.338815 | 0.293876 | 0.119916 |
    |     1 | 1168108 | rs11260566 |  0.191589 |   0.128568 |          0.618855 |        -0.146819 | 0.590189 | 0.262791 |
    |     1 | 1375074 |  rs1312568 |  0.441358 | 2.5376e-19 |         -0.456003 |         0.492582 | 0.574737 | 0.250501 |



```julia
# clean up
rm("trajgwas.null.txt", force=true)
rm("gxe_snp.pval.txt", force=true)
```

## SNP-set testing

In many applications, we want to test a SNP-set. The function with keyword `analysistype = "snpset"` can be used to do this. To specify the type of snpset test, use the `snpset` keyword argument. 

The snpset can be specified as either:
- a window (test every X snps) => `snpset = X`
- an annotated file.  This requires `snpset = filename` where filename is an input file, with no header and two columns separated by a space. The first column must contain the snpset ID and the second column must contain the snpid's (rsid) identical to the bimfile, or in the case of a BGEN format, the order the data will be read (see BGEN above for more details).
- a joint test on only a specific set of snps. `snpset = AbstractVector` where the vector specifies the snps you want to perform one joint snpset test for. The vector can either be a vector of integers where the elements are the indicies of the SNPs to test, a vector of booleans, where true represents that you want to select that SNP index, or a range indicating the indicies of the SNPs to test. 

SNPset testing is currently only implemented for the score test (not Wald). 

In the following example, we perform a SNP-set test on the 50th to 55th snps. 


```julia
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = "snpset.pval.txt", 
        analysistype = "snpset",
        snpset = 50:55)
```

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043069
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.040745





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 324
    Total observations: 3240
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  13.2282     0.146459    90.32    <1e-99
    β2: sex          -3.29295    0.2101     -15.67    <1e-54
    β3: onMeds        0.459585   0.0596002    7.71    <1e-13
    τ1: (Intercept)   0.792508   0.0850728    9.32    <1e-19
    τ2: sex          -0.2865     0.0970732   -2.95    0.0032
    τ3: onMeds        0.422303   0.063825     6.62    <1e-10
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  3.32057
    





```julia
run(`head snpset.pval.txt`);
```

    The pvalue of snps indexed at 50:55 is betapval: 2.315617121196509e-5, taupval: 3.972993706371914e-10



```julia
# clean up
rm("snpset.pval.txt", force=true)
rm("trajgwas.null.txt", force=true)
```

In the following example we run a SNP-set test on the annotated SNP-set file.


```julia
run(`head ../data/hapmap_snpsetfile.txt`);
```

    gene1 rs10458597
    gene1 rs12562034
    gene1 rs2710875
    gene1 rs11260566
    gene1 rs1312568
    gene1 rs35154105
    gene1 rs16824508
    gene1 rs2678939
    gene1 rs7553178
    gene1 rs13376356



```julia
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = "snpset.pval.txt", 
        analysistype = "snpset",
        snpset = datadir * "/hapmap_snpsetfile.txt")
```

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044896
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.040974





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 324
    Total observations: 3240
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  13.2282     0.146459    90.32    <1e-99
    β2: sex          -3.29295    0.2101     -15.67    <1e-54
    β3: onMeds        0.459585   0.0596002    7.71    <1e-13
    τ1: (Intercept)   0.792508   0.0850728    9.32    <1e-19
    τ2: sex          -0.2865     0.0970732   -2.95    0.0032
    τ3: onMeds        0.422303   0.063825     6.62    <1e-10
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  3.32057
    





```julia
pretty_table(first(CSV.read("snpset.pval.txt", DataFrame; delim="\t"), 5), tf = tf_markdown)
```

    |[1m snpsetid [0m|[1m nsnps [0m|[1m  betapval [0m|[1m   taupval [0m|
    |[90m  String7 [0m|[90m Int64 [0m|[90m   Float64 [0m|[90m   Float64 [0m|
    |----------|-------|-----------|-----------|
    |    gene1 |    93 |  0.111864 |  0.011595 |
    |    gene2 |    93 | 0.0249929 | 0.0648265 |
    |    gene3 |    93 |  0.131741 | 0.0298884 |
    |    gene4 |    92 |  0.010327 | 0.0584591 |
    |    gene5 |    93 | 0.0303905 | 0.0924954 |



```julia
# clean up
rm("snpset.pval.txt", force=true)
rm("trajgwas.null.txt", force=true)
```

In the following example we run a SNP-set test on every 15 SNPs.


```julia
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = "snpset.pval.txt", 
        analysistype = "snpset",
        snpset = 15)
```

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043705
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.039385





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 324
    Total observations: 3240
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  13.2282     0.146459    90.32    <1e-99
    β2: sex          -3.29295    0.2101     -15.67    <1e-54
    β3: onMeds        0.459585   0.0596002    7.71    <1e-13
    τ1: (Intercept)   0.792508   0.0850728    9.32    <1e-19
    τ2: sex          -0.2865     0.0970732   -2.95    0.0032
    τ3: onMeds        0.422303   0.063825     6.62    <1e-10
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  3.32057
    





```julia
pretty_table(first(CSV.read("snpset.pval.txt", DataFrame), 5), tf = tf_markdown)
```

    |[1m startchr [0m|[1m startpos [0m|[1m startsnpid [0m|[1m endchr [0m|[1m   endpos [0m|[1m   endsnpid [0m|[1m    betapval [0m|[1m     taupval [0m|
    |[90m    Int64 [0m|[90m    Int64 [0m|[90m   String15 [0m|[90m  Int64 [0m|[90m    Int64 [0m|[90m   String15 [0m|[90m     Float64 [0m|[90m     Float64 [0m|
    |----------|----------|------------|--------|----------|------------|-------------|-------------|
    |        1 |   554484 | rs10458597 |      1 |  3431124 | rs12093117 |  3.83046e-6 | 2.24303e-10 |
    |        1 |  3633945 | rs10910017 |      1 |  6514524 |   rs932112 | 0.000127922 |  7.94764e-6 |
    |        1 |  6715827 |   rs441515 |      1 |  9534606 |  rs4926480 |  8.10353e-5 |  6.11848e-7 |
    |        1 |  9737551 | rs12047054 |      1 | 12559747 |  rs4845907 | 0.000425108 |  1.15311e-8 |
    |        1 | 12760427 |   rs848577 |      1 | 16021797 |  rs6679870 |  2.75832e-5 |  0.00028553 |



```julia
# clean up
rm("snpset.pval.txt", force=true)
rm("trajgwas.null.txt", force=true)
```

## Matching Indicies

In some cases, there are only a subset of individuals with both genetic data and covariate information available. The null model must be fit on a subset of the individuals with the genetic data. The rows can be specified with the argument `covrowinds` if you pass in a covariate file. The genetic indicies can be specified with the `geneticrowinds` argument. 

For simplicity, we have implemented a function `matchindices(meanformula, reformula, wsvarformula, idvar, df, geneticsamples)` which can be used to do this.
Input the mean, random effects, and within-subject variance formulas, the grouping (id) variable,
the dataframe (or table), and a vector of the IDs in the order of the genetic file and it returns `covrowmask, geneticrowmask` for matching indicies in a covariate file and geneticfile.

Note: the `idvar` in the dataframe and the `geneticsamples` vector must have the same element type. To convert a vector of integers to strings `string.(intvector)` can be used. To go from a vector of strings to integers you can use"  `map(x -> parse(Int, x), stringvectoparse)`.


```julia
famfileids = CSV.read(datadir * "hapmap3.fam", DataFrame, header = false)[!, 1] # famfile contains the sample IDs for PLINK files
```




    324-element Vector{String7}:
     "A1"
     "2"
     "3"
     "4"
     "5"
     "6"
     "7"
     "8"
     "9"
     "10"
     "11"
     "12"
     "13"
     ⋮
     "313"
     "314"
     "315"
     "316"
     "317"
     "318"
     "319"
     "320"
     "321"
     "322"
     "323"
     "Z324"




```julia
covdf = CSV.read(datadir * "trajgwas_plinkex.csv", DataFrame)
pretty_table(first(covdf, 11), tf = tf_markdown)
```

    |[1m     sex [0m|[1m  onMeds [0m|[1m    snp1 [0m|[1m    snp2 [0m|[1m    snp3 [0m|[1m    snp4 [0m|[1m       y [0m|[1m      id [0m|
    |[90m Float64 [0m|[90m Float64 [0m|[90m Float64 [0m|[90m Float64 [0m|[90m Float64 [0m|[90m Float64 [0m|[90m Float64 [0m|[90m String7 [0m|
    |---------|---------|---------|---------|---------|---------|---------|---------|
    |     0.0 |     1.0 |     0.0 |     1.0 |     2.0 |     0.0 | 12.2667 |      A1 |
    |     0.0 |     0.0 |     0.0 |     1.0 |     2.0 |     0.0 | 10.2681 |      A1 |
    |     0.0 |     0.0 |     0.0 |     1.0 |     2.0 |     0.0 |  12.166 |      A1 |
    |     0.0 |     0.0 |     0.0 |     1.0 |     2.0 |     0.0 | 11.8797 |      A1 |
    |     0.0 |     0.0 |     0.0 |     1.0 |     2.0 |     0.0 | 12.8127 |      A1 |
    |     0.0 |     0.0 |     0.0 |     1.0 |     2.0 |     0.0 | 9.98766 |      A1 |
    |     0.0 |     0.0 |     0.0 |     1.0 |     2.0 |     0.0 | 12.1408 |      A1 |
    |     0.0 |     1.0 |     0.0 |     1.0 |     2.0 |     0.0 | 13.2058 |      A1 |
    |     0.0 |     0.0 |     0.0 |     1.0 |     2.0 |     0.0 | 11.3631 |      A1 |
    |     0.0 |     1.0 |     0.0 |     1.0 |     2.0 |     0.0 | 15.2511 |      A1 |
    |     0.0 |     0.0 |     0.0 |     0.0 |     2.0 |     2.0 |  12.746 |       2 |



```julia
covrowmask, geneticrowmask = matchindices(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id, covdf, famfileids);
```


```julia
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        datadir * "hapmap3",
        pvalfile = pvalpath,
        nullfile = nullpath,
        covrowinds = covrowmask,
        geneticrowinds = geneticrowmask)
```

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.043990
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.040579





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 324
    Total observations: 3240
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  13.2282     0.146459    90.32    <1e-99
    β2: sex          -3.29295    0.2101     -15.67    <1e-54
    β3: onMeds        0.459585   0.0596002    7.71    <1e-13
    τ1: (Intercept)   0.792508   0.0850728    9.32    <1e-19
    τ2: sex          -0.2865     0.0970732   -2.95    0.0032
    τ3: onMeds        0.422303   0.063825     6.62    <1e-10
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  3.32057
    





```julia
# clean up
rm("trajgwas.null.txt", force=true)
rm("trajgwas.pval.txt", force=true)
```

## Plotting Results

To plot the GWAS results, we recommend using the [MendelPlots.jl package](https://openmendel.github.io/MendelPlots.jl/latest/).

## Multiple Plink file sets

In large scale studies, genotypes data are split into multiple Plink files, e.g., by chromosome. Then GWAS analysis can be done in parallel. This can be achieved by two steps.

Let's first create demo data by splitting hapmap3 according to chromosome:


```julia
# split example hapmap3 data according to chromosome
SnpArrays.split_plink(datadir * "hapmap3", :chromosome; prefix=datadir * "hapmap3.chr.")
readdir(glob"hapmap3.chr.*", datadir)
```




    75-element Vector{String}:
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.1.bed"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.1.bim"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.1.fam"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.10.bed"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.10.bim"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.10.fam"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.11.bed"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.11.bim"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.11.fam"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.12.bed"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.12.bim"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.12.fam"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.13.bed"
     ⋮
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.6.bed"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.6.bim"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.6.fam"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.7.bed"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.7.bim"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.7.fam"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.8.bed"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.8.bim"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.8.fam"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.9.bed"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.9.bim"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.9.fam"



Step 1: Fit the null model. Setting third argument `geneticfile` to `nothing` instructs `trajgwas` function to fit the null model only.


```julia
nm = trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        nothing)
```

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = LOCALLY_SOLVED, time(s) = 0.044113
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = LOCALLY_SOLVED, time(s) = 0.038771





    
    Within-subject variance estimation by robust regression (WiSER)
    
    Mean Formula:
    y ~ 1 + sex + onMeds
    Random Effects Formula:
    y ~ 1
    Within-Subject Variance Formula:
    y ~ 1 + sex + onMeds
    
    Number of individuals/clusters: 324
    Total observations: 3240
    
    Fixed-effects parameters:
    ────────────────────────────────────────────────────────
                      Estimate  Std. Error       Z  Pr(>|Z|)
    ────────────────────────────────────────────────────────
    β1: (Intercept)  13.2282     0.146459    90.32    <1e-99
    β2: sex          -3.29295    0.2101     -15.67    <1e-54
    β3: onMeds        0.459585   0.0596002    7.71    <1e-13
    τ1: (Intercept)   0.792508   0.0850728    9.32    <1e-19
    τ2: sex          -0.2865     0.0970732   -2.95    0.0032
    τ3: onMeds        0.422303   0.063825     6.62    <1e-10
    ────────────────────────────────────────────────────────
    Random effects covariance matrix Σγ:
     "γ1: (Intercept)"  3.32057
    




Step 2: GWAS for each chromosome.


```julia
# this part can be submitted as separate jobs
for chr in 1:23
    plinkfile = datadir * "hapmap3.chr." * string(chr)
    pvalfile = plinkfile * ".pval.txt" 
    trajgwas(nm, plinkfile, pvalfile=pvalfile)
end
```


```julia
# show the result files
readdir(glob"*.pval.txt", datadir)
```




    23-element Vector{String}:
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.1.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.10.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.11.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.12.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.13.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.14.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.15.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.16.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.17.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.18.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.19.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.2.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.20.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.21.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.22.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.23.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.3.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.4.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.5.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.6.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.7.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.8.pval.txt"
     "/Users/xyz/.julia/dev/TrajGWAS/data/hapmap3.chr.9.pval.txt"



In the rare situations where the multiple sets of Plink files lack the `fam` file or the corresponding bed and bim files have different filenames, users can explicitly supply bed filename, bim file name, and number of individuals. Replace Step 2 by 

Step 2': GWAS for each chromosome.


```julia
# this part can be submitted as separate jobs
for chr in 1:23
    bedfile = datadir * "hapmap3.chr." * string(chr) * ".bed"
    bimfile = datadir * "hapmap3.chr." * string(chr) * ".bim"
    pvalfile = datadir * "hapmap3.chr." * string(chr) * ".pval.txt"
    trajgwas(nm, bedfile, bimfile, 324; pvalfile=pvalfile)
end
```


```julia
# clean up
rm("trajgwas.null.txt", force=true)
isfile(datadir * "fittednullmodel.jld2") && rm(datadir * "fittednullmodel.jld2")
for chr in 1:23
    pvalfile = datadir * "hapmap3.chr." * string(chr) * ".pval.txt"
    rm(pvalfile, force=true)
end
for chr in 1:26
    plinkfile = datadir * "hapmap3.chr." * string(chr)
    rm(plinkfile * ".bed", force=true)
    rm(plinkfile * ".fam", force=true)
    rm(plinkfile * ".bim", force=true)
end
```

## Multiple file sets on cluster

For running the score tests on a cluster, it would be desirable to fit a null model on a single machine first, and save the model as a serialized Julia object (`.jls`). For example:

```julia
using DataFrames, CSV
using Statistics
using TrajGWAS
using Ipopt, WiSER
using LinearAlgebra
using BGEN
# fit the null model
BLAS.set_num_threads(1)
solver = Ipopt.Optimizer()
solver_config = Dict("print_level" => 1, 
    "watchdog_shortened_iter_trigger" => 5,
    "max_iter" => 120)

bp_data = CSV.read("bp.csv", DataFrame)
@time nm = trajgwas(@formula(dbp ~ 1 + SEX + age + age_sq +
        PC1 + PC2 + PC3 + PC4 + PC5 + bmi),
    @formula(dbp ~ 1 + age),
    @formula(dbp ~ 1 + SEX + age + age_sq +
        PC1 + PC2 + PC3 + PC4 + PC5 +
        bmi),
    :FID, # subject ID
    bp_data,
    nothing;
    nullfile="dbp.null.txt",
    solver=solver,
    solver_config=solver_config,
    init = nothing, # may change to `x -> WiSER.init_ls!(x; gniters=0)`
    runs=10
)

println(nm)
using Serialization
open("null.model.jls", "w") do io
    Serialization.serialize(io, nm)
end
```

Then, the fitted null model would be used for the score test. It could be desirable to configure each job to run on a slice of SNPs on a BGEN file for higher throughput. The julia script for score test (`scoretest_bp.jl`) for a slice of BGEN file would look like:

```julia
using DataFrames, CSV
using Statistics
using TrajGWAS
using WiSER
using LinearAlgebra
using BGEN
# fit the null model
BLAS.set_num_threads(1)


using Serialization

bgendir = ARGS[1] 
chr = ARGS[2] # 1 to 22
fitted_null = ARGS[3] # "null.model.jls"
pvalfile = ARGS[4] # "sbp.test.diabetics.chr$(chr).txt"
chunkidx = parse(Int, ARGS[5])
nchunks  = parse(Int, ARGS[6])

nm = open(deserialize, fitted_null)
genetic_iids_subsample = nm.ids

bgenfilename = bgendir * "/Chr$(chr)" # to analyze, for example, /path/to/bgen/Chr5.bgen
samplefilename = bgendir * "/Samples.sample" # .sample file compatible with BGEN
mfifilename = bgendir * "/mfi_chr$(chr).txt" # "MFI" file, an external file with MAF + info score	
ukb_data = Bgen(bgenfilename * ".bgen"; sample_path = samplefilename)
genetic_iids = map(x -> parse(Int, split(x, " ")[1]), samples(ukb_data))

order_dict = Dict{Int, Int}()
for (i, iid) in enumerate(genetic_iids)
    order_dict[iid] = i
end

sample_indicator = falses(length(genetic_iids))
for v in genetic_iids_subsample
    sample_indicator[order_dict[v]] = true # extract only the samples being used for the analysis
end

# GWAS for each chromosome

## pre-filtering SNPs not passing the criteria (MAF > 0.002, info score > 0.3)
min_maf = 0.002
min_info_score = 0.3
min_hwe_pval = 1e-10
mfi = CSV.read(mfifilename, DataFrame; header=false)
mfi.Column8 = map(x -> x == "NA" ? NaN : parse(Float64, x), mfi.Column8) # Column8: info score
snpmask = (mfi.Column6 .> min_maf) .& (mfi.Column8 .> 0.3) # Column6: MAF

## compute range to run the analysis
chunksize = n_variants(ukb_data) ÷ nchunks + (n_variants(ukb_data) % nchunks > 0 ? 1 : 0)
startidx = chunksize * (chunkidx - 1) + 1
endidx = min(chunksize * chunkidx, n_variants(ukb_data))
snpmask = snpmask[startidx:endidx]

println("running for variants $startidx to $endidx")

# rearrange data in nm so that it matches bgen data
nullinds = indexin(genetic_iids[sample_indicator], nm.ids)
nm.obswts .= isempty(nm.obswts) ? nm.obswts : nm.obswts[nullinds]
nm.ids .= nm.ids[nullinds]
nm.nis .= nm.nis[nullinds]
nm.data .= nm.data[nullinds]
@assert genetic_iids[sample_indicator] == nm.ids "there is some issue -- sampleids not matching"
    
trajgwas(nm, bgenfilename * ".bgen", count(sample_indicator);
    samplepath=samplefilename,
    pvalfile=pvalfile,
    snpinds=snpmask,
    min_hwe_pval = min_hwe_pval,
    bgenrowinds = sample_indicator,
    startidx = startidx,
    endidx = endidx,
    usespa=true)

```
Note that an index file (`.bgen.bgi`) is required for this slicing of BGEN file. See [this link](https://enkre.net/cgi-bin/code/bgen/doc/trunk/doc/wiki/bgenix.md) to see how to create one. 


- Command-line arguments
    - Argument 1: directory for the BGEN files. BGEN files (.bgen), BGEN index files (.bgen.bgi), and MFI files (.txt) should be included there.
    - Argument 2: chromosome
    - Argument 3: fitted null model (.jls)
    - Argument 4: path for the result p-value file
    - Argument 5: chunk index (1-based)
    - Argument 6: number of chunks
    
The code above runs the analysis on `ARGS[5]`-th slice out of `ARGS[6]` slices of chromosome `ARGS[2]`.

Then, the following script could be used for a cluster managed by Sun Grid Engine: (`sbp_diabetes.sh`)

```bash
#!/bin/bash
#$ -cwd
#$ -o joblog.$JOB_ID.$TASK_ID
#$ -j y
#$ -pe shared 2
#$ -l h_rt=8:00:00,h_data=8G,arch=intel*
# Email address to notify
##$ -M $USER@mail
# Notify when
#$ -m a
#  Job array indexes
#$ -t 1-352:1

NCHUNKS=16
CHUNKIDX=$(( (${SGE_TASK_ID} - 1) % ${NCHUNKS} + 1 ))
CHR=$(( (${SGE_TASK_ID} - 1) / ${NCHUNKS} + 1))

PROJECTDIR=/user/dir/jobscripts
BGENDIR=/user/dir/imputed
FITTED_NULL=/user/dir/null.model.jls
PVALFILE=/user/dir/pvalfiles/sbp.test.diabetes.chr${CHR}.${CHUNKIDX}of${NCHUNKS}.txt

module load julia
time julia --project=${PROJECTDIR} ${PROJECTDIR}/scoretest_bp.jl ${BGENDIR} ${CHR} ${FITTED_NULL} ${PVALFILE} ${CHUNKIDX} ${NCHUNKS}
```

and this could be submitted using the `qsub` command.

## Troubleshooting

If there are issues you're encountering with running TrajGWAS, the following are possible remedies.

- Null Model
    - Issues with null model convergence may be solved by choosing alternate starting values for parameters, using a different solver, transforming variables, and increasing the number of runs (WiSER runs). These are detailed in the [WiSER documentation here](https://github.com/OpenMendel/WiSER.jl/blob/master/docs/src/model_fitting.md#tips-for-improving-estimation). `init = x -> WiSER.init_ls!(x; gniters=0)` works fine in general. 
    
- GWAS Results
    - If you use the score test instead of the SPA-score test (SPA is default for single-SNP analyses), then there can be inflation in type I error and decreased power when (a) the sample size is small, (b) the number of repeated measures is low, or (c) the variants analyzed are rare with low minor allele frequencies. In these cases, the score test is not optimal and it is suggested to use the SPA version (`usespa=true`). SPA is only implemented for single-SNP analyses. These issues can occur in both Wald and score tests. 
    
If you notice any problems with your output or results, [file an issue](https://github.com/OpenMendel/TrajGWAS.jl/issues). 
