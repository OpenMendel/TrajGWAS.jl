# TrajGWAS.jl

TrajGWAS.jl is a Julia package for performing genome-wide association studies (GWAS) for continuous longitudinal phenotypes using a modified linear mixed effects model. It builds upon the [within-subject variance estimation by robust regression (WiSER) method](https://github.com/OpenMendel/WiSER.jl) and can be used to identify variants associated with changes in the mean and within-subject variability of the longitduinal trait. The estimation procedure is robust to distributional misspecifications of both the random effects and the response. A saddlepoint approximation (SPA) option is implemented to provide improved power and calibrated type I error for rare variants. 


TrajGWAS.jl currently supports [PLINK](https://zzz.bwh.harvard.edu/plink/), [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) (both dosage and genotype data), and [BGEN](https://www.well.ox.ac.uk/~gav/bgen_format/) file formats. We plan to add [PGEN](https://www.cog-genomics.org/plink/2.0/formats#pgen) support in the future. 

## Installation

This package requires Julia v1.5 or later and four other unregistered packages SnpArrays.jl, VCFTools.jl, BGEN.jl, and WiSER.jl. The package has not yet been registered and must be installed using the repository location. Execute the following code to install the package:

```julia
using Pkg
pkg"add https://github.com/OpenMendel/SnpArrays.jl"
pkg"add https://github.com/OpenMendel/VCFTools.jl"
pkg"add https://github.com/OpenMendel/BGEN.jl"
pkg"add https://github.com/OpenMendel/WiSER.jl"
pkg"add https://github.com/OpenMendel/TrajGWAS.jl" 
```

To run the code in this document, the packages installed by the following command are also necessary:
```julia
pkg"add BenchmarkTools CSV Glob"
```


```julia
# machine information for this tutorial
versioninfo()
```

    Julia Version 1.6.0
    Commit f9720dc2eb (2021-03-24 12:55 UTC)
    Platform Info:
      OS: macOS (x86_64-apple-darwin19.6.0)
      CPU: Intel(R) Core(TM) i7-6920HQ CPU @ 2.90GHz
      WORD_SIZE: 64
      LIBM: libopenlibm
      LLVM: libLLVM-11.0.1 (ORCJIT, skylake)
    Environment:
      JULIA_EDITOR = code
      JULIA_NUM_THREADS = 4



```julia
# for use in this tutorial
ENV["COLUMNS"] = 250
using BenchmarkTools, CSV, Glob, SnpArrays, TrajGWAS
```

## Example data sets

The `data` folder of the package contains the example data sets for use with PLINK and VCF Files. In general, the user can locate this folder by command:


```julia
pvalpath = "trajgwas.pval.txt"
nullpath = "trajgwas.null.txt"
const datadir = normpath(joinpath(dirname(pathof(TrajGWAS)), "../data/"))
```




    "/Users/huazhou/.julia/dev/TrajGWAS/data/"




```julia
# content of the data folder
readdir(glob"*.*", datadir)
```




    14-element Vector{String}:
     "/Users/huazhou/.julia/dev/TrajGWAS/data/bgen_snpsetfile.txt"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/covariate.txt"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/example.8bits.bgen"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/example.8bits.bgen.bgi"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/hapmap3.bed"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/hapmap3.bim"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/hapmap3.fam"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/hapmap_snpsetfile.txt"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/sim_data.jl"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/snpsetfile_vcf.txt"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/test_vcf.vcf.gz"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/trajgwas_bgen_ex.csv"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/trajgwas_plinkex.csv"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/trajgwas_vcfex.csv"



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
    
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.259617
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.034867





    
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
     "/Users/huazhou/.julia/dev/TrajGWAS/data/hapmap3.bed"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/hapmap3.bim"
     "/Users/huazhou/.julia/dev/TrajGWAS/data/hapmap3.fam"



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
first(CSV.read("trajgwas.pval.txt", DataFrame), 8)
```




<div class="data-frame"><p>8 rows × 10 columns</p><table class="data-frame"><thead><tr><th></th><th>chr</th><th>pos</th><th>snpid</th><th>maf</th><th>hwepval</th><th>betapval</th><th>betadir</th><th>taupval</th><th>taudir</th><th>jointpval</th></tr><tr><th></th><th title="Int64">Int64</th><th title="Int64">Int64</th><th title="String">String</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Int64">Int64</th><th title="Float64">Float64</th><th title="Int64">Int64</th><th title="Float64">Float64</th></tr></thead><tbody><tr><th>1</th><td>1</td><td>554484</td><td>rs10458597</td><td>0.0</td><td>1.0</td><td>1.0</td><td>0</td><td>1.0</td><td>0</td><td>1.0</td></tr><tr><th>2</th><td>1</td><td>758311</td><td>rs12562034</td><td>0.0776398</td><td>0.409876</td><td>0.38939</td><td>-1</td><td>0.922822</td><td>-1</td><td>0.547682</td></tr><tr><th>3</th><td>1</td><td>967643</td><td>rs2710875</td><td>0.324074</td><td>4.07625e-7</td><td>5.83856e-6</td><td>1</td><td>3.35408e-6</td><td>-1</td><td>4.93598e-7</td></tr><tr><th>4</th><td>1</td><td>1168108</td><td>rs11260566</td><td>0.191589</td><td>0.128568</td><td>0.000850889</td><td>1</td><td>0.0559055</td><td>1</td><td>0.00101827</td></tr><tr><th>5</th><td>1</td><td>1375074</td><td>rs1312568</td><td>0.441358</td><td>2.5376e-19</td><td>0.000171683</td><td>-1</td><td>1.84811e-13</td><td>1</td><td>2.07091e-15</td></tr><tr><th>6</th><td>1</td><td>1588771</td><td>rs35154105</td><td>0.0</td><td>1.0</td><td>1.0</td><td>0</td><td>1.0</td><td>0</td><td>1.0</td></tr><tr><th>7</th><td>1</td><td>1789051</td><td>rs16824508</td><td>0.00462963</td><td>0.933278</td><td>0.295035</td><td>1</td><td>0.304109</td><td>1</td><td>0.299504</td></tr><tr><th>8</th><td>1</td><td>1990452</td><td>rs2678939</td><td>0.453704</td><td>5.07696e-11</td><td>1.27871e-7</td><td>1</td><td>7.73809e-9</td><td>-1</td><td>3.3986e-9</td></tr></tbody></table></div>



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
covdf
```




<div class="data-frame"><p>324 rows × 6 columns</p><table class="data-frame"><thead><tr><th></th><th>famid</th><th>perid</th><th>faid</th><th>moid</th><th>sex</th><th>trait</th></tr><tr><th></th><th title="String">String</th><th title="String">String</th><th title="Int64">Int64</th><th title="Int64">Int64</th><th title="Int64">Int64</th><th title="Int64">Int64</th></tr></thead><tbody><tr><th>1</th><td>2431</td><td>NA19916</td><td>0</td><td>0</td><td>1</td><td>4</td></tr><tr><th>2</th><td>2424</td><td>NA19835</td><td>0</td><td>0</td><td>2</td><td>4</td></tr><tr><th>3</th><td>2469</td><td>NA20282</td><td>0</td><td>0</td><td>2</td><td>4</td></tr><tr><th>4</th><td>2368</td><td>NA19703</td><td>0</td><td>0</td><td>1</td><td>3</td></tr><tr><th>5</th><td>2425</td><td>NA19901</td><td>0</td><td>0</td><td>2</td><td>3</td></tr><tr><th>6</th><td>2427</td><td>NA19908</td><td>0</td><td>0</td><td>1</td><td>4</td></tr><tr><th>7</th><td>2430</td><td>NA19914</td><td>0</td><td>0</td><td>2</td><td>4</td></tr><tr><th>8</th><td>2470</td><td>NA20287</td><td>0</td><td>0</td><td>2</td><td>1</td></tr><tr><th>9</th><td>2436</td><td>NA19713</td><td>0</td><td>0</td><td>2</td><td>3</td></tr><tr><th>10</th><td>2426</td><td>NA19904</td><td>0</td><td>0</td><td>1</td><td>1</td></tr><tr><th>11</th><td>2431</td><td>NA19917</td><td>0</td><td>0</td><td>2</td><td>1</td></tr><tr><th>12</th><td>2436</td><td>NA19982</td><td>0</td><td>0</td><td>1</td><td>2</td></tr><tr><th>13</th><td>2487</td><td>NA20340</td><td>0</td><td>0</td><td>1</td><td>4</td></tr><tr><th>14</th><td>2427</td><td>NA19909</td><td>0</td><td>0</td><td>2</td><td>4</td></tr><tr><th>15</th><td>2424</td><td>NA19834</td><td>0</td><td>0</td><td>1</td><td>4</td></tr><tr><th>16</th><td>2480</td><td>NA20317</td><td>0</td><td>0</td><td>2</td><td>4</td></tr><tr><th>17</th><td>2418</td><td>NA19818</td><td>0</td><td>0</td><td>1</td><td>1</td></tr><tr><th>18</th><td>2490</td><td>NA20346</td><td>0</td><td>0</td><td>1</td><td>2</td></tr><tr><th>19</th><td>2433</td><td>NA19921</td><td>0</td><td>0</td><td>2</td><td>4</td></tr><tr><th>20</th><td>2469</td><td>NA20281</td><td>0</td><td>0</td><td>1</td><td>4</td></tr><tr><th>21</th><td>2495</td><td>NA20359</td><td>0</td><td>0</td><td>2</td><td>4</td></tr><tr><th>22</th><td>2477</td><td>NA20301</td><td>0</td><td>0</td><td>2</td><td>2</td></tr><tr><th>23</th><td>2492</td><td>NA20349</td><td>0</td><td>0</td><td>1</td><td>3</td></tr><tr><th>24</th><td>2474</td><td>NA20294</td><td>0</td><td>0</td><td>2</td><td>4</td></tr><tr><th>25</th><td>2494</td><td>NA20357</td><td>0</td><td>0</td><td>2</td><td>3</td></tr><tr><th>26</th><td>2425</td><td>NA19900</td><td>0</td><td>0</td><td>1</td><td>4</td></tr><tr><th>27</th><td>2491</td><td>NA20348</td><td>0</td><td>0</td><td>1</td><td>4</td></tr><tr><th>28</th><td>2471</td><td>NA20289</td><td>0</td><td>0</td><td>2</td><td>4</td></tr><tr><th>29</th><td>2489</td><td>NA20344</td><td>0</td><td>0</td><td>2</td><td>4</td></tr><tr><th>30</th><td>2418</td><td>NA19819</td><td>0</td><td>0</td><td>2</td><td>4</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table></div>




```julia
plkfam
```




<div class="data-frame"><p>324 rows × 6 columns</p><table class="data-frame"><thead><tr><th></th><th>Column1</th><th>Column2</th><th>Column3</th><th>Column4</th><th>Column5</th><th>Column6</th></tr><tr><th></th><th title="String">String</th><th title="String">String</th><th title="Int64">Int64</th><th title="Int64">Int64</th><th title="Int64">Int64</th><th title="Int64">Int64</th></tr></thead><tbody><tr><th>1</th><td>A1</td><td>NA19916</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr><tr><th>2</th><td>2</td><td>NA19835</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>3</th><td>3</td><td>NA20282</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>4</th><td>4</td><td>NA19703</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr><tr><th>5</th><td>5</td><td>NA19901</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>6</th><td>6</td><td>NA19908</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr><tr><th>7</th><td>7</td><td>NA19914</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>8</th><td>8</td><td>NA20287</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>9</th><td>9</td><td>NA19713</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>10</th><td>10</td><td>NA19904</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr><tr><th>11</th><td>11</td><td>NA19917</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>12</th><td>12</td><td>NA19982</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr><tr><th>13</th><td>13</td><td>NA20340</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr><tr><th>14</th><td>14</td><td>NA19909</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>15</th><td>15</td><td>NA19834</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr><tr><th>16</th><td>16</td><td>NA20317</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>17</th><td>17</td><td>NA19818</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr><tr><th>18</th><td>18</td><td>NA20346</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr><tr><th>19</th><td>19</td><td>NA19921</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>20</th><td>20</td><td>NA20281</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr><tr><th>21</th><td>21</td><td>NA20359</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>22</th><td>22</td><td>NA20301</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>23</th><td>23</td><td>NA20349</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr><tr><th>24</th><td>24</td><td>NA20294</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>25</th><td>25</td><td>NA20357</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>26</th><td>26</td><td>NA19900</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr><tr><th>27</th><td>27</td><td>NA20348</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr><tr><th>28</th><td>28</td><td>NA20289</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>29</th><td>29</td><td>NA20344</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>30</th><td>30</td><td>NA19819</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table></div>



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

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.042119
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.036129
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.050494
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.051440
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.046484
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.034452
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.048439
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.050991
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.037938
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.034008
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.050020
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.040776
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.049833
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.043173
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.037266
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.032848
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.062879
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.038666
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.037709
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.034644
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.047753
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.051048
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.036596
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.035788
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.057062
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.041374
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.040622
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.035942
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.053378
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.039402
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.041604
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.035013
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.055528
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.039201
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.048910
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.036422
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.050033
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.044553
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.040895
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.041951
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.060043
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.049582
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.041068
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.038690
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.055143
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.041776
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.036662
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.034396
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.056156
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.038946
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.037163
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.032339
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.058409
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.038301
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.041117
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.036138
    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.054768
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.037206
      406.061 ms (1450847 allocations: 136.02 MiB)



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

    run = 1, ‖Δβ‖ = 0.003459, ‖Δτ‖ = 0.421747, ‖ΔL‖ = 0.002492, status = Optimal, time(s) = 0.017842
    run = 2, ‖Δβ‖ = 0.001369, ‖Δτ‖ = 0.015998, ‖ΔL‖ = 0.003244, status = Optimal, time(s) = 0.014044





    
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
first(CSV.read("trajgwas.pval.txt", DataFrame), 8)
```




<div class="data-frame"><p>8 rows × 8 columns</p><table class="data-frame"><thead><tr><th></th><th>chr</th><th>pos</th><th>snpid</th><th>betapval</th><th>betadir</th><th>taupval</th><th>taudir</th><th>jointpval</th></tr><tr><th></th><th title="Int64">Int64</th><th title="Int64">Int64</th><th title="String">String</th><th title="Float64">Float64</th><th title="Int64">Int64</th><th title="Float64">Float64</th><th title="Int64">Int64</th><th title="Float64">Float64</th></tr></thead><tbody><tr><th>1</th><td>22</td><td>20000086</td><td>rs138720731</td><td>0.322684</td><td>-1</td><td>0.61956</td><td>-1</td><td>0.424353</td></tr><tr><th>2</th><td>22</td><td>20000146</td><td>rs73387790</td><td>1.0</td><td>0</td><td>1.0</td><td>0</td><td>1.0</td></tr><tr><th>3</th><td>22</td><td>20000199</td><td>rs183293480</td><td>0.214207</td><td>1</td><td>0.190785</td><td>1</td><td>0.201819</td></tr><tr><th>4</th><td>22</td><td>20000291</td><td>rs185807825</td><td>0.292231</td><td>1</td><td>0.200029</td><td>1</td><td>0.237495</td></tr><tr><th>5</th><td>22</td><td>20000428</td><td>rs55902548</td><td>0.0140114</td><td>1</td><td>0.00369005</td><td>1</td><td>0.00576184</td></tr><tr><th>6</th><td>22</td><td>20000683</td><td>rs142720028</td><td>1.0</td><td>0</td><td>1.0</td><td>0</td><td>1.0</td></tr><tr><th>7</th><td>22</td><td>20000771</td><td>rs114690707</td><td>0.536406</td><td>-1</td><td>0.145113</td><td>1</td><td>0.22843</td></tr><tr><th>8</th><td>22</td><td>20000793</td><td>rs189842693</td><td>0.161414</td><td>-1</td><td>0.757017</td><td>1</td><td>0.266091</td></tr></tbody></table></div>



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

    run = 1, ‖Δβ‖ = 0.096578, ‖Δτ‖ = 0.162563, ‖ΔL‖ = 0.007625, status = Optimal, time(s) = 0.058796
    run = 2, ‖Δβ‖ = 0.003173, ‖Δτ‖ = 0.005584, ‖ΔL‖ = 0.001468, status = Optimal, time(s) = 0.046212





    
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
first(CSV.read("trajgwas.pval.txt", DataFrame), 8)
```




<div class="data-frame"><p>8 rows × 12 columns</p><table class="data-frame"><thead><tr><th></th><th>chr</th><th>pos</th><th>snpid</th><th>varid</th><th>hwepval</th><th>maf</th><th>infoscore</th><th>betapval</th><th>betadir</th><th>taupval</th><th>taudir</th><th>jointpval</th></tr><tr><th></th><th title="Int64">Int64</th><th title="Int64">Int64</th><th title="String">String</th><th title="String">String</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Int64">Int64</th><th title="Float64">Float64</th><th title="Int64">Int64</th><th title="Float64">Float64</th></tr></thead><tbody><tr><th>1</th><td>1</td><td>1001</td><td>RSID_101</td><td>SNPID_101</td><td>0.86274</td><td>0.416977</td><td>0.984639</td><td>0.510629</td><td>-1</td><td>0.425112</td><td>-1</td><td>0.463963</td></tr><tr><th>2</th><td>1</td><td>2000</td><td>RSID_2</td><td>SNPID_2</td><td>0.192181</td><td>0.19751</td><td>9.0</td><td>4.5062e-10</td><td>1</td><td>1.40173e-21</td><td>-1</td><td>6.67473e-21</td></tr><tr><th>3</th><td>1</td><td>2001</td><td>RSID_102</td><td>SNPID_102</td><td>0.1844</td><td>0.197667</td><td>0.727309</td><td>4.47847e-10</td><td>-1</td><td>1.79121e-21</td><td>1</td><td>8.45783e-21</td></tr><tr><th>4</th><td>1</td><td>3000</td><td>RSID_3</td><td>SNPID_3</td><td>0.965354</td><td>0.483396</td><td>0.955355</td><td>0.0721822</td><td>-1</td><td>0.600949</td><td>-1</td><td>0.128884</td></tr><tr><th>5</th><td>1</td><td>3001</td><td>RSID_103</td><td>SNPID_103</td><td>0.965354</td><td>0.483396</td><td>0.955355</td><td>0.0721821</td><td>1</td><td>0.600949</td><td>1</td><td>0.128884</td></tr><tr><th>6</th><td>1</td><td>4000</td><td>RSID_4</td><td>SNPID_4</td><td>0.371927</td><td>0.21671</td><td>0.991768</td><td>0.112472</td><td>-1</td><td>0.133752</td><td>-1</td><td>0.122192</td></tr><tr><th>7</th><td>1</td><td>4001</td><td>RSID_104</td><td>SNPID_104</td><td>0.371928</td><td>0.21671</td><td>0.991768</td><td>0.112471</td><td>1</td><td>0.133752</td><td>1</td><td>0.122192</td></tr><tr><th>8</th><td>1</td><td>5000</td><td>RSID_5</td><td>SNPID_5</td><td>0.587013</td><td>0.388082</td><td>0.968258</td><td>0.526325</td><td>-1</td><td>0.542332</td><td>-1</td><td>0.534209</td></tr></tbody></table></div>



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

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.037584
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.032523
      1.810614 seconds (6.52 M allocations: 446.105 MiB, 7.47% gc time, 55.30% compilation time)





    
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
    1	758311	rs12562034	0.07763975155279501	0.4098763332666681	0.38938987159522687	-1	0.9228220976232072	-1	0.5476822137398479
    1	967643	rs2710875	0.32407407407407407	4.076249100705747e-7	5.838562452864651e-6	1	3.354075049151666e-6	-1	4.935983727867444e-7
    1	1168108	rs11260566	0.19158878504672894	0.1285682279446898	0.0008508894259509042	1	0.055905499774689345	1	0.0010182674469333265
    1	1375074	rs1312568	0.441358024691358	2.5376019650614977e-19	0.00017168324418138014	-1	1.8481146312336014e-13	1	2.0709090755361968e-15
    1	1990452	rs2678939	0.4537037037037037	5.07695957708431e-11	1.278709529589261e-7	1	7.7380891475822e-9	-1	3.3986012947082215e-9
    1	2194615	rs7553178	0.22685185185185186	0.17056143157457776	0.07404001756940146	-1	0.00023109048716816823	1	0.00010021413726889936
    1	2396747	rs13376356	0.1448598130841121	0.9053079215078139	0.48972088680069503	-1	0.6595593862223496	1	0.5620909278620404
    1	2823603	rs1563468	0.4830246913580247	4.23065537243926e-9	2.581476982676523e-7	-1	1.1700736777037956e-8	1	6.367979793593832e-8
    1	3025087	rs6690373	0.2538699690402477	9.238641887192776e-8	3.111965106232366e-5	1	2.0479935907554532e-7	-1	2.1230949377453417e-5



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

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.040074
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.037214
    run = 1, ‖Δβ‖ = 0.033243, ‖Δτ‖ = 0.146511, ‖ΔL‖ = 0.005476, status = Optimal, time(s) = 0.050248
    run = 2, ‖Δβ‖ = 0.005774, ‖Δτ‖ = 0.042246, ‖ΔL‖ = 0.001784, status = Optimal, time(s) = 0.042017
    run = 1, ‖Δβ‖ = 0.013090, ‖Δτ‖ = 0.130781, ‖ΔL‖ = 0.005011, status = Optimal, time(s) = 0.051252
    run = 2, ‖Δβ‖ = 0.003913, ‖Δτ‖ = 0.037309, ‖ΔL‖ = 0.001516, status = Optimal, time(s) = 0.061290
    run = 1, ‖Δβ‖ = 0.022159, ‖Δτ‖ = 0.141135, ‖ΔL‖ = 0.005554, status = Optimal, time(s) = 0.046682
    run = 2, ‖Δβ‖ = 0.001482, ‖Δτ‖ = 0.021700, ‖ΔL‖ = 0.001435, status = Optimal, time(s) = 0.050569
    run = 1, ‖Δβ‖ = 0.026764, ‖Δτ‖ = 0.368620, ‖ΔL‖ = 0.000317, status = Optimal, time(s) = 0.043348
    run = 2, ‖Δβ‖ = 0.003023, ‖Δτ‖ = 0.030938, ‖ΔL‖ = 0.003568, status = Optimal, time(s) = 0.049007
      3.558469 seconds (6.55 M allocations: 436.689 MiB, 3.70% gc time, 79.41% compilation time)





    
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

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.043004
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.035886
      1.051142 seconds (4.00 M allocations: 290.773 MiB, 16.08% gc time)





    
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
first(CSV.read("score.pval.txt", DataFrame), 8)
```




<div class="data-frame"><p>8 rows × 10 columns</p><table class="data-frame"><thead><tr><th></th><th>chr</th><th>pos</th><th>snpid</th><th>maf</th><th>hwepval</th><th>betapval</th><th>betadir</th><th>taupval</th><th>taudir</th><th>jointpval</th></tr><tr><th></th><th title="Int64">Int64</th><th title="Int64">Int64</th><th title="String">String</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Int64">Int64</th><th title="Float64">Float64</th><th title="Int64">Int64</th><th title="Float64">Float64</th></tr></thead><tbody><tr><th>1</th><td>1</td><td>554484</td><td>rs10458597</td><td>0.0</td><td>1.0</td><td>1.0</td><td>0</td><td>1.0</td><td>0</td><td>1.0</td></tr><tr><th>2</th><td>1</td><td>758311</td><td>rs12562034</td><td>0.0776398</td><td>0.409876</td><td>0.38939</td><td>-1</td><td>0.922822</td><td>-1</td><td>0.547682</td></tr><tr><th>3</th><td>1</td><td>967643</td><td>rs2710875</td><td>0.324074</td><td>4.07625e-7</td><td>5.83856e-6</td><td>1</td><td>3.35408e-6</td><td>-1</td><td>4.93598e-7</td></tr><tr><th>4</th><td>1</td><td>1168108</td><td>rs11260566</td><td>0.191589</td><td>0.128568</td><td>0.000850889</td><td>1</td><td>0.0559055</td><td>1</td><td>0.00101827</td></tr><tr><th>5</th><td>1</td><td>1375074</td><td>rs1312568</td><td>0.441358</td><td>2.5376e-19</td><td>0.000171683</td><td>-1</td><td>1.84811e-13</td><td>1</td><td>2.07091e-15</td></tr><tr><th>6</th><td>1</td><td>1588771</td><td>rs35154105</td><td>0.0</td><td>1.0</td><td>1.0</td><td>0</td><td>1.0</td><td>0</td><td>1.0</td></tr><tr><th>7</th><td>1</td><td>1789051</td><td>rs16824508</td><td>0.00462963</td><td>0.933278</td><td>0.295035</td><td>1</td><td>0.304109</td><td>1</td><td>0.299504</td></tr><tr><th>8</th><td>1</td><td>1990452</td><td>rs2678939</td><td>0.453704</td><td>5.07696e-11</td><td>1.27871e-7</td><td>1</td><td>7.73809e-9</td><td>-1</td><td>3.3986e-9</td></tr></tbody></table></div>



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

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.038970
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.037406
    run = 1, ‖Δβ‖ = 0.026764, ‖Δτ‖ = 0.368620, ‖ΔL‖ = 0.000317, status = Optimal, time(s) = 0.060669
    run = 2, ‖Δβ‖ = 0.003023, ‖Δτ‖ = 0.030938, ‖ΔL‖ = 0.003568, status = Optimal, time(s) = 0.046794
    run = 1, ‖Δβ‖ = 0.039946, ‖Δτ‖ = 0.185049, ‖ΔL‖ = 0.005690, status = Optimal, time(s) = 0.053884
    run = 2, ‖Δβ‖ = 0.001640, ‖Δτ‖ = 0.029639, ‖ΔL‖ = 0.002028, status = Optimal, time(s) = 0.043935
    run = 1, ‖Δβ‖ = 0.040207, ‖Δτ‖ = 0.407863, ‖ΔL‖ = 0.001924, status = Optimal, time(s) = 0.052120
    run = 2, ‖Δβ‖ = 0.002632, ‖Δτ‖ = 0.039609, ‖ΔL‖ = 0.005131, status = Optimal, time(s) = 0.043479
    run = 1, ‖Δβ‖ = 0.052372, ‖Δτ‖ = 0.828946, ‖ΔL‖ = 0.001277, status = Optimal, time(s) = 0.064348
    run = 2, ‖Δβ‖ = 0.013694, ‖Δτ‖ = 0.185493, ‖ΔL‖ = 0.007118, status = Optimal, time(s) = 0.046386
    run = 1, ‖Δβ‖ = 0.032435, ‖Δτ‖ = 0.199659, ‖ΔL‖ = 0.004973, status = Optimal, time(s) = 0.053088
    run = 2, ‖Δβ‖ = 0.002847, ‖Δτ‖ = 0.032438, ‖ΔL‖ = 0.002941, status = Optimal, time(s) = 0.043273
    run = 1, ‖Δβ‖ = 0.009362, ‖Δτ‖ = 0.625969, ‖ΔL‖ = 0.002237, status = Optimal, time(s) = 0.064251
    run = 2, ‖Δβ‖ = 0.005893, ‖Δτ‖ = 0.176417, ‖ΔL‖ = 0.005017, status = Optimal, time(s) = 0.047283
    run = 1, ‖Δβ‖ = 0.012805, ‖Δτ‖ = 0.243313, ‖ΔL‖ = 0.006289, status = Optimal, time(s) = 0.056666
    run = 2, ‖Δβ‖ = 0.002636, ‖Δτ‖ = 0.044210, ‖ΔL‖ = 0.003312, status = Optimal, time(s) = 0.046489
    run = 1, ‖Δβ‖ = 0.014062, ‖Δτ‖ = 0.225841, ‖ΔL‖ = 0.003899, status = Optimal, time(s) = 0.057465
    run = 2, ‖Δβ‖ = 0.002112, ‖Δτ‖ = 0.030955, ‖ΔL‖ = 0.003860, status = Optimal, time(s) = 0.047704
    run = 1, ‖Δβ‖ = 0.026425, ‖Δτ‖ = 0.258805, ‖ΔL‖ = 0.004204, status = Optimal, time(s) = 0.056906
    run = 2, ‖Δβ‖ = 0.001090, ‖Δτ‖ = 0.051281, ‖ΔL‖ = 0.003822, status = Optimal, time(s) = 0.047899
    run = 1, ‖Δβ‖ = 0.035737, ‖Δτ‖ = 0.174702, ‖ΔL‖ = 0.001439, status = Optimal, time(s) = 0.063558
    run = 2, ‖Δβ‖ = 0.004776, ‖Δτ‖ = 0.019790, ‖ΔL‖ = 0.001867, status = Optimal, time(s) = 0.047032
      2.290992 seconds (3.89 M allocations: 297.829 MiB, 2.45% gc time, 45.79% compilation time)





    
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
CSV.read("wald.pval.txt", DataFrame)
```




<div class="data-frame"><p>10 rows × 9 columns</p><table class="data-frame"><thead><tr><th></th><th>chr</th><th>pos</th><th>snpid</th><th>maf</th><th>hwepval</th><th>betaeffect</th><th>betapval</th><th>taueffect</th><th>taupval</th></tr><tr><th></th><th title="Int64">Int64</th><th title="Int64">Int64</th><th title="String">String</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th></tr></thead><tbody><tr><th>1</th><td>1</td><td>1375074</td><td>rs1312568</td><td>0.441358</td><td>2.5376e-19</td><td>-0.456032</td><td>0.000165921</td><td>0.496452</td><td>1.17735e-36</td></tr><tr><th>2</th><td>1</td><td>11552817</td><td>rs2745282</td><td>0.421053</td><td>1.3537e-14</td><td>-0.67282</td><td>7.6492e-9</td><td>0.425653</td><td>1.45395e-19</td></tr><tr><th>3</th><td>1</td><td>120276030</td><td>rs6688004</td><td>0.469136</td><td>1.75473e-27</td><td>-0.637817</td><td>5.13062e-9</td><td>0.430562</td><td>6.07457e-21</td></tr><tr><th>4</th><td>2</td><td>135623558</td><td>rs6730157</td><td>0.313272</td><td>7.76753e-16</td><td>-0.556166</td><td>1.10163e-5</td><td>0.462346</td><td>2.88628e-20</td></tr><tr><th>5</th><td>9</td><td>126307510</td><td>rs3814134</td><td>0.427469</td><td>4.28276e-33</td><td>0.723472</td><td>1.41601e-11</td><td>-0.518216</td><td>4.75338e-39</td></tr><tr><th>6</th><td>15</td><td>40603025</td><td>rs2617236</td><td>0.390966</td><td>1.16466e-14</td><td>-0.644375</td><td>1.44948e-7</td><td>0.445647</td><td>9.43451e-19</td></tr><tr><th>7</th><td>15</td><td>40803767</td><td>rs3742988</td><td>0.42284</td><td>4.26181e-30</td><td>-0.796743</td><td>2.24785e-13</td><td>0.389234</td><td>7.09349e-17</td></tr><tr><th>8</th><td>17</td><td>56509992</td><td>rs8064681</td><td>0.481481</td><td>1.37485e-21</td><td>-0.58961</td><td>1.9256e-7</td><td>0.425417</td><td>5.60171e-25</td></tr><tr><th>9</th><td>17</td><td>71293786</td><td>rs2125345</td><td>0.339009</td><td>1.75527e-14</td><td>-0.635669</td><td>1.82739e-7</td><td>0.449612</td><td>1.7146e-21</td></tr><tr><th>10</th><td>23</td><td>64815688</td><td>rs5964999</td><td>0.475078</td><td>3.63212e-56</td><td>0.75886</td><td>5.04109e-14</td><td>-0.38849</td><td>5.39098e-23</td></tr></tbody></table></div>




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

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.041253
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.035768





    
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
first(CSV.read("GxE.pval.txt", DataFrame), 5)
```




<div class="data-frame"><p>5 rows × 8 columns</p><table class="data-frame"><thead><tr><th></th><th>chr</th><th>pos</th><th>snpid</th><th>maf</th><th>hwepval</th><th>betapval</th><th>taupval</th><th>jointpval</th></tr><tr><th></th><th title="Int64">Int64</th><th title="Int64">Int64</th><th title="String">String</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th></tr></thead><tbody><tr><th>1</th><td>1</td><td>554484</td><td>rs10458597</td><td>0.0</td><td>1.0</td><td>1.0</td><td>1.0</td><td>1.0</td></tr><tr><th>2</th><td>1</td><td>758311</td><td>rs12562034</td><td>0.0776398</td><td>0.409876</td><td>0.383</td><td>0.33511</td><td>0.357458</td></tr><tr><th>3</th><td>1</td><td>967643</td><td>rs2710875</td><td>0.324074</td><td>4.07625e-7</td><td>9.98217e-6</td><td>5.73202e-7</td><td>1.08415e-6</td></tr><tr><th>4</th><td>1</td><td>1168108</td><td>rs11260566</td><td>0.191589</td><td>0.128568</td><td>0.00225789</td><td>0.105482</td><td>0.00442114</td></tr><tr><th>5</th><td>1</td><td>1375074</td><td>rs1312568</td><td>0.441358</td><td>2.5376e-19</td><td>0.000928728</td><td>3.48807e-15</td><td>6.97615e-15</td></tr></tbody></table></div>




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

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.042947
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.051226
    run = 1, ‖Δβ‖ = 0.480377, ‖Δτ‖ = 0.023510, ‖ΔL‖ = 0.001822, status = Optimal, time(s) = 0.045999
    run = 2, ‖Δβ‖ = 0.000323, ‖Δτ‖ = 0.002136, ‖ΔL‖ = 0.000025, status = Optimal, time(s) = 0.044744
    run = 1, ‖Δβ‖ = 1.059989, ‖Δτ‖ = 0.458239, ‖ΔL‖ = 0.063093, status = Optimal, time(s) = 0.065265
    run = 2, ‖Δβ‖ = 0.009118, ‖Δτ‖ = 0.105352, ‖ΔL‖ = 0.000056, status = Optimal, time(s) = 0.061506
    run = 1, ‖Δβ‖ = 1.144781, ‖Δτ‖ = 0.251459, ‖ΔL‖ = 0.034318, status = Optimal, time(s) = 0.060325
    run = 2, ‖Δβ‖ = 0.000830, ‖Δτ‖ = 0.025631, ‖ΔL‖ = 0.000042, status = Optimal, time(s) = 0.043208
    run = 1, ‖Δβ‖ = 0.683108, ‖Δτ‖ = 0.773919, ‖ΔL‖ = 0.042397, status = Optimal, time(s) = 0.052987
    run = 2, ‖Δβ‖ = 0.008159, ‖Δτ‖ = 0.069178, ‖ΔL‖ = 0.003574, status = Optimal, time(s) = 0.045152





    
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
CSV.read("gxe_snp.pval.txt", DataFrame)
```




<div class="data-frame"><p>5 rows × 9 columns</p><table class="data-frame"><thead><tr><th></th><th>chr</th><th>pos</th><th>snpid</th><th>maf</th><th>hwepval</th><th>snpeffectnullbeta</th><th>snpeffectnulltau</th><th>betapval</th><th>taupval</th></tr><tr><th></th><th title="Int64">Int64</th><th title="Int64">Int64</th><th title="String">String</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th></tr></thead><tbody><tr><th>1</th><td>1</td><td>554484</td><td>rs10458597</td><td>0.0</td><td>1.0</td><td>0.0</td><td>0.0</td><td>1.0</td><td>1.0</td></tr><tr><th>2</th><td>1</td><td>758311</td><td>rs12562034</td><td>0.0776398</td><td>0.409876</td><td>-0.230593</td><td>0.0102355</td><td>0.448376</td><td>0.149797</td></tr><tr><th>3</th><td>1</td><td>967643</td><td>rs2710875</td><td>0.324074</td><td>4.07625e-7</td><td>0.639871</td><td>-0.338815</td><td>0.293876</td><td>0.119916</td></tr><tr><th>4</th><td>1</td><td>1168108</td><td>rs11260566</td><td>0.191589</td><td>0.128568</td><td>0.618855</td><td>-0.146819</td><td>0.590189</td><td>0.262791</td></tr><tr><th>5</th><td>1</td><td>1375074</td><td>rs1312568</td><td>0.441358</td><td>2.5376e-19</td><td>-0.456003</td><td>0.492582</td><td>0.574737</td><td>0.250501</td></tr></tbody></table></div>




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

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.059821
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.037852





    
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

    The pvalue of snps indexed at 50:55 is betapval: 2.3156171211965083e-5, taupval: 3.9729937063719077e-10



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

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.040899
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.036685





    
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
first(CSV.read("snpset.pval.txt", DataFrame; delim="\t"), 5)
```




<div class="data-frame"><p>5 rows × 4 columns</p><table class="data-frame"><thead><tr><th></th><th>snpsetid</th><th>nsnps</th><th>betapval</th><th>taupval</th></tr><tr><th></th><th title="String">String</th><th title="Int64">Int64</th><th title="Float64">Float64</th><th title="Float64">Float64</th></tr></thead><tbody><tr><th>1</th><td>gene1</td><td>93</td><td>0.111864</td><td>0.011595</td></tr><tr><th>2</th><td>gene2</td><td>93</td><td>0.0249929</td><td>0.0648265</td></tr><tr><th>3</th><td>gene3</td><td>93</td><td>0.131741</td><td>0.0298884</td></tr><tr><th>4</th><td>gene4</td><td>92</td><td>0.010327</td><td>0.0584591</td></tr><tr><th>5</th><td>gene5</td><td>93</td><td>0.0303905</td><td>0.0924954</td></tr></tbody></table></div>




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

    run = 1, ‖Δβ‖ = 0.037090, ‖Δτ‖ = 0.136339, ‖ΔL‖ = 0.005441, status = Optimal, time(s) = 0.039755
    run = 2, ‖Δβ‖ = 0.000913, ‖Δτ‖ = 0.019810, ‖ΔL‖ = 0.001582, status = Optimal, time(s) = 0.043883





    
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
first(CSV.read("snpset.pval.txt", DataFrame), 5)
```




<div class="data-frame"><p>5 rows × 8 columns</p><table class="data-frame"><thead><tr><th></th><th>startchr</th><th>startpos</th><th>startsnpid</th><th>endchr</th><th>endpos</th><th>endsnpid</th><th>betapval</th><th>taupval</th></tr><tr><th></th><th title="Int64">Int64</th><th title="Int64">Int64</th><th title="String">String</th><th title="Int64">Int64</th><th title="Int64">Int64</th><th title="String">String</th><th title="Float64">Float64</th><th title="Float64">Float64</th></tr></thead><tbody><tr><th>1</th><td>1</td><td>554484</td><td>rs10458597</td><td>1</td><td>3431124</td><td>rs12093117</td><td>3.83046e-6</td><td>2.24303e-10</td></tr><tr><th>2</th><td>1</td><td>3633945</td><td>rs10910017</td><td>1</td><td>6514524</td><td>rs932112</td><td>0.000127922</td><td>7.94764e-6</td></tr><tr><th>3</th><td>1</td><td>6715827</td><td>rs441515</td><td>1</td><td>9534606</td><td>rs4926480</td><td>8.10353e-5</td><td>6.11848e-7</td></tr><tr><th>4</th><td>1</td><td>9737551</td><td>rs12047054</td><td>1</td><td>12559747</td><td>rs4845907</td><td>0.000425108</td><td>1.15311e-8</td></tr><tr><th>5</th><td>1</td><td>12760427</td><td>rs848577</td><td>1</td><td>16021797</td><td>rs6679870</td><td>2.75832e-5</td><td>0.00028553</td></tr></tbody></table></div>




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




    324-element Vector{String}:
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
first(covdf, 11)
```




<div class="data-frame"><p>11 rows × 8 columns</p><table class="data-frame"><thead><tr><th></th><th>sex</th><th>onMeds</th><th>snp1</th><th>snp2</th><th>snp3</th><th>snp4</th><th>y</th><th>id</th></tr><tr><th></th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="String">String</th></tr></thead><tbody><tr><th>1</th><td>0.0</td><td>1.0</td><td>0.0</td><td>1.0</td><td>2.0</td><td>0.0</td><td>12.2667</td><td>A1</td></tr><tr><th>2</th><td>0.0</td><td>0.0</td><td>0.0</td><td>1.0</td><td>2.0</td><td>0.0</td><td>10.2681</td><td>A1</td></tr><tr><th>3</th><td>0.0</td><td>0.0</td><td>0.0</td><td>1.0</td><td>2.0</td><td>0.0</td><td>12.166</td><td>A1</td></tr><tr><th>4</th><td>0.0</td><td>0.0</td><td>0.0</td><td>1.0</td><td>2.0</td><td>0.0</td><td>11.8797</td><td>A1</td></tr><tr><th>5</th><td>0.0</td><td>0.0</td><td>0.0</td><td>1.0</td><td>2.0</td><td>0.0</td><td>12.8127</td><td>A1</td></tr><tr><th>6</th><td>0.0</td><td>0.0</td><td>0.0</td><td>1.0</td><td>2.0</td><td>0.0</td><td>9.98766</td><td>A1</td></tr><tr><th>7</th><td>0.0</td><td>0.0</td><td>0.0</td><td>1.0</td><td>2.0</td><td>0.0</td><td>12.1408</td><td>A1</td></tr><tr><th>8</th><td>0.0</td><td>1.0</td><td>0.0</td><td>1.0</td><td>2.0</td><td>0.0</td><td>13.2058</td><td>A1</td></tr><tr><th>9</th><td>0.0</td><td>0.0</td><td>0.0</td><td>1.0</td><td>2.0</td><td>0.0</td><td>11.3631</td><td>A1</td></tr><tr><th>10</th><td>0.0</td><td>1.0</td><td>0.0</td><td>1.0</td><td>2.0</td><td>0.0</td><td>15.2511</td><td>A1</td></tr><tr><th>11</th><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>2.0</td><td>2.0</td><td>12.746</td><td>2</td></tr></tbody></table></div>




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

Step 1: Fit the null model. Setting third argument `geneticfile` to `nothing` instructs `trajgwas` function to fit the null model only.


```julia
nm = trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        datadir * "trajgwas_plinkex.csv",
        nothing)
```

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
solver = Ipopt.IpoptSolver(print_level=1, watchdog_shortened_iter_trigger=5, max_iter=120)

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
    - Issues with null model convergence may be solved by choosing alternate starting values for parameters, using a different solver, transforming variables, and increasing the number of runs (WiSER runs). These are detailed in the [WiSER documentation here](https://github.com/OpenMendel/WiSER.jl/blob/master/docs/src/model_fitting.md#tips-for-improving-estimation).
    
- GWAS Results
    - If you use the score test instead of the SPA-score test (SPA is default for single-SNP analyses), then there can be inflation in type I error and decreased power when (a) the sample size is small, (b) the number of repeated measures is low, or (c) the variants analyzed are rare with low minor allele frequencies. In these cases, the score test is not optimal and it is suggested to use the SPA version (`usespa=true`). SPA is only implemented for single-SNP analyses. These issues can occur in both Wald and score tests. 
    
If you notice any problems with your output or results, [file an issue](https://github.com/OpenMendel/TrajGWAS.jl/issues). 
