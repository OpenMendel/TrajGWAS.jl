# TrajGWAS

<!--[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/TrajGWAS.jl/stable)-->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://OpenMendel.github.io/TrajGWAS.jl/dev)
[![build Actions Status](https://github.com/OpenMendel/TrajGWAS.jl/workflows/CI/badge.svg)](https://github.com/OpenMendel/TrajGWAS.jl/actions)
[![codecov](https://codecov.io/gh/OpenMendel/TrajGWAS.jl/branch/main/graph/badge.svg?token=MABDDASURF)](https://codecov.io/gh/OpenMendel/TrajGWAS.jl)

TrajGWAS.jl is a Julia package for performing genome-wide association studies (GWAS) for continuous longitudinal phenotypes using a modified linear mixed effects model. It builds upon the [within-subject variance estimation by robust regression (WiSER) method](https://github.com/OpenMendel/WiSER.jl) and can be used to identify variants associated with changes in the mean and within-subject variability of the longitduinal trait. The estimation procedure is robust to distributional misspecifications of both the random effects and the response. A saddlepoint approximation (SPA) option is implemented to provide improved power and calibrated type I error for rare variants. It runs efficiently and scales well to very large datasets. The package currently supports [PLINK](https://zzz.bwh.harvard.edu/plink/), [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) (both dosage and genotype data), and [BGEN](https://www.well.ox.ac.uk/~gav/bgen_format/) file formats. We plan to add [PGEN](https://www.cog-genomics.org/plink/2.0/formats#pgen) support in the future. 

TrajGWAS.jl supports Julia v1.5 or later. See the [documentation](https://openmendel.github.io/TrajGWAS.jl/dev/) for usage.  
<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/TrajGWAS.jl/stable) [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/TrajGWAS.jl/latest) -->

TrajGWAS.jl is a registered package, and it will require running the following code to install. 

```{julia}
using Pkg
pkg"add TrajGWAS"
```

## Citation

The methods and applications of this software package are detailed in the following publication:

*Ko S, German CA, Jensen A, Shen J, Wang A, Mehrotra DV, Sun YV, Sinsheimer JS, Zhou H, Zhou JJ. GWAS of longitudinal trajectories at biobank scale. Am J Hum Genet. 2022 Mar 3;109(3):433-445. doi: 10.1016/j.ajhg.2022.01.018. Epub 2022 Feb 22. PMID: 35196515; PMCID: PMC8948167.*

If you use [OpenMendel](https://openmendel.github.io) analysis packages in your research, please cite the following reference in the resulting publications:

*OPENMENDEL: a cooperative programming project for statistical genetics. Zhou H, Sinsheimer JS, Bates DM, Chu BB, German CA, Ji SS, Keys KL, Kim J, Ko S, Mosher GD, Papp JC, Sobel EM, Zhai J, Zhou JJ, Lange K. Hum Genet. 139, 61–71 (2020). [doi:10.1007/s00439-019-02001-z](10.1007/s00439-019-02001-z). PMCID:[PMC6763373](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6763373/).*

## Acknowledgments

This project has been supported by the National Institutes of Health under awards R01GM053275, R01HG006139, R25GM103774, and 1R25HG011845.
