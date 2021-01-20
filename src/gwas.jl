"""
    vgwas(nullmeanformula, reformula, nullwsvarformula, idvar, covfile, geneticfile; kwargs...)
    vgwas(nullmeanformula, reformula, nullwsvarformula, idvar, df, geneticfile; kwargs...)
    vgwas(fittednullmodel, geneticfile; kwargs...)
    vgwas(fittednullmodel, bedfile, bimfile, bedn; kwargs...)
    vgwas(fittednullmodel, vcffile, nsamples, vcftype; kwargs...)

# Positional arguments
- `nullmeanformula::FormulaTerm`: mean formula (β) for the null model.
- `reformula::FormulaTerm`: random effects formula (γ) for the model.
- `nullmeanformula::FormulaTerm`: within-subject variance formula (τ) for the null model.
- `idvar::Union{String, Symbol}`: id variable for groupings.
- `covfile::AbstractString`: covariate file (csv) with one header line. One column
    should be the longitudinal phenotype.
- `df::DataFrame`: DataFrame containing response and regressors for null model.
- `geneticfile::Union{Nothing, AbstractString}`: File containing genetic information for GWAS.
    This includes a PLINK file name without the .bed, .fam, or .bim
    extensions or a VCF file without the .vcf extension. If `geneticfile==nothing`,
    only null model is fitted. If `geneticfile` is provided, bed, bim, and fam file (or vcf) with
    the same `geneticfile` prefix need to exist. Compressed file formats such as gz and bz2
    are allowed. Check all allowed formats by `SnpArrays.ALLOWED_FORMAT`. If you're using a VCF file,
    make sure to use the `geneticformat = "VCF"` keyword option, and specificy dosage (:DS) or
    genotype (:GT) data with the `vcftype` command.
- `fittednullmodel::StatsModels.TableRegressionModel`: the fitted null model
    output from `vgwas(nullformula, covfile)` or `vgwas(nullformula, df)`.
    **NOTE** If the nullmodel is passed in with the `bedfile, bimfile, bedn` or 
    `vcffile, nsamples, vcftype` arguments, the IDs/data in the null model must match the
    order of the IDs in the PLINK/VCF file.  
- `bedfile::Union{AbstractString,IOStream}`: path to Plink bed file with full file name.
- `bimfile::Union{AbstractString,IOStream}`: path to Plink bim file with full file name.
- `bedn::Integer`: number of samples in bed/vcf file.

# Keyword arguments
- `analysistype`::AbstractString: Type of analysis to conduct. Default is `singlesnp`. Other options are `snpset` and `gxe`.
- `geneticformat`::AbstractString: Type of file used for the genetic analysis. `"PLINK"` and `"VCF"` are currently supported. Default is PLINK.
- `vcftype`::Union{Symbol, Nothing}: Data to extract from the VCF file for the GWAS analysis. `:DS` for dosage or `:GT` for genotypes. Default is nothing.
- `nullfile::Union{AbstractString, IOStream}`: output file for the fitted null model;
    default is `vgwas.null.txt`.
- `pvalfile::Union{AbstractString, IOStream}`: output file for the gwas p-values; default is
    `vgwas.pval.txt`.
- `covtype::Vector{DataType}`: type information for `covfile`. This is useful
    when `CSV.read(covarfile)` has parsing errors.
- `covrowinds::Union{Nothing,AbstractVector{<:Integer}}`: sample indices for covariate file.
- `testformula::FormulaTerm`: formula for test unit. Default is `@formula(trait ~ 0 + snp)`.
- `test::Symbol`: `:score` (default) or `:wald`.
- `link::GLM.Link`: `LogitLink()` (default), `ProbitLink()`, `CauchitLink()`,
    or `CloglogLink()`.
- `snpmodel`: `ADDITIVE_MODEL` (default), `DOMINANT_MODEL`, or `RECESSIVE_MODEL`.
- `snpinds::Union{Nothing,AbstractVector{<:Integer}}`: SNP indices for bed/vcf file.
- `geneticrowinds::Union{Nothing,AbstractVector{<:Integer}}`: sample indices for bed/vcf file.
- `solver`: an optimization solver supported by MathProgBase. Default is
    `NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000)`. Another common choice is
    `IpoptSolver(print_level=0)`.
- `runs::Integer`: Number of weighted NLS runs for the null model; default is 2.
    Each run will use the newest `m.τ` and `m.Lγ` to update the weight matrices
    `Vi` and solve the new weighted NLS.
- `parallel::Bool`: Multi-threading or not for fitting the null model. Default is `false`.
- `verbose::Bool`: default is `false`.
- `snpset::Union{Nothing, Integer, AbstractString, AbstractVector{<:Integer}}`: Only include
    if you are conducting a snpset analysis. An integer indicates a window of SNPs
    (i.e. every 500 snps). An abstract string allows you to specify an input file,
    with no header and two columns separated by a space. The first column must contain the snpset ID
    and the second column must contain the snpid's identical to the bimfile. An AbstractVector
    allows you to specify the snps you want to perform one joint snpset test for.
- `e::Union{AbstractString,Symbol}`: Only include if you are conducting a GxE analysis.
    Enviromental variable to be used to test the GxE interaction.
    For instance, for testing `sex & snp` interaction, use `:sex` or `"sex"`.


# Examples
The following is an example of basic GWAS with PLINK files:
```julia
plinkfile = "plinkexample"
covfile = "covexample"
vgwas(@formula(trait ~ sex), covfile, plkfile)
```

The following is an example of basic GWAS with a VCF file using dosages then genotypes:
```julia
vcffile = "vcfexample"
covfile = "covexample"
vgwas(@formula(trait ~ sex), covfile, vcffile;
    geneticfile = "VCF", vcftype = :DS)

vgwas(@formula(trait ~ sex), covfile, vcffile;
    geneticfile = "VCF", vcftype = :GT)
```

The following is an example of snpset GWAS (every 50 snps). For more types of snpset analyses see documentation:
```julia
vgwas(@formula(trait ~ sex), covfile, plkfile;
    analysistype = "snpset", snpset = 50)
```

The following is an example of GxE GWAS testing the interaction effect:
```julia
vgwas(@formula(trait ~ sex), covfile, plkfile;
    analysistype = "gxe", e = :sex)
```

"""
function vgwas(
    # positional arguments
    nullmeanformula::FormulaTerm,
    reformula::FormulaTerm,
    nullwsvarformula::FormulaTerm,
    idvar::Union{String, Symbol},
    covfile::AbstractString,
    geneticfile::Union{Nothing, AbstractString} = nothing;
    # keyword arguments
    covtype::Union{Nothing, Vector{DataType}} = nothing,
    covrowinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    kwargs...
    )
    covdf = SnpArrays.makestream(covfile) do io
        CSV.read(io, DataFrame; types=covtype)
    end
    vgwas(nullmeanformula, reformula, nullwsvarformula, idvar,
     covrowinds == nothing ? covdf : covdf[covrowinds, :],
        geneticfile; kwargs...)
end

function vgwas(
    nullmeanformula::FormulaTerm,
    reformula::FormulaTerm,
    nullwsvarformula::FormulaTerm,
    idvar::Union{String, Symbol},
    nulldf::DataFrame,
    geneticfile::Union{Nothing, AbstractString} = nothing;
    nullfile::Union{AbstractString, IOStream} = "vgwas.null.txt",
    solver = Ipopt.IpoptSolver(print_level=0, mehrotra_algorithm = "yes",
    warm_start_init_point="yes", max_iter=100),
    parallel::Bool = false,
    runs::Int = 2,
    verbose::Bool = false,
    kwargs...
    )
    # fit and output null model
    nm = WSVarLmmModel(nullmeanformula, reformula, nullwsvarformula,
     idvar, nulldf)
    WiSER.fit!(nm, solver, parallel = parallel, runs = runs)
    verbose && show(nm)
    SnpArrays.makestream(nullfile, "w") do io
        show(io, nm)
    end
    geneticfile == nothing && (return nm)
    vgwas(nm, geneticfile; solver = solver, verbose = verbose, kwargs...)
end

function vgwas(
    # positional arguments
    fittednullmodel::WSVarLmmModel,
    geneticfile::AbstractString;
    # keyword arguments
    analysistype::AbstractString = "singlesnp",
    geneticformat::AbstractString = "PLINK",
    vcftype::Union{Symbol, Nothing} = nothing,
    testformula::FormulaTerm = fittednullmodel.meanformula.lhs ~ Term(:snp),
    test::Symbol = :score,
    pvalfile::Union{AbstractString, IOStream} = "vgwas.pval.txt",
    snpmodel::Union{Val{1}, Val{2}, Val{3}} = ADDITIVE_MODEL,
    snpinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    usespa::Bool = true,
    geneticrowinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    solver = Ipopt.IpoptSolver(print_level=0, mehrotra_algorithm = "yes",
    warm_start_init_point="yes", max_iter=100),
    parallel::Bool = false,
    runs::Int = 2,
    verbose::Bool = false,
    snpset::Union{Nothing, Integer, AbstractString, #for snpset analysis
        AbstractVector{<:Integer}} = nothing,
    e::Union{Nothing, AbstractString, Symbol} = nothing # for GxE analysis
    )

    # locate plink bed, fam, bim files or VCF file
    lowercase(geneticformat) in ["plink", "vcf"] || error("`geneticformat` $geneticformat not valid. Please use 'VCF' or 'PLINK'.")
    isplink = "plink" == lowercase(geneticformat)
    if isplink
        if isfile(geneticfile * ".bed")
            bedfile = geneticfile * ".bed"
        else
            fmt = findfirst(isfile, geneticfile * ".bed." .* SnpArrays.ALLOWED_FORMAT)
            fmt == nothing && throw(ArgumentError("bed file not found"))
            bedfile = geneticfile * ".bed." * SnpArrays.ALLOWED_FORMAT[fmt]
        end
        famfile = replace(bedfile, ".bed" => ".fam")
        isfile(famfile) || throw(ArgumentError("fam file not found"))
        bimfile = replace(bedfile, ".bed" => ".bim")
        isfile(bimfile) || throw(ArgumentError("bim file not found"))
        # selected rows should match nobs in null model
        bedn = SnpArrays.makestream(countlines, famfile)
        bedids = SnpArrays.makestream(famfile) do io
            CSV.read(io, DataFrame, header=false)[!, 1] # assuming first column is IDs column
        end
    else
        vcftype in [:GT, :DS] || throw(ArgumentError("vcftype not specified. Allowable types are :GT for genotypes and :DS for dosages."))
        if isfile(geneticfile * ".vcf")
            vcffile = geneticfile * ".vcf"
        else
            fmt = findfirst(isfile, geneticfile * ".vcf." .* SnpArrays.ALLOWED_FORMAT)
            fmt == nothing && throw(ArgumentError("VCF file not found"))
            vcffile = geneticfile * ".vcf." * SnpArrays.ALLOWED_FORMAT[fmt]
        end
        bedn = VCFTools.nsamples(vcffile)
        reader = VCF.Reader(openvcf(vcffile, "r"))
        bedids = VCF.header(reader).sampleID
        close(reader)
    end
    if geneticrowinds == nothing
        nbedrows = bedn
        rowinds = 1:bedn
    else
        nbedrows = eltype(geneticrowinds) == Bool ? count(geneticrowinds) : length(geneticrowinds)
        rowinds = geneticrowinds
    end


    # Match null model ids to genetic file
    nullinds = indexin(bedids[rowinds], fittednullmodel.ids)
    # make sure all ids in null model are in bedids[geneticrowinds]
    @assert all([id in bedids[rowinds] for id in fittednullmodel.ids]) "Null model contains IDs not in the $geneticfile file applying `geneticrowinds` mask."
    # make sure there aren't more rows in genetic file after applying mask
    @assert count(isnothing.(nullinds)) == 0 "$geneticfile file IDs with `geneticrowinds` contains IDs not in null model."
    # make sure the sets of IDs are equal
    @assert bedids[rowinds] == fittednullmodel.ids[nullinds] "$geneticfile file IDs with `geneticrowinds` contains IDs not in null model."

    nbedrows == fittednullmodel.m ||
        throw(ArgumentError("number of samples in geneticrowinds does not match null model"))

    # rearrange null model to match bedids[rowinds] order
    inorder = fittednullmodel.ids == fittednullmodel.ids[nullinds]
    if !inorder
        fittednullmodel.obswts .= isempty(fittednullmodel.obswts) ?
            fittednullmodel.obswts : fittednullmodel.obswts[nullinds]
        fittednullmodel.ids .= fittednullmodel.ids[nullinds]
        fittednullmodel.nis .= fittednullmodel.nis[nullinds]
        fittednullmodel.data .= fittednullmodel.data[nullinds]
        warn("Null model data was rearranged to match genetic file order")
    end

    # validate testing method
    test = Symbol(lowercase(string(test)))
    test == :score || test == :wald || throw(ArgumentError("unrecognized test $test"))

    # gwas
    if isplink #plink
        vgwas(fittednullmodel, bedfile, bimfile, bedn;
            analysistype = analysistype,
            testformula = testformula,
            test = test,
            pvalfile = pvalfile,
            snpmodel = snpmodel,
            snpinds = snpinds,
            usespa = usespa,
            bedrowinds = rowinds,
            solver = solver,
            parallel = parallel,
            runs = runs,
            verbose = verbose,
            snpset = snpset,
            e = e)
    else #vcf
        vgwas(fittednullmodel, vcffile, bedn, vcftype;
            analysistype = analysistype,
            testformula = testformula,
            test = test,
            pvalfile = pvalfile,
            snpmodel = snpmodel,
            snpinds = snpinds,
            usespa = usespa,
            vcfrowinds = rowinds,
            solver = solver,
            parallel = parallel,
            runs = runs,
            verbose = verbose,
            snpset = snpset,
            e = e)
    end
end

# For PLINK Analysis
function vgwas(
    fittednullmodel::WSVarLmmModel,
    bedfile::Union{AbstractString, IOStream}, # full path and bed file name
    bimfile::Union{AbstractString, IOStream}, # full path and bim file name
    bedn::Integer;           # number of samples in bed file
    analysistype::AbstractString = "singlesnp",
    testformula::FormulaTerm = fittednullmodel.meanformula.lhs ~ Term(:snp),
    test::Symbol = :score,
    pvalfile::Union{AbstractString, IOStream} = "vgwas.pval.txt",
    snpmodel::Union{Val{1}, Val{2}, Val{3}} = ADDITIVE_MODEL,
    snpinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    usespa::Bool = true,
    bedrowinds::AbstractVector{<:Integer} = 1:bedn, # row indices for SnpArray
    solver = Ipopt.IpoptSolver(print_level=0, mehrotra_algorithm = "yes",
    warm_start_init_point="yes", max_iter=100),
    parallel::Bool = false,
    runs::Int = 2,
    verbose::Bool = false,
    snpset::Union{Nothing, Integer, AbstractString, #for snpset analysis
        AbstractVector{<:Integer}} = nothing,
    e::Union{Nothing, AbstractString, Symbol} = nothing # for GxE analysis
    )
    # create SnpArray
    genomat = SnpArrays.SnpArray(bedfile, bedn)
    cc = SnpArrays.counts(genomat, dims=1) # column counts of genomat
    mafs = SnpArrays.maf(genomat)
    nbedrows = eltype(bedrowinds) == Bool ? count(bedrowinds) : length(bedrowinds)

    # storage vectors for SPA if it is set to true
    if (usespa == true) & (analysistype == "singlesnp")
        tmp = Array{Float64}(undef, fittednullmodel.m)
        tmp2 = Array{Float64}(undef, 3)
    end

    # create SNP mask vector
    if snpinds == nothing
        snpmask = trues(SnpArrays.makestream(countlines, bimfile))
    elseif eltype(snpinds) == Bool
        snpmask = snpinds
    else
        snpmask = falses(SnpArrays.makestream(countlines, bimfile))
        snpmask[snpinds] .= true
    end

    analysistype = lowercase(analysistype)
    analysistype in ["singlesnp", "snpset", "gxe"] ||
        error("Analysis type $analysis invalid option.
        Available options are 'singlesnp', 'snpset' and 'gxe'.")

    # determine analysis type
    if analysistype == "singlesnp"

        # carry out score or wald test SNP by SNP
        snponly = testformula.rhs == Term(:snp)
        # extra columns in design matrix to be tested
        testdf = DataFrame(fittednullmodel)
        testdf[!, :snp] = zeros(fittednullmodel.nsum)
        snpholder = zeros(fittednullmodel.m)
        if snponly
            q = 1
        else
            Z = similar(modelmatrix(testformula, testdf))
            q = size(Z, 2)
        end
        SnpArrays.makestream(pvalfile, "w") do io
            if test == :score
                println(io, "chr,pos,snpid,maf,hwepval,betapval,taupval,jointpval")
                if snponly
                    ts = WSVarScoreTestInvariant(fittednullmodel, 1, 1)
                    if usespa
                        Ks = ecgf(ts)
                    end
                else
                    ts = WSVarScoreTest(fittednullmodel, q, q)
                    testvec = [Matrix{Float64}(undef, ni, q) for
                    ni in fittednullmodel.nis]
                end
            else # wald
                γ̂β = Vector{Float64}(undef, q) # effect size for columns being tested
                pvalsβ = Vector{Float64}(undef, q) # effect size for columns being tested
                γ̂τ = Vector{Float64}(undef, q) # effect size for columns being tested
                pvalsτ = Vector{Float64}(undef, q) # effect size for columns being tested
                if snponly
                    println(io, "chr,pos,snpid,maf,hwepval,betaeffect,betapval,",
                    "taueffect,taupval")
                    fullmeanformula = FormulaTerm(fittednullmodel.meanformula.lhs,
                    sum(union(fittednullmodel.meanformula.rhs, [testformula.rhs])))
                    fullwsvarformula = FormulaTerm(fittednullmodel.meanformula.lhs,
                    sum(union(fittednullmodel.wsvarformula.rhs, [testformula.rhs])))
                else
                    fullmeanformula = FormulaTerm(fittednullmodel.meanformula.lhs,
                    sum(union(fittednullmodel.meanformula.rhs, testformula.rhs)))
                    fullwsvarformula = FormulaTerm(fittednullmodel.meanformula.lhs,
                    sum(union(fittednullmodel.wsvarformula.rhs, testformula.rhs)))
                    print(io, "chr,pos,snpid,maf,hwepval,")
                    for j in 1:q
                        print(io, "betaeffect$j,betapval$j,")
                    end
                    for j in 1:q
                        print(io, "taueffect$j,taupval$j")
                        if j != q
                            print(io, ",")
                        end
                    end
                    println(io, "")
                end
            end
            SnpArrays.makestream(bimfile) do bimio
                for j in eachindex(snpmask)
                    row = readline(bimio)
                    snpmask[j] || continue
                    snpj = split(row)
                    hwepval = SnpArrays.hwe(cc[1, j], cc[3, j], cc[4, j])
                    maf = mafs[j]
                    if test == :score
                        if maf == 0 # mono-allelic
                            betapval = 1.0
                            taupval = 1.0
                            jointpval = 1.0
                        else
                            if snponly
                                copyto!(snpholder, @view(genomat[bedrowinds, j]),
                                impute=true, model=snpmodel)
                                betapval, taupval, jointpval = test!(ts, snpholder, snpholder)
                                ps = test!(ts, snpholder, snpholder)
                                betapval, taupval, jointpval = ps
                                if usespa
                                    betapval, taupval, jointpval = spa(snpholder, ts, 
                                        ps, Ks; tmp_g = tmp, tmp_g2 = tmp2)
                                end
                            else # snp + other terms
                                copyto!(snpholder, @view(genomat[bedrowinds, j]),
                                    impute=true, model=snpmodel)
                                snptodf!(testdf[!, :snp], snpholder, fittednullmodel)
                                copyto!(Z, modelmatrix(testformula, testdf))
                                loadtimevar!(testvec, Z, fittednullmodel)
                                betapval, taupval, jointpval = test!(ts, testvec, testvec)
                            end
                        end
                        println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),",
                        "$maf,$hwepval,$betapval, $taupval, $jointpval")
                    elseif test == :wald
                        if maf == 0 # mono-allelic
                            fill!(γ̂β, 0)
                            fill!(pvalsβ, 1.0)
                            fill!(γ̂τ, 0)
                            fill!(pvalsτ, 1.0)
                        else
                            copyto!(snpholder, @view(genomat[bedrowinds, j]),
                                    impute=true, model=snpmodel)
                            snptodf!(testdf[!, :snp], snpholder, fittednullmodel)
                            altmodel = WSVarLmmModel(fullmeanformula,
                                fittednullmodel.reformula, fullwsvarformula,
                                :id, testdf)
                            altmodel.obswts .= fittednullmodel.obswts
                            WiSER.fit!(altmodel, solver, parallel = parallel, runs = runs)
                            copyto!(γ̂β, 1, altmodel.β, fittednullmodel.p + 1, q)
                            copyto!(γ̂τ, 1, altmodel.τ, fittednullmodel.l + 1, q)
                            copyto!(pvalsβ, 1, coeftable(altmodel).cols[4], fittednullmodel.p + 1, q)
                            copyto!(pvalsτ, 1, coeftable(altmodel).cols[4], 
                                altmodel.p + fittednullmodel.l + 1, q)
                        end
                        if snponly
                            println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$maf,$hwepval,",
                            "$(γ̂β[1]),$(pvalsβ[1]),$(γ̂τ[1]),$(pvalsτ[1])")
                        else
                            print(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$maf,$hwepval,")
                            for j in 1:q
                                print(io, "$(γ̂β[j]),$(pvalsβ[j]),")
                            end
                            for j in 1:q
                                print(io, "$(γ̂τ[j]),$(pvalsτ[j])")
                                if j != q
                                    print(io, ",")
                                end
                            end
                            println(io, "") #end line
                        end
                    end
                end
            end
        end
    elseif analysistype == "snpset"
        # determine type of snpset analysis
        if isa(snpset, Nothing)
            @warn("Nothing set for `snpset`.
            This will default to `singlesnp` analysis (windowsize = 1).")
            setlength = 1
        elseif isa(snpset, AbstractString)
            isfile(snpset) || throw(ArgumentError("snpset file not found,
            to specify a window replace snpset string with a window size"))
            #first column SNPset ID, second column SNP ID
            snpsetFile = CSV.read(snpset, DataFrame, header = [:snpset_id, :snp_id], delim = " ")
            #make sure it matches bim file
            biminfo = CSV.read(bimfile, DataFrame, header = [:chr, :snp_id, :c3, :bp, :c5, :c6], delim = "\t")
            snpsetFile[!, :snp_id] == biminfo[!, :snp_id] || throw(ArgumentError("snp order in snpset file
            must match (in the same order) bimfile"))
            snpset_ids = unique(snpsetFile[!, :snpset_id])
            nSets = length(snpset_ids)
            setlength = 0
        elseif isa(snpset, Integer)
            if snpset == 1
                println("SNP-set length `snpset = 1`, a single-SNP analysis will be run.")
                # single-snp analysis
                vgwas(fittednullmodel, bedfile, bimfile, bedn;
                analysistype = "singlesnp",
                test = test,
                pvalfile = pvalfile,
                snpmodel = snpmodel,
                snpinds = snpinds,
                usespa = usespa,
                bedrowinds = bedrowinds,
                solver = solver,
                parallel = parallel,
                runs = runs,
                verbose = verbose)
                return fittednullmodel
            end
            setlength = snpset
        else #abstract vector (boolean of true at indicies or range or indicies)
            setlength = -1
        end

        # conduct analysis based on type
        if setlength > 0 #single snp analysis or window
            Z = zeros(fittednullmodel.m, setlength) # column counts of genomat
            totalsnps = SnpArrays.makestream(countlines, bimfile)
            SnpArrays.makestream(pvalfile, "w") do io
                if test == :score
                    println(io, "startchr,startpos,startsnpid,endchr,endpos,",
                        "endsnpid,betapval,taupval,jointpval")
                    ts = WSVarScoreTestInvariant(fittednullmodel, setlength, setlength)
                else # wald
                    # #TODO
                    # println(io, "startchr,startpos,startsnpid,endchr,",
                    # "endpos,endsnpid,l2normeffect,betapval,taupval")
                    # nulldev = deviance(fittednullmodel.model)
                    # Xaug = [fittednullmodel.model.X Z]
                    # γ̂ = Vector{Float64}(undef, setlength) # effect size for columns being tested
                end
                SnpArrays.makestream(bimfile) do bimio
                    q = setlength
                    for j in 1:q:totalsnps
                        endj = j + q - 1
                        rowj = readline(bimio)
                        if endj >= totalsnps
                            endj = totalsnps
                            q = totalsnps - j + 1
                            #length of Z will be different
                            if test == :score
                                Z = zeros(fittednullmodel.m, q)
                                ts = WSVarScoreTestInvariant(fittednullmodel, q, q)
                            else
                                # #TODO
                                # Xaug = [fittednullmodel.model.X zeros(size(
                                #     fittednullmodel.mm, 1), q)]
                            end
                        end
                        for i in 1:(q - 2) #
                            readline(bimio)
                        end
                        endj == totalsnps ? rowj_s = rowj : rowj_s = readline(bimio)
                        snpj = split(rowj)
                        snpj_s = split(rowj_s)
                        if test == :score
                            if all(@view(mafs[j:endj]) .== 0) # all mono-allelic, unlikely but just in case
                                betapval = 1.0
                                taupval = 1.0
                                jointpval = 1.0
                            else
                                copyto!(Z, @view(genomat[bedrowinds, j:endj]), impute=true, model=snpmodel)
                                betapval, taupval, jointpval = test!(ts, Z, Z)
                            end
                            println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$(snpj_s[1]),",
                                "$(snpj_s[4]),$(snpj_s[2]),$betapval,$taupval,$jointpval")
                        elseif test == :wald
                            # # TODO
                            # if all(@view(mafs[j:endj]) .== 0) # all mono-allelic, unlikely but just in case
                            #     fill!(γ̂, 0)
                            #     pval = 1.0
                            # else
                            #     copyto!(@view(Xaug[:, (fittednullmodel.model.p+1):end]),
                            #             @view(genomat[bedrowinds, j:endj]),
                            #             impute=true, model=snpmodel)
                            #     altmodel = polr(Xaug, fittednullmodel.model.Y,
                            #         fittednullmodel.model.link, solver,
                            #         wts = fittednullmodel.model.wts)
                            #     copyto!(γ̂, @view(altmodel.β[(fittednullmodel.model.p+1):end]))#, fittednullmodel.model.p + 1, setlength)
                            #     l2normeffect = norm(γ̂)
                            #     pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                            # end
                            # println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),",
                            #     "$(snpj_s[1]),$(snpj_s[4]),$(snpj_s[2]),$l2normeffect,$pval")
                        end
                    end
                end
            end
        elseif setlength == 0 #snpset is defined by snpset file
            # Report effects and pvals in beta and tau for each SNP
            SnpArrays.makestream(pvalfile, "w") do io
                test == :score ? println(io, "snpsetid,nsnps,betapval,taupval,",
                "jointpval") : println(io, "snpsetid,nsnps,l2normeffect,betapval,",
                "taupval,jointpval")
                for j in eachindex(snpset_ids)
                    snpset_id = snpset_ids[j]
                    snpinds = findall(snpsetFile[!, :snpset_id] .== snpset_id)
                    q = length(snpinds)
                    Z = zeros(fittednullmodel.m, q) ### m or nsum?
                    γ̂ = Vector{Float64}(undef, q)

                    if all(@view(mafs[snpinds]) .== 0) # all mono-allelic, unlikely but just in case
                        betapval = 1.0
                        taupval = 1.0
                        jointpval = 1.0
                        l2normeffect = 0.0
                        test == :score ? println(io, "$(snpset_id),$q,$betapval,$taupval,$jointpval") :
                        println(io, "$(snpset_id),$q,$l2normeffect,$betapval,$taupval,$jointpval")
                    elseif test == :score
                        ts = WSVarScoreTestInvariant(fittednullmodel, q, q)
                        copyto!(Z, @view(genomat[bedrowinds, snpinds]), impute=true,
                            model=snpmodel)
                        betapval, taupval, jointpval = test!(ts, Z, Z)
                        ### CHANGE TO WSVAR score test
                        println(io, "$(snpset_id),$q,$betapval,$taupval,$jointpval")
                    elseif test == :wald
                        # # TODO
                        # nulldev = deviance(fittednullmodel.model)
                        # copyto!(@view(Xaug[:, fittednullmodel.model.p+1:end]),
                        #         @view(genomat[bedrowinds, snpinds]),
                        #         impute=true, model=snpmodel)
                        # altmodel = polr(Xaug, fittednullmodel.model.Y,
                        #     fittednullmodel.model.link, solver,
                        #     wts = fittednullmodel.model.wts)
                        # copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                        # l2normeffect = norm(γ̂)
                        # pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                        # println(io, "$(snpset_id),$q,$l2normeffect,$pval")
                    end
                end
            end
        else #setlength == -1 (testing just one set with specified snps in snpset)
            SnpArrays.makestream(pvalfile, "w") do io
                if all(@view(mafs[snpset]) .== 0) # all mono-allelic, unlikely but just in case
                    betapval = 1.0
                    taupval = 1.0
                    jointpval = 1.0
                    l2normeffect = 0.0
                    test == :score ? println(io, "The joint pvalue of snps indexed",
                    " at $(snpset) is $pval") : println(io, "The l2norm of the effect size vector",
                    " is $l2normeffect and joint pvalue of snps indexed",
                    " at $(snpset) is $pval")
                else
                    q = length(snpset)
                    γ̂ = Vector{Float64}(undef, q)
                    Z = zeros(fittednullmodel.m, q) # m or nsum?
                    if test == :score
                        ts = WSVarScoreTestInvariant(fittednullmodel, q, q)
                        copyto!(Z, @view(genomat[bedrowinds, snpset]), impute=true, model=snpmodel)
                        betapval, taupval, jointpval = test!(ts, Z, Z)
                        ### CHANGE TO WSVAR score test
                        println(io, "The joint pvalue of snps indexed",
                         " at $(snpset) is betapval: $betapval, taupval: $taupval, jointpval: $jointpval")
                    elseif test == :wald
                        # TODO
                        # nulldev = deviance(fittednullmodel.model)
                        # Xaug = [fittednullmodel.model.X Z]
                        # copyto!(@view(Xaug[:, fittednullmodel.model.p+1:end]),
                        #         @view(genomat[bedrowinds, snpset]),
                        #         impute=true, model=snpmodel)
                        # altmodel = polr(Xaug, fittednullmodel.model.Y,
                        #     fittednullmodel.model.link, solver,
                        #     wts = fittednullmodel.model.wts)
                        # copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                        # l2normeffect = norm(γ̂)
                        # pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                        # println(io, "The l2norm of the effect size vector",
                        # " is $l2normeffect and joint pvalue of snps indexed",
                        # " at $(snpset) is $pval")
                    end
                end
            end
        end
    else #analysistype == "gxe"
        isnothing(e) &&
            @error("GxE analysis indicated but no environmental variable keyword argument: `e` set.")

        # extra columns in design matrix to be tested
        testdf = DataFrame(fittednullmodel)
        @assert Symbol(e) in propertynames(testdf) "environmental variable $e not found in dataframe"
        testdf[!, :snp] = zeros(fittednullmodel.nsum)

        snpholder = zeros(fittednullmodel.m)
        snpeffectnullbeta = 0.0
        snpeffectnulltau = 0.0
        gxeformula = FormulaTerm(fittednullmodel.meanformula.lhs,
        InteractionTerm(term.((:snp, Symbol(e)))))
        nullmeanformula = FormulaTerm(fittednullmodel.meanformula.lhs,
            sum(union(fittednullmodel.meanformula.rhs, [term(:snp)])))
        nullwsvarformula = FormulaTerm(fittednullmodel.meanformula.lhs,
            sum(union(fittednullmodel.wsvarformula.rhs, [term(:snp)])))
        fullmeanformula = FormulaTerm(fittednullmodel.meanformula.lhs,
            sum(union(fittednullmodel.meanformula.rhs, [term(Symbol(e))],
            [term(:snp)], [InteractionTerm(term.((:snp, Symbol(e))))] )))
        fullwsvarformula = FormulaTerm(fittednullmodel.meanformula.lhs,
            sum(union(fittednullmodel.wsvarformula.rhs, [term(Symbol(e))],
            [term(:snp)], [InteractionTerm(term.((:snp, Symbol(e))))] )))
        SnpArrays.makestream(pvalfile, "w") do io
            if test == :score
                println(io, "chr,pos,snpid,maf,hwepval,snpeffectnullbeta,",
                "snpeffectnulltau,betapval,taupval,jointpval")
                # e may be factor - Z should match dimension
                Z = similar(modelmatrix(FormulaTerm(term(:y), term(Symbol(e))), testdf))
                # create vector of arrays for score test
                q = size(Z, 2)
                testvec = [Matrix{Float64}(undef, ni, q) for
                ni in fittednullmodel.nis]
            else
                γ̂β = 0.0 # effect size for beta gxe effect
                γ̂τ = 0.0 # effect size for tau gxe effect
                snpeffectbeta = 0.0
                snpeffecttau = 0.0
                println(io, "chr,pos,snpid,maf,hwepval,snpeffectbeta,snpeffecttau,",
                "GxEeffectbeta,GxEeffecttau,betapval,taupval")
            end
            SnpArrays.makestream(bimfile) do bimio
                for j in eachindex(snpmask)
                    row = readline(bimio)
                    snpmask[j] || continue
                    hwepval = SnpArrays.hwe(cc[1, j], cc[3, j], cc[4, j])
                    maf = mafs[j]
                    snpj = split(row)
                    if maf == 0 # mono-allelic
                        γ̂β = 0.0
                        γ̂τ = 0.0
                        betapval = 1.0
                        taupval = 1.0
                        jointpval = 1.0
                        snpeffectfullbeta = 0.0
                        snpeffectfulltau = 0.0
                        snpeffectnullbeta = 0.0
                        snpeffectnulltau = 0.0
                    else
                        copyto!(snpholder, @view(genomat[bedrowinds, j]),
                                impute=true, model=snpmodel)
                        snptodf!(testdf[!, :snp], snpholder, fittednullmodel)
                        if test == :score
                            nm = WSVarLmmModel(nullmeanformula,
                                fittednullmodel.reformula, nullwsvarformula,
                                :id, testdf)
                            nm.obswts .= fittednullmodel.obswts
                            @assert nm.ids == fittednullmodel.ids "IDs not matching for GxE."
                            WiSER.fit!(nm, solver, parallel = parallel, runs = runs)
                            snpeffectnullbeta = nm.β[end]
                            snpeffectnulltau = nm.τ[end]
                            copyto!(Z, modelmatrix(gxeformula, testdf))
                            loadtimevar!(testvec, Z, nm)
                            ts = WSVarScoreTest(nm, q, q)
                            betapval, taupval, jointpval = test!(ts, testvec, testvec)
                        elseif test == :wald
                            fullmod = WSVarLmmModel(fullmeanformula,
                                fittednullmodel.reformula, fullwsvarformula,
                                :id, testdf)
                            fullmod.obswts .= fittednullmodel.obswts
                            WiSER.fit!(fullmod, solver, parallel = parallel, runs = runs)
                            γ̂β = fullmod.β[end]
                            γ̂τ = fullmod.β[end]
                            snpeffectbeta = fullmod.β[end-1]
                            snpeffecttau = fullmod.τ[end-1]
                            betapval = coeftable(fullmod).cols[4][fullmod.p]
                            taupval = coeftable(fullmod).cols[4][end]
                        end
                    end
                    if test == :score
                        println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$maf,",
                        "$hwepval,$snpeffectnullbeta,$snpeffectnulltau,",
                        "$betapval,$taupval,$jointpval")
                    else
                        println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$maf,",
                        "$hwepval,$snpeffectbeta,$snpeffecttau,$γ̂β,$γ̂τ,",
                        "$betapval,$taupval")
                    end
                end
            end
        end
    end
    return fittednullmodel
end


# For VCF Analysis
function vgwas(
    fittednullmodel::WSVarLmmModel,
    vcffile::Union{AbstractString, IOStream}, # full path and vcf file name
    nsamples::Integer,          # number of samples in bed file
    vcftype::Symbol;           # :GT = genotype, :DS = dosage
    analysistype::AbstractString = "singlesnp",
    testformula::FormulaTerm = fittednullmodel.meanformula.lhs ~ Term(:snp),
    test::Symbol = :score,
    pvalfile::Union{AbstractString, IOStream} = "vgwas.pval.txt",
    snpmodel::Union{Val{1}, Val{2}, Val{3}} = ADDITIVE_MODEL,
    snpinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    usespa::Bool = true,
    vcfrowinds::AbstractVector{<:Integer} = 1:nsamples, # row indices for VCF array
    solver = Ipopt.IpoptSolver(print_level=0, mehrotra_algorithm = "yes",
    warm_start_init_point="yes", max_iter=100),
    parallel::Bool = false,
    runs::Int = 2,
    verbose::Bool = false,
    snpset::Union{Nothing, Integer, AbstractString, #for snpset analysis
        AbstractVector{<:Integer}} = nothing,
    e::Union{Nothing, AbstractString, Symbol} = nothing # for GxE analysis
    )

    # get number of SNPs in file
    nsnps = nrecords(vcffile)

    # for VCFTools, snpmodel is coded differently
    snpmodel = modelingdict[snpmodel]

    # create SNP mask vector
    if snpinds == nothing
        snpmask = trues(nsnps)
    elseif eltype(snpinds) == Bool
        snpmask = snpinds
    else
        snpmask = falses(nsnps)
        snpmask[snpinds] .= true
    end

    # storage vectors for SPA if it is set to true
    if (usespa == true) & (analysistype == "singlesnp")
        tmp = Array{Float64}(undef, fittednullmodel.m)
        tmp2 = Array{Float64}(undef, 3)
    end
    

    analysistype = lowercase(analysistype)
    analysistype in ["singlesnp", "snpset", "gxe"] || error("Analysis type $analysis invalid option.
    Available options are 'singlesnp', 'snpset' and 'gxe'.")
    # open VCF File
    reader = VCF.Reader(openvcf(vcffile))

    # determine analysis type
    if analysistype == "singlesnp"
        # carry out score or wald test SNP by SNP
        snponly = testformula.rhs == Term(:snp)
        # extra columns in design matrix to be tested
        testdf = DataFrame(fittednullmodel)
        testdf[!, :snp] = zeros(fittednullmodel.nsum)
        if snponly
            q = 1
        else
            Z = similar(modelmatrix(testformula, testdf))
            q = size(Z, 2)
        end


        # create holders for chromome, position, id, dosage/snps
        rec_chr = Array{Any, 1}(undef, 1)
        rec_pos = Array{Any, 1}(undef, 1)
        rec_ids = Array{Any, 1}(undef, 1)
        gholder = zeros(Union{Missing, Float64}, nsamples)
        snpholder = zeros(fittednullmodel.m)

        # carry out score or LRT test SNP by SNP

        SnpArrays.makestream(pvalfile, "w") do io
            if test == :score
                println(io, "chr,pos,snpid,betapval,taupval,jointpval")
                if snponly
                    ts = WSVarScoreTestInvariant(fittednullmodel, 1, 1)
                    if usespa
                        Ks = ecgf(ts)
                    end
                else
                    ts = WSVarScoreTest(fittednullmodel, q, q)
                    testvec = [Matrix{Float64}(undef, ni, q) for
                    ni in fittednullmodel.nis]
                end
            else # wald
                γ̂β = Vector{Float64}(undef, q) # effect size for columns being tested
                pvalsβ = Vector{Float64}(undef, q) # effect size for columns being tested
                γ̂τ = Vector{Float64}(undef, q) # effect size for columns being tested
                pvalsτ = Vector{Float64}(undef, q) # effect size for columns being tested
                if snponly
                    println(io, "chr,pos,snpid,betaeffect,betapval,",
                    "taueffect,taupval")
                    fullmeanformula = FormulaTerm(fittednullmodel.meanformula.lhs,
                    sum(union(fittednullmodel.meanformula.rhs, [testformula.rhs])))
                    fullwsvarformula = FormulaTerm(fittednullmodel.meanformula.lhs,
                    sum(union(fittednullmodel.wsvarformula.rhs, [testformula.rhs])))
                else
                    fullmeanformula = FormulaTerm(fittednullmodel.meanformula.lhs,
                    sum(union(fittednullmodel.meanformula.rhs, testformula.rhs)))
                    fullwsvarformula = FormulaTerm(fittednullmodel.meanformula.lhs,
                    sum(union(fittednullmodel.wsvarformula.rhs, testformula.rhs)))
                    print(io, "chr,pos,snpid,")
                    for j in 1:q
                        print(io, "betaeffect$j,betapval$j,")
                    end
                    for j in 1:q
                        print(io, "taueffect$j,taupval$j")
                        if j != q
                            print(io, ",")
                        end
                    end
                    println(io, "")
                end
            end
            for j in eachindex(snpmask)
                if vcftype == :GT #genotype
                    copy_gt!(gholder, reader; model = snpmodel, impute = true,
                    record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                else #dosage
                    copy_ds!(gholder, reader; model = snpmodel, impute = true,
                    record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                end
                if !snpmask[j] #skip snp, must have read marker still.
                    continue
                end
                if test == :score
                    if snponly
                        copyto!(snpholder, @view(gholder[vcfrowinds]))
                        ps = test!(ts, snpholder, snpholder)
                        betapval, taupval, jointpval = ps
                        if usespa
                            betapval, taupval, jointpval = spa(snpholder, ts, 
                                ps, Ks; tmp_g = tmp, tmp_g2 = tmp2)
                        end
                    else # snp + other terms
                        snptodf!(testdf[!, :snp], @view(gholder[vcfrowinds]), fittednullmodel)
                        copyto!(Z, modelmatrix(testformula, testdf))
                        loadtimevar!(testvec, Z, fittednullmodel)
                        betapval, taupval, jointpval = test!(ts, testvec, testvec)
                    end
                    println(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),",
                    "$betapval,$taupval,$jointpval")
                elseif test == :wald
                    snptodf!(testdf[!, :snp], @view(gholder[vcfrowinds]), fittednullmodel)
                    altmodel = WSVarLmmModel(fullmeanformula,
                        fittednullmodel.reformula, fullwsvarformula,
                        :id, testdf)
                    altmodel.obswts .= fittednullmodel.obswts
                    WiSER.fit!(altmodel, solver, parallel = parallel, runs = runs)
                    copyto!(γ̂β, 1, altmodel.β, fittednullmodel.p + 1, q)
                    copyto!(γ̂τ, 1, altmodel.τ, fittednullmodel.l + 1, q)
                    copyto!(pvalsβ, 1, coeftable(altmodel).cols[4], fittednullmodel.p + 1, q)
                    copyto!(pvalsτ, 1, coeftable(altmodel).cols[4], 
                    altmodel.p + fittednullmodel.l + 1, q)
                    if snponly
                        println(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),",
                        "$(γ̂β[1]),$(pvalsβ[1]),$(γ̂τ[1]),$(pvalsτ[1])")
                    else
                        print(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),")
                        for j in 1:q
                            print(io, "$(γ̂β[j]),$(pvalsβ[j]),")
                        end
                        for j in 1:q
                            print(io, "$(γ̂τ[j]),$(pvalsτ[j])")
                            if j != q
                                print(io, ",")
                            end
                        end
                        println(io, "") #end line
                    end
                end
            end
        end
    elseif analysistype == "snpset"
        # max size of a snpset length
        maxsnpset = 1

        #determine snpset
        if isa(snpset, Nothing)
            setlength = 1
            maxsnpset = 1
        elseif isa(snpset, AbstractString)
            isfile(snpset) || throw(ArgumentError("snpset file not found,
            to specify a window replace snpset string with a window size"))
            #first column SNPset ID, second column SNP ID
            snpsetFile = CSV.read(snpset, DataFrame, header = [:snpset_id, :snp_id], delim = " ")
            maxsnpset = combine(groupby(snpsetFile, :snpset_id), :snp_id => length => :snpset_length) |>
                x -> maximum(x.snpset_length)
            snpset_ids = unique(snpsetFile[!, :snpset_id])
            nSets = length(snpset_ids)
            setlength = 0
        elseif isa(snpset, Integer)
            if snpset == 1
                println("SNP-set length `snpset = 1`, a single-SNP analysis will be run.")
                # single-snp analysis
                vgwas(fittednullmodel, vcffile, nsamples, vcftype;
                analysistype = "singlesnp",
                test = test,
                pvalfile = pvalfile,
                snpmodel = snpmodel,
                snpinds = snpinds,
                vcfrowinds = vcfrowinds,
                solver = solver,
                parallel = parallel,
                runs = runs,
                verbose = verbose)
                return fittednullmodel
            end
            setlength = snpset
            maxsnpset = snpset
        else #abstract vector (boolean of true at indicies or range or indicies)
            setlength = -1
            maxsnpset = count(snpset .!= 0)
        end

        # create holders for chromome, position, id, dosage/gt
        rec_chr = Array{Any, 1}(undef, maxsnpset)
        rec_pos = Array{Any, 1}(undef, maxsnpset)
        rec_ids = Array{Any, 1}(undef, maxsnpset)
        gholder = zeros(Union{Missing, Float64}, nsamples, maxsnpset)

        if setlength > 0 #single snp analysis or window
            Z = zeros(fittednullmodel.m, setlength) #
            q = setlength
            SnpArrays.makestream(pvalfile, "w") do io
                if test == :score
                    println(io, "startchr,startpos,startsnpid,endchr,endpos,",
                    "endsnpid,betapval,taupval,jointpval")
                    ts = WSVarScoreTestInvariant(fittednullmodel, q, q)
                elseif test == :wald
                    # TODO
                    # println(io, "startchr,startpos,startsnpid,endchr,",
                    # "endpos,endsnpid,l2normeffect,betapval,taupval,jointpval")
                    # nulldev = deviance(fittednullmodel.model)
                    # Xaug = [fittednullmodel.model.X Z]
                    # γ̂ = Vector{Float64}(undef, setlength) # effect size for columns being tested
                end
                for j in 1:q:nsnps
                    endj = j + q - 1
                    if endj >= nsnps
                        endj = nsnps
                        q = nsnps - j + 1
                        #length of Z will be different
                        gholder = zeros(Union{Missing, Float64}, nsamples, q)
                        rec_chr = Array{Any, 1}(undef, q)
                        rec_pos = Array{Any, 1}(undef, q)
                        rec_ids = Array{Any, 1}(undef, q)
                        if test == :score
                            Z = zeros(fittednullmodel.m, q)
                            ts = WSVarScoreTestInvariant(fittednullmodel, q, q)
                        elseif test == :wald
                            ## TODO
                            # Xaug = [fittednullmodel.model.X zeros(size(
                            # fittednullmodel.mm, 1), q)]
                        end
                    end
                    if vcftype == :GT #genotype
                        copy_gt!(gholder, reader; model = snpmodel, impute = true,
                        record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                    else #dosage
                        copy_ds!(gholder, reader; model = snpmodel, impute = true,
                        record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                    end
                    if test == :score
                        copyto!(Z, @view(gholder[vcfrowinds, :]))
                        betapval, taupval, jointpval = test!(ts, Z, Z)
                        println(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),",
                        "$(rec_chr[end]),$(rec_pos[end]),$(rec_ids[end][end]),$betapval,$taupval,$jointpval")
                    elseif test == :wald
                        # # TODO
                        # copyto!(@view(Xaug[:, (fittednullmodel.model.p+1):end]),
                        # @view(gholder[vcfrowinds]))
                        # altmodel = polr(Xaug, fittednullmodel.model.Y,
                        # fittednullmodel.model.link, solver,
                        # wts = fittednullmodel.model.wts)
                        # copyto!(γ̂, @view(altmodel.β[(fittednullmodel.model.p+1):end]))#, fittednullmodel.model.p + 1, setlength)
                        # l2normeffect = norm(γ̂)
                        # pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                        # println(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),",
                        # "$(rec_chr[end]),$(rec_pos[end]),$(rec_ids[end][end]),$l2normeffect,$betapval,taupval,jointpval")
                    end
                end
            end
        elseif setlength == 0 #snpset is defined by snpset file
            @warn("This method requires reading in the entire VCF File.
             This can take a lot of memory for large files, as they must be brought into memory.")
            if vcftype == :GT #genotype
                genomat = convert_gt(Float64, vcffile;
                model = snpmodel, impute = true,
                center = false, scale = false)
            else #dosage
                genomat = convert_ds(Float64, vcffile; model = snpmodel,
                key="DS", impute = true, center = false, scale = false)
            end
            SnpArrays.makestream(pvalfile, "w") do io
                test == :score ? println(io, "snpsetid,nsnps,betapval,taupval,jointpval") :
                    println(io, "snpsetid,nsnps,l2normeffect,betapval,taupval,jointpval")
                for j in eachindex(snpset_ids)
                    snpset_id = snpset_ids[j]
                    snpinds = findall(snpsetFile[!, :snpset_id] .== snpset_id)
                    q = length(snpinds)
                    Z = zeros(fittednullmodel.m, q)
                    if test == :score
                        ts = WSVarScoreTestInvariant(fittednullmodel, q, q)
                        copyto!(Z, @view(genomat[vcfrowinds, snpinds]))
                        betapval, taupval, jointpval = test!(ts, Z, Z)
                        println(io, "$(snpset_id),$q,$betapval,$taupval,$jointpval")
                    elseif test == :wald
                        # # TODO
                        # γ̂ = Vector{Float64}(undef, q)
                        # Xaug = [fittednullmodel.model.X Z]
                        # nulldev = deviance(fittednullmodel.model)
                        # copyto!(@view(Xaug[:, fittednullmodel.model.p+1:end]),
                        #         @view(genomat[vcfrowinds, snpinds]))
                        # altmodel = polr(Xaug, fittednullmodel.model.Y,
                        #     fittednullmodel.model.link, solver,
                        #     wts = fittednullmodel.model.wts)
                        # copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                        # l2normeffect = norm(γ̂)
                        # pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                        # println(io, "$(snpset_id),$q,$l2normeffect,$betapval,taupval,jointpval")
                    end
                end
            end
        else #setlength == -1 (testing just one set with specified snps in snpset)
            @warn("This method requires reading in the entire VCF File.
            This can take a lot of memory for large files, as they must be brought into memory.")
            if vcftype == :GT #genotype
                genomat = convert_gt(Float64, vcffile;
                model = snpmodel, impute = true,
                center = false, scale = false)
            else #dosage
                genomat = convert_ds(Float64, vcffile; model = snpmodel,
                key="DS", impute=true, center=false, scale=false)
            end
            SnpArrays.makestream(pvalfile, "w") do io
                q = length(snpset)
                γ̂ = Vector{Float64}(undef, q)
                Z = zeros(fittednullmodel.m, q)
                if test == :score
                    ts = WSVarScoreTestInvariant(fittednullmodel, q, q)
                    copyto!(Z, @view(genomat[vcfrowinds, snpset]))
                    betapval, taupval, jointpval = test!(ts, Z, Z)
                    println(io, "The joint pvalue of snps indexed",
                     " at $(snpset) is betapval: $betapval, taupval: $taupval,",
                     " jointpval: $jointpval")
                elseif test == :wald
                    # #TODO
                    # nulldev = deviance(fittednullmodel.model)
                    # Xaug = [fittednullmodel.model.X Z]
                    # copyto!(@view(Xaug[:, fittednullmodel.model.p+1:end]),
                    #         @view(genomat[vcfrowinds, snpset]))
                    # altmodel = polr(Xaug, fittednullmodel.model.Y,
                    #     fittednullmodel.model.link, solver,
                    #     wts = fittednullmodel.model.wts)
                    # copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                    # l2normeffect = norm(γ̂)
                    # pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                    # println(io, "The l2norm of the effect size vector",
                    # " is $l2normeffect and joint pvalue of snps indexed",
                    # " at $(snpset) is $pval")
                end
            end
        end
    else #analysistype == "gxe"
        isnothing(e) &&
            @error("GxE analysis indicated but no environmental variable keyword argument: `e` set.")

        # extra columns in design matrix to be tested
        testdf = DataFrame(fittednullmodel)
        @assert Symbol(e) in propertynames(testdf) "environmental variable $e not found in dataframe"
        testdf[!, :snp] = zeros(fittednullmodel.nsum)

        # create holders for chromome, position, id
        rec_chr = Array{Any, 1}(undef, 1)
        rec_pos = Array{Any, 1}(undef, 1)
        rec_ids = Array{Any, 1}(undef, 1)
        gholder = zeros(Union{Missing, Float64}, nsamples)

        snpeffectnullbeta = 0.0
        snpeffectnulltau = 0.0
        gxeformula = FormulaTerm(fittednullmodel.meanformula.lhs,
        InteractionTerm(term.((:snp, Symbol(e)))))
        nullmeanformula = FormulaTerm(fittednullmodel.meanformula.lhs,
            sum(union(fittednullmodel.meanformula.rhs, [term(:snp)])))
        nullwsvarformula = FormulaTerm(fittednullmodel.meanformula.lhs,
            sum(union(fittednullmodel.wsvarformula.rhs, [term(:snp)])))
        fullmeanformula = FormulaTerm(fittednullmodel.meanformula.lhs,
            sum(union(fittednullmodel.meanformula.rhs, [term(Symbol(e))],
            [term(:snp)], [InteractionTerm(term.((:snp, Symbol(e))))] )))
        fullwsvarformula = FormulaTerm(fittednullmodel.meanformula.lhs,
            sum(union(fittednullmodel.wsvarformula.rhs, [term(Symbol(e))],
            [term(:snp)], [InteractionTerm(term.((:snp, Symbol(e))))] )))

        SnpArrays.makestream(pvalfile, "w") do io
            if test == :score
                println(io, "chr,pos,snpid,snpeffectnullbeta,",
                "snpeffectnulltau,betapval,taupval,jointpval")
                # e may be factor - Z should match dimension
                Z = similar(modelmatrix(FormulaTerm(term(:y), term(Symbol(e))), testdf))
                # create vector of arrays for score test
                q = size(Z, 2)
                testvec = [Matrix{Float64}(undef, ni, q) for
                ni in fittednullmodel.nis]
            else
                γ̂β = 0.0 # effect size for beta gxe effect
                γ̂τ = 0.0 # effect size for tau gxe effect
                snpeffectbeta = 0.0
                snpeffecttau = 0.0
                println(io, "chr,pos,snpid,snpeffectbeta,snpeffecttau,",
                "GxEeffectbeta,GxEeffecttau,betapval,taupval")
            end
            for j in eachindex(snpmask)
                if vcftype == :GT #genotype
                    copy_gt!(gholder, reader; model = snpmodel, impute = true,
                    record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                else #dosage
                    copy_ds!(gholder, reader; model = snpmodel, impute = true,
                    record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                end
                if !snpmask[j] #skip snp, must read marker still.
                    continue
                end
                # add SNP values to testdf
                snptodf!(testdf[!, :snp], @view(gholder[vcfrowinds]), fittednullmodel)
                #copyto!(testvec, modelmatrix(@formula(trait ~ snp & e), testdf))
                if test == :score
                    nm = WSVarLmmModel(nullmeanformula,
                    fittednullmodel.reformula, nullwsvarformula,
                    :id, testdf)
                    nm.obswts .= fittednullmodel.obswts
                    @assert nm.ids == fittednullmodel.ids "IDs not matching for GxE."
                    WiSER.fit!(nm, solver, parallel = parallel, runs = runs)
                    snpeffectnullbeta = nm.β[end]
                    snpeffectnulltau = nm.τ[end]
                    copyto!(Z, modelmatrix(gxeformula, testdf))
                    loadtimevar!(testvec, Z, nm)
                    ts = WSVarScoreTest(nm, q, q)
                    betapval, taupval, jointpval = test!(ts, testvec, testvec)
                    ### CHANGE to WSVAR score test
                    println(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),",
                        "$snpeffectnullbeta,$snpeffectnulltau,",
                        "$betapval,$taupval,$jointpval")
                elseif test == :wald
                    fullmod = WSVarLmmModel(fullmeanformula,
                        fittednullmodel.reformula, fullwsvarformula,
                        :id, testdf)
                    fullmod.obswts .= fittednullmodel.obswts
                    WiSER.fit!(fullmod, solver, parallel = parallel, runs = runs)
                    γ̂β = fullmod.β[end]
                    γ̂τ = fullmod.β[end]
                    snpeffectbeta = fullmod.β[end-1]
                    snpeffecttau = fullmod.τ[end-1]
                    betapval = coeftable(fullmod).cols[4][fullmod.p]
                    taupval = coeftable(fullmod).cols[4][end]
                    println(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),",
                        "$snpeffectbeta,$snpeffecttau,$γ̂β,$γ̂τ,",
                        "$betapval,$taupval")
                end
            end
        end
    end
    close(reader)
    return fittednullmodel
end

# VCFTools uses different coding for additive, dominant, recessive models than SnpArrays
modelingdict = Dict(
    Val{1}() => :additive,
    Val{2}() => :dominant,
    Val{3}() => :recessive
    )


"""
    snptodf!(dfvec, snpholder, nullmodel)

Takes SNPs in `snpholder` and puts them into the dfvec (df col) based on nullmodel observation counts.
"""
function snptodf!(dfvec::AbstractVector,
                 snpholder::AbstractArray,
                 nullmodel::WSVarLmmModel)
    copyto!(dfvec, vcat(map((x, y) -> fill(x, y), snpholder, nullmodel.nis)...))
end

"""
    loadtimevar!(testvec, longmat, nullmodel)

Takes entries from a matrix/vector of observations across all individuals
and loads them into testvec which is a Vector{Array} based on nullmodel.nis.
"""

function loadtimevar!(testvec::AbstractVector,
                    longmat::AbstractArray,
                    nullmodel::WSVarLmmModel)
    offset = 1
    for i in 1:length(nullmodel.nis)
        ni            = nullmodel.nis[i]
        rangei        = offset:(offset + ni - 1)
        testvec[i]    = @view(longmat[rangei, :])
        offset       += ni
    end
end
