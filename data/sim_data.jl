using SnpArrays, WiSER, vGWAS, Random, DataFrames, CSV, Statistics, Distributions, LinearAlgebra
using GeneticVariation, VCFTools

# Simulate PLINK data 
hapmap = SnpArray("hapmap3.bed")
famids = CSV.read("hapmap3.fam", header=false)[!, 1]
hapmap = convert(Matrix{Float64}, hapmap, model=ADDITIVE_MODEL, impute=true)
toplist = sortperm(vec(var(hapmap, dims=1))[1:12000], rev=true)[1:30]
# pick a few SNPs with not too much correlation 
lastcor = 1
for i in toplist, j in toplist #find most independent SNPs 
    curcor = abs(Statistics.cor(@view(hapmap[:, i]), @view(hapmap[:, j])))
    if curcor < lastcor
        lastcor = curcor
        println("$i and $j cor is $curcor")
    end
end

m, _ = size(hapmap)
ni = 10
Random.seed!(310)
sex = rand(0:1, m)
timeinvarX = [sex hapmap[:, 7828] hapmap[:, 8250] hapmap[:, 5] hapmap[:, 11793]]
p = 7
l = 7
q = 1

βtrue = [10.0; -3.2; 0.5; 0.5; 1.5; 0.25; 0.75]
τtrue = [0.2; -0.3; 0.5; -0.5; 0.5; 0.25; 0.15]
Σγ    = Matrix(Diagonal([2.0]))
δγω   = [0.0]
σω    = [0.0]

Σγω   = [Σγ δγω; δγω' σω]
Lγω   = cholesky(Symmetric(Σγω), check = false).L
Lγ    = Lγω[1:q, 1:q]
lγω   = Lγω[q + 1, 1:q]
lω    = Lγω[q + 1, q + 1]

γω = Vector{Float64}(undef, q + 1)
z  = similar(γω) # hold vector of iid std normal

df = DataFrame()
obsvec = Vector{WSVarLmmObs{Float64}}(undef, m)
for i in 1:m
    sex_i = timeinvarX[i, 1]
    snp1_i = timeinvarX[i, 2]
    snp2_i = timeinvarX[i, 3]
    snp3_i = timeinvarX[i, 4]
    snp4_i = timeinvarX[i, 5]
    # Xu from a normal dist 0 mean 2 std 
    OnMeds = rand([0, 1], ni)
    # first column intercept, remaining entries iid std normal
    X = Matrix{Float64}(undef, ni, p)
    X[:, 1] .= 1
    X[:, 2] .= sex_i
    X[:, 3] = OnMeds
    X[:, 4] .= snp1_i
    X[:, 5] .= snp2_i
    X[:, 6] .= snp3_i
    X[:, 7] .= snp4_i

    # first column intercept, remaining entries iid std normal
    Z = Matrix{Float64}(undef, ni, q)
    Z[:, 1] .= 1
    # first column intercept, remaining entries iid std normal
    W = Matrix{Float64}(undef, ni, l)
    W[:, 1] .= 1
    W[:, 2] .= sex_i
    W[:, 3] = OnMeds
    W[:, 4] .= snp1_i
    W[:, 5] .= snp2_i
    W[:, 6] .= snp3_i
    W[:, 7] .= snp4_i
    # generate random effects: γω = Lγω * z
    mul!(γω, Lγω, randn!(z))
    # generate y
    μy = X * βtrue + Z * γω[1:q]
    @views ysd = exp.(0.5 .* (W * τtrue .+ dot(γω[1:q], lγω) .+ γω[end]))
    y = ysd .* randn(ni) .+ μy
    # form a VarLmmObs instance
    obsvec[i] = WSVarLmmObs(y, X, Z, W)
    tempdf = hcat(DataFrame(X[:, 2:end], 
        [:sex, :onMeds, :snp1, :snp2, :snp3, :snp4]),
        DataFrame(y = y, id = fill(famids[i], ni))
    )
    df = vcat(df, tempdf)
end

fullmodel = WSVarLmmModel(@formula(y ~ 1 + sex + onMeds + snp1 + snp2 + snp3 + snp4),
              @formula(y ~ 1),
              @formula(y ~ 1 + sex + onMeds + snp1 + snp2 + snp3 + snp4),
              :id,
              df)
nullmodel = WSVarLmmModel(@formula(y ~ 1 + sex + onMeds),
                @formula(y ~ 1),
                @formula(y ~ 1 + sex + onMeds),
                :id,
                df)
WiSER.fit!(fullmodel)
WiSER.fit!(nullmodel)

CSV.write("vgwas_plinkex.csv", df)

vgwas(@formula(y ~ 1 + sex + onMeds),
@formula(y ~ 1),
@formula(y ~ 1 + sex + onMeds),
:id, "vgwas_plinkex.csv", "hapmap3")


vgwas(@formula(y ~ 1 + sex + onMeds),
       @formula(y ~ 1),
       @formula(y ~ 1 + sex + onMeds),
       :id, 
       "vgwas_plinkex.csv", 
       "hapmap3", 
       pvalfile = "snpset.window5.pval.txt", 
       analysistype="snpset", 
       snpset = 10)


# Simulate VCF data



D = convert_ds(Float64, "test_vcf.vcf.gz"; key="DS", impute=true, center=false, scale=false)

# pick a few SNPs with not too much correlation 
using Statistics
toplist = sortperm(vec(var(D, dims=[1])), rev=true)[1:30]
lastcor = 1
for i in toplist, j in toplist #find most independent SNPs 
    curcor = abs(Statistics.cor(D[:, i], D[:, j]))
    if curcor < lastcor
        lastcor = curcor
        println("$i and $j cor is $curcor")
    end
end

people, _ = VCFTools.nsamples("test_vcf.vcf.gz"), nrecords("test_vcf.vcf.gz")
reader = VCF.Reader(openvcf("test_vcf.vcf.gz", "r"))
famids = VCF.header(reader).sampleID
close(reader)
Random.seed!(310)
sex = rand(0:1, people)
m = people

ni = 10
Random.seed!(310)
sex = rand(0:1, m)
timeinvarX = [sex D[:, 86] D[:, 316] D[:, 656]]
p = 6
l = 6
q = 1

βtrue = [10.0; -3.2; 0.5; 0.5; 1.5; 0.25]
τtrue = [0.2; -0.3; 0.5; -0.5; 0.5; 0.25]
Σγ    = Matrix(Diagonal([2.0]))
δγω   = [0.0]
σω    = [0.0]

Σγω   = [Σγ δγω; δγω' σω]
Lγω   = cholesky(Symmetric(Σγω), check = false).L
Lγ    = Lγω[1:q, 1:q]
lγω   = Lγω[q + 1, 1:q]
lω    = Lγω[q + 1, q + 1]

γω = Vector{Float64}(undef, q + 1)
z  = similar(γω) # hold vector of iid std normal

df = DataFrame()
obsvec = Vector{WSVarLmmObs{Float64}}(undef, m)
for i in 1:m
    sex_i = timeinvarX[i, 1]
    snp1_i = timeinvarX[i, 2]
    snp2_i = timeinvarX[i, 3]
    snp3_i = timeinvarX[i, 4]
    # Xu from a normal dist 0 mean 2 std 
    OnMeds = rand([0, 1], ni)
    # first column intercept, remaining entries iid std normal
    X = Matrix{Float64}(undef, ni, p)
    X[:, 1] .= 1
    X[:, 2] .= sex_i
    X[:, 3] = OnMeds
    X[:, 4] .= snp1_i
    X[:, 5] .= snp2_i
    X[:, 6] .= snp3_i

    # first column intercept, remaining entries iid std normal
    Z = Matrix{Float64}(undef, ni, q)
    Z[:, 1] .= 1
    # first column intercept, remaining entries iid std normal
    W = Matrix{Float64}(undef, ni, l)
    W[:, 1] .= 1
    W[:, 2] .= sex_i
    W[:, 3] = OnMeds
    W[:, 4] .= snp1_i
    W[:, 5] .= snp2_i
    W[:, 6] .= snp3_i
    # generate random effects: γω = Lγω * z
    mul!(γω, Lγω, randn!(z))
    # generate y
    μy = X * βtrue + Z * γω[1:q]
    @views ysd = exp.(0.5 .* (W * τtrue .+ dot(γω[1:q], lγω) .+ γω[end]))
    y = ysd .* randn(ni) .+ μy
    # form a VarLmmObs instance
    obsvec[i] = WSVarLmmObs(y, X, Z, W)
    tempdf = hcat(DataFrame(X[:, 2:end], 
        [:sex, :onMeds, :snp1, :snp2, :snp3]),
        DataFrame(y = y, id = fill(famids[i], ni))
    )
    df = vcat(df, tempdf)
end

fullmodel = WSVarLmmModel(@formula(y ~ 1 + sex + onMeds + snp1 + snp2 + snp3),
              @formula(y ~ 1),
              @formula(y ~ 1 + sex + onMeds + snp1 + snp2 + snp3),
              :id,
              df)
nullmodel = WSVarLmmModel(@formula(y ~ 1 + sex + onMeds),
                @formula(y ~ 1),
                @formula(y ~ 1 + sex + onMeds),
                :id,
                df)
WiSER.fit!(fullmodel)
WiSER.fit!(nullmodel)


CSV.write("vgwas_vcf_ex.csv", df)

vgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id, 
        "vgwas_vcf_ex.csv", 
        "test_vcf", 
        geneticformat = "VCF", 
        vcftype = :DS,
        pvalfile = "vgwas_vcfex.pval.txt")



