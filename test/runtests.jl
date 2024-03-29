using TrajGWAS
using Test
using CSV, DataFrames, WiSER

filepath = normpath(joinpath(dirname(pathof(WiSER)), "../data/"))
df = DataFrame(CSV.File(filepath * "sbp.csv"), copycols = false)
vlmm = WSVarLmmModel(
    @formula(sbp ~ 1 + agegroup + bmi_std + meds), # removed obswt, gender
    @formula(sbp ~ 1 + bmi_std),
    @formula(sbp ~ 1 + agegroup + bmi_std), # removed obswt, gender, meds
    :id, df
);
WiSER.fit!(vlmm)
df_cut = df[!, [:obswt, :gender, :meds]]
df_cut[!, :gender_bool] = map(x -> x == "Male", df_cut[!, :gender])
df_cut[!, :meds_bool] = map(x -> x == "OnMeds", df_cut[!, :meds])
df_cut = df_cut[!, setdiff(names(df_cut), ["meds", "gender"])];
X1vec = Matrix{Float64}[]
W1vec = Matrix{Float64}[]
for i in unique(df[!, :id])
    push!(X1vec, Matrix{Float64}(df_cut[df[!, :id] .== i, :]))
    push!(W1vec, Matrix{Float64}(df_cut[df[!, :id] .== i, :]))
end

@testset "generic" begin
st = TrajGWAS.WSVarScoreTest(vlmm, 3, 3)
_, v1, _, _ = TrajGWAS.test!(st, X1vec, W1vec)
@test isapprox(v1, 6.236443105947029e-11)

st = TrajGWAS.WSVarScoreTest(vlmm, 1, 1);

# testing for obswt
X1vec_obswt = map(x -> x[:, 1], X1vec)
v1, v2, _, _ = TrajGWAS.test!(st, X1vec_obswt, X1vec_obswt)
@test all(isapprox.((v1, v2),
    (0.4409730648525017, 0.7156750808059059)))

# testing for gender
X1vec_gender = map(x -> x[:, 2], X1vec)
v1, v2, _, _ = TrajGWAS.test!(st, X1vec_gender, X1vec_gender)
@test isapprox(v2, 0.8535798925691855)

# testing for meds. meds is already included in the null model.
st = TrajGWAS.WSVarScoreTest(vlmm, 0, 1)
X1vec_meds = map(x -> x[:, 3], X1vec)
v1, v2, _, _ = TrajGWAS.test!(st, nothing, X1vec_meds)
@test all(isapprox.((v1, v2),
    (-1.0, 1.3935485985991808e-12)))

end

@testset "generic_invariant" begin
for i in 1:length(X1vec)
    for j in 1:size(X1vec[i], 1)
        @assert all(X1vec[i][j, :] .== X1vec[i][1, :]) # in fact, we are dealing with time-invariant features
    end
end
X1 = Matrix{Float64}(undef, length(X1vec), 3)
for i in 1:length(X1vec)
    X1[i, :] .= X1vec[i][1, :]
end
W1 = copy(X1);

st_i = TrajGWAS.WSVarScoreTestInvariant(vlmm, 3, 3)
_, v1, _, _ = TrajGWAS.test!(st_i, X1, W1)
@test isapprox(v1, 6.236443105947029e-11)

st_i = TrajGWAS.WSVarScoreTestInvariant(vlmm, 1, 1);
# obswt
v1, v2, _, _ = TrajGWAS.test!(st_i, X1[:, 1], W1[:, 1])
@test all(isapprox.((v1, v2),
    (0.4409730648525017, 0.7156750808059059)))

# testing for gender
v1, v2, _, _ = TrajGWAS.test!(st_i, X1[:, 2], W1[:, 2])
@test isapprox(v2, 0.8535798925691855)

# testing for meds. meds is already included in the null model.
# meds
st_i = TrajGWAS.WSVarScoreTestInvariant(vlmm, 0, 1);
v1, v2, _, _ = TrajGWAS.test!(st_i, nothing, W1[:, 3])
@test all(isapprox.((v1, v2),
    (-1.0, 1.3935485985991808e-12)))
end

filepath = joinpath(dirname(pathof(TrajGWAS)), "../data/")

pvalpath = "testpval"
@testset "trajgwas_singlesnp_plink" begin
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        filepath * "trajgwas_plinkex.csv",
        filepath * "hapmap3",
        pvalfile = pvalpath,
        usespa = false)
results = CSV.read(pvalpath, DataFrame)
println(mean(results.betapval))
println(mean(results.taupval))
@test all(isapprox.((mean(results.betapval),
    mean(results.taupval)),
    (0.20892290777343045, 0.1905578803635017); atol=1e-3))

#spa 
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        filepath * "trajgwas_plinkex.csv",
        filepath * "hapmap3",
        pvalfile = pvalpath,
        usespa = true)
results = CSV.read(pvalpath, DataFrame)
println(mean(results.betapval))
println(mean(results.taupval))
@test all(isapprox.((mean(results.betapval),
    mean(results.taupval)),
    (0.20841900621550116, 0.18913803849088423); atol=1e-3))
end

@testset "trajgwas_singlesnp_vcf" begin
    trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        filepath * "trajgwas_vcfex.csv",
        filepath * "test_vcf",
        pvalfile = pvalpath,
        geneticformat = "VCF",
        vcftype = :DS,
        usespa = false)
    results = CSV.read(pvalpath, DataFrame)
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval)),
        (0.47580382024875084, 0.4752535666464088); atol=1e-3))

    #spa version
    trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        filepath * "trajgwas_vcfex.csv",
        filepath * "test_vcf",
        pvalfile = pvalpath,
        geneticformat = "VCF",
        vcftype = :DS,
        usespa = true)
    results = CSV.read(pvalpath, DataFrame)

    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval)),
        (0.4709957279235137, 0.4691934843663542); atol=1e-3))
end

@testset "trajgwas_singlesnp_bgen" begin
    trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id, 
        filepath * "trajgwas_bgen_ex.csv", 
        filepath * "example.8bits", 
        usespa = false,
        geneticformat = "BGEN", 
        pvalfile = pvalpath,
        first_dosage=false)
    results = CSV.read(pvalpath, DataFrame)
    println(mean(results.betapval), mean(results.taupval))
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval)),
        (0.46429181963000665, 0.4931759050587741); atol=1e-3))

    #spa version
    trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id, 
        filepath * "trajgwas_bgen_ex.csv", 
        filepath * "example.8bits", 
        usespa = true,
        geneticformat = "BGEN", 
        pvalfile = pvalpath,
        reportchisq=true,
        first_dosage=false)#Adjustor(ones(500, 1)))
    results = CSV.read(pvalpath, DataFrame)
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval)),
        (0.46428062017147764, 0.49230409589144597); atol=1e-3))
end

@testset "trajgwas_snpset_plink" begin
    trajgwas(@formula(y ~ 1 + sex + onMeds),
            @formula(y ~ 1),
            @formula(y ~ 1 + sex),
            :id,
            filepath * "trajgwas_plinkex.csv",
            filepath * "hapmap3",
            pvalfile = pvalpath,
            analysistype = "snpset",
            snpset = 2)
    results = CSV.read(pvalpath, DataFrame)
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval)),
        (0.10426027515117509, 0.08863103713602266); atol=1e-3))

    trajgwas(@formula(y ~ 1 + sex + onMeds),
            @formula(y ~ 1),
            @formula(y ~ 1 + sex),
            :id,
            filepath * "trajgwas_plinkex.csv",
            filepath * "hapmap3",
            pvalfile = pvalpath,
            analysistype = "snpset",
            snpset = filepath * "hapmap_snpsetfile.txt")
    results = CSV.read(pvalpath, DataFrame)
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval)),
        (0.08547206177854592, 0.07865194590141143); atol=1e-3))

    trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex),
        :id,
        filepath * "trajgwas_plinkex.csv",
        filepath * "hapmap3",
        pvalfile = pvalpath,
        analysistype = "snpset",
        snpset = 1:5)
    pvalfile = open(pvalpath)
    pvals = split(readline(pvalfile), r"[, ]")[[end-3, end]]
    close(pvalfile)
    @test all(isapprox.(parse.(Float64, pvals), (9.948443168811026e-6,
    1.2760805146705195e-13), atol=1e-3))
end

@testset "trajgwas_snpset_vcf" begin
    trajgwas(@formula(y ~ 1 + sex + onMeds),
            @formula(y ~ 1),
            @formula(y ~ 1 + sex + onMeds),
            :id,
            filepath * "trajgwas_vcfex.csv",
            filepath * "test_vcf",
            pvalfile = pvalpath,
            geneticformat = "VCF",
            vcftype = :DS,
            analysistype = "snpset",
            snpset = 20)
    results = CSV.read(pvalpath, DataFrame)
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval)),
        (0.25347423101832717, 0.19024751040999272); atol=1e-3))

    trajgwas(@formula(y ~ 1 + sex + onMeds),
            @formula(y ~ 1),
            @formula(y ~ 1 + sex),
            :id,
            filepath * "trajgwas_vcfex.csv",
            filepath * "test_vcf",
            pvalfile = pvalpath,
            geneticformat = "VCF",
            vcftype = :DS,
            analysistype = "snpset",
            snpset = filepath * "snpsetfile_vcf.txt")
    results = CSV.read(pvalpath, DataFrame)
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval)),
        (0.24940609907211, 0.22513112148296); atol=1e-3))

    trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex),
        :id,
        filepath * "trajgwas_vcfex.csv",
        filepath * "test_vcf",
        pvalfile = pvalpath,
        geneticformat = "VCF",
        vcftype = :DS,
        analysistype = "snpset",
        snpset = 1:5)
    pvalfile = open(pvalpath)
    pvals = split(readline(pvalfile), r"[, ]")[[end-3, end]]
    close(pvalfile)
    @test all(isapprox.(parse.(Float64, pvals), (0.079501633468,
    0.01825642314521), atol=1e-3))
end


@testset "trajgwas_snpset_bgen" begin
    trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        filepath * "trajgwas_bgen_ex.csv", 
        filepath * "example.8bits", 
        geneticformat = "BGEN", 
        pvalfile = pvalpath,
        analysistype = "snpset",
        snpset = 20,
        first_dosage=false)
    results = CSV.read(pvalpath, DataFrame)
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval)),
        (0.4444189887865223, 0.43837358376815266); atol=1e-3))

    trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex),
        :id,
        filepath * "trajgwas_bgen_ex.csv", 
        filepath * "example.8bits", 
        geneticformat = "BGEN", 
        pvalfile = pvalpath,
        analysistype = "snpset",
        snpset = filepath * "bgen_snpsetfile.txt",
        first_dosage=false,
        r=2.0,
        adjustor=Adjustor(ones(500, 1)))
        #adjustor=Adjustor(ones(500, 1)))
    results = CSV.read(pvalpath, DataFrame)
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval)),
        (0.43077865459266, 0.38136082362710); atol=1e-3))

    trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex),
        :id,
        filepath * "trajgwas_bgen_ex.csv", 
        filepath * "example.8bits", 
        geneticformat = "BGEN", 
        pvalfile = pvalpath,
        analysistype = "snpset",
        snpset = 1:5,
        first_dosage=false)
    pvalfile = open(pvalpath)
    pvals = split(readline(pvalfile), r"[, ]")[[end-3, end]]
    close(pvalfile)
    @test all(isapprox.(parse.(Float64, pvals), (7.975567153147366e-7,
    7.722715354932478e-20), atol=1e-3))
end

@testset "trajgwas_gxe_plink" begin
trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id,
        filepath * "trajgwas_plinkex.csv",
        filepath * "hapmap3",
        pvalfile = pvalpath,
        analysistype = "gxe",
        e = :sex,
        snpinds = 1:10)
    results = CSV.read(pvalpath, DataFrame)
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval)),
        (0.5917252588216351, 0.4047112266291169); atol=1e-3))
end

@testset "trajgwas_gxe_vcf" begin
    trajgwas(@formula(y ~ 1 + sex + onMeds),
            @formula(y ~ 1),
            @formula(y ~ 1 + sex + onMeds),
            :id,
            filepath * "trajgwas_vcfex.csv",
            filepath * "test_vcf",
            pvalfile = pvalpath,
            geneticformat = "VCF",
            vcftype = :DS,
            analysistype = "gxe",
            e = :sex,
            snpinds = 1:10)
    results = CSV.read(pvalpath, DataFrame)
    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval)),
        (0.6296884957250368, 0.711909935890595); atol=1e-3))
    end

@testset "trajgwas_gxe_bgen" begin
    trajgwas(@formula(y ~ 1 + sex + onMeds),
        @formula(y ~ 1),
        @formula(y ~ 1 + sex + onMeds),
        :id, 
        filepath * "trajgwas_bgen_ex.csv", 
        filepath * "example.8bits", 
        geneticformat = "BGEN", 
        pvalfile = pvalpath,
        analysistype = "gxe",
        e = :sex,
        snpinds = 1:5)
    results = CSV.read(pvalpath, DataFrame)

    @test all(isapprox.((mean(results.betapval),
        mean(results.taupval)),
        (0.6339056748172345, 0.34364220119286115); atol=1e-3))
    @test all(isapprox.((mean(results.snpeffectnullbeta),
        mean(results.snpeffectnulltau)),
        (-0.0144768768325708, 0.005777868071713646); atol=1e-3))
end

rm(pvalpath)
rm("trajgwas.null.txt")
