using TestItems, TestItemRunner
@run_package_tests

@testitem "ExtendPdf" begin
    d1 = ExtendPdf((x, _)->x, (1,3))
    d2 = ExtendPdf((x, _)->abs2(x), (1,3))
    d3 = ExtendPdf((x, _) -> xor((x << 2), x) % (x^2 + 1), (1,3))
    
    @test BinnedDistributionFit.scalar_eval(d1, 3) == 3
    @test BinnedDistributionFit.vector_eval(d1, [3, 5]) == [3, 5]

    @test BinnedDistributionFit.scalar_eval(d2, 3) == 9
    @test BinnedDistributionFit.vector_eval(d2, [3, 5]) == [9, 25]

    @test BinnedDistributionFit.scalar_eval(d3, 3) == 5
    @test BinnedDistributionFit.vector_eval(d3, [3, 5]) == [5, 17]
    end

@testitem "SumOfPdfs" begin
    d1 = ExtendPdf((x, _)->x, (1,3))
    d2 = ExtendPdf((x, _)->abs2(x), (1,3))
    d3 = ExtendPdf((x, _) -> xor((x << 2), x) % (x^2 + 1), (1,3))

    @test BinnedDistributionFit.scalar_eval(d1+d2, 3) == 12
    @test BinnedDistributionFit.vector_eval(d1+d2, [3, 5]) == [12, 30]

    @test BinnedDistributionFit.scalar_eval(d1+d3, 3) == 8
    @test BinnedDistributionFit.vector_eval(d1+d3, [3, 5]) == [8, 22]

    @test BinnedDistributionFit.scalar_eval(d2+d3, 3) == 14
    @test BinnedDistributionFit.vector_eval(d2+d3, [3, 5]) == [14, 42]
    end

@testitem "RooFitNLL ExtendPdf" begin
    using FHist
    # https://www.desmos.com/calculator/fhmq5rnwql finds manual solutions to RooFitNLL
    d1 = ExtendPdf((x, _)->x, (1,3))
    h1 = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
    NF1 = RooFitNLL_functor(d1, h1; num_integrator=BinnedDistributionFit.QuadGKIntegrator())
    @test NF1([2.0]) ≈ 1.6827899396467232

    d2 = ExtendPdf((x, _)->abs2(x), (1,3))
    h2 = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
    NF2 = RooFitNLL_functor(d2, h2; num_integrator=BinnedDistributionFit.QuadGKIntegrator())
    @test NF2([2.0]) ≈ 1.84583612533

    d3 = ExtendPdf((x, _)-> exp(x) - abs2(x) + 3, (1,3))
    h3 = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
    NF3 = RooFitNLL_functor(d3, h3; num_integrator=BinnedDistributionFit.QuadGKIntegrator())
    @test NF3([2.0]) ≈ 1.90019114214
end


@testitem "RooFit SumOfPdfs" begin
    using FHist
    d1 = ExtendPdf((x, _)->x, (1,3))
    d2 = ExtendPdf((x, _)->x^2, (1,3))
    sd1 = d1+d2
    h1 = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
    NF1 = RooFitNLL_functor(sd1, h1; num_integrator=BinnedDistributionFit.QuadGKIntegrator())
    @test NF1([[2.0, 0.5], [], []]) ≈ 0.8497340452248006
end


@testitem "_" begin
    using Aqua
    Aqua.test_all(BinnedDistributionFit)
end
