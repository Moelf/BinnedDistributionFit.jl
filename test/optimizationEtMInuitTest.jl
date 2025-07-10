@testitem "optimization.jl optimization" begin
    using FHist, ComponentArrays, Optimization, ForwardDiff, Distributions, CairoMakie, BinnedDistributionFit
    N = Normal(100.0, 5.0)
    N_gaus_truth = 50000
    dat1 = rand(N, N_gaus_truth)

    E = Exponential(50)
    N_exp_truth = 5 * 10^6
    dat2 = rand(E, N_exp_truth)

    N_total = N_exp_truth + N_gaus_truth
    h1 = Hist1D(append!(dat1, dat2); binedges = 0:200)
    f1(x, ps) = pdf(Exponential(ps[1]), x)
    f2(x, ps) = pdf(Normal(ps[1], ps[2]), x)
    pdfs = BinnedDistributionFit.ExtendPdf(f1, (0, 200)) + BinnedDistributionFit.ExtendPdf(f2, (0, 200))
    NLL = BinnedDistributionFit.LikelihoodSpec(pdfs, h1; num_int = BinnedDistributionFit.SimpleSumIntegrator())

    #para_guess = ComponentArray(norms = [N_total * 0.99, N_total * 0.01], p1 = [30.0], p2 = [70.0, 6])
    para_guess = ComponentArray(norms = [N_total*0.99, N_total*0.01], p1 = [40.0], p2 = [150.0, 3.])

    optf2 = OptimizationFunction(NLL, AutoForwardDiff())
    prob2 = OptimizationProblem(
        optf2, para_guess;
        lb = ComponentArray(norms = [eps(), eps()], p1 = [eps()], p2 = [eps(), 1]),
        ub = ComponentArray(norms = [N_total, N_total], p1 = [100.0], p2 = [200.0, 100.0])
    )
    sol = solve(prob2, Optimization.LBFGS(); maxiters = 500000)

    para_truth = ComponentArray(norms = [N_exp_truth, N_gaus_truth], p1 = [50.0], p2 = [100.0, 5.0])
    @test all(isapprox.(sol.u, para_truth; rtol=0.1))
    @show para_truth
    @show sol.u
end

@testitem "minuit optimization" begin
    using FHist, ComponentArrays, Distributions, CairoMakie, BinnedDistributionFit, Minuit2
    N = Normal(100.0, 5.0)
    N_gaus_truth = 5000
    dat1 = rand(N, N_gaus_truth)

    E = Exponential(50)
    N_exp_truth = 5 * 10^5
    dat2 = rand(E, N_exp_truth)

    N_total = N_exp_truth + N_gaus_truth
    h1 = Hist1D(append!(dat1, dat2); binedges = 0:5:200)
    f1(x, ps) = pdf(Exponential(ps[1]), x)
    f2(x, ps) = pdf(Normal(ps[1], ps[2]), x)
    pdfs = BinnedDistributionFit.ExtendPdf(f1, (0, 200)) + BinnedDistributionFit.ExtendPdf(f2, (0, 200))
    NLL = BinnedDistributionFit.LikelihoodSpec(pdfs, h1; num_int = BinnedDistributionFit.SimpleSumIntegrator())

    para_guess = ComponentArray(norms = [N_total * 0.99, N_total * 0.01], p1 = [20.0], p2 = [70.0, 1.1])
    para_axes = getaxes(para_guess)
    NLL_minuit_wrapper(x) = NLL(ComponentArray(x, para_axes))

    mi = Minuit(
        NLL_minuit_wrapper, getdata(para_guess); strategy = 2, tol = 0.0001,
        limits = [(eps(), N_total), (eps(), N_total), (eps(), 100.0), (eps(), 200.0), (1.0, 10.0)], arraycall = true
    )
    migrad!(mi)
    para_opt = ComponentArray(collect(mi.values), para_axes)
    para_truth = ComponentArray(norms = [N_exp_truth, N_gaus_truth], p1 = [50.0], p2 = [100.0, 5.0])
    @test all(isapprox.(para_opt, para_truth; rtol = 0.1))
end

begin
    using FHist, ComponentArrays, Optimization, ForwardDiff, Distributions, CairoMakie, BinnedDistributionFit
    # Creating our sample data and sampling 5,000 points from it
    L_dist = Laplace(50, 20)
    L_data = rand(L_dist, 50000)

    # Creating our pdf with two unknown parameters, both stored in the vector ps
    f(x, ps) = pdf(Laplace(ps[1], ps[2]), x)

    # Wrapping our data and pdf into the proper form for LikelihoodSpec
    hist = Hist1D(L_data; binedges = 0:100)
    pdf_input = BinnedDistributionFit.ExtendPdf(f, (0, 100))

    # Creating the actual object that will be optimized
    NLL = BinnedDistributionFit.LikelihoodSpec(pdf_input, hist)

    # This will act as our initial guesses for the values of the parameters and an overall norm.
    # Note that each vector element of parameter ComponentArrays must be of type Vector{float}
    para_guess = ComponentArray(norm = [47000.0], p1 = [70.0, 30.0])

    # Defines the optimization function and problem within the Optimization.jl format
    opt_f = OptimizationFunction(NLL, AutoForwardDiff())
    opt_p = OptimizationProblem(opt_f, para_guess)

    sol = solve(opt_p, Optimization.LBFGS())

    # Graphing the data and the optimized function
    fig = BinnedDistributionFit.plotthing(NLL, sol.u)
    @show sol
    @show fig
end
