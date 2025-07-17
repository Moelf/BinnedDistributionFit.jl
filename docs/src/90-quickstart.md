```@contents
Pages = ["showcase.md"]
```

# [Quick Start](@id quickstart)
To optimize the fit of one function over a data set, 
```julia
using BinnedDistributionFit, FHist, Optimize, ForwardDiff
# User inputs
pdf(x, ps...) = x^ps[1] + ps[2]                       # <-- Enter function here
hist = Hist1D(; bincounts = 1:100, binedges = 1:100)   # <-- Enter histogram here

support = extrema(binedges(hist))

NLL = LikelihoodSpec(ExtendPdf(pdf, support), hist)

# Guess of the values of ps
ps = [1, 2]
# Optimization implementation
optf = OptimizationFunction(NLL_wrapper, AutoForwardDiff())
prob = OptimizationProblem(optf, ComponentArray(norms = [sum(bincounts(hist))], p1 = ps))
sol = solve(prob, Optimization.LBFGS())

@show sol.u
```

# Introduction
This page will act as a quick start guide to the BinnedDistributionFit package. To create a fit for one function, first create an `ExtendPdf` which contains the desired function. We define both the value of the function and its support. It is expected that the function input has some number of free parameters, labeled as `ps` throughout the documentation, though they are not necessary in the creation of an `ExtendPdf`.

Once an `ExtendPdf` is created, we can then evaluate it at a single point using the [`scalar_eval`](@ref BinnedDistributionFit.scalar_eval) function. To evaluate a pdf at multiple points, we can use [`vector_eval`](@ref BinnedDistributionFit.vector_eval) instead.

```julia
f(x) = x^2 + x*2
d1 = ExtendPdf(f, (1,2))
scalar_eval(d1, 2)
# Output: 8
vector_eval(d1, [1, 2])
# Output: [3, 8]

g(x, ps) = x^ps[1] + x*ps[2]
d2 = ExtendPdf(g, (1, 2))
scalar_eval(d1, 2, [2, 3])
# Output: 10
vector_eval(d1, [1, 2], [2, 3])
# Output: [4, 10]
```

## NLL Loss Function
`BinnedDistributionFit` uses the FHist package to create histograms. Using. `Hist1D`, we score a certain function on its fit to the data using several methods. The simplest method is using negative log likelihood scoring, as it does not require the standard deviation of each histogram bin. 

We create a [`LikelihoodSpec`](@ref BinnedDistributionFit.LikelihoodSpec) with our data and pdf, which we can then call with specific values of the parameters of pdf.
```julia
pdf = ExtendPdf((x, _) -> x, (1, 3)) 
hist = Hist1D(; binedges = 1:3, bincounts = [2.0, 4.0])
NLL = LikelihoodSpec(
    pdf, hist; 
    num_integrator = BinnedDistributionFit.QuadGKIntegrator()
)
```
The method of integration can be specified. By default, it is set to `SimpleSumIntegrator()`.

In this case, as there are no parameters `ps` specified when `pdf` was created, we call `NLL` with just the norm of `pdf`. Either a vector of vectors or a `ComponentArray` may be used to call `NLL`.
```julia
NLL([[2.0]])
# Output: 1.6827899396467232

NLL(ComponentArray(norm = [2.0]))
# Output: 1.6827899396467232
```

## Using Multiple Functions

In addition to the `ExtendPdf` struct, `BinnedDistributionFit` also uses the `SumOfPdfs` struct to handle multiple input to be evaluated over the same dataset. The same `ExtendPdf` syntax is used. Internally, whenever two `ExtendPdf` are summed the output is a `SumOfPdfs`.

The order of arguments for both `NLL` methods (the one taking a vector of vectors and a `ComponentArray`) is the same: `norms, pN...`, where `pN` is the Nth set of parameters.
```julia
pdf1 = ExtendPdf((x, _) -> x, (1,3))
pdf2 = ExtendPdf((x, _) -> x^2, (1,3))

hist = Hist1D(; binedges = 1:3, bincounts = [2.0, 4.0])
NLL = LikelihoodSpec(
    pdf1 + pdf2, hist; 
    num_integrator = BinnedDistributionFit.QuadGKIntegrator()
)

NLL([[2.0, 0.5], [], []])
# Output: 0.8497340452248006

NLL(ComponentArray(norms = [2.0, 0.5], p1 = [0], p2 = [0]))
# Output: 0.8497340452248006
```

## Chi Squared Loss Function
Another loss function implemented into the `BinnedDistributionFit` is the chi squared function (referenced as `chi2`). The NLL and chi squared loss functions are not commensurable. Only one should be used for analysis of a given data set. 

We use the keyword argument `loss_type` when defining the `LikelihoodSpec` to specify chi squared as the loss type. The default of `loss_type` is `NLL`.
```julia
pdf = ExtendPdf((x, _) -> x-0.5, (1,3))
hist = Hist1D(; binedges = 1:3, bincounts = [1, 2], sumw2 = [1, 1])
x = LikelihoodSpec(
    pdf, hist; 
    loss_type = BinnedDistributionFit.CSQ(), 
    num_integrator = BinnedDistributionFit.QuadGKIntegrator()
)

x([3.0])
# Output: 0

x(ComponentArray(norm = [3.0]))
# Output: 0
```
Note that the `sumw2` field must contain only nonzero values or else the value of the `LikelihoodSpec` will be `Inf`. 

Example with multiple pdfs:
```julia
pdf1 = ExtendPdf((x, _) -> x-0.5, (1, 4))
pdf2 = ExtendPdf((x, _) -> 3x^2, (1, 4))

hist = Hist1D(; binedges = 1:4, bincounts = [1, 2, 3], sumw2 = [1, 2, 1])
x = BinnedDistributionFit.chi2_functor(
    pdf1 + pdf2, hist; 
    loss_type = BinnedDistributionFit.CSQ(), 
    num_integrator = BinnedDistributionFit.QuadGKIntegrator()
)

x([[3.0, 4.0], [], []])
# Output: 1.73774816049

x(ComponentArray(norms = [3.0, 4.0], p1 = [0], p2 = [0]))
# Output: 1.73774816049
```

## Usage with Optimization.jl
`BinnedDistributionFit` can integrate with various optimization packages. This example shows how to fit and graph one pdf using `Optimization` and `CairoMakie `. First, we can create sample data data and choose 5,0000 random points from it.
```@setup Plotting_Example
using Optimization, ForwardDiff, Distributions, CairoMakie, BinnedDistributionFit
using ComponentArrays, FHist
CairoMakie.activate!(; type="svg")
```
```@example Plotting_Example; continued = true
L_dist = Laplace(50, 20)
L_data = rand(L_dist, 50000)
```
Since we know the shape of the data, we can use a Laplace distribution as our function. We then create a pdf with two unknown parameters, both stored in the vector ps.
```@example Plotting_Example; continued = true
f(x, ps) = pdf(Laplace(ps[1], ps[2]), x)
```
We then use the core components of the `BinnedDistributionFit` package: `ExtendPdf` and `LikelihoodSpec`. We create a histogram to wrap the data we generated and then wrap that and the `ExtendPdf` into `LikelihoodSpec`, our loss function.
```@example Plotting_Example; continued = true
hist = BinnedDistributionFit.Hist1D(L_data; binedges = 0:100)
pdf_input = BinnedDistributionFit.ExtendPdf(f, (0,100))

NLL = BinnedDistributionFit.LikelihoodSpec(pdf_input, hist) 
```
We then generate an initial guesses for the values of the parameters and an overall norm. Note that each element of the 'ComponentArray' must be of type `Vector{float}`.
```@example Plotting_Example; continued = true
para_guess = ComponentArray(norm = [47000.], p1 = [70., 30.])
```
Using the `Optimization` package we define the function and problem to optimize.
```@example Plotting_Example; continued = true
optf = OptimizationFunction(NLL, AutoForwardDiff())
prob = OptimizationProblem(
    optf, para_guess;
    # lower and upper bounds on the parameters
    lb = ComponentArray(norm = [eps()], p1 = [eps(), eps()]), 
    ub = ComponentArray(norm = [100000], p1 = [100, 50])
)

sol = solve(prob, Optimization.LBFGS())
```
Finally, we can use the `BinnedDistributionFit` extension of `Makie` to plot the fitted pdf and the histogram.
```@example Plotting_Example
fig = BinnedDistributionFit.plotthing(NLL, sol.u)

fig # hide
```