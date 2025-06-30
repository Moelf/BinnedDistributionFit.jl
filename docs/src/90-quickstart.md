# Quicker Start

```julia
using FHist, Optimize, ForwardDiff

d(x, ps...) =   # <-- Enter function here with some number of parameters in the list ps
h =             # <-- Enter histogram here as an FHist histogram

support = extrema(binedges(h1))

NLL = BinnedDistributionFit.RooFitNLL_functor(ExtendPdf(d, support), h1)
NLL_wrapper(x,_) = NLL(x)

optf = OptimizationFunction(NLL_wrapper, AutoForwardDiff())
prob = OptimizationProblem(optf, [integral(h); ps])
sol = solve(prob, Optimization.LBFGS())
sol
# Answer comes out here
```

# [Quick Start](@id quickstart)

## Contents

```@contents
Pages = ["90-quickstart.md"]
```

## Introduction
This page will act as a quick start guide to the BinnedDistributionFit package. To create a fit for one function, first create an `ExtendPdf` which contains the desired function. We define both the value of the function and its support. \[Functions can only be evaluated on their support.\] It is expected that the function input has some number of free parameters, labeled as `ps` throughout the documentation, though they are not necessary in the creation of an `ExtendPdf`.

Once an `ExtendPdf` is created, we can then evaluate it at a single point using the [`scalar_eval`](@ref BinnedDistributionFit.scalar_eval) function. To evaluate a pdf at multiple points, we can use [`vector_eval`](@ref BinnedDistributionFit.vector_eval) instead.

```julia
f(x) = x^2 + x*2
d1 = BinnedDistributionFit.ExtendPdf(f, (1,2))
scalar_eval(d1, 2)
# Output: 8
vector_eval(d1, [1, 2])
# Output: [3, 8]

g(x, ps) = x^ps[1] + x*ps[2]
d2 = BinnedDistributionFit.ExtendPdf(g, (1, 2))
scalar_eval(d1, 2, [2, 3])
# Output: 10
vector_eval(d1, [1, 2], [2, 3])
# Output: [4, 10]
```

## NLL Loss Function
This package uses the FHist package to create histograms. For more information regarding the FHist package, visit [link]. Using. `Hist1D`, we score a certain function on its fit to the data using several methods. The simplest method is using negative log likelihood scoring, as it does not require the standard deviation of each histogram bin. 

We use `RooFitNLL_functor` to create a functor which can then be instantiated with particular parameter values. We call instances of the functor with the syntax `[norm, ps...]`, where `norm` is the norm of the function, and `ps...` contains all of the free parameters of the function.
```julia
d = ExtendPdf((x, _) -> x, (1, 3)) 
h = Hist1D(; binedges = 1:3, bincounts = [2.0, 4.0])
NF = RooFitNLL_functor(d, h; num_integrator=BinnedDistributionFit.QuadGKIntegrator())
```
The method of integration can be chosen. By default, it is set to `SimpleSumIntegrator()`.

In this case, as there are no parameters `ps` specified when `d` was created, our input is simply the norm of `d`.
```julia
NF([2.0])
# Output: 1.6827899396467232
```
For more details on the math behind `Roo_FitNLL`, see [link].

## Using Multiple Functions

In addition to the `ExtendPdf` struct, `BinnedDistributionFit` also uses the `SumOfPdfs` struct to handle multiple input to be evaluated over the same dataset. The same `ExtendPdf` syntax is used, however internally whenever two `ExtendPdf` are summed the output is a `SumOfPdfs`.

```julia
d1 = ExtendPdf((x, _) -> x, (1,3))
d2 = ExtendPdf((x, _) -> x^2, (1,3))

h = Hist1D(; binedges = 1:3, bincounts = [2.0, 4.0], sumw2 = [2.0, 4.0])
NF = RooFitNLL_functor(d1 + d2, h; num_integrator = BinnedDistributionFit.QuadGKIntegrator())

NF([[2.0, 0.5], [], []])
# Output: 0.8497340452248006
```

Additionally, when instantiating the functor the input syntax is different from the previous section. A `SumOfPdfs` functor `NF` expects an array of arrays as an input, the first element of which contains the norm of each function, and each subsequent array contains each parameter of a function, like so:
```julia
NF([[norm1, norm2], [ps1...], [ps2...]])
```

## Chi Squared Loss Function
Another loss function implemented into the `BinnedDistributionFit` is the chi squared function (referenced as `chi2`). The NLL and chi squared loss functions are not commensurable. Only one should be used for analysis of a given data set. 

Similarly to the `RooFitNLL_functor`, we can call the `chi2_functor` like so:
```julia
d = ExtendPdf((x, _) -> x-0.5, (1,3))
h = Hist1D(; binedges = 1:3, bincounts = [1, 2], sumw2 = [1, 1])
x = BinnedDistributionFit.chi2_functor(d, h; num_integrator = BinnedDistributionFit.QuadGKIntegrator())

x([3])
# Output: 0
```
Note that the `sumw2` field must contain only nonzero values or else the value of the functor will be `Inf`. 

Just as with the `RooFitNLL_functor`, we can also apply the `chi2_functor` to multiple pdfs.
```julia
d1 = ExtendPdf((x, _) -> x-0.5, (1, 4))
d2 = ExtendPdf((x, _) -> 3x^2, (1, 4))

h = Hist1D(; binedges = 1:4, bincounts = [1, 2, 3], sumw2 = [1, 2, 1])
x = BinnedDistributionFit.chi2_functor(d1 + d2, h; num_integrator = BinnedDistributionFit.QuadGKIntegrator())

x([[3, 4], [], []])
# Output: 1.73774816049
```