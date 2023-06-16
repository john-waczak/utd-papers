# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.6
#   kernelspec:
#     display_name: Julia 1.9.1
#     language: julia
#     name: julia-1.9
# ---

using Pkg
Pkg.activate(".") 
Pkg.status()

Pkg.add(url="https://github.com/john-waczak/TimeSeriesTools.jl")

using CSV, DataFrames 
using Plots 
using TimeSeriesTools
using ParameterHandling
using Dates, TimeZones
using Unitful
using Markdown, LaTeXStrings
using Statistics, StatsBase
using BenchmarkTools

# # Download Data from OSN

# +
df1 = CSV.File(download("https://ncsa.osn.xsede.org/ees230012-bucket01/AirQualityNetwork/data/raw/Central_Hub_1/2023/03/04/MINTS_001e06318c91_IPS7100_2023_03_04.csv")) |> DataFrame
CSV.write("data/central_hub_1_ips7100__1.csv", df1)

df2 = CSV.File(download("https://ncsa.osn.xsede.org/ees230012-bucket01/AirQualityNetwork/data/raw/Central_Hub_1/2023/03/05/MINTS_001e06318c91_IPS7100_2023_03_05.csv")) |> DataFrame
CSV.write("data/central_hub_1_ips7100__2.csv", df2)

df3 = CSV.File(download("https://ncsa.osn.xsede.org/ees230012-bucket01/AirQualityNetwork/data/raw/Central_Hub_1/2023/03/06/MINTS_001e06318c91_IPS7100_2023_03_06.csv")) |> DataFrame
CSV.write("data/central_hub_1_ips7100__3.csv", df3)

df = vcat(df1, df2, df3)
# -

nrow(df)

# verify smooth transition between files"
println(df1.dateTime[end], "\t", df2.dateTime[1])
println(df2.dateTime[end], "\t", df3.dateTime[1])

# # Exploratory data analysis

# First, let's take our donloaded data and construct some timeseries to analyze: 

describe(df)

md""" 
The starting datetime for this dataset is $(df.dateTime[1]).
"""

pm2_5 = RegularTimeSeries(
    df.pm2_5, 
    1.0, 
    u"μg/m^3", 
    u"s",
    ZonedDateTime(2023, 3, 4, 0, 0, 0, 035, tz"UTC")
)

# Here I demonstrate a variety of plot customizations. Also note that we can automatically convert the units by using the broadcasted pipe `|>` operator

# +
myfont = Plots.font("sans-serif", :black)

default(
    linewidth=2,
    grid=true,
    tickdirection=:out,
    minorgrid=true,
    gridwidth=2,
    minorgridwidth=1,
    gridalpha=0.3,
    minorgridalpha=0.25,
    margin=5*Plots.mm,
    topmargin=2*Plots.mm,
    framestyle=:box,
    titlefont=myfont,
    guidefont=myfont,
    legendfont=myfont,
    tickfont=myfont,
    titlefontsize=16,
    labelfontsize=13,
    tickfontsize=11,
    colorbar_titlefontsize=13,
)

#scalefontsizes(0.7)  # scales all font sizes by 30 %
# -

p1 = plot(
    times(pm2_5) * pm2_5.t_units .|> u"hr", 
    pm2_5.z * pm2_5.z_units,
    xlabel="time",
    ylabel="Concentration",
    title="Time series for PM 10 starting on $(Date(pm2_5.start_time))",
    label="",
    xlims=(0, Inf),  # lims only considers numbers and will autoselect for any Inf supplied.
    ylims=(0, Inf),
)

savefig("figures/time-series.png")
savefig("figures/time-series.pdf")
savefig("figures/time-series.svg")

# # Uncertainty Estimation

# ## Method 1: Representativeness uncertainty via deviation in a roling window
#
# The first way we seek to estimate the uncertainty for our time series is by evaluating the standard or mean deviation of values within a rolling window centered at each point. If our time series is denoted by $Z(t)$, then we compute
#
# \begin{equation}
#     \Delta Z_{\text{std}} \approx \sqrt{\dfrac{\sum\limits_{(t-\Delta t/2) \leq t\_i \leq (t + \Delta t/2}^{N} \big(Z(t_i)-\mu_i\big)^2}{N}}
# \end{equation}
#
# for the standard deviation or
#
# \begin{equation}
#     \Delta Z_{\text{mean}} \approx \dfrac{\sum\limits_{(t-\Delta t/2) \leq t\_i \leq (t + \Delta t/2}^{N} \big\lvert Z(t_i)-\mu_i \big\rvert }{N}
# \end{equation}
#
#
# where a key hyperparameter for this method is the window size $\Delta t$. 
#
# Let's first begin by creating a function to return the appropriate indices for our centered window. For the edge cases, we can either shift the window to maintain it's length, or we can can drop points from the time series which don't have enough neighbors to compute the deviation.

# for simplicity, let N be an odd number: 
N = 15 
@assert N < length(pm2_5)
println("Window duration is $(N*pm2_5.Δt*pm2_5.t_units)")

# +
function window_idxs(i, N, Z::AbstractVector)
    if N > length(Z)
        thow(ArgumentError("N must be less than length of data"))
    end
    
    n = Int((N-1)/2)

    idx_out = i-n:i+n
    
    # handle edge cases
    if idx_out[1] < 1
        offset = 1-idx_out[1]
        idx_out = idx_out .+ offset
    end
    
    if idx_out[end] > length(Z)
        offset = idx_out[end] - length(Z)
        idx_out = idx_out .- offset
    end
    
    return idx_out
end

window_idxs(i, N, Z::AbstractRegularTimeSeries) = window_idxs(i, N, Z.z)


idxs1 = window_idxs(1, N, pm2_5)
idxs2 = window_idxs(length(pm2_5), N, pm2_5)
pm2_5.z[idxs1]
pm2_5.z[idxs2]
# -

# Now we should be able to easily compute the representativeness 

# +
function rolling_deviation(Z::AbstractRegularTimeSeries, Nwindow::Int; func=std)
    Δz = zeros(length(Z))
    for i ∈ 1:length(Z)
        Δz[i] = func(Z.z[window_idxs(i,Nwindow,Z)])
    end
    return Δz
end

@benchmark rolling_deviation(pm2_5, N)

# +
function mean_dev(Z::AbstractVector)
    μ = mean(Z)
    return sum(abs.(Z .- μ))/length(Z)
end


@benchmark rolling_deviation(pm2_5, N, func=mean_dev)
# -

rolling_deviation(pm2_5, 61)

# Now let's plot the original data with a ribon for the new uncertainty. Another question to address is the number of deviations to use for our uncertainty window... 

p1 = plot(
    times(pm2_5)* pm2_5.t_units .|> u"hr",
    pm2_5.z * pm2_5.t_units,
    ribbon = 2*rolling_deviation(pm2_5, 5*60+1; func=mean_dev),  # 5 minute window
    xlabel="time",
    ylabel="Concentration",
    title="Time series for PM 10 on $(Date(pm2_5.start_time))",
    label="",
    xlims=(0, 0.5),  # lims only considers numbers and will autoselect for any Inf supplied.
    ylims=(0, 0.5), 
    lw=3
)

# ## Model 2: Uncertainty via Semivariogram
#
# Alternatively, we can compute the semivariogram $\gamma$ which measures the expected variance in our data as a function of time lag $h$. Extrapolating to a lag of $h=0$ then provides a measure of the intrinsic uncertainty in our measurements.

γ, h = semivariogram(pm2_5; lag_max=60*60)

# +
γ_params = (
    nugget=positive(0.01),
    sill=positive(0.1),
    range=positive(100.0)
)

γ_fit = fit_spherical_γ(h, γ, γ_params)

# +
scatter(
    h * pm2_5.t_units .|> u"minute",
    γ,
    xlabel=L"\Delta t",
    ylabel=L"\gamma(\Delta t)",
    title="IPS7100 semivariogram for $(Date(pm2_5.start_time))",
    label=L"Empirical $\gamma$",
    ms=2,
    msw=0,
    xlims=(-0.075, Inf),
    ylims=(0, Inf)
)

h_fit = 0.0:1.0:h[end]
γ_out = γ_fit.(h_fit)

plot!(
    h_fit * pm2_5.t_units .|> u"minute",
    γ_out,
    label="Spherical Model",
    lw=3,
)

scatter!([0], [nugget(γ_fit)], label="nugget", color=:grey)
vline!([γ_range(γ_fit) * pm2_5.t_units |> u"minute"], ls=:dash, color=:grey, label="range")
hline!([sill(γ_fit) + nugget(γ_fit)], ls=:dash, color=:brown, label="sill")

# -

#| echo: false
md"""
So what do we do with this information? The nugget provides an estimate of the expected variance of our data for a time lag $\Delta t = 0$. 
Consequently, the square root of the nugget (in this case $(round(sqrt(nugget(γ_fit)), digits=3))) provides an intrinsic estimation of the uncertainty for
our measurements. Alternatively, we can use the range as an estimate for the relevant timescale in our dataset and construct our rolling
window deviation using that discovered duration.
"""

Nwindow = round(Int, γ_range(γ_fit)/pm2_5.Δt)
if Nwindow % 2 == 0 
    Nwindow += 1
end

p1 = plot(
    times(pm2_5)* pm2_5.t_units .|> u"hr",
    pm2_5.z * pm2_5.t_units,
    ribbon = 2*rolling_deviation(pm2_5, Nwindow; func=mean_dev),  # 5 minute window
    xlabel="time",
    ylabel="Concentration",
    title="Time series for PM 10 on $(Date(pm2_5.start_time))",
    label="",
    xlims=(0, 0.5),  # lims only considers numbers and will autoselect for any Inf supplied.
    ylims=(0, 0.5), 
    lw=3
)

# # Chaotic Systems: HAVOK Analysis
#
# Now that we have established methods for evaluating the uncertainty of naked time series data, let's continue by examining an
# interesting method for evaluating chaotic external forcing via a natural embedding of our time series into an $N$-dimensional 
# manifold. We use [this paper](https://arxiv.org/pdf/1608.05306.pdf) as the model for this approach.

#H = TimeDelayEmbedding(pm2_5; method=:backward)
H = TimeDelayEmbedding(pm2_5; nrow=Nwindow)

Hdemo = TimeDelayEmbedding(1:10; nrow=5, method=:forward)
Hdemo = Int.(Hdemo)
Hprint = ["z_$(Hdemo[i,j])" for i ∈ axes(Hdemo,1), j ∈ axes(Hdemo,2) ]

# With my choice of `method=:forward`, we have defined our embedding to construct our a matrix of the following form: 

Hprint

# Therefore, we see that our embedding maps the time series point 
#
# \begin{equation}
#     Z_i \mapsto [Z_i, Z_{i+1}, ..., Z_{i+n}]
# \end{equation}
#
# where $n$ is the desired output dimension (controlled by the `nrows` argument). We therefore interpret each column of our embedding matrix (aka *Hankel* matrix) as a vector in $\mathbb{R}^n$. Similarly, each row of our embedding matrix $H$ can be thought of as an *offset* or *delayed* timeseries giving the trajectory of the i$^{th}$ delay coordinate.
#
# With these embeddings collected into a matrix, a natural question is what should the eigenvetors and eigenvalues represent? Since the $H$ is not square, the next best thing is to look at the singular value decomposition: 
#
# \begin{equation}
#     H = U \Sigma V^\dagger
# \end{equation}

using LinearAlgebra

# +
U, σ, V = svd(H)

println("U has dimensions ", size(U))
println("V has dimensions ", size(V))
println("we have $(length(σ)) singular values")
# -

# we can verify that this decompostion is accurate: 

all(H .≈ U*Diagonal(σ)*V')

# The columns of $V$ represent timeseries of the *eigenmodes* defined by the columns of $U\Sigma$ in our data. We can visualize the first 3 of these to get a sense for the shape of our learned attractor: 

# +
plot_times =  collect(1:size(V,1)) ./ (60.0)^2
p1 = scatter(
    V[:,1], V[:,1], V[:,3],
    ms=2,
    msw=0,
    msa=0,
    marker_z = plot_times,
    frame=:none,
    ticks=nothing,
    xlabel="",
    ylabel="",
    zlabel="",
    label="",
    cbar=false,
    margins=0*Plots.mm,
    background_color=:transparent
)

#p2 = plot([NaN], lims=(0,1), framestyle=:none, label="")
#p3 = heatmap([NaN;;], framestyle=:none, clims=extrema(plot_times), cbar=true, colorbar_fontsize=8)
#l = @layout [a{0.95w} b]
#plot(p1,p2,p3, layout=l, margins=0Plots.mm)
# -

# Excellent! So we see that in this 3 dimensional view of the embeddeing a lot of time is spent in the center region with ocasional devations likely representing switching between attractor fixed points. 
#
# NOTE: Is this strictly true? I think there is a threorem about the number of attractors of a time-delay embedding...

# + active=""
# Key to HAVOK analysis is the idea that there is a cutoff value $r$ below which the dynamics seen in the columns of $V$ can be well modeled as a linear syste with external forcing. Let's look at the singular values and see if we can observe a natural cuttoff point: 
# -

plot(σ./sum(σ), xlabel="index, i", ylabel="σᵢ/(Σᵢσᵢ)", lw=3, label="")

# Great, with this in mind, let's define a simple function to identify a cutoff point once this ratio falls below a certain threshold. This is similar to the idea of truncating the output of PCA by the explained variance. 

# +
function r_cutoff(σ; ratio=0.01, rmax=15)
    return min(sum(σ ./ sum(σ) .> ratio) + 1, rmax)
end

r = r_cutoff(σ,rmax=20)
# -

# Now we can truncate our matrices

# truncate the matrices
Vr = @view V[:,1:r]
Ur = @view U[:,1:r]
σr = @view σ[1:r]

# ## Fitting the HAVOK model
#
# The HAVOK model seeks to identify matrices $A$ and $B$ such that 
# \begin{equation}
#     \dfrac{dV_{[1:r-1]}}{dt} = AV_{[1:r-1]} + BV_{r}
# \end{equation}
#
# To perform the fit, we need data for $\dot{V}$ which we obtain via fourth order central difference, discarding values at the edge to keep things simple. Alternatively, we could construct an interpolant from [Interpolations.jl](http://juliamath.github.io/Interpolations.jl/latest/convenience-construction/) and obtain derivative estimates by calling the inbuilt function 
# ```julia
# g = Interpolations.gradient(itp, x, y, ...)
# ```

# +
dt = pm2_5.Δt
size(Vr)

dVr = zeros(size(Vr,1)-5, r-1)
for k ∈ 1:r-1, i ∈ 3:size(Vr,1)-3
    dVr[i-2,k] = (1/(12*dt)) * (-Vr[i+2, k] + 8*Vr[i+1, k] - 8*Vr[i-1,k] + Vr[i-2,k])
end

@assert size(dVr,2) == r-1

# +
# chop off edges to size of data matches size of derivative
X = @view Vr[3:end-3, :]
dX = @view dVr[:,:]

@assert size(X,2) == size(dX,2)  + 1
# -

# thinking of data as rows of X, we now solve the regression problem given by $dX = X\Xi$ (where $\Xi$ is composed of our A and B matrices) to obtain $\Xi$ via ordinary least squares (i.e. via normal equations)

# +
Ξ = (X\dX)'  # now Ξx = dx for a single column vector view

A = Ξ[:, 1:end-1]   # State matrix A
B = Ξ[:, end]'      # Control matrix B
# -

# It is interesting to observe the structure of the linear component of our model, `A`: 

p1 = heatmap(A, yflip=true, xlabel="A", ylabel="", showaxis=false,link=:y, cbar=false,clims=(minimum(A), maximum(A)), leftmargin=0*Plots.mm,rightmargin=0*Plots.mm)
p2 = heatmap(reshape(B', length(B), 1), yflip=true, xlabel="B", ylabel="", showaxis=false, link=:y, clims=(minimum(A),maximum(A)), color=:inferno,rightmargin=10*Plots.mm, leftmargin=0*Plots.mm)
plot(p1, p2, layout = @layout([a{0.8w} b{0.2w} ]))
#plot(p1, p2, layout=grid(1,2, widths=[0.8, 0.2]))

# We might also be interested in visualizing the "eigenmodes" discovered by the SVD:

p = plot(legend=:outerright)
for i ∈ 1:r
    plot!(Ur[:,i], label=L"u_{%$i}")
end
p

# Interestingly, these look like a sort of polynomial basis for $\mathbb{R}^n$.

using DifferentialEquations
using Interpolations
using StaticArrays

size(X)

# X data starts at index 3 in V, i.e. at 2Δt
ts = range(0, step=pm2_5.Δt, length=size(X,1))
tspan = (ts[1], ts[end])

# let's construct an interpolation function to return the value of our forcing function given by $X_r$ 

# +
xᵣ = linear_interpolation(ts, X[:,end], extrapolation_bc=Line())

xᵣ(ts[1])

# +
xs = xᵣ.(ts)

p1 = plot(
    ts ./3600, 
    X[:,1],
    xlabel="",
    ylabel="v₁",
    label=""
    ) 

p2 = plot(
    ts ./3600,
    xs.^2,
    ylabel="vᵣ²",
    xlabel="time (hr)",
    label="",
    color=:red,
    link=:x,
    grid=false,
    minorgrid=false,
    yticks=[0.0, 0.0025]
)

l = @layout [ 
    a{0.8h}
    b
]
plot(p1, p2, layout=l)
# -

# With our function in hand, we can now simulate the dynamics using our HAVOK model: 
# \begin{equation}
#     \dot{X} = AX + Bx_r(t)
# \end{equation}
# This will let us evaluate how well we have reconstructed the time series attractor: 

# Generate our A and B matrices as static arrays
sA = @SMatrix[A[i,j] for i ∈ axes(A,1), j ∈ axes(A,2)]
sB = @SVector[b for b ∈ B'] # j ∈ axes(B,2)]

# +
function f!(du, u, p, t)
    A,B = p
    du .= A*u + B*xᵣ(t)
end

params = (sA, sB)
x₀ = X[1,1:end-1]
dx = copy(x₀)

@assert size(x₀) == size(dx)
@benchmark f!(dx, x₀, params, ts[1])
# -

prob = ODEProblem(f!, x₀, tspan, params)
sol = solve(prob, saveat=ts);
size(sol)

# +
p1 = plot(
    ts ./3600, 
    X[:,1],
    xlabel="",
    ylabel="coordinate",
    label="v₁",
) 

scatter!(
    sol.t[1:20:end]./3600,
    sol[1,1:20:end],
    label="v₁ prediction"
)
# -












