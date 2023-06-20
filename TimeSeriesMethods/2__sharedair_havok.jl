using CSV, DataFrames
using Plots, StatsPlots
using TimeSeriesTools
using ParameterHandling
using Dates, TimeZones
using Unitful
using Markdown, LaTeXStrings
using Statistics, StatsBase, Distributions, KernelDensity
using BenchmarkTools
using LinearAlgebra
using DifferentialEquations
using DataInterpolations
using StaticArrays

include("plot_defaults.jl")
add_mints_theme()
theme(:mints)

include("plot_recipes.jl")
include("utils.jl")

if !ispath("figures/sharedair/havok")
    mkpath("figures/sharedair/havok")
end


# 0. Load in data and set window size (for ncol argument)
println("loading in the data...")

function get_data(path, cols_to_use)
    df = CSV.File("data/sharedair/data.csv") |> DataFrame
    ts = df.t
    return ts, Matrix(df[:, cols_to_use])
end

#ts, Data = get_data("data/sharedair/data.csv", [:pm1_0, :pm2_5, :pm10_0])
ts, Data = get_data("data/sharedair/data.csv", [:pm2_5]) #, :pm2_5, :pm10_0])
size(ts)
size(Data)

Nwindow = 81 # determined in step 1 by range of variogram fit

# 1. Generate Hankel matrix for each individual time series and concatenate together
println("computing time delay embedding Hankel-Takens matrix...")
H = vcat([TimeDelayEmbedding(Data[:,i]; nrow=Nwindow, method=:backward) for i∈axes(Data,2)]...)

size(H)  # so we have 10 columns of data with 81 times each giving 810 elements per embedding vector


# 2. Decompose via SVD
println("computing SVD... this could take a while")
U, σ, V = svd(H)

println("U has dimensions ", size(U))
println("V has dimensions ", size(V))
println("we have $(length(σ)) singular values")

@assert all(H .≈ U*Diagonal(σ)*V')  # verify that decomposition works

# 3. visualize the attractor:
#Nmax = 200000
Nmax = 400000
plot_times =  collect(1:size(V,1)) ./ (60.0)^2
p1 = scatter(
    V[1:Nmax,1], V[1:Nmax,2], V[1:Nmax,3],
    ms=2,
    msw=0,
    msa=0,
    marker_z = plot_times[1:Nmax],
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

savefig("figures/sharedair/havok/attractor.png")
savefig("figures/sharedair/havok/attractor.pdf")

# 4. visualize singular values
plot(σ./sum(σ), xlabel=L"index, $i$", ylabel=L"\sigma_{i}/(\Sigma_{i}\sigma_{i})"*"\n", lw=3, label="")
savefig("figures/sharedair/havok/singular_values.png")
savefig("figures/sharedair/havok/singular_values.pdf")


# 5. compute cut-off value for singular values based on magnitude
r = r_cutoff(σ,ratio=0.005, rmax=30)

size(V)

# r = r_optimal_approx(σ, size(V,2), size(V,1))
# r = 46  # too high
# r = 30 # too high
# r = 20 # slightly too big
# r = 15 # a smidge too small
r = 18  # really close but slightly too big
# r = 17 # too small
# r = 19


# add extra dimensions to r for > 1 control variable
# n_control = 3
# n_control = 2  # not enough
# n_control = 5  
n_control = 10  # not enough

r = r + n_control - 1


# 6. truncate the matrices
Vr = @view V[:,1:r]
Ur = @view U[:,1:r]
σr = @view σ[1:r]


# 7. compute derivative using fourth order central difference scheme
dt = mean(ts[2:end] .- ts[1:end-1])

dVr = zeros(size(Vr,1)-5, r-n_control)

Threads.@threads for k ∈ 1:r-n_control
    for i ∈ 3:size(Vr,1)-3
        @inbounds dVr[i-2,k] = (1/(12*dt)) * (-Vr[i+2, k] + 8*Vr[i+1, k] - 8*Vr[i-1,k] + Vr[i-2,k])
    end
end

@assert size(dVr,2) == r-n_control


# 8. chop off edges to size of data matches size of derivative
X = @view Vr[3:end-3, :]
dX = @view dVr[:,:]

# Xtrain = @view X[1:end-1000,:]
# Xtest = @view X[end-1000+1:end,:]

# dXtrain = @view dX[1:end-1000, :]
# dXtest = @view dX[end-1000+1:end, :]


Ntrain = size(X,1) - 1000

@assert size(dX,2) == size(X,2)  - n_control


# 9. Compute model matrix via least squares
Ξ = (X\dX)'  # now Ξx = dx for a single column vector view

A = Ξ[:, 1:r-n_control]   # State matrix A
B = Ξ[:, r-n_control+1:end]      # Control matrix B



# 10. visualize matrices
p1 = heatmap(A, yflip=true, xlabel="A", ylabel="", showaxis=false,link=:y, cbar=false,clims=(minimum(A), maximum(A)), leftmargin=0*Plots.mm,rightmargin=0*Plots.mm)
p2 = heatmap(B, yflip=true, xlabel="B", ylabel="", showaxis=false, link=:y, clims=(minimum(A),maximum(A)), color=:inferno,rightmargin=10*Plots.mm, leftmargin=0*Plots.mm)
plot(p1, p2, layout = @layout([a{0.7w} b{0.3w} ]))

savefig("figures/sharedair/havok/heatmap.png")
savefig("figures/sharedair/havok/heatmap.pdf")


# 11. visualize eigenmodes
p = plot([], yticks=[0.0,], legend=:outerright, label="")
for i ∈ 1:r
    if i ≤ 3
        plot!(Ur[:,i], label=L"u_{%$i}", color=:blue)
    elseif i > 3 && i ≤ r-n_control
        plot!(Ur[:,i], label="", color=:grey, alpha=0.5)
    elseif i == r - n_control + 1
        plot!(Ur[:,i], label=L"u_{r}", color=:red)
    else
        plot!(Ur[:,i], label="", color=:red)
    end
end

for i ∈ 1:Int(size(Ur, 1) / Nwindow) - 1
    vline!([Nwindow*i], color = :black, alpha=1, label="", lw=2)
end

display(p)

savefig("figures/sharedair/havok/eigenmodes.png")
savefig("figures/sharedair/havok/eigenmodes.pdf")



# 12. define interpolation function for forcing coordinate(s)
# +3 because of offset from derivative...
ts = range(dt*(Nwindow+3), step=dt, length=size(X,1))

#itps = [CubicSpline(X[:,j], ts) for j ∈ r-n_control+1:r]
itps = [DataInterpolations.LinearInterpolation(X[:,j], ts) for j ∈ r-n_control+1:r]
u(t) = [itp(t) for itp ∈ itps]

#xᵣ = CubicSpline(X[:,end], ts)
#xᵣ = linear_interpolation(Float64.(ts), X[:, r-n_control+1:end])

# xᵣ(ts[1])


# 13. visualize first embedding coordinate + the forcing term
xs = u.(ts)
p1 = plot(
    ts[1:Nmax] ./ (60^2),
    X[1:Nmax,1],
    xlabel="",
    ylabel="v₁",
    label=""
)

p2  = plot()
for j ∈ 1:1
    plot!(
        ts[1:Nmax] ./ (60^2),
        map(x->x[j]^2, xs[1:Nmax]),
        ylabel="vᵣ²",
        xlabel="time",
        label="",
        color=:red,
        link=:x,
        grid=false,
        minorgrid=false,
        yticks=[0.0]
    )
end

l = @layout [
    a{0.8h}
    b
]
p = plot(p1, p2, layout=l)
display(p)

savefig("figures/sharedair/havok/v1_with_forcing.png")
savefig("figures/sharedair/havok/v1_with_forcing.pdf")





# 14. Integrate model forward in time
sA = @SMatrix[A[i,j] for i ∈ axes(A,1), j ∈ axes(A,2)]
sB = @SMatrix[B[i,j] for i ∈ axes(B,1), j ∈ axes(B,2)]
#sB = @SVector[b for b ∈ B'] # j ∈ axes(B,2)]


# define function and integrate to get model predictions
function f!(dx, x, p, t)
    A,B = p
    dx .= A*x + B*u(t)
end

params = (sA, sB)
x₀ = X[1,1:r-n_control]
dx = copy(x₀)

# A*x₀ + sB*u(ts[1])

@assert size(x₀) == size(dx)
@benchmark f!(dx, x₀, params, ts[1])



prob = ODEProblem(f!, x₀, (ts[1], ts[end]), params)
sol = solve(prob, saveat=ts)# , abstol=1e-12, reltol=1e-12);
size(sol)


# 15. visualize results
offset = 1
#offset = Ntrain - 10000
Nmax = 10000
p1 = plot(
    ts[offset:offset+Nmax] ./ (60^2),
    X[offset:offset+Nmax, 1],
    # xlabel="time (hours)",
    ylabel="v₁",
    label="embedding",
    lw=2,
    legend=:topright
)


plot!(
    sol.t[offset:offset+Nmax] ./ (60^2),
    sol[1,offset:offset+Nmax],
    label="fit",
    ls=:dot,
    lw=2,
)


p2 = plot(
    ts[offset:offset+Nmax] ./ (60^2),
    X[offset:offset+Nmax, 2],
    # xlabel="time (hours)",
    ylabel="v₂",
    label="embedding",
    lw=2,
    legend=:topright,
    link=:x,
)


plot!(
    sol.t[offset:offset+Nmax] ./ (60^2),
    sol[2,offset:offset+Nmax],
    label="fit",
    ls=:dot,
    lw=2
)


p3 = plot(
    ts[offset:offset+Nmax] ./ (60^2),
    X[offset:offset+Nmax, 3],
    xlabel="time (hours)",
    ylabel="v₃",
    label="embedding",
    lw=2,
    legend=:topright,
    link=:x
)


plot!(
    sol.t[offset:offset+Nmax] ./ (60^2),
    sol[3,offset:offset+Nmax],
    label="fit",
    ls=:dot,
    lw=2
)

plot(p1, p2, p3, layout=(3,1), size=(1600, 1000), margin=5*Plots.mm)


savefig("figures/sharedair/havok/timeseries_reconstructed.png")
savefig("figures/sharedair/havok/timeseries_reconstructed.pdf")



# 16. visualize the fitted attractor:
Nmax = 15000
p1 = scatter(
    sol[1,1:Nmax], sol[2, 1:Nmax], sol[3, 1:Nmax],
    ms=2,
    msw=0,
    msa=0,
    marker_z = sol.t[1:Nmax],
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

savefig("figures/sharedair/havok/attractor_reconstructed.png")
savefig("figures/sharedair/havok/attractor_reconstructed.pdf")


train_shift = 30000
# 17. scatter plot and quantile quantile of fit
p1 = scatterresult(
    X[1:Nmax,1], sol[1,1:Nmax],
    X[train_shift+1:train_shift+Nmax,1], sol[1, train_shift+1:train_shift+Nmax],
    xlabel="True v₁",
    ylabel="Predicted v₁",
    plot_title="HAVOK Fit for v₁",
)

savefig("figures/sharedair/havok/scatterplot.png")
savefig("figures/sharedair/havok/scatterplot.pdf")


p1 = quantilequantile(
    X[1:Nmax,1], sol[1,1:Nmax],
    X[train_shift+1:train_shift+Nmax,1], sol[1, train_shift+1:train_shift+Nmax],
    xlabel="True v₁",
    ylabel="Predicted v₁",
    title="HAVOK Fit for v₁",
)

savefig("figures/sharedair/havok/quantile-quantile.png")
savefig("figures/sharedair/havok/quantile-quantile.pdf")


pdf = kde(X[:, r-n_control + 1])
idxs_nozero = pdf.density .> 0
gauss = fit(Normal, X[:, r-n_control+1])

plot(gauss, label="gaussian fit", yaxis=:log, ls=:dash)
plot!(pdf.x[idxs_nozero], pdf.density[idxs_nozero], label="pdf")
ylims!(1e-1, 1e3)
xlims!(-0.01, 0.01)
title!("Forcing Statistics")

savefig("figures/sharedair/havok/forcing-stats.png")
savefig("figures/sharedair/havok/forcing-stats.pdf")


