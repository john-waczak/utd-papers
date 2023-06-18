using CSV, DataFrames
using Plots
using TimeSeriesTools
using ParameterHandling
using Dates, TimeZones
using Unitful
using Markdown, LaTeXStrings
using Statistics, StatsBase
using BenchmarkTools
using LinearAlgebra, StaticArrays
using DifferentialEquations
using DataInterpolations

include("utils.jl")

if !ispath("figures/lorenz")
    mkpath("figures/lorenz")
end


# load data
Data = Matrix(CSV.File("data/lorenz/data.csv") |> DataFrame)


# compute time delay embedding
H = TimeDelayEmbedding(Data[:,2], method=:backward)

# compute singular value decomposition
U, σ, V = svd(H)

size(V)

# visualize attractor
dt = 0.001
tspan = range(dt, step=dt, length=size(Data,1))
Nmax = 50000  # max value for plotting

p1 = scatter(
    Data[1:Nmax,2], Data[1:Nmax,3], Data[1:Nmax,4],
    ms=2,
    msw=0,
    msa=0,
    marker_z = tspan[1:Nmax],
    frame=:semi,
    ticks=nothing,
    xlabel="x",
    ylabel="y",
    zlabel="z",
    label="",
    cbar=false,
    margins=0*Plots.mm,
    background_color=:transparent,
    title="attractor"
)


p2 = scatter(
    V[1:Nmax,1], V[1:Nmax,2], V[1:Nmax,3],
    ms=2,
    msw=0,
    msa=0,
    marker_z = tspan[1:Nmax],
    frame=:semi,
    ticks=nothing,
    xlabel="v₁",
    ylabel="v₂",
    zlabel="v₃",
    label="",
    cbar=false,
    margins=0*Plots.mm,
    background_color=:transparent,
    title="embedded attractor"
)

plot(p1, p2)

savefig("figures/lorenz/attractors1.png")
savefig("figures/lorenz/attractors1.pdf")


# set r value to 15 as in paper
r = 15
# truncate the matrices
Vr = @view V[:,1:r]
Ur = @view U[:,1:r]
σr = @view σ[1:r]


# compute derivatives with fourth order finite difference scheme
dVr = zeros(size(Vr,1)-5, r-1)
for k ∈ 1:r-1, i ∈ 3:size(Vr,1)-3
    dVr[i-2,k] = (1/(12*dt)) * (-Vr[i+2, k] + 8*Vr[i+1, k] - 8*Vr[i-1,k] + Vr[i-2,k])
end

@assert size(dVr,2) == r-1


# chop off edges so size of data matches size of derivative
X = @view Vr[3:end-3, :]
dX = @view dVr[:,:]

@assert size(X,2) == size(dX,2)  + 1

# compute matrix such that dX = XΞ'
Ξ = (X\dX)'  # now Ξx = dx for a single column vector view

A = Ξ[:, 1:end-1]   # State matrix A
B = Ξ[:, end]'      # Control matrix B


# visualize matrices
p1 = heatmap(A, yflip=true, xlabel="A", ylabel="", showaxis=false,link=:y, cbar=false,clims=(minimum(A), maximum(A)), leftmargin=0*Plots.mm,rightmargin=0*Plots.mm)
p2 = heatmap(reshape(B', length(B), 1), yflip=true, xlabel="B", ylabel="", showaxis=false, link=:y, clims=(minimum(A),maximum(A)), color=:inferno,rightmargin=10*Plots.mm, leftmargin=0*Plots.mm)
plot(p1, p2, layout = @layout([a{0.8w} b{0.2w} ]))

savefig("figures/lorenz/heatmap.png")
savefig("figures/lorenz/heatmap.pdf")

# visualize eigenmodes
p = plot([], yticks=[-0.3, 0.0, 0.3], legend=:outerright, label="")
for i ∈ 1:r
    if i ≤ 3
        plot!(Ur[:,i], label=L"u_{%$i}", color=:blue)
    elseif i > 3 && i < r
        plot!(Ur[:,i], label=L"u_{%$i}", color=:grey, alpha=0.5)
    else
        plot!(Ur[:,i], label=L"u_{%$i}", color=:red)
    end
end
display(p)

savefig("figures/lorenz/eigenmodes.png")
savefig("figures/lorenz/eigenmodes.pdf")


# define interpolation function for forcing coordinate
ts = range(dt, step=dt, length=size(X,1))
xᵣ = CubicSpline(X[:,end], ts)
xᵣ(ts[1])

# visualize first embedding coordinate that we want to fit:
xs = xᵣ.(ts)

p1 = plot(
    ts[1:Nmax],
    X[1:Nmax,1],
    xlabel="",
    ylabel="v₁",
    label=""
)

p2 = plot(
    ts[1:Nmax],
    xs[1:Nmax].^2,
    ylabel="vᵣ²",
    xlabel="time",
    label="",
    color=:red,
    link=:x,
    grid=false,
    minorgrid=false,
    yticks=[0.0, 0.0005]
)

l = @layout [
    a{0.8h}
    b
]
plot(p1, p2, layout=l)

savefig("figures/lorenz/v1_with_forcing.png")
savefig("figures/lorenz/v1_with_forcing.pdf")


# Generate our A and B matrices as static arrays
sA = @SMatrix[A[i,j] for i ∈ axes(A,1), j ∈ axes(A,2)]
sB = @SVector[b for b ∈ B'] # j ∈ axes(B,2)]


# define function and integrate to get model predictions
function f!(du, u, p, t)
    A,B = p
    du .= A*u + B*xᵣ(t)
end

params = (sA, sB)
x₀ = X[1,1:end-1]
dx = copy(x₀)

@assert size(x₀) == size(dx)
@benchmark f!(dx, x₀, params, ts[1])

prob = ODEProblem(f!, x₀, (ts[1], ts[end]), params)
sol = solve(prob, saveat=ts);
size(sol)

# visualize results
L = 300:25000

p1 = plot(
    tspan[L],
    X[L],
    xlabel="time",
    ylabel="v₁",
    label="embedding",
    lw=2
)


plot!(
    sol.t[L],
    sol[1,L],
    label="fit",
    ls=:dash,
    lw=2
)

savefig("figures/lorenz/timeseries_reconstructed.png")
savefig("figures/lorenz/timeseries_reconstructed.pdf")
