using CSV, DataFrames
using Plots
using TimeSeriesTools
using ParameterHandling
using Dates, TimeZones
using Unitful
using Markdown, LaTeXStrings
using Statistics, StatsBase
using BenchmarkTools
using LinearAlgebra
using DifferentialEquations
using DataInterpolations
using StaticArrays

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

ts, Data = get_data("data/sharedair/data.csv", [:pm1_0, :pm2_5, :pm10_0])
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
Nmax = 200000
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


# 6. truncate the matrices
Vr = @view V[:,1:r]
Ur = @view U[:,1:r]
σr = @view σ[1:r]


# 7. compute derivative using fourth order central difference scheme
dt = mean(ts[2:end] .- ts[1:end-1])

dVr = zeros(size(Vr,1)-5, r-1)
for k ∈ 1:r-1, i ∈ 3:size(Vr,1)-3
    dVr[i-2,k] = (1/(12*dt)) * (-Vr[i+2, k] + 8*Vr[i+1, k] - 8*Vr[i-1,k] + Vr[i-2,k])
end

@assert size(dVr,2) == r-1


# 8. chop off edges to size of data matches size of derivative
X = @view Vr[3:end-3, :]
dX = @view dVr[:,:]

@assert size(X,2) == size(dX,2)  + 1


# 9. Compute model matrix via least squares
Ξ = (X\dX)'  # now Ξx = dx for a single column vector view

A = Ξ[:, 1:end-1]   # State matrix A
B = Ξ[:, end]'      # Control matrix B


# 10. visualize matrices
p1 = heatmap(A, yflip=true, xlabel="A", ylabel="", showaxis=false,link=:y, cbar=false,clims=(minimum(A), maximum(A)), leftmargin=0*Plots.mm,rightmargin=0*Plots.mm)
p2 = heatmap(reshape(B', length(B), 1), yflip=true, xlabel="B", ylabel="", showaxis=false, link=:y, clims=(minimum(A),maximum(A)), color=:inferno,rightmargin=10*Plots.mm, leftmargin=0*Plots.mm)
plot(p1, p2, layout = @layout([a{0.8w} b{0.2w} ]))

savefig("figures/sharedair/havok/heatmap.png")
savefig("figures/sharedair/havok/heatmap.pdf")


# 11. visualize eigenmodes
p = plot([], yticks=[0.0,], legend=:outerright, label="")
for i ∈ 1:r
    if i ≤ 3
        plot!(Ur[:,i], label=L"u_{%$i}", color=:blue)
    elseif i > 3 && i < r
        plot!(Ur[:,i], label="", color=:grey, alpha=0.5)
    else
        plot!(Ur[:,i], label=L"u_{r}", color=:red)
    end
end

for i ∈ 1:Int(size(Ur, 1) / Nwindow) - 1
    vline!([Nwindow*i], color = :black, alpha=1, label="", lw=2)
end

display(p)

savefig("figures/sharedair/havok/eigenmodes.png")
savefig("figures/sharedair/havok/eigenmodes.pdf")



# 12. define interpolation function for forcing coordinate
# +3 because of offset from derivative...
ts = range(dt*(Nwindow+3), step=dt, length=size(X,1))
xᵣ = CubicSpline(X[:,end], ts)
xᵣ(ts[1])


# 13. visualize first embedding coordinate + the forcing term
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
    yticks=[0.0]
)

l = @layout [
    a{0.8h}
    b
]
plot(p1, p2, layout=l)

savefig("figures/sharedair/havok/v1_with_forcing.png")
savefig("figures/sharedair/havok/v1_with_forcing.pdf")





# 14. Integrate model forward in time
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


# 15. visualize results
p1 = plot(
    ts[1:Nmax] ./ (60^2),
    X[1:Nmax, 1],
    xlabel="time",
    ylabel="v₁",
    label="embedding",
    lw=2
)


plot!(
    sol.t[1:Nmax] ./ (60^2),
    sol[1,1:Nmax],
    label="fit",
    ls=:dot,
    lw=2
)

savefig("figures/lorenz/timeseries_reconstructed.png")
savefig("figures/lorenz/timeseries_reconstructed.pdf")

