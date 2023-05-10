#Coarse-grained simulations of elastic filaments, implements method outlined in https://www.biorxiv.org/content/10.1101/2023.01.17.524366v1
#Currently no event condition implemnted. If any of ||X[4:6]||, ..., ||X[end-2:end]|| approach pi, system becomes very stiff and X[4:6]... need rescaling
#I run code with N = 3 to compile before setting N = 20

using LinearAlgebra, StaticArrays, DifferentialEquations, Sundials, Plots
using BenchmarkTools, ODEInterfaceDiffEq, DiffEqDevTools, OrdinaryDiffEq
using Logging: global_logger
using TerminalLoggers: TerminalLogger
using Profile, ProfileSVG

global_logger(TerminalLogger())
gr() # gr(fmt=:png)

#Run required_functions.jl first
N = 20#segments per filament, must be > 2
n = 1 #spheres per segment
Nbody = 0 #must be 0, not fully implemented yet
Ntot = N * n + Nbody
Nfil = 1 #number of filaments (must be 0 currently)
S = 4 #filament stiffness parameter 
a = 1 / (2 * N * n) * ones(Ntot) #sphere radii
b = [0.0 0.0 0.0] #body-filament separation in body frame
swimming_strength = 10.0

include("GetSimulationParameters.jl")

X = zeros(6 + combined_params.params.N * 3);

X[4:3:end] .= 0.0
X[5:3:end] .= 0.0
X[6:3:end] .= 0.0
# X[8:3:end] .= LinRange(0, pi/2, params.N) #set swimming strength to 0 and uncomment this line to see a relaxing filament

prob = ODEProblem(dynamics_combined, X, (0.0, 10.0), combined_params)
sol = solve(prob, CVODE_BDF(),abstol=1e-6, reltol=1e-3, progress = true, progress_steps=1)

# Plot solution over time:
d1 = matrices.d1
d2 = matrices.d2
d3 = matrices.d3
X3 = matrices.X3
@gif for t = 0:0.2:10.0
X = sol(t)
# calc_directors!(d1, d2, d3, X, params) #director basis along filament
calc_directors_n!(d1, d2, d3, X, params)
calc_X3!(X3, X, params.N, params.Nbody, params.b, params.Nfil, params.n, d1, d2, d3) #sphere positions in lab frame
plot(X3[1:3:end], X3[2:3:end],X3[3:3:end], xlims = (-1., 1.), ylims = (-1., 1.),zlims = (-1., 1.), aspect_ratio = 1, legend = false)
end

# Tests:

# prob = ODEProblem(dynamics_combined, X, (0.0, 10.0), combined_params)
# sol = solve(prob,FBDF(autodiff=false),abstol=1/10^14,reltol=1/10^14,progress = true,progress_steps=1)

# test_sol = TestSolution(sol)
# abstols = 1.0 ./ 10.0 .^ (4:13)
# reltols = 1.0 ./ 10.0 .^ (1:10)
# setups = [
#           Dict(:alg=>CVODE_BDF()),
#           Dict(:alg=>FBDF(autodiff=false)),
#           Dict(:alg=>QNDF(autodiff=false)),
#           Dict(:alg=>RadauIIA5(autodiff=false)),
#           ]

# wp = WorkPrecisionSet(prob,abstols,reltols,setups;verbose=false,
#                       save_everystep=false,appxsol=test_sol,maxiters=Int(1e5),numruns=2)
# h = plot(wp)
# savefig(h, "WorkPrecisionSet_2.png")