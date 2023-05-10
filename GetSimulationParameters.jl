# Create the SimulationParameters instance

params = SimulationParameters(N, n, Nbody, N * n + Nbody, Nfil, S, a, b,swimming_strength)

# Create the PreallocatedMatrices instance
matrices = PreallocatedMatrices(
    zeros(6 * params.Ntot, 6 * params.Ntot),
    zeros(3 + 3 * params.N, 6 * params.Ntot),
    zeros(6 * (params.Nbody + (params.Nfil * params.N * params.n)), 3 + 3 * params.N * params.Nfil),
    zeros(3 + 3 * params.N, 3 + 3 * params.N),
    zeros(1 + params.N, 3),
    zeros(1 + params.N, 3),
    zeros(1 + params.N, 3),
    zeros(3 * params.Ntot),
    zeros(3+3*N)
)

# Create the CombinedParameters instance
combined_params = CombinedParameters(params, matrices)

