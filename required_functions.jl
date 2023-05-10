function dynamics_combined(X, combined_params, t)
    # Access the preallocated matrices and parameters
    matrices = combined_params.matrices
    params = combined_params.params

    Mh = matrices.Mh
    Mf = matrices.Mf
    Q = matrices.Q
    Qomega = matrices.Qomega
    d1 = matrices.d1
    d2 = matrices.d2
    d3 = matrices.d3
    X3 = matrices.X3
    K = matrices.K

    # calc_directors!(d1, d2, d3, X, params) #director basis along filament
    # calc_RHS!(K, t, X, params, d1, d2, d3) #right-hand-side bending moments and any internal forcing
    # calc_Qomega!(Qomega, X, params, d1, d2, d3) #quaternion exponential map -> angular velocity

    calc_Qomega_n!(Qomega, X, params)
    calc_directors_n!(d1, d2, d3, X, params)
    calc_RHS_n!(K, t, X, params, d1, d2, d3)

    calc_X3!(X3, X, params.N, params.Nbody, params.b, params.Nfil, params.n, d1, d2, d3) #sphere positions in lab frame
    calc_Mf!(Mf, X, X3, params.N, params.n, params.Ntot, d3, params.Nbody) #encodes force/torque balance
    calc_Mh!(Mh, X3, params.Ntot, params.a) #calculates RPY matrix (need to implement overlapping sphere case)
    calc_Q!(Q, X, params, d1, d2, d3, params.Nbody, params.Ntot, params.Nfil, params.N, params.n, params.b) #dimensionality reduction matrix         

    A = params.S^4 * Mf
    B = Mh
    C = Q * Qomega

    dX = (A * (B \ C)) \ K

    #force body frame and filament segment 1 frame to evolve together, implement novel conditions later
    insert!(dX, 4, dX[4])
    insert!(dX, 5, dX[6])
    insert!(dX, 6, dX[8])
    return dX
end

using StaticArrays

function skew!(s, a)
    s[1, 2] = -a[3]
    s[1, 3] = a[2]
    s[2, 1] = a[3]
    s[2, 3] = -a[1]
    s[3, 1] = -a[2]
    s[3, 2] = a[1]
end

skewS(a) = @SMatrix [0 -a[3] a[2]; a[3] 0 -a[1]; -a[2] a[1] 0]

function calc_Mh!(Mh, X3, Ntot, a)
    I3 = Diagonal(SVector(true, true, true))
    for i = 1:Ntot
        xi = @views SVector{3}(X3[3(i-1)+1:3*i])
        for j = i+1:Ntot
            xj = @views SVector{3}(X3[3*(j-1)+1:3*j])

            rij = xi - xj
            r = norm(rij)
            rhat = rij / r
            rhat_skew = skewS(-rhat)

            sqrd = (a[i]^2 + a[j]^2) / r^2
            rhat_op = rhat * rhat'

            Mtt = 1 / (2 * r) * ((1 + sqrd / 3) * I3 + (1 - sqrd) * rhat_op)
            Mtr = 1 / (2 * r^2) * rhat_skew
            Mrt = Mtr
            Mrr = 1 / (4 * r^3) * (3 * rhat_op - I3)

            Mh[3*(i-1)+1:3*i, 3*(j-1)+1:3*j] = Mtt
            Mh[3*(i-1)+1:3*i, 3*(Ntot)+3*(j-1)+1:3*(Ntot)+3*j] = Mtr
            Mh[3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i, 3*(j-1)+1:3*j] = Mrt
            Mh[3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i, 3*(Ntot)+3*(j-1)+1:3*(Ntot)+3*j] = Mrr

            Mh[3*(j-1)+1:3*j, 3*(i-1)+1:3*i] = Mtt'
            Mh[3*(j-1)+1:3*j, 3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i] = Mtr'
            Mh[3*(Ntot)+3*(j-1)+1:3*(Ntot)+3*j, 3*(i-1)+1:3*i] = Mrt'
            Mh[3*(Ntot)+3*(j-1)+1:3*(Ntot)+3*j, 3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i] = Mrr'
        end

        Mtt = 2 / (3 * a[i]) * I3
        Mrr = 1 / (2 * a[i]^3) * I3
        Mh[3*(i-1)+1:3*i, 3*(i-1)+1:3*i] = Mtt
        Mh[3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i, 3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i] = Mrr
    end
    Mh .= Mh .* (1 / (4 * pi))
end

function calc_Mf!(Mf, X, X3, N, n, Ntot, d3, Nbody)
    I3 = Diagonal(SVector(true, true, true))
    ds = 1 / (2 * N * n)
    # display(N*n)
    for i = 1:N*n
        Mf[1:3, 3*(i-1)+1:3*i] = I3
    end

    x0 = @views SVector{3}(X[1:3])

    for i = 1:Ntot
        xi = @views SVector{3}(X3[3*(i-1)+1:3*i])
        dx = xi - 0 * x0
        # skew!(Bji,dx)

        Bji = skewS(dx)

        Mf[3+3*(1-1)+1:3+3*1, 3*(i-1)+1:3*i] = Bji
        Mf[3+3*(1-1)+1:3+3*1, 3*Ntot+3*(i-1)+1:3*Ntot+3*i] = I3
    end

    for j = 2:N
        sn = Nbody + (j - 1) * n + 1
        t1 = @views SVector{3}(X3[3*(sn-1)+1:3*sn])
        t2 = @views SVector{3}(d3[j+1, :])
        # xj = SVector{3}(X3[3*(sn-1)+1:3*sn] .- ds .* d3[j+1, :])
        xj = SVector{3}(t1 - ds * t2)
        for i = Nbody+(j-1)*n+1:Nbody+N*n
            xi = @views SVector{3}(X3[3*(i-1)+1:3*i])
            dx = xi - xj

            Bji = skewS(dx)
            Mf[3*j+1:3*(j+1), 3*(i-1)+1:3*i] = Bji
            Mf[3*j+1:3*(j+1), 3*Ntot+3*(i-1)+1:3*Ntot+3*i] = I3
        end
    end
end




function calc_Q!(Q, X, params, d1, d2, d3, Nbody, Ntot, Nfil, N, n, b)
    Δ₁ = zeros(3, 3)
    skew!(Δ₁, b[end, :])
    a::Float64 = 1 / (2 * N * n)

    rcross = zeros(3, 3)
    iden = Diagonal(SVector(true, true, true))
    x0 = X[1:3]
    for i = 1:Nbody
        r = reshape(x0, (1, 3)) + reshape(b[i, :], (1, 3)) * [d1[1, :] d2[1, :] d3[1, :]]
        skew!(rcross, r)
        Q[3*(i-1)+1:3*i, 1:3] = iden
        Q[3*(i-1)+1:3*i, 4:6] = -rcross
        Q[3*Ntot+3*(i-1)+1:3*Ntot+3*i, 4:6] = iden
    end

    Q[3*Nbody+3*(1-1)+1:3*Nbody+3*N*n, 1:3] = repeat(iden, N * n)
    Q[3*Nbody+3*(1-1)+1:3*Nbody+3*N*n, 4:6] = repeat(Δ₁, N * n)

    #Segment 1
    i = 1
    d3i = @views SVector{3}(d3[i+1, :])
    s1 = Nbody + (i - 1) * n + 1
    sn = Nbody + (i - 1) * n + n
    dist_mat = zeros(3 * n, 3)
    d3skew = skewS(d3i)
    for i = 1:n
        dist_mat[3*(i-1)+1:3*i, 1:3] = -(a + (i - 1) * a * 2) * d3skew
    end

    # for i = 1:n
    Q[3*(s1-1)+1:3*(sn), 3+3*(i-1)+1:3+3*i] .= dist_mat#cat([dist[i] * d3skew for i in 1:n]...,dims=1)

    #Subsequent segments
    for i = 2:N

        d3i1 = @views SVector{3}(d3[i+1, :])
        d3i = @views SVector{3}(d3[i, :])
        # skew!(d3skew,d3[i+1,:])
        d3skew = skewS(d3i1)
        s1 = Nbody + (i - 1) * n + 1
        sn = Nbody + (i - 1) * n + n

        #Get all terms from previous segment
        s_prev = Nbody + (i - 1 - 1) * n

        Q[3*(s1-1)+1:3*(sn), :] .= Q[3*(s_prev)+1:3*(s_prev+n), :]

        for i = 1:n
            dist_mat[3*(i-1)+1:3*i, 1:3] .= -(a + (i - 1) * a * 2) .* d3skew
        end

        #Contribution from current segment
        Q[3*(s1-1)+1:3*(sn), 3+3*(i-1)+1:3+3*i] .= dist_mat#cat([dist[i] * d3skew for i in 1:n]...,dims=1)

        d3skew = skewS(d3i)
        for i = 1:n
            dist_mat[3*(i-1)+1:3*i, 1:3] .= -2.0 .* n .* a .* d3skew
        end
        Q[3*(s1-1)+1:3*(sn), 3+3*(i-1-1)+1:3+3*(i-1)] .= dist_mat#repeat(-2.0 .*n.*a.*d3skew,n)
    end

    for i = 1:params.N*params.n
        k = ceil(Int, i / params.n) #seg_num
        Q[3*Ntot+3*Nbody+3(i-1)+1:3*Ntot+3*Nbody+3*i, 3+(k-1)*3+1:3+k*3] .= iden
    end
    return Q
end

function calc_X3!(X3, X, N, Nbody, b, Nfil, n, d1, d2, d3)
    x0 = @views SVector{3}(X[1:3])
    Δs = 1 / N
    a = Δs / 2

    for i = 1:Nbody
        X3[3*(i-1)+1:3*i] .= x0 + reshape(b[i, :], (1, 3)) * [d1[1, :] d2[2, :] d3[3, :]]
    end

    for i = 1:Nfil

        xfil0 = reshape(x0, 1, 3) + reshape(b[Nbody+i, :], (1, 3)) * [d1[1, :] d2[2, :] d3[3, :]]
        X3_seg = zeros(3 * (N + 1))

        X3_seg[1:3:end] .= xfil0[1] .+ [0; cumsum(Δs * d3[2:end, 1], dims=1)]
        X3_seg[2:3:end] .= xfil0[2] .+ [0; cumsum(Δs * d3[2:end, 2], dims=1)]
        X3_seg[3:3:end] .= xfil0[3] .+ [0; cumsum(Δs * d3[2:end, 3], dims=1)]

        for j = 1:N
            s1 = Nbody + (j - 1) * n + 1
            sn = Nbody + j * n

            Xj = @views SVector{3}(X3_seg[3*(j-1)+1:3*j])
            Xj1 = @views SVector{3}(X3_seg[3*j+1:3*(j+1)])

            dx = (Xj1 - Xj) ./ (2 * n)

            for k = 1:n
                X3[3*Nbody+3*(s1-1)+3*(k-1)+1:3*Nbody+3*(s1-1)+3*k] .= Xj .+ (2 * k - 1) .* dx
            end

        end
    end

end



function generatortoomega!(r)
    r = SVector{3}([r[1] + 1e-16, r[2], r[3]])
    I3 = Diagonal(SVector(true, true, true))
    rmag = norm(r)
    cross_term = sinc(rmag / pi)^2 * skewS(r)

    term3 = @SMatrix [r[2]^2+r[3]^2 -r[1]*r[2] -r[1]*r[3]
        -r[1]*r[2] r[1]^2+r[3]^2 -r[2]*r[3]
        -r[1]*r[3] -r[2]*r[3] r[1]^2+r[2]^2]
    term3 = (rmag - sin(rmag)cos(rmag)) / (rmag^3) * term3
    term3[isnan.(term3)] .= 0 #NaNs should be avoided by + 1e-16 above, can probably remove this
    D = 2 * (I3 + cross_term + -term3)
    return D
end

function calc_RHS_n!(K, t, X, params, d1, d2, d3)
    s = collect(1/params.N:1/params.N:1-1/params.N)

    d1s = d1[3:end, :] .- d1[2:end-1, :]
    d2s = d2[3:end, :] .- d2[2:end-1, :]
    d3s = d3[3:end, :] .- d3[2:end-1, :]

    k1 = N * [dot(d2s[i, :], d3[i+2, :]) for i in 1:size(d2s, 1)]  .- params.strength .* sin.(2 * pi .* s .- t)
    k2 = N * [dot(d3s[i, :], d1[i+2, :]) for i in 1:size(d3s, 1)]
    k3 = N * [dot(d1s[i, :], d2[i+2, :]) for i in 1:size(d1s, 1)]

    # k1 = zeros(params.N - 1)
    # k2 = zeros(params.N - 1)
    # k3 = zeros(params.N - 1)

    kx = zeros(params.N - 1)
    ky = zeros(params.N - 1)
    kz = zeros(params.N - 1)

    kx .= k1 .* d1[3:end, 1] .+ k2 .* d2[3:end, 1] .+ k3 .* d3[3:end, 1]
    ky .= k1 .* d1[3:end, 2] .+ k2 .* d2[3:end, 2] .+ k3 .* d3[3:end, 2]
    kz .= k1 .* d1[3:end, 3] .+ k2 .* d2[3:end, 3] .+ k3 .* d3[3:end, 3]

    K[10-3:3:end] .= -kx
    K[11-3:3:end] .= -ky
    K[12-3:3:end] .= -kz
end

function calc_directors_n!(d1, d2, d3, X, params)
    for i = 0:params.N
        r = @views SVector{3}(X[7+3*(i-1):6+3*i])
        
        q = [cos(norm(r)); sinc(norm(r) / pi) * r]
        q0 = q[1]
        q1 = q[2]
        q2 = q[3]
        q3 = q[4]

        d1[1+i, 1] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3
        d1[1+i, 2] = 2 * (q2 * q1 + q0 * q3)
        d1[1+i, 3] = 2 * (q3 * q1 - q0 * q2)

        d2[1+i, 1] = 2 * (q1 .* q2 - q0 .* q3)
        d2[1+i, 2] = q0 .* q0 - q1 .* q1 + q2 .* q2 - q3 .* q3
        d2[1+i, 3] = 2 * (q3 .* q2 + q0 .* q1)

        d3[1+i, 1] = 2 * (q1 .* q3 + q0 .* q2)
        d3[1+i, 2] = 2 * (q2 .* q3 - q0 .* q1)
        d3[1+i, 3] = q0 .* q0 - q1 .* q1 - q2 .* q2 + q3 .* q3
    end
end

function calc_Qomega_n!(Q, X, params)
    I3 = Diagonal(SVector(true, true, true))
    r = @views SVector{3}(X[4:6])
    Db = generatortoomega!(r)

    Q[1:3, 1:3] = I3
    Q[4:6, 4:6] = Db

    for j = 2:params.N
        rows = 3+3*(j-1)+1:3+3*j
        cols = 3+3*(j-1)+1:3+3*j
        r = @views SVector{3}(X[6+3*(j-1)+1:6+3*j])
        Q[rows, cols] = generatortoomega!(r)
    end
end


struct SimulationParameters
    N::Int
    n::Int
    Nbody::Int
    Ntot::Int
    Nfil::Int
    S::Int
    a::Vector{Float64}
    b::Matrix{Float64}
    strength::Float64
end

struct PreallocatedMatrices
    Mh::Matrix{Float64}
    Mf::Matrix{Float64}
    Q::Matrix{Float64}
    Qomega::Matrix{Float64}
    d1::Matrix{Float64}
    d2::Matrix{Float64}
    d3::Matrix{Float64}
    X3::Vector{Float64}
    K::Vector{Float64}
end

struct CombinedParameters
    params::SimulationParameters
    matrices::PreallocatedMatrices
end