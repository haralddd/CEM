using RayleighSolver

"""
    kspp(below::Material)
Calculates the parallel wavenumber for the upper medium in units of k0

# Arguments:
- `below`: [`Material`](@ref) containing information about the lower medium

# Returns:
- Parallel wavenumber for the upper medium
"""
function kspp(below::Material)
    if below isa Uniaxial
        A = get_A(below)
        eps_para = below.eps_para
        mu_para = below.mu_para
    elseif below isa Isotropic
        A = 1.0
        eps_para = below.eps
        mu_para = below.mu
    end
    return sqrt(eps_para * (mu_para - eps_para) / (A - eps_para^2))
end

"""
    ksmp(below::Material)
Calculates the perpendicular wavenumber for the upper medium in units of k0

# Arguments:
- `below`: [`Material`](@ref) containing information about the lower medium

# Returns:
- Perpendicular wavenumber for the upper medium
"""
function ksmp(below::Material)
    if below isa Uniaxial
        eps_para = below.eps_para
        mu_para = below.mu_para
    elseif below isa Isotropic
        eps_para = below.eps
        mu_para = below.mu
    end
    return sqrt(mu_para * (eps_para - mu_para) / (1 - mu_para^2))
end

LSP(k) = 1.0 / (2imag(k))
θmax(km, kspp) = asind(real(km - kspp))

eps = -7.5 + 0.24im   # permittivity of metallic film

metallic1 = Uniaxial(eps, 0.5 * eps, 1.0, 1.0)
metallic2 = Uniaxial(eps, 1.0 * eps, 1.0, 1.0)
metallic3 = Uniaxial(eps, 1.5 * eps, 1.0, 1.0)

hm11 = Uniaxial(eps, -0.5 * eps, 1.0, 1.0)
hm12 = Uniaxial(eps, -1.0 * eps, 1.0, 1.0)
hm13 = Uniaxial(eps, -1.5 * eps, 1.0, 1.0)

hm21 = Uniaxial(-0.5 * eps, eps, 1.0, 1.0)
hm22 = Uniaxial(-1.0 * eps, eps, 1.0, 1.0)
hm23 = Uniaxial(-1.5 * eps, eps, 1.0, 1.0)

kmm1 = kspp(metallic1)
kmm2 = kspp(metallic2)
kmm3 = kspp(metallic3)

khm21 = kspp(hm21)
khm22 = kspp(hm22)
khm23 = kspp(hm23)

Lmm1 = LSP(kmm1) / 2π
Lmm2 = LSP(kmm2) / 2π
Lmm3 = LSP(kmm3) / 2π

Lhm21 = LSP(khm21) / 2π
Lhm22 = LSP(khm22) / 2π
Lhm23 = LSP(khm23) / 2π

θmm1 = θmax(kmm1, 0.782)
θmm2 = θmax(kmm2, 0.782)
θmm3 = θmax(kmm3, 0.782)

θhm21 = θmax(khm21, 0.782)
θhm22 = θmax(khm22, 0.782)
θhm23 = θmax(khm23, 0.782)


# Create uniaxial medium with anisotropic properties
bianisotropy = Uniaxial(-1.5+0.0im, -3.0+0.1im, -0.4+0.0im, -0.4+0.0im)

kspp_bi = kspp(bianisotropy)
ksmp_bi = ksmp(bianisotropy)

Lbi1 = LSP(kspp_bi) / 2π
Lbi2 = LSP(ksmp_bi) / 2π

θbi1 = θmax(kspp_bi, 0.782)
θbi2 = θmax(ksmp_bi, 0.782)


sg1 = Uniaxial(3.688 + 0.017im, -0.675 + 0.072im, 1.0+0.0im, 1.0+0.0im)
sg2 = Uniaxial(6.425 + 0.088im, -2.625 + 0.120im, 1.0 + 0.0im, 1.0 + 0.0im)
sg3 = Uniaxial(24.803 + 1.846im, -4.575 + 0.168im, 1.0 + 0.0im, 1.0 + 0.0im)
kspp_sg1 = kspp(sg1)
kspp_sg2 = kspp(sg2)
kspp_sg3 = kspp(sg3)
Lsg1 = LSP(kspp_sg1) / 2π
Lsg2 = LSP(kspp_sg2) / 2π
Lsg3 = LSP(kspp_sg3) / 2π
θsg1 = θmax(kspp_sg1, 0.782)
θsg2 = θmax(kspp_sg2, 0.782)
θsg3 = θmax(kspp_sg3, 0.782)
