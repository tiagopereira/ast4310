using Unitful

struct Atmosphere1D{}
    nz::Int64
    z
    temperature
    velocity_z
    velocity_turb
    electron_density
    hydrogen1_density  # neutral hydrogen across all levels
    proton_density
end

struct AtomicLine{}
    λ0
    loggf
    γrad
    χ
    σ_ABO
    α_ABO
    const_ABO
    mass
    abundance
end

function blackbody_λ(λ, temp) 
    radiation = 2h * c_0^2 * λ^-5 / (exp(h * c_0 / k_B / (λ * temp)) - 1)
    return radiation |> u"kW / (m^2 * nm * sr)"
end

"""
Formal solver, single wavelength
"""
function formal_simple(α, S, z)
    nz, ny, nx = size(α)
    intensity = zeros(eltype(S), ny, nx)
    for i in 1:nx, j in 1:ny
        tau = cumul_integrate(-z, α[:, j, i])
        intensity[j, i] = integrate(tau, S[:, j, i] .* exp.(-tau))
    end
    return intensity
end


const Ar_H = 1.007975  # Atomic weight of hydrogen
const hc_k = h * c_0 / k_B
const hminusχ = 0.754u"eV"
const saha_const = h^2 / (2 * π * m_e * k_B)  # Constant for Saha equation
const σ_thomson = (e^4 / (6π * ε_0^2 * m_e^2 * c_0^4)) |> u"m^2"  
const αl_const = (e^2 / (4 * ε_0 * m_e * c_0^2)) |> u"m"   # constant for line extinction

@derived_dimension NumberDensity Unitful.𝐋^-3
@derived_dimension PerLength Unitful.𝐋^-1

#=----------------------------------------------------------------------------
        Recipes from Wishart (1979) and Broad and Reinhardt (1976)
----------------------------------------------------------------------------=#
const wbr_λ = [     18, 19.6, 21.4, 23.6, 26.4, 29.8, 34.3, 40.4, 49.1, 62.6,  121, 139,
                   164,  175,  200,  225,  250,  275,  300,  325,  350,  375,  400, 425,
                   450,  475,  500,  525,  550,  575,  600,  625,  650,  675,  700, 725,
                   750,  775,  800,  825,  850,  875,  900,  925,  950,  975, 1000, 1025,
                  1050, 1075, 1100, 1125, 1150, 1175, 1200, 1225, 1250, 1275, 1300, 1325,
                  1350, 1375, 1400, 1425, 1450, 1475, 1500, 1525, 1550, 1575, 1600, 1610,
                  1620, 1630]   # in nm
const wbr_σ = [0.067, 0.088, 0.117, 0.155, 0.206, 0.283, 0.414, 0.703,  1.24,  2.33,
                5.43,  5.91,  7.29, 7.918, 9.453, 11.08, 12.75, 14.46, 16.19, 17.92,
               19.65, 21.35, 23.02, 24.65, 26.24, 27.77, 29.23, 30.62, 31.94, 33.17,
               34.32, 35.37, 36.32, 37.17, 37.91, 38.54, 39.07, 39.48, 39.77, 39.95,
               40.01, 39.95, 39.77, 39.48, 39.06, 38.53, 37.89, 37.13, 36.25, 35.28,
               34.19, 33.01, 31.72, 30.34, 28.87, 27.33, 25.71, 24.02, 22.26, 20.46,
               18.62, 16.74, 14.85, 12.95, 11.07, 9.211, 7.407, 5.677, 4.052, 2.575,
               1.302, 0.8697, 0.4974, 0.1989]  # in 1e-22 m^2
const wbr_bf_interp = linear_interpolation(wbr_λ, wbr_σ, extrapolation_bc=0)

#=----------------------------------------------------------------------------
                            Recipes from Stilley
----------------------------------------------------------------------------=#
const stilley_ff_λ = [0.0, 303.8, 455.6, 506.3, 569.5, 650.9, 759.4, 911.3,
                      1013.0, 1139.0, 1302.0, 1519.0, 1823.0, 2278.0, 3038.0,
                      4556.0, 9113.0]  # in nm
const stilley_ff_t = 5040.0 ./ [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
                                1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]  # in K
const stilley_ff_table =
   [[0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00  #=    0.0 nm
=#   0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00];
    [3.44e-02 4.18e-02 4.91e-02 5.65e-02 6.39e-02 7.13e-02 7.87e-02 8.62e-02  #=  303.8 nm
=#   9.36e-02 1.01e-01 1.08e-01 1.16e-01 1.23e-01 1.30e-01 1.38e-01 1.45e-01];
    [7.80e-02 9.41e-02 1.10e-01 1.25e-01 1.40e-01 1.56e-01 1.71e-01 1.86e-01  #=  455.6 nm
=#   2.01e-01 2.16e-01 2.31e-01 2.45e-01 2.60e-01 2.75e-01 2.89e-01 3.03e-01];
    [9.59e-02 1.16e-01 1.35e-01 1.53e-01 1.72e-01 1.90e-01 2.08e-01 2.25e-01  #=  506.3 nm
=#   2.43e-01 2.61e-01 2.78e-01 2.96e-01 3.13e-01 3.30e-01 3.47e-01 3.64e-01];
    [1.21e-01 1.45e-01 1.69e-01 1.92e-01 2.14e-01 2.36e-01 2.58e-01 2.80e-01  #=  569.5 nm
=#   3.01e-01 3.22e-01 3.43e-01 3.64e-01 3.85e-01 4.06e-01 4.26e-01 4.46e-01];
    [1.56e-01 1.88e-01 2.18e-01 2.47e-01 2.76e-01 3.03e-01 3.31e-01 3.57e-01  #=  650.9 nm
=#   3.84e-01 4.10e-01 4.36e-01 4.62e-01 4.87e-01 5.12e-01 5.37e-01 5.62e-01];
    [2.10e-01 2.53e-01 2.93e-01 3.32e-01 3.69e-01 4.06e-01 4.41e-01 4.75e-01  #=  759.4 nm
=#   5.09e-01 5.43e-01 5.76e-01 6.08e-01 6.40e-01 6.72e-01 7.03e-01 7.34e-01];
    [2.98e-01 3.59e-01 4.16e-01 4.70e-01 5.22e-01 5.73e-01 6.21e-01 6.68e-01  #=  911.3 nm
=#   7.15e-01 7.60e-01 8.04e-01 8.47e-01 8.90e-01 9.32e-01 9.73e-01 1.01e+00];
    [3.65e-01 4.39e-01 5.09e-01 5.75e-01 6.39e-01 7.00e-01 7.58e-01 8.15e-01  #= 1013.0 nm
=#   8.71e-01 9.25e-01 9.77e-01 1.03e+00 1.08e+00 1.13e+00 1.18e+00 1.23e+00];
    [4.58e-01 5.50e-01 6.37e-01 7.21e-01 8.00e-01 8.76e-01 9.49e-01 1.02e+00  #= 1139.0 nm
=#   1.09e+00 1.15e+00 1.22e+00 1.28e+00 1.34e+00 1.40e+00 1.46e+00 1.52e+00];
    [5.92e-01 7.11e-01 8.24e-01 9.31e-01 1.03e+00 1.13e+00 1.23e+00 1.32e+00  #= 1302.0 nm
=#   1.40e+00 1.49e+00 1.57e+00 1.65e+00 1.73e+00 1.80e+00 1.88e+00 1.95e+00];
    [7.98e-01 9.58e-01 1.11e+00 1.25e+00 1.39e+00 1.52e+00 1.65e+00 1.77e+00  #= 1519.0 nm
=#   1.89e+00 2.00e+00 2.11e+00 2.21e+00 2.32e+00 2.42e+00 2.51e+00 2.61e+00];
    [1.14e+00 1.36e+00 1.58e+00 1.78e+00 1.98e+00 2.17e+00 2.34e+00 2.52e+00  #= 1823.0 nm
=#   2.68e+00 2.84e+00 3.00e+00 3.15e+00 3.29e+00 3.43e+00 3.57e+00 3.70e+00];
    [1.77e+00 2.11e+00 2.44e+00 2.75e+00 3.05e+00 3.34e+00 3.62e+00 3.89e+00  #= 2278.0 nm
=#   4.14e+00 4.39e+00 4.63e+00 4.86e+00 5.08e+00 5.30e+00 5.51e+00 5.71e+00];
    [3.10e+00 3.71e+00 4.29e+00 4.84e+00 5.37e+00 5.87e+00 6.36e+00 6.83e+00  #= 3038.0 nm
=#   7.28e+00 7.72e+00 8.14e+00 8.55e+00 8.95e+00 9.33e+00 9.71e+00 1.01e+01];
    [6.92e+00 8.27e+00 9.56e+00 1.08e+01 1.19e+01 1.31e+01 1.42e+01 1.52e+01  #= 4556.0 nm
=#   1.62e+01 1.72e+01 1.82e+01 1.91e+01 2.00e+01 2.09e+01 2.17e+01 2.25e+01];
    [2.75e+01 3.29e+01 3.80e+01 4.28e+01 4.75e+01 5.19e+01 5.62e+01 6.04e+01  #= 9133.0 nm
=#   6.45e+01 6.84e+01 7.23e+01 7.60e+01 7.97e+01 8.32e+01 8.67e+01 9.01e+01]]
const stilley_ff_interp = linear_interpolation((stilley_ff_λ, stilley_ff_t[end:-1:1]),
                             stilley_ff_table[:, end:-1:1], extrapolation_bc=Line());

"""
    calc_hminus_density(
        h_neutral_density::NumberDensity,
        temperature::Unitful.Temperature,
        electron_density::NumberDensity
    )

Compute H minus populations based on electron density, temperature, and
density of neutral hydrogen atoms.
"""
function calc_hminus_density(
    h_neutral_density::NumberDensity,
    temperature::Unitful.Temperature,
    electron_density::NumberDensity
)
    tmp = (ustrip(saha_const) / ustrip(temperature |> u"K"))^(3/2) * u"m^3"
    ϕ = tmp * exp(hminusχ / (k_B * temperature)) / 4
    return (ϕ * h_neutral_density) * electron_density
end


"""
    α_thomson(electron_density::NumberDensity)

Compute the Thomson extinction as a function of electron density.
"""
α_thomson(electron_density::NumberDensity) = σ_thomson * electron_density


"""
    α_hminus_bf(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        h_neutral_density::NumberDensity,
        electron_density::NumberDensity
    )

Compute extinction from H minus ion, from input H minus populations. Uses recipe from
[Wishart (1979)](https://ui.adsabs.harvard.edu/abs/1979MNRAS.187P..59W) for λ down to 175 nm,
and recipe from [Broad and Reinhardt (1976)](https://ui.adsabs.harvard.edu/abs/1976PhRvA..14.2159B)
for λ=164 nm and below, following the recommendation from Mathisen (1984, MSC thesis).
"""
function α_hminus_bf(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
    electron_density::NumberDensity
)
    λi = ustrip(λ |> u"nm")   # convert to units of table
    κ = wbr_bf_interp(λi)::Float64 * 1e-22u"m^2"
    stimulated_emission = exp(-hc_k / (λ * temperature))
    h_minus_density = calc_hminus_density(h_neutral_density, temperature, electron_density)
    return (κ * h_minus_density * (1 - stimulated_emission)) |> u"m^-1"
end


"""
    σ_hminus_ff(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        h_neutral_density::NumberDensity,
        electron_density::NumberDensity
    )

Compute free-free cross section coefficient from H minus ion, for a given wavelength
and temperature. Units are m^5, needs to be multiplied by electron density and
density of neutral hydrogen atoms to obtain linear extinction. Interpolates table from
[Stilley & Callaway (1970)](https://ui.adsabs.harvard.edu/abs/1970ApJ...160..245S/abstract),
page 255, which is valid for λ up to 9113 nm.
"""
function α_hminus_ff(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
    electron_density::NumberDensity
)
    λi = ustrip(λ |> u"nm")   # convert to units of table
    temp = ustrip(temperature |> u"K")
    kappa = max(0.0, stilley_ff_interp(λi, temp)::Float64) * 1e-29u"m^4/N"
    σ = k_B * temperature * kappa 
    return (σ * h_neutral_density * electron_density) |> u"m^-1"
end


function α_cont(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
    electron_density::NumberDensity
)
    return (α_hminus_bf(λ, temperature, h_neutral_density, electron_density) + 
            α_hminus_ff(λ, temperature, h_neutral_density, electron_density) +
            α_thomson(electron_density))
end


const invSqrtPi = 1. / sqrt(π)

function humlicek(z::Complex)
    s = abs(real(z)) + imag(z)
    if s > 15.0
        # region I
        w = im * invSqrtPi * z / (z * z - 0.5)
    elseif s > 5.5
        # region II
        zz = z * z
        w = im * (z * (zz * invSqrtPi - 1.4104739589)) / (0.75 + zz * (zz - 3.0))
    else
        x, y = real(z), imag(z)
        t = y - im * x
        if y >= 0.195 * abs(x) - 0.176
            # region III
            w = ((16.4955 + t * (20.20933 + t * (11.96482 + t * (3.778987 + 0.5642236 * t))))
               / (16.4955 + t * (38.82363 + t * (39.27121 + t * (21.69274 + t * (6.699398 + t))))))
        else
            # region IV
            u = t * t
            nom = t * (36183.31 - u * (3321.99 - u * (1540.787 -  u *
                   (219.031 - u * (35.7668 - u * (1.320522 - u * .56419))))))
            den = 32066.6 - u * (24322.8 - u * (9022.23 - u * (2186.18 -
                    u * (364.219 - u * (61.5704 - u * (1.84144 - u))))))
            w = exp(u) - nom / den
        end
    end
    return w
end

"""
About twice as fast as voigt()
"""
function voigt(a::T, u::AbstractFloat)::T where T <: AbstractFloat
    z = u + a * im
    return real(humlicek(z))
end

"""
Damping constant.
"""
function damping(γ::Unitful.Frequency, λ::Unitful.Length, ΔλD::Unitful.Length)
    c1 = 1 / (4 * π * c_0)
    return ustrip((c1 * γ * λ^2 / ΔλD) |> u"m/m")
end


"""
    function γ_barklem(
        α::AbstractFloat,
        barklem_const::Unitful.VolumeFlow,
        temperature::Unitful.Temperature,
        h_neutral_density::NumberDensity,
    )

Compute van der Waals broadening from collisions with neutral hydrogen atoms following the
theory from Barklem/O'Mara/Anstee.

# Arguments
- `α::AbstractFloat`: velocity exponent from Barklem/O'Mara/Anstee tables
- `barklem_const::Unitful.VolumeFlow`: atmosphere-independent constant computed from
   `const_barklem()`.
- `temperature::Unitful.Temperature`
- `h_neutral_density::NumberDensity`: number density of neutral hydrogen atoms

# Returns
- `γ::Unitful.Frequency`: broadening in units of rad / s.

# Examples
```
julia> bconst = const_barklem(my_atom_weight, 0.275, 291)
7.495208174533257e-16 m³ rad s⁻¹

julia> γ = γ_barklem(0.275, bconst, 6000u"K", 1e23u"m^-3")
1.3596876505340942e11 rad s⁻¹
```
"""
function γ_barklem(
    α::AbstractFloat,
    barklem_const::Unitful.VolumeFlow,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
)
    return barklem_const * ustrip(temperature |> u"K")^((1 - α)/2) * h_neutral_density
end


function const_barklem(atomic_mass::Unitful.Mass, α::Real, σ::Real)
    α < 0 && throw(DomainError("α must be non-negative"))
    σ < 0 && throw(DomainError("σ must be non-negative"))
    μ = m_u / (1 / Ar_H + 1 / (atomic_mass / m_u))
    # Using 1 K to keep units right for later multiplication by correct temperature
    v_bar = sqrt(8 * k_B * u"K"/ (π * μ)) |> u"m/s"
    v_ratio = (1e4u"m/s" / v_bar) |> u"m/m"
    # Squared Bohr radius is to convert from atomic units to m^2, factor of 2 from HW to FW
    return (a_0^2 * 2 * (4 / π)^(α / 2) * gamma((4 - α) / 2) * v_bar * σ *
            v_ratio^α) |> u"m^3 * rad / s"
end