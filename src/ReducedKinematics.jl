
# --- Volumetric deformation ---

"""
Return the isotropic (volumetric) deformation gradient corresponding to a
given Jacobian `J`, i.e. `F = J^(1/3) * I`.
"""
function F_volumetric(J::Real)
  λ = J^(1/3)
  TensorValue(λ, 0, 0, 0, λ, 0, 0, 0, λ)
end

"""
Return the thermally-induced volumetric Jacobian `J` for a given physical
model and temperature. Defaults to `1.0` (no thermal expansion) for models
that do not couple volume to temperature.
"""
function J_thermal(::PhysicalModel, ::Float64)
  1.0
end

"""
Return the thermally-induced volumetric Jacobian `J` for a
[`ThermalVolumetric`](@ref) model at a given temperature `θ`, solved by a
Newton iteration on the zero-pressure condition `pressure(J) = 0`.
"""
function J_thermal(m::ThermalVolumetric, θ::Float64)
  _, ∂Ψ∂F, _, ∂∂Ψ∂FF, _, _ = m()
  pressure(J) = 1/3 * tr(∂Ψ∂F(F_volumetric(J), θ)) * J^(-2/3)
  ∂pressure∂J(J) = 1/9 * tr(∂∂Ψ∂FF(F_volumetric(J), θ) ⊙ I3) * J^(-4/3) - 2/9 * tr(∂Ψ∂F(F_volumetric(J), θ)) * J^(-5/3)

  J0 = 1.0
  p0 = pressure(J0)

  tol = abs(p0) * 1e-10
  maxiter = 20

  for _ in 1:maxiter
    if abs(p0) < tol
      return J0
    end
    dp = ∂pressure∂J(J0)
    J0 -= p0 / dp
    p0 = pressure(J0)
  end
  @debug "Jacobian not converged after $(maxiter) iterations, with J=$(J0) and p=$(p0) at θ=$(θ)"
  J0
end

"""
Return the thermally-induced volumetric Jacobian for a given physical model
and [`AbstractCondition`](@ref), by extracting the temperature from the
condition. Defaults to `1.0` for models that do not couple volume to
temperature.
"""
function J_thermal(m::PhysicalModel, c::AbstractCondition)
  1.0
end

"""
Return the thermally-induced volumetric Jacobian for a `ThermoMechano` model
under the given [`AbstractCondition`](@ref).
"""
function J_thermal(m::ThermoMechano, c::AbstractCondition)
  θ = temperature(c)
  J_thermal(m, θ)
end


# --- Deformation gradient ---
 
"""
Return the purely volumetric deformation gradient for a given physical model
at temperature `θ`, with the Jacobian solved by [`J_thermal`](@ref).
"""
function calculate_F(m::PhysicalModel, θ::Float64)
  J = J_thermal(m, θ)
  F_volumetric(J)
end

"""
Return the deformation gradient for a [`Uniaxial`](@ref) kinematics test at
stretch `λ` under the given [`AbstractCondition`](@ref), including any
thermally-induced volumetric contribution.
"""
function calculate_F(m::PhysicalModel, ::Type{Uniaxial}, λ::Real, c::AbstractCondition)
  J = J_thermal(m, c)
  TensorValue(λ, 0, 0, 0, λ^(-1/2), 0, 0, 0, λ^(-1/2)) .* J^(1/3)
end

"""
Return the deformation gradient for a [`Biaxial`](@ref) kinematics test at
stretch `λ` under the given [`AbstractCondition`](@ref), including any
thermally-induced volumetric contribution.
"""
function calculate_F(m::PhysicalModel, ::Type{Biaxial}, λ::Real, c::AbstractCondition)
  J = J_thermal(m, c)
  TensorValue(λ, 0, 0, 0, λ, 0, 0, 0, λ^(-2)) .* J^(1/3)
end

"""
Return the deformation gradient for an independently biaxial test with
stretches `λ1` and `λ2` under the given [`AbstractCondition`](@ref),
including any thermally-induced volumetric contribution.
"""
function calculate_F(m::PhysicalModel, λ1::Real, λ2::Real, c::AbstractCondition)
  J = J_thermal(m, c)
  TensorValue(λ1, 0, 0, 0, λ2, 0, 0, 0, (λ1*λ2)^(-1)) .* J^(1/3)
end
