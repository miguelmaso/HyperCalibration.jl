
function F_volumetric(J::Real)
  λ = J^(-1/3)
  TensorValue(λ, 0, 0, 0, λ, 0, 0, 0, λ)
end

function J_thermal(::PhysicalModel, ::AbstractCondition)
  1.0
end

function J_thermal(m::ThermalVolumetric, c::AbstractCondition)
  _, ∂Ψ∂F, _, ∂∂Ψ∂FF, _, _ = m()
  θ = temperature(c)
  pressure(J) = 1/3 * tr(∂Ψ∂F(F_volumetric(J), θ)) * J^(-2/3)
  ∂pressure∂J(J) = 1/9 * tr(∂∂Ψ∂FF(F_volumetric(J), θ) * I3) * J^(-4/3) - 2/9 * tr(∂Ψ∂F(F_volumetric(J), θ)) * J^(-5/3)

  J0 = 1.0
  p0 = pressure(J0)

  tol = 1e-9
  maxiter = 20

  for _ in 1:maxiter
    if abs(p0) < tol
      return J0
    end
    dp = ∂pressure∂J(J0)
    if dp == 0
      break
    end
    J0 -= p0 / dp
  end
  J0
end

function calculate_F(m::PhysicalModel, c::AbstractCondition)
  J = J_thermal(m, c)
  F_volumetric(J)
end

function calculate_F(m::PhysicalModel, ::Type{Uniaxial}, λ::Real, c::AbstractCondition)
  J = J_thermal(m, c)
  TensorValue(λ, 0, 0, 0, λ^(-1/2), 0, 0, 0, λ^(-1/2)) .* J^(1/3)
end

function calculate_F(m::PhysicalModel, ::Type{Biaxial}, λ::Real, c::AbstractCondition)
  J = J_thermal(m, c)
  TensorValue(λ, 0, 0, 0, λ, 0, 0, 0, λ^(-2)) .* J^(1/3)
end

function calculate_F(m::PhysicalModel, λ1::Real, λ2::Real, c::AbstractCondition)
  J = J_thermal(m, c)
  TensorValue(λ1, 0, 0, 0, λ2, 0, 0, 0, (λ1*λ2)^(-1)) .* J^(1/3)
end
