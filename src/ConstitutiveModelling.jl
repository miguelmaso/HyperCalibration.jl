
# --- Tensile tests: elastic material ---

"""
Evaluate the stress predicted by a constitutive model
under the given mechanical protocol, conditions, and geometry of an experimental test.
"""
function evaluate_stress(model::Elasto, protocol::MechanicalProtocol{K}, cond::AbstractCondition, ::AbstractGeometry)
  P_func = model()[2]
  map(stretches(protocol)) do λ
    F = calculate_F(model, K, λ, cond)
    P = P_func(F)
    p = -P[3,3] * F[3,3]  # Volumetric pressure term
    return P[1,1] + p / F[1,1]
  end
end

function evaluate_stress(model::ThermoMechano{<:Any,<:Elasto}, protocol::MechanicalProtocol{K}, cond::AbstractCondition, ::AbstractGeometry)
  P_func = model()[2]
  θ = temperature(cond)
  map(stretches(protocol)) do λ
    F = calculate_F(model, K, λ, cond)
    P = P_func(F,θ)
    p = -P[3,3] * F[3,3]  # Volumetric pressure term
    return P[1,1] + p / F[1,1]
  end
end

# --- Tensile tests: viscoelastic material ---

function new_state(model::ViscoElastic, F, Fn, A...)
  map(model.branches, A) do b, Ai
    _, Se, ∂Se∂Ce = SecondPiola(b.elasto)
    HyperFEM.PhysicalModels.ReturnMapping(b, Se, ∂Se∂Ce, F, Fn, Ai)[2]
  end
end

function evaluate_stress(model::ViscoElastic, protocol::MechanicalProtocol{K}, cond::AbstractCondition, ::AbstractGeometry)
  update_time_step!(model, time_step(protocol))
  P_func = model()[2]
  n  = length(model.branches)
  A  = ntuple(_ -> VectorValue(I3..., 0.0), Val(n))
  Fn = calculate_F(model, K, 1.0, cond)
  map(stretches(protocol)) do λ
    F = calculate_F(model, K, λ, cond)
    P = try P_func(F, Fn, A...) catch; zeros(3,3) end
    A = try new_state(model, F, Fn, A...) catch; A end
    Fn = F               # Update the previous deformation gradient for the next iteration
    p = -P[3,3] * F[3,3]  # Volumetric pressure term
    return P[1,1] + p / F[1,1]
  end
end

function evaluate_stress(model::ThermoMechano{<:Any,<:ViscoElastic}, protocol::MechanicalProtocol{K}, cond::AbstractCondition, ::AbstractGeometry)
  update_time_step!(model, time_step(protocol))
  θ = temperature(cond)
  P_func = model()[2]
  n  = length(model.branches)
  A  = ntuple(_ -> VectorValue(I3..., 0.0), Val(n))
  Fn = calculate_F(model, K, 1.0, cond)
  map(stretches(protocol)) do λ
    F = calculate_F(model, K, λ, cond)
    P = try P_func(F, θ, Fn, A...) catch; zeros(3,3) end
    A = try new_state(model, F, θ, Fn, A...) catch; A end
    Fn = F               # Update the previous deformation gradient for the next iteration
    p = -P[3,3] * F[3,3]  # Volumetric pressure term
    return P[1,1] + p / F[1,1]
  end
end

# --- Tensile tests: electro-mechanical material ---

function evaluate_stress(model::ThermoElectroMechano{<:Any,<:Electro,<:ViscoElastic}, protocol::MechanicalProtocol{K}, cond::AbstractCondition, geom::AbstractGeometry)
  update_time_step!(model, Δt)
  θ = temperature(cond)
  E0 = electric_field(cond, geom)
  P_func, ∂P_func = model()[[2,5]]
  n  = length(model.mechano.branches)
  A  = ntuple(_ -> VectorValue(I3..., 0.0), Val(n))
  λ2 = 1.0
  Fn = calculate_F(model, 1.0, λ2, cond)
  
  function evaluate_P(λ1, λ2, E)
    Fi = calculate_F(model, λ1, λ2, cond)
    Pi = P_func(Fi, E, θ, Fn, A...)
    p_ext = -Pi[3,3] * Fi[3,3]
    P_tot = Pi + p_ext*inv(Fi)
    return P_tot, Fi
  end

  function evaluate_∂P22_∂λ2(λ1, λ2, E)
    Fi = calculate_F(model, λ1, λ2, cond)
    Pi = P_func(Fi, E, θ, Fn, A...)
    ∂Pi = ∂P_func(Fi, E, θ, Fn, A...)
    ∂Piso22_∂λ22 = ∂Pi[5,5] - ∂Pi[5,9]*Fi[3,3]/λ2
    ∂Piso33_∂λ22 = ∂Pi[9,5] - ∂Pi[9,9]*Fi[3,3]/λ2
    P22 = Pi[2,2] -Pi[3,3]*Fi[3,3]/Fi[2,2]
    ∂P22_∂λ2 = ∂Piso22_∂λ22 - ∂Piso33_∂λ22*Fi[3,3]/λ2 + 2.0*Pi[3,3]*Fi[3,3]/λ2^2
    return P22, ∂P22_∂λ2
  end

  function evaluate_P11_impl(λ1, E, λ2_guess)
    P22, dP22_dλ2 = evaluate_∂P22_∂λ2(λ1, λ2_guess, E)
    tol = 1e-6
    iter = 0
    maxiter = 10
    while abs(P22) > tol && iter < maxiter
      λ2_guess -= P22 / dP22_dλ2  # Update λ2
      P22, dP22_dλ2 = evaluate_∂P22_∂λ2(λ1, λ2_guess, E)
      iter += 1
    end
    if iter == maxiter
      @warn "Not converged, V=$V, θ=$θ, λ=$λ"
    end
    P, F = evaluate_P(λ1, λ2, E)
    A = new_state(model.mechano, F, Fn, A...)
    Fn = F
    return P[1,1], λ2_guess
  end

  for Ei in range(0.0, E0, length=10)  # Incrementally apply initial voltage
    _, λ2 = evaluate_P11_impl(λ1, Ei, λ2)
  end
  
  P_values = zeros(length(stretches(protocol)))
  for (i, λ) in enumerate(stretches(protocol))
    P11, λ2 = evaluate_P11_impl(λ, E0, λ2)
    P_values[i] = P11
  end
  return P_values
end

# --- Thermal tests ---

"""
Evaluate the specific heat capacity predicted by a constitutive model
under the given thermal protocol, conditions and geometry of an experimental test.
"""
function evaluate_cv(model::ThermoMechano, protocol::TemperatureSweepProtocol, ::AbstractCondition, ::AbstractGeometry)
  J(θ) = J_temp(model, θ)
  ∂∂Ψ = model()[5]
  if model.mechano isa Elasto
    return map(θ -> -θ*∂∂Ψ(calculate_F(model, θ), θ), temperatures(protocol))
  else
    update_time_step!(model, 1.0)
    n = length(model.mechano.branches)
    A = ntuple(_ -> VectorValue(I3..., 0.0), Val(n))
    return map(θ -> -θ*∂∂Ψ(calculate_F(model, θ), θ, calculate_F(model, θ), A...), temperatures(protocol))
  end
end

# --- Dielectric tests ---

"""
Evaluate the dielectric permittivity predicted by a constitutive model
under the given thermal protocol, conditions and geometry of an experimental test.
"""
function evaluate_epsilon(model::ThermoElectro, protocol::TemperatureSweepProtocol, ::AbstractCondition, ::AbstractGeometry)
  ∂∂Ψ∂EE = model()[6]
  F1 = F_volumetric(1.0)
  E0 = 0.0
  map(θi -> -1/ϵ0*∂∂Ψ∂EE(F1, E0, θi)[1], temperatures(protocol))
end
