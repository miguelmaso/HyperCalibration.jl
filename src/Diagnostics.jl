
using LinearAlgebra, Statistics

# --- Storage & Wrapper Structs ---

"""
A struct that holds the results of a calibration process.
"""
struct CalibrationResult{B, P, D}
    builder::B
    params::P
    data::D
    model::PhysicalModel
    CalibrationResult(builder::B, params::P, data::D) where {B,P,D} = new{B,P,D}(builder, params, data, builder(params...))
end

"""
A struct returned by [`parameter_stats`](@ref), gathering the
standard errors, confidence intervals and normalized sensitivity
for the calibrated parameters.
"""
struct ParameterStats{P, E, C, S}
    names::Vector{String}
    params::P
    std_errs::E
    ci_bounds::C       # Matrix/Tuple of (lower, upper)
    sensitivities::S
    r2::Float64
end


# --- Residuals & Jacobian Engine ---

function residuals(model::PhysicalModel, data)
    y_true, y_pred = experiment_prediction(model, data)
    return y_true .- y_pred
end

function finite_difference_jacobian(model_builder, params, data; h_rel=1e-5)
    r0 = residuals(model_builder(params...), data)
    n_res = length(r0)
    n_params = length(params)
    J = Matrix{Float64}(undef, n_res, n_params)
    
    p_step = copy(params)
    for j in 1:n_params
        h = h_rel * max(abs(params[j]), 1.0)
        
        p_step[j] = params[j] + h
        r_plus = residuals(model_builder(p_step...), data)
        
        p_step[j] = params[j] - h
        r_minus = residuals(model_builder(p_step...), data)
        
        p_step[j] = params[j] # reset
        @. J[:, j] = (r_plus - r_minus) / (2 * h)
    end
    return J
end


# --- Covariance Matrix via SVD ---

"""
Compute the covariance matrix of a [`CalibrationResult`](@ref) via
the Jacobian of the residual ``J_{ij} = \\frac{\\partial r_i}{\\partial p_j}`` evaluated with central differences.
"""
function covariance_matrix(res::CalibrationResult)
    J = finite_difference_jacobian(res.builder, res.params, res.data)
    r = residuals(res.model, res.data)
    
    n_data = length(r)
    n_params = length(res.params)
    n_dof = max(1, n_data - n_params)
    
    sse = sum(abs2, r)
    res_var = sse / n_dof
    
    # Stable pseudoinverse of J'J using SVD of J directly
    # J = U * S * V' => (J'J)^+ = V * S^-2 * V'
    F = svd(J)
    tol = maximum(F.S) * max(n_data, n_params) * eps(Float64)
    inv_S2 = [s > tol ? 1.0 / (s^2) : 0.0 for s in F.S]
    
    cov_mat = res_var .* (F.V * Diagonal(inv_S2) * F.V')
    JtJ = J' * J
    
    return Symmetric(cov_mat), JtJ, sse
end


# --- Predictions & Goodness of Fit ---

"""Calculate the R² score for the model predictions against the experimental data."""
function r2_score(model::PhysicalModel, data)
    y_true, y_pred = experiment_prediction(model, data)
    ss_res = sum(abs2, y_true .- y_pred)
    ss_tot = sum(abs2, y_true .- mean(y_true))
    return 1.0 - (ss_res / ss_tot)
end

r2_score(res::CalibrationResult) = r2_score(res.model, res.data)


# --- Parameter Statistics & Pretty Printing ---

function t_critical_975(dof::Int)
  dof >= 30 && return 1.96
  # Rational approximation for two-tailed 95% t-critical values without Distributions.jl
  return 1.96 + (2.378 / dof) + (2.64 / (dof^2))
end

"""
Evaluate the standard errors, confidence intervals and
normalized sensitivity for the calibrated parameters.
Returns a [`ParameterStats`](@ref) struct.
"""
function parameter_stats(res::CalibrationResult; names=map(i -> "p$i", 1:length(res.params)))
  n_params = length(res.params)
  n_data = npoints(res.data)
  n_dof = max(1, n_data - n_params)

  cov_mat, JtJ, sse = covariance_matrix(res)
  t_crit = t_critical_975(n_dof)

  std_errs = sqrt.(abs.(diag(cov_mat)))
  ci_lower = res.params .- t_crit .* std_errs
  ci_upper = res.params .+ t_crit .* std_errs

  # Normalized sensitivity coefficient: (∂²E/∂p²) * p² / SSE ≈ 2*(J'J)_ii * p² / SSE
  sensitivities = [2JtJ[i,i] * res.params[i]^2 / sse for i in 1:n_params]

  return ParameterStats(names, res.params, std_errs, zip(ci_lower, ci_upper) |> collect, sensitivities, r2_score(res))
end

Base.show(io::IO, stats::ParameterStats) = show(io, MIME("text/plain"), stats)

function Base.show(io::IO, ::MIME"text/plain", stats::ParameterStats)
  println(io, "Model Calibration Summary (R² = $(round(stats.r2, digits=4))):")
  println(io, "--------------------------------------------------------")
  @printf(io, "%-8s | %-18s | %-12s | %-10s\n", "Param", "Estimate ± Margin", "Rel. Err (%)", "Sensitivity")
  println(io, "--------------------------------------------------------")
  for i in eachindex(stats.params)
    abs_e = stats.params[i] - stats.ci_bounds[i][1]
    rel_e = abs(abs_e / stats.params[i]) * 100
    @printf(io, "%-8s | %8.3g ± %-7.2g | %-12.1f | %-10.1f\n", 
            stats.names[i], stats.params[i], abs_e, rel_e, stats.sensitivities[i])
  end
end

function Base.show(io::IO, ::MIME"text/latex", stats::ParameterStats)
    println(io, "Model Calibration Summary (\$R^2 = $(round(stats.r2, digits=4))\$) \\\\")
    println(io, "\\begin{tabular}{l c r r}")
    println(io, "\\hline")
    println(io, "Param & Estimate \$\\pm\$ Margin & Rel. Err (\\%) & Sensitivity \\\\")
    println(io, "\\hline")
    for i in eachindex(stats.params)
        abs_e = stats.params[i] - stats.ci_bounds[i][1]
        rel_e = abs(abs_e / stats.params[i]) * 100
        @printf(io, "\$%-8s\$ & \$%8.3g \\pm %-7.2g\$ & %12.1f & %10.1f \\\\\n", 
                stats.names[i], stats.params[i], abs_e, rel_e, stats.sensitivities[i])
    end
    println(io, "\\hline")
    println(io, "\\end{tabular}")
end

"""
    latex_string(param_stats)

Generate a string representing a summary with the [`ParameterStats`](@ref) to be copy-pasted into a ``\\LaTeX`` document.
"""
latex_string(x) = sprint(show, MIME("text/latex"), x)

"""
    latex_print(param_stats)

Print a string representing a summary with the [`ParameterStats`](@ref) to be copy-pasted into a ``\\LaTeX`` document.
"""
latex_print(x) = print(latex_string(x))


# --- Uncertainty Ensemble Generation ---

"""
    sample_parameters(res::CalibrationResult, n_samples=100)

Generates *n* parameter samples according to a normal multivariate distribution `N(params, Cov)`.
"""
function sample_parameters(res::CalibrationResult, n_samples::Int=100)
  cov_mat, _, _ = covariance_matrix(res)
  
  # Eigen-decomposition for stable Gaussian sampling: X = μ + V * √Δ * Z
  E = eigen(Symmetric(cov_mat))
  vals_clean = max.(E.values, 0.0) # Enforce positive semi-definiteness
  transform = E.vectors * Diagonal(sqrt.(vals_clean))
  
  n_params = length(res.params)
  samples = Matrix{Float64}(undef, n_params, n_samples)
  
  for k in 1:n_samples
    z = randn(n_params)
    samples[:, k] = res.params .+ (transform * z)
  end
  return samples
end
