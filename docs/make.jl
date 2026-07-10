using Pkg

Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using HyperCalibration

# The documentation home page is always generated from the repository README.
cp(joinpath(@__DIR__, "..", "README.md"), joinpath(@__DIR__, "src", "index.md"), force=true)

makedocs(
  modules=[HyperCalibration],
  sitename="HyperCalibration.jl",
  pages=[
    "Home" => "index.md",
    "API reference" => "api.md",
  ],
)
