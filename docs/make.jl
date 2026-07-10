using Pkg

Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using HyperCalibration

# The documentation home page is always generated from the repository README.
cp(joinpath(@__DIR__, "..", "README.md"), joinpath(@__DIR__, "src", "index.md"), force=true)

makedocs(
  sitename="HyperCalibration.jl",
  modules=[HyperCalibration],
  pages=[
    "Home" => "index.md",
    "API reference" => "api.md",
  ],
)

deploydocs(
  repo = "github.com/miguelmaso/HyperCalibration.jl.git",
  devbranch = "main"
)
