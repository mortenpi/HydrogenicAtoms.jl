using Documenter, Plots
using HydrogenicAtoms

makedocs(
    sitename="HydrogenicAtoms",
    modules = [HydrogenicAtoms],
    pages = [
        "Home" => "index.md",
        "dirac.md"
    ]
)
