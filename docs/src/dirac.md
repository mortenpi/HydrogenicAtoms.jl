# Hydrogenic Dirac equation

```@example plots
using Plots
using AtomicLevels: RelativisticOrbital, @ro_str
using HydrogenicAtoms.DiracAnalytic: DiracSolution, PQ, WikiSolution

function plot_solution!(s::DiracSolution; rmax=15.0, length=500, kwargs...)
    rs = range(0, rmax, length=length)
    pqs = PQ.(s, rs)
    pq = [pqs[i][j] for i=1:length, j=1:2]
    # plot(layout=(2,1),
    #     title=["Large component" "Small component"],
    #     ylabel = ["P(r)" "Q(r)"],
    #     xlabel = "r",
    # )
    display(pq)

    plot!([rs rs], pq; layout=(2,1), kwargs...)
    # plot_fg!(1, -1; label="1s(1/2)")
    # plot_fg!(2, -1; label="2s(1/2)")
    # #plot_fg!(2,  1; label="2s(3/2)")
    # plot_fg!(3, -1; label="3s(1/2)")
    # xlims!(0, 15)
end
```

```@example plots
let Z = 1, rmax=20
    plot(layout=(2,1),
        title=["Large component" "Small component"],
        ylabel = ["P(r)" "Q(r)"],
        xlabel = "r",
    )
    for orb in [ro"1s", ro"2s", ro"3s", ro"4s"]
        plot_solution!(WikiSolution(Z, orb.n, orb.κ); rmax=rmax, label="$orb")
    end
    xlims!(0, rmax)
end
```

```@example plots
let Z = 100, rmax=0.5
    plot(layout=(2,1),
        title=["Large component" "Small component"],
        ylabel = ["P(r)" "Q(r)"],
        xlabel = "r",
    )
    for orb in [ro"2p-", ro"2p", ro"3p-", ro"3p"]
        plot_solution!(WikiSolution(Z, orb.n, orb.κ); rmax=rmax, label="$orb")
    end
    xlims!(0, rmax)
end
```
