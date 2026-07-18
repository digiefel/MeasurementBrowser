# Post-fix remasure: fresh process, warmed, no prior make_io_thread call.
# Relies on open_browser -> _ensure_dedicated_io_thread!.
module DataBrowserPostFixRemeasure

using LinearAlgebra

const DIP = Main.DataBrowserInteractiveProfile

function main(; budget_h=10.0, budget_g=20.0)
    println("BLAS=", BLAS.get_num_threads())
    println("interactive tids=", collect(Threads.threadpooltids(:interactive)))
    println("default tids=", collect(Threads.threadpooltids(:default)))

    println("warmup...")
    DIP.profile_scan!(mode=:headless, budget=5, profiler=:none, cache_mode=:fresh)

    h = DIP.profile_scan!(mode=:headless, budget=budget_h, profiler=:none, cache_mode=:fresh)
    println("HEADLESS ", h.outdir)
    println(read(joinpath(h.outdir, "summary.txt"), String))

    g = DIP.profile_scan!(mode=:gui, budget=budget_g, profiler=:none, cache_mode=:fresh)
    println("GUI ", g.outdir)
    println(read(joinpath(g.outdir, "summary.txt"), String))

    return (; headless=h.outdir, gui=g.outdir)
end

end
