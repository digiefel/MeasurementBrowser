# Controlled GUI vs headless throughput experiments in one warmed session.
# Requires DataBrowserInteractiveProfile already loaded in Main.
module DataBrowserGuiThroughputDiag

using LinearAlgebra

const DIP = Main.DataBrowserInteractiveProfile

function _print_summary(label, result)
    println("==== ", label, " ====")
    println("outdir=", result.outdir)
    summary = joinpath(result.outdir, "summary.txt")
    println(read(summary, String))
end

function _warmup!()
    println("warmup headless 5s...")
    DIP.profile_scan!(mode=:headless, budget=5, profiler=:none, cache_mode=:fresh)
    return nothing
end

"""Run baseline then optional mitigations. Call once in a clean jmux session."""
function main(;
    warmup::Bool=true,
    budget_h::Real=10,
    budget_g::Real=20,
)
    println("nthreads default=", Threads.nthreads(:default),
            " interactive=", Threads.nthreads(:interactive))
    println("default tids=", collect(Threads.threadpooltids(:default)))
    println("interactive tids=", collect(Threads.threadpooltids(:interactive)))
    println("BLAS=", BLAS.get_num_threads())

    warmup && _warmup!()

    h0 = DIP.profile_scan!(mode=:headless, budget=budget_h, profiler=:none, cache_mode=:fresh)
    _print_summary("headless baseline", h0)

    g0 = DIP.profile_scan!(mode=:gui, budget=budget_g, profiler=:none, cache_mode=:fresh)
    _print_summary("gui baseline", g0)

    # Mitigation A: BLAS single-threaded to cut oversubscription with OpenGL.
    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    println("BLAS now=", BLAS.get_num_threads())
    g_blas = DIP.profile_scan!(mode=:gui, budget=budget_g, profiler=:none, cache_mode=:fresh)
    _print_summary("gui BLAS=1", g_blas)
    BLAS.set_num_threads(old_blas)

    # Mitigation B: dedicated libuv IO thread (Julia 1.12 Experimental).
    made = Base.Experimental.make_io_thread()
    println("make_io_thread -> ", made)
    g_io = DIP.profile_scan!(mode=:gui, budget=budget_g, profiler=:none, cache_mode=:fresh)
    _print_summary("gui after make_io_thread", g_io)

    h_io = DIP.profile_scan!(mode=:headless, budget=budget_h, profiler=:none, cache_mode=:fresh)
    _print_summary("headless after make_io_thread", h_io)

    return (;
        headless_baseline=h0.outdir,
        gui_baseline=g0.outdir,
        gui_blas1=g_blas.outdir,
        gui_make_io=g_io.outdir,
        headless_make_io=h_io.outdir,
    )
end

end
