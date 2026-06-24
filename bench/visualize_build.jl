# Visualize a build's event log + profile CSV (produced by open_workspace(...; event_log=, build_profile=)).
# Renders a multi-panel PNG and prints a numeric summary. Headless via CairoMakie.
#
#   julia --project=<env-with-CairoMakie,CSV,DataFrames> bench/visualize_build.jl \
#       [events.csv] [profile.csv] [out.png]
#
# The headline question: does per-operation cost grow with cache size (the O(N^2) signature)?

using CairoMakie, CSV, DataFrames, Statistics, Printf

events_path  = length(ARGS) >= 1 ? ARGS[1] : "/tmp/build_events.csv"
profile_path = length(ARGS) >= 2 ? ARGS[2] : "/tmp/build_profile.csv"
out_path     = length(ARGS) >= 3 ? ARGS[3] : "/tmp/build_analysis.png"

# ---- load -----------------------------------------------------------------------------------------
ev = CSV.read(events_path, DataFrame)
done = ev[ev.phase .== "end", :]                 # end rows carry dur_ms; begin rows are for crash detection
done.dur_ms = Float64.(done.dur_ms)
pr = CSV.read(profile_path, DataFrame)

ops = ["delete", "insert", "read"]
opcolor = Dict("delete" => :crimson, "insert" => :royalblue, "read" => :seagreen)
sub(op) = done[done.op .== op, :]

# Mean of y within equal-width bins of x — turns a noisy scatter into the trend line.
function binned(x, y; nbins=40)
    isempty(x) && return (Float64[], Float64[])
    lo, hi = extrema(x); hi == lo && return ([Float64(lo)], [mean(y)])
    edges = range(lo, hi; length=nbins + 1)
    cx, cy = Float64[], Float64[]
    for i in 1:nbins
        m = (x .>= edges[i]) .& (i == nbins ? x .<= edges[i+1] : x .< edges[i+1])
        any(m) || continue
        push!(cx, (edges[i] + edges[i+1]) / 2); push!(cy, mean(y[m]))
    end
    return (cx, cy)
end

# Least-squares slope/intercept of y ~ a + b*x.
function linfit(x, y)
    length(x) < 2 && return (0.0, 0.0)
    x̄, ȳ = mean(x), mean(y)
    vx = sum((x .- x̄).^2)
    vx == 0 && return (ȳ, 0.0)
    b = sum((x .- x̄) .* (y .- ȳ)) / vx
    return (ȳ - b * x̄, b)
end

# ---- figure ---------------------------------------------------------------------------------------
fig = Figure(size=(1700, 1900), fontsize=15)
Label(fig[0, 1:3], "Cache build analysis — $(nrow(done)) ops, $(round(pr.elapsed_s[end]; digits=0))s, " *
    "$(pr.analysis_done[end])/$(pr.analysis_total[end]) items analysed",
    fontsize=20, font=:bold)

# Row 1: the O(N^2) test — latency vs how full the cache is, per op.
for (col, op) in enumerate(["delete", "insert"])
    d = sub(op)
    ax = Axis(fig[1, col], title="$op latency vs cache size", xlabel="items already in cache",
        ylabel="dur (ms)")
    scatter!(ax, d.prior_items, d.dur_ms; color=(opcolor[op], 0.12), markersize=3)
    bx, by = binned(d.prior_items, d.dur_ms)
    lines!(ax, bx, by; color=opcolor[op], linewidth=3)
    a, b = linfit(Float64.(d.prior_items), d.dur_ms)
    lines!(ax, bx, a .+ b .* bx; color=:black, linestyle=:dash, linewidth=2)
    text!(ax, 0.03, 0.95; text=@sprintf("slope = %.4f ms / 1000 items", 1000b), space=:relative,
        align=(:left, :top), fontsize=13)
end
let op = "read", d = sub("read")
    ax = Axis(fig[1, 3], title="read latency vs cache size", xlabel="items already in cache",
        ylabel="dur (ms)")
    scatter!(ax, d.prior_items, d.dur_ms; color=(opcolor[op], 0.12), markersize=3)
    bx, by = binned(d.prior_items, d.dur_ms)
    lines!(ax, bx, by; color=opcolor[op], linewidth=3)
    a, b = linfit(Float64.(d.prior_items), d.dur_ms)
    text!(ax, 0.03, 0.95; text=@sprintf("slope = %.4f ms / 1000 items", 1000b), space=:relative,
        align=(:left, :top), fontsize=13)
end

# Row 2: where the time goes.
let ax = Axis(fig[2, 1], title="binned mean latency vs cache size (all ops)",
        xlabel="items already in cache", ylabel="mean dur (ms)")
    for op in ops
        d = sub(op); bx, by = binned(d.prior_items, d.dur_ms)
        lines!(ax, bx, by; color=opcolor[op], linewidth=3, label=op)
    end
    axislegend(ax; position=:lt)
end
let ax = Axis(fig[2, 2], title="cumulative time spent, by op", xlabel="build time (s)",
        ylabel="cumulative dur (s)")
    for op in ops
        d = sort(sub(op), :t_s)
        lines!(ax, d.t_s, cumsum(d.dur_ms) ./ 1000; color=opcolor[op], linewidth=3, label=op)
    end
    axislegend(ax; position=:lt)
end
let ax = Axis(fig[2, 3], title="total time by op (s)", xticks=(1:3, ops), ylabel="seconds")
    barplot!(ax, 1:3, [sum(sub(op).dur_ms) / 1000 for op in ops];
        color=[opcolor[op] for op in ops])
end

# Row 3: the profile time series.
let ax = Axis(fig[3, 1], title="cumulative write time vs elapsed", xlabel="elapsed (s)",
        ylabel="cumulative (s)")
    lines!(ax, pr.elapsed_s, pr.interp_write_s; label="interp", linewidth=3)
    lines!(ax, pr.elapsed_s, pr.processed_write_s; label="processed", linewidth=3)
    lines!(ax, pr.elapsed_s, pr.stats_write_s; label="stats", linewidth=3)
    lines!(ax, pr.elapsed_s, pr.writer_busy_s; label="writer busy", linewidth=3, color=:black)
    axislegend(ax; position=:lt)
end
let ax = Axis(fig[3, 2], title="progress vs elapsed", xlabel="elapsed (s)", ylabel="items / source items")
    lines!(ax, pr.elapsed_s, pr.scan_done; label="scan_done", linewidth=3)
    lines!(ax, pr.elapsed_s, pr.analysis_done; label="analysis_done", linewidth=3)
    hlines!(ax, [pr.analysis_total[end]]; color=:gray, linestyle=:dash)
    axislegend(ax; position=:lt)
end
let ax = Axis(fig[3, 3], title="analysis throughput vs cache fill", xlabel="items analysed",
        ylabel="items / s (instantaneous)")
    de = diff(pr.elapsed_s); da = diff(Float64.(pr.analysis_done))
    rate = da ./ de; x = pr.analysis_done[2:end]
    keep = isfinite.(rate)
    scatter!(ax, x[keep], rate[keep]; markersize=5, color=(:purple, 0.5))
    bx, by = binned(Float64.(x[keep]), rate[keep]; nbins=30)
    lines!(ax, bx, by; color=:purple, linewidth=3)
end

# Row 4: memory, writer concurrency, per-item write cost.
let ax = Axis(fig[4, 1], title="RSS vs elapsed", xlabel="elapsed (s)", ylabel="RSS (GiB)")
    lines!(ax, pr.elapsed_s, pr.rss_bytes ./ 1024^3; linewidth=3, color=:darkorange)
end
let ax = Axis(fig[4, 2], title="writers busy concurrently vs elapsed", xlabel="elapsed (s)",
        ylabel="d(writer_busy)/d(elapsed)")
    de = diff(pr.elapsed_s); db = diff(pr.writer_busy_s)
    lines!(ax, pr.elapsed_s[2:end], db ./ de; linewidth=2, color=:black)
end
let ax = Axis(fig[4, 3], title="per-item write cost vs cache fill", xlabel="items analysed",
        ylabel="processed write s / item (instantaneous)")
    de = diff(Float64.(pr.analysis_done)); dw = diff(pr.processed_write_s)
    keep = de .> 0; x = pr.analysis_done[2:end][keep]
    lines!(ax, x, (dw[keep] ./ de[keep]) .* 1000; linewidth=2, color=:teal)
    # ylabel says s but we scaled to ms:
    ax.ylabel = "processed write ms / item"
end

save(out_path, fig)

# ---- numeric summary ------------------------------------------------------------------------------
println("\n==================== BUILD ANALYSIS ====================")
@printf("ops: %d delete, %d insert, %d read   ·   build %.0fs\n",
    nrow(sub("delete")), nrow(sub("insert")), nrow(sub("read")), pr.elapsed_s[end])
for op in ops
    d = sub(op)
    a, b = linfit(Float64.(d.prior_items), d.dur_ms)
    e0, e1 = mean(d.dur_ms[d.prior_items .<= quantile(d.prior_items, 0.1)]),
             mean(d.dur_ms[d.prior_items .>= quantile(d.prior_items, 0.9)])
    @printf("  %-7s  median %.2fms  total %.1fs  slope %.4f ms/1k items  (early %.2fms → late %.2fms = %.1fx)\n",
        op, median(d.dur_ms), sum(d.dur_ms)/1000, 1000b, e0, e1, e1/max(e0, 1e-9))
end
@printf("writer busy %.0fs over %.0fs elapsed = %.1f writers avg · wait %.2fs\n",
    pr.writer_busy_s[end], pr.elapsed_s[end], pr.writer_busy_s[end]/pr.elapsed_s[end], pr.writer_wait_s[end])
@printf("analysis: %d/%d done (%.1f%%) · RSS %.1f GiB\n",
    pr.analysis_done[end], pr.analysis_total[end], 100*pr.analysis_done[end]/pr.analysis_total[end],
    pr.rss_bytes[end]/1024^3)
println("saved: $out_path")
