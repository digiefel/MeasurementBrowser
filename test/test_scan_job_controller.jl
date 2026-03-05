using MeasurementBrowser
using Test

@testset "scan job controller guards" begin
    ui = Dict{Symbol,Any}()
    MeasurementBrowser._init_scan_state!(ui)

    # Stale scan_id events must be ignored.
    ui[:active_scan_id] = 2
    ui[:scan_events] = Channel{NamedTuple}(8)
    ui[:scan_cancel_token] = Base.Threads.Atomic{Bool}(false)
    put!(ui[:scan_events], (
        kind=:progress,
        scan_id=1,
        progress=(
            phase=:scanning,
            total_csv=10,
            processed_csv=9,
            loaded_measurements=9,
            skipped_csv=0,
            current_path="stale",
        ),
    ))
    put!(ui[:scan_events], (
        kind=:progress,
        scan_id=2,
        progress=(
            phase=:scanning,
            total_csv=10,
            processed_csv=2,
            loaded_measurements=2,
            skipped_csv=1,
            current_path="fresh",
        ),
    ))
    MeasurementBrowser._poll_scan_job!(ui)
    @test ui[:scan_state] == :scanning
    @test ui[:scan_progress][:processed_csv] == 2
    @test ui[:scan_progress][:current_path] == "fresh"

    # Request while active should queue path and cancel current scan.
    token = Base.Threads.Atomic{Bool}(false)
    current_events = Channel{NamedTuple}(8)
    ui[:scan_state] = :scanning
    ui[:scan_cancel_token] = token
    ui[:scan_events] = current_events
    MeasurementBrowser._request_scan!(ui, "/tmp/queued_path"; persist_on_success=true)
    @test ui[:pending_scan_path] == "/tmp/queued_path"
    @test ui[:pending_scan_persist_on_success] == true
    @test token[] == true
    @test ui[:scan_state] == :canceling
    @test ui[:scan_events] === current_events
end
