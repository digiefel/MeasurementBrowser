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
    MeasurementBrowser._poll_scan_events!(ui)
    @test ui[:scan_state] == :scanning
    @test ui[:scan_progress][:processed_csv] == 2
    @test ui[:scan_progress][:current_path] == "fresh"

    # Request while active should queue path and cancel current scan.
    token = Base.Threads.Atomic{Bool}(false)
    current_events = Channel{NamedTuple}(8)
    ui[:scan_state] = :scanning
    ui[:scan_cancel_token] = token
    ui[:scan_events] = current_events
    MeasurementBrowser._queue_scan!(ui, "/tmp/queued_path"; persist_on_success=true)
    @test ui[:pending_scan_path] == "/tmp/queued_path"
    @test ui[:pending_scan_persist_on_success] == true
    @test token[] == true
    @test ui[:scan_state] == :canceling
    @test ui[:scan_events] === current_events

    # Same-path request while active should not queue/restart.
    ui2 = Dict{Symbol,Any}()
    MeasurementBrowser._init_scan_state!(ui2)
    ui2[:scan_state] = :scanning
    ui2[:scan_path] = "/tmp/same_path"
    ui2[:scan_persist_on_success] = false
    ui2[:scan_events] = Channel{NamedTuple}(8)
    token2 = Base.Threads.Atomic{Bool}(false)
    ui2[:scan_cancel_token] = token2
    MeasurementBrowser._queue_scan!(ui2, "/tmp/same_path"; persist_on_success=true)
    @test ui2[:pending_scan_path] === nothing
    @test token2[] == false
    @test ui2[:scan_state] == :scanning
    @test ui2[:scan_persist_on_success] == true

    # Same-path request should clear stale queued path.
    ui3 = Dict{Symbol,Any}()
    MeasurementBrowser._init_scan_state!(ui3)
    ui3[:scan_state] = :scanning
    ui3[:scan_path] = "/tmp/active"
    ui3[:pending_scan_path] = "/tmp/old_pending"
    ui3[:pending_scan_persist_on_success] = true
    token3 = Base.Threads.Atomic{Bool}(false)
    ui3[:scan_cancel_token] = token3
    MeasurementBrowser._queue_scan!(ui3, "/tmp/active")
    @test ui3[:pending_scan_path] === nothing
    @test ui3[:pending_scan_persist_on_success] == false
    @test token3[] == false

    # Force restart should queue and cancel even for same path.
    ui4 = Dict{Symbol,Any}()
    MeasurementBrowser._init_scan_state!(ui4)
    ui4[:scan_state] = :scanning
    ui4[:scan_path] = "/tmp/active"
    token4 = Base.Threads.Atomic{Bool}(false)
    ui4[:scan_cancel_token] = token4
    MeasurementBrowser._queue_scan!(ui4, "/tmp/active"; force_restart=true)
    @test ui4[:pending_scan_path] == "/tmp/active"
    @test ui4[:scan_state] == :canceling
    @test token4[] == true
end
