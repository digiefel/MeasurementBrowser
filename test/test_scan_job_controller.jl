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

    # Incremental scan events should populate live state before completion.
    ui5 = Dict{Symbol,Any}()
    MeasurementBrowser._init_scan_state!(ui5)
    MeasurementBrowser._init_plot_state!(ui5)
    ui5[:active_scan_id] = 7
    ui5[:scan_events] = Channel{NamedTuple}(8)
    ui5[:scan_persist_on_success] = false

    indexed = IndexedCsvFile(
        "/tmp/file.csv",
        "/tmp/file.csv",
        "file.csv",
        nothing,
        Dict{Symbol,Any}(),
    )
    item = MeasurementItem(
        "/tmp/file.csv",
        "/tmp/file.csv",
        "/tmp/file.csv",
        :pund,
        ["Chip", "Block", "D1"],
        nothing,
        Dict{Symbol,Any}(:area_um2 => 10.0),
        Dict{Symbol,Any}(:voltage_V => 3.0),
        "FE PUND Chip_Block_D1",
    )

    put!(ui5[:scan_events], (
        kind=:scan_start,
        scan_id=7,
        path="/tmp",
        project=MeasurementBrowser.RUO2_PROJECT,
        has_device_metadata=true,
    ))
    put!(ui5[:scan_events], (kind=:items, scan_id=7, items=MeasurementItem[item]))
    put!(ui5[:scan_events], (
        kind=:progress,
        scan_id=7,
        progress=(
            phase=:scanning,
            total_csv=1,
            processed_csv=1,
            loaded_measurements=1,
            skipped_csv=0,
            current_path="/tmp/file.csv",
        ),
    ))
    put!(ui5[:scan_events], (kind=:result, scan_id=7, path="/tmp"))
    MeasurementBrowser._poll_scan_events!(ui5)
    @test ui5[:scan_state] == :done
    @test length(ui5[:all_measurements]) == 1
    @test ui5[:all_measurements][1].filepath == "/tmp/file.csv"
    @test ui5[:device_metadata_keys] == [:area_um2]

    # Canceling an incremental rescan should restore the previous visible state.
    ui6 = Dict{Symbol,Any}()
    MeasurementBrowser._init_scan_state!(ui6)
    MeasurementBrowser._init_plot_state!(ui6)
    old_measurement = MeasurementInfo(
        "old.csv",
        "/old/old.csv",
        "old",
        :pund,
        nothing,
        DeviceInfo(["Chip", "Block", "D1"], Dict{Symbol,Any}()),
        Dict{Symbol,Any}(:voltage_V => 2.0),
        nothing,
    )
    old_hierarchy = MeasurementHierarchy([old_measurement], "/old", false, MeasurementBrowser.RUO2_PROJECT, 0)
    ui6[:scan_hierarchy] = old_hierarchy
    ui6[:hierarchy_root] = old_hierarchy.root
    ui6[:all_measurements] = old_hierarchy.all_measurements
    ui6[:root_path] = "/old"
    ui6[:project] = MeasurementBrowser.RUO2_PROJECT
    ui6[:device_metadata_keys] = Symbol[]
    ui6[:active_scan_id] = 9
    ui6[:scan_events] = Channel{NamedTuple}(8)
    ui6[:scan_persist_on_success] = false
    MeasurementBrowser._capture_scan_restore_state!(ui6)

    put!(ui6[:scan_events], (
        kind=:scan_start,
        scan_id=9,
        path="/new",
        project=MeasurementBrowser.RUO2_PROJECT,
        has_device_metadata=false,
    ))
    put!(ui6[:scan_events], (kind=:items, scan_id=9, items=MeasurementItem[item]))
    put!(ui6[:scan_events], (kind=:canceled, scan_id=9))
    MeasurementBrowser._poll_scan_events!(ui6)
    @test ui6[:root_path] == "/old"
    @test length(ui6[:all_measurements]) == 1
    @test ui6[:all_measurements][1].filepath == "/old/old.csv"
end
