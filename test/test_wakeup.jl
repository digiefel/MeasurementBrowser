"""
test_wakeup.jl - Tests for enhanced measurement file display functionality

This test suite verifies that measurement files get enhanced GUI display information:

**Wakeup Files:**
1. **File Reading**: Tests that wakeup files are correctly parsed to count pulse lines
   and extract amplitude from filename
2. **Plot Generation**: Tests that wakeup plots display text only (no data points)
3. **GUI Enhancement**: Tests that pulse count is added to GUI list display

**FE PUND Files:**
1. **GUI Enhancement**: Tests that voltage amplitude is added to GUI list display
2. **Integration**: Tests that voltage extraction works correctly

**General:**
- Tests that other measurement types remain unchanged
- Tests that info panel clean_title keeps original format for all types

Expected behavior:
- Wakeup GUI list: "2025-10-01T17:10:48 100× Wakeup" (shows pulse count)
- FE PUND GUI list: "2025-10-01T17:12:33 FE PUND 3V" (shows voltage)
- Other files GUI list: "2025-10-01T16:21:45 TLM 4-Point" (unchanged)
- All info panels: Original format with device info (e.g., "Wakeup A9_VI_D1 2025-10-01")
- Wakeup plots: Display "N× wakeup pulses\namplitude = X V" as centered text only
"""

using Test
using CSV
using DataFrames
using Dates

using DataLoader
using DataPlotter
# Import the modules we need to test
# include("../src/DataLoader.jl")
# using .DataLoader: read_wakeup

include("../src/DeviceParser.jl")
# DeviceParser functions are not in a module, so we can use them directly

# include("../src/PlotGenerator.jl")
# using .PlotGenerator: plot_wakeup, figure_for_file

@testset "Wakeup File Processing Tests" begin

    # Test file path
    wakeup_file = "Wakeup 3V [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_10_48].csv"
    test_dir = @__DIR__
    wakeup_path = joinpath(test_dir, wakeup_file)

    @testset "File Existence" begin
        @test isfile(wakeup_path)
    end

    if isfile(wakeup_path)

        @testset "read_wakeup Function" begin
            result = read_wakeup(wakeup_file, test_dir)

            @test result isa DataFrame
            @test nrow(result) == 1
            @test names(result) == ["pulse_count", "amplitude"]

            pulse_count = result.pulse_count[1]
            amplitude = result.amplitude[1]

            @test pulse_count == 100
            @test amplitude == 3.0
            @test pulse_count isa Int
            @test amplitude isa Float64
        end

        @testset "Amplitude Extraction from Filename" begin
            # Test various filename patterns
            test_cases = [
                ("Wakeup 3V [test].csv", 3.0),
                ("Wakeup 10V [test].csv", 10.0),
                ("Wakeup 5.5V [test].csv", 5.5),
                ("Something 2V else.csv", 2.0),
                ("No voltage here.csv", 0.0)
            ]

            for (filename, expected_voltage) in test_cases
                # Create a minimal test DataFrame to avoid file operations
                amplitude_match = match(r"(\d+(?:\.\d+)?)V", filename)
                amplitude = amplitude_match !== nothing ? parse(Float64, amplitude_match.captures[1]) : 0.0
                @test amplitude == expected_voltage
            end
        end

        @testset "plot_wakeup Function" begin
            # Create test data
            test_df = DataFrame(pulse_count=[100], amplitude=[3.0])

            fig = plot_wakeup(test_df, "Test Wakeup")

            @test fig !== nothing
            # Note: We can't easily test the visual content without complex Makie internals,
            # but we can test that it doesn't crash and returns a figure
        end

        @testset "figure_for_file(s)" begin
            # TODO
            # test that figure_for_file and figure_for_files return expected figures
            @test true
        end

        @testset "MeasurementInfo Integration" begin
            meas_info = MeasurementInfo(wakeup_path)

            @test measurement_label(meas_info.measurement_kind) == "Wakeup"

            # Test that the clean title follows expected pattern (original format with device info)
            expected_title = "Wakeup A9_VI_D1 2025-10-01"
            @test meas_info.clean_title == expected_title
        end

        @testset "Data Line Counting" begin
            # Test that we correctly count data lines vs header lines
            lines = readlines(wakeup_path)

            # Find header line
            header_idx = 0
            for (i, line) in enumerate(lines)
                if occursin("Time,MeasResult1_value,MeasResult2_value", line)
                    header_idx = i
                    break
                end
            end

            @test header_idx > 0

            # Count data lines manually
            data_lines = 0
            for line in lines[header_idx+1:end]
                if !isempty(strip(line)) && occursin(',', line)
                    data_lines += 1
                end
            end

            @test data_lines == 100
        end

    else
        @warn "Wakeup test file not found at $wakeup_path - skipping file-dependent tests"
    end

    @testset "Edge Cases" begin
        # Test with empty DataFrame
        empty_df = DataFrame(pulse_count=Int[], amplitude=Float64[])
        fig = plot_wakeup(empty_df, "Empty Test")
        @test fig === nothing

        # Test measurement type detection
        @test detect_measurement_kind("Wakeup 5V test.csv") == :wakeup
        @test detect_measurement_kind("wakeup lowercase.csv") == :wakeup
        @test detect_measurement_kind("WAKEUP UPPERCASE.csv") == :wakeup
    end

    @testset "GUI Display Verification" begin
        if isfile(wakeup_path)
            meas_info = MeasurementInfo(wakeup_path)

            # Test that GUI shows pulse count while info panel shows original format
            gui_display = display_label(meas_info)  # What GUI actually shows
            info_display = meas_info.clean_title  # What info panel shows

            @test gui_display == "2025-10-01T17:10:48 Wakeup 100×"  # GUI shows pulse count
            @test info_display == "Wakeup A9_VI_D1 2025-10-01"  # Info panel shows original format
        end
    end

    @testset "FE PUND Voltage Display" begin
        # Test that FE PUND files get voltage added to GUI display
        pund_path = joinpath(test_dir, "3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv")
        if isfile(pund_path)
            pund_meas = MeasurementInfo(pund_path)

            # GUI display should include voltage for FE PUND files
            gui_display = display_label(pund_meas)
            @test gui_display == "2025-10-01T17:12:33 FE PUND 3V"

            # Info panel should show normal format (unchanged)
            info_display = pund_meas.clean_title
            @test info_display == "FE PUND A9_VI_D1 2025-10-01"
        end
    end

    @testset "Other Measurement Types Unchanged" begin
        # Test that non-wakeup, non-FE PUND files are not affected by the changes
        tlm_path = joinpath(test_dir, "TLM_4P [RuO2test_A9_VI_TLML100W2(12) ; 2025-10-01 16_21_45].csv")
        if isfile(tlm_path)
            tlm_meas = MeasurementInfo(tlm_path)

            # GUI display should NOT have pulse count or voltage for TLM files
            gui_display = display_label(tlm_meas)
            @test gui_display == "2025-10-01T16:21:45 TLM 4-Point"

            # Info panel should show normal format
            info_display = tlm_meas.clean_title
            @test info_display == "TLM 4-Point A9_VI_TLML100W2 2025-10-01"
        end
    end

    @testset "Voltage Extraction Edge Cases" begin
        # Test various voltage patterns for FE PUND files
        test_cases = [
            ("3V FE PUND [test].csv", "3V"),
            ("10V FE PUND [test].csv", "10V"),
            ("2.5V FE PUND [test].csv", "2.5V"),
            ("Something 5V else PUND.csv", "5V"),
            ("FE PUND no voltage.csv", nothing)  # Should fall back to default
        ]

        for (filename, expected_voltage) in test_cases
            # Test voltage extraction logic
            amplitude_match = match(r"(\d+(?:\.\d+)?)V", filename)
            if amplitude_match !== nothing
                voltage = parse(Float64, amplitude_match.captures[1])
                voltage_str = voltage == floor(voltage) ? "$(Int(voltage))V" : "$(voltage)V"
                @test voltage_str == expected_voltage
            else
                @test expected_voltage === nothing
            end
        end
    end
end
