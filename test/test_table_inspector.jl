using MeasurementBrowser
using MeasurementBrowser: inspect_table
using Test

@testset "table inspector" begin
    mktempdir() do dir
        comma_path = joinpath(dir, "with_preamble.csv")
        write(
            comma_path,
            "# instrument note\n" *
            "sweep id: 3\n" *
            "time_s,current_A,voltage_V\n" *
            "0.0,1e-9,0.0\n" *
            "1.0,2e-9,0.1\n",
        )

        preview = inspect_table(comma_path)
        @test preview.delimiter == ','
        @test preview.header_row == 3
        @test preview.data_start_row == 4
        @test preview.columns == ["time_s", "current_A", "voltage_V"]
        @test preview.row_count == 2
        @test preview.preview_rows == 2

        tab_path = joinpath(dir, "no_header.tsv")
        write(tab_path, "1\t2\t3\n4\t5\t6\n")
        preview = inspect_table(tab_path; max_rows=1)
        @test preview.delimiter == '\t'
        @test preview.header_row === nothing
        @test preview.data_start_row == 1
        @test preview.row_count == 2
        @test preview.preview_rows == 1
        @test !isempty(preview.warnings)
    end
end
