using MeasurementBrowser
using Test

@testset "bad registry" begin
    mktempdir() do dir
        registry = MeasurementBrowser.BadRegistry(
            Set(["RuO2test/A9/VI/D1", "RuO2test/A9/VI/D2"]),
            Set(["/tmp/m1.csv", "/tmp/m2.csv#cycle=2"]),
        )

        MeasurementBrowser.save_bad_registry(dir, registry)
        path = MeasurementBrowser.bad_registry_path(dir)
        @test isfile(path)
        @test read(path, String) == join([
            "device RuO2test/A9/VI/D1",
            "device RuO2test/A9/VI/D2",
            "measurement /tmp/m1.csv",
            "measurement /tmp/m2.csv#cycle=2",
            "",
        ], "\n")

        loaded = MeasurementBrowser.load_bad_registry(dir)
        @test loaded.devices == registry.devices
        @test loaded.measurements == registry.measurements
    end

    mktempdir() do dir
        path = joinpath(dir, "bad_measurements")
        write(path, "# comment\n\ndevice RuO2test/A9/VI/D1\nmeasurement /tmp/x.csv#split=A\n")
        loaded = MeasurementBrowser.load_bad_registry(dir)
        @test loaded.devices == Set(["RuO2test/A9/VI/D1"])
        @test loaded.measurements == Set(["/tmp/x.csv#split=A"])
    end

    mktempdir() do dir
        MeasurementBrowser.save_bad_registry(dir, MeasurementBrowser.BadRegistry())
        @test !isfile(joinpath(dir, "bad_measurements"))
    end

    mktempdir() do dir
        write(joinpath(dir, "bad_measurements"), "bogus\n")
        err = try
            MeasurementBrowser.load_bad_registry(dir)
            nothing
        catch thrown
            thrown
        end
        @test err isa MeasurementBrowser.BadRegistryParseError
        @test occursin("expected 'device <path>' or 'measurement <id>'", sprint(showerror, err))
    end
end
